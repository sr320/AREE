"""One chat turn: Claude answers a visitor's question by calling AREE's query tools.

`run_turn` streams events as they happen, for the server to relay to the browser:

* {"type": "text", "text": ...}           a fragment of Claude's answer
* {"type": "tool_call", ...}              a query Claude is running (name, input)
* {"type": "tool_result", ...}            its result, exactly as Claude received it
* {"type": "round_reset"}                 discard text since the last tool call (a retried round)
* {"type": "notice", "text": ...}         something the visitor should know (a cut-off answer)
* {"type": "error", "code", "message"}    the turn could not be answered
* {"type": "done", ...}                   the turn's usage and estimated cost

Conversations are held in memory for an hour; Claude's responses are appended to
the history unchanged, as the API requires for its thinking blocks.
"""
from __future__ import annotations

import threading
import time
import uuid
from collections.abc import Iterator
from dataclasses import dataclass, field

import query
from chat.limits import cost_usd
from chat.tools import TOOLS, run_tool

FALLBACK_BETA = "server-side-fallback-2026-07-01"


@dataclass
class ChatConfig:
    model: str = "claude-opus-5-5"
    effort: str = "medium"
    # Per response. Answers are short; this leaves room for thinking and tool calls.
    max_tokens: int = 16000
    # Tool-calling rounds per question; the round after the last must answer.
    max_tool_rounds: int = 6
    max_turns: int = 10
    # Server-side refusal fallbacks (not available for Haiku models).
    fallbacks: bool = True


@dataclass
class Conversation:
    id: str
    messages: list = field(default_factory=list)
    turns: int = 0
    last_used: float = field(default_factory=time.monotonic)
    lock: threading.Lock = field(default_factory=threading.Lock)


class Conversations:
    """In-memory conversations, dropped an hour after their last turn."""

    def __init__(self, ttl_seconds: float = 3600, max_conversations: int = 2000, clock=time.monotonic):
        self.ttl, self.max, self.clock = ttl_seconds, max_conversations, clock
        self.items: dict[str, Conversation] = {}
        self.lock = threading.Lock()

    def get_or_create(self, conversation_id: str | None) -> Conversation:
        now = self.clock()
        with self.lock:
            for cid in [c for c, conv in self.items.items() if now - conv.last_used > self.ttl]:
                del self.items[cid]
            conv = self.items.get(conversation_id or "")
            if conv is None:
                if len(self.items) >= self.max:  # evict the least recently used
                    del self.items[min(self.items, key=lambda c: self.items[c].last_used)]
                conv = Conversation(id=uuid.uuid4().hex)
                self.items[conv.id] = conv
            conv.last_used = now
            return conv


def _request(config: ChatConfig, system: str, messages: list, final_round: bool) -> dict:
    params = {
        "model": config.model,
        "max_tokens": config.max_tokens,
        # Tools render before system, so this one breakpoint caches both; the
        # top-level marker caches the growing conversation.
        "system": [{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}],
        "tools": TOOLS,
        "messages": messages,
        "cache_control": {"type": "ephemeral"},
        "output_config": {"effort": config.effort},
    }
    if final_round:
        params["tool_choice"] = {"type": "none"}
    if config.fallbacks:
        params["betas"] = [FALLBACK_BETA]
        params["fallbacks"] = "default"
    return params


def _tool_result(block_id: str, content: str, is_error: bool) -> dict:
    result = {"type": "tool_result", "tool_use_id": block_id, "content": content}
    if is_error:
        result["is_error"] = True
    return result


def run_turn(client, conversation: Conversation, question: str, sources: query.Sources,
             config: ChatConfig, system: str, spend=lambda usd: None) -> Iterator[dict]:
    """Answer one question, streaming events.

    `spend` receives each response's estimated cost as soon as it is known, so a
    turn abandoned midway (the visitor closed the page) still counts. On an
    exception, or when the caller stops reading, the conversation is rolled back to
    before the question, so a failed turn leaves no half-built history.
    """
    if conversation.turns >= config.max_turns:
        yield {"type": "error", "code": "conversation_limit",
               "message": f"This conversation has reached {config.max_turns} questions. Start a new one."}
        return

    start = len(conversation.messages)
    conversation.messages.append({"role": "user", "content": question})
    usage = {"input_tokens": 0, "output_tokens": 0, "cache_creation_input_tokens": 0,
             "cache_read_input_tokens": 0}
    cost, tools_used, json_retries, model = 0.0, [], 0, config.model
    try:
        rounds = 0
        while True:
            final_round = rounds >= config.max_tool_rounds
            try:
                with client.beta.messages.stream(
                        **_request(config, system, conversation.messages, final_round)) as stream:
                    for event in stream:
                        if event.type == "text":
                            yield {"type": "text", "text": event.text}
                    response = stream.get_final_message()
            except ValueError:
                # A tool input the SDK could not parse at all: there is no tool_use
                # block to answer, so ask again (bounded).
                json_retries += 1
                if json_retries > 2:
                    raise
                yield {"type": "round_reset"}
                continue
            json_retries = 0
            model = getattr(response, "model", None) or config.model
            response_cost = cost_usd(model, response.usage)
            spend(response_cost)
            cost += response_cost
            for key in usage:
                usage[key] += getattr(response.usage, key, 0) or 0
            conversation.messages.append({"role": "assistant", "content": response.content})

            if response.stop_reason == "refusal":
                yield {"type": "error", "code": "refusal",
                       "message": "The assistant could not answer this question. Try rephrasing it, "
                                  "or use the Explore page."}
                break
            tool_uses = [b for b in response.content if b.type == "tool_use"]
            if not tool_uses:
                if response.stop_reason == "max_tokens":
                    yield {"type": "notice", "text": "The answer was cut off. Ask for a shorter answer "
                                                     "or a narrower question."}
                break

            results = []
            for block in tool_uses:
                tools_used.append(block.name)
                yield {"type": "tool_call", "id": block.id, "name": block.name, "input": block.input}
                if response.stop_reason == "max_tokens":
                    # A truncated input still parses as a partial object; don't run it.
                    content, is_error = ('{"error": "The tool input was cut off; ask for less at once."}',
                                         True)
                else:
                    content, is_error = run_tool(block.name, block.input, sources)
                yield {"type": "tool_result", "id": block.id, "name": block.name,
                       "is_error": is_error, "content": content}
                results.append(_tool_result(block.id, content, is_error))
            conversation.messages.append({"role": "user", "content": results})
            rounds += 1
    except BaseException:
        del conversation.messages[start:]
        raise

    conversation.turns += 1
    yield {"type": "done", "conversation_id": conversation.id, "model": model, "tools": tools_used,
           "usage": usage, "cost_usd": round(cost, 6)}
