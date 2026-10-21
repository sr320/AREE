"""Cost and abuse limits for the public chat, held in memory.

Limits are per server process. Run the service with a single instance (Cloud Run
--max-instances=1) for them to be exact; with N instances each limit is up to N
times looser. The API workspace's monthly spend limit is the hard backstop.
"""
from __future__ import annotations

import datetime as _dt
import threading
import time
from collections import defaultdict, deque

# USD per million tokens: input, output, 5-minute cache write, cache read.
PRICES = {
    "claude-opus-5-5": (4.00, 20.00, 5.00, 0.20),
    "claude-sonnet-5-5": (2.00, 10.00, 2.50, 0.20),
    "claude-haiku-5-5": (0.10, 0.50, 0.125, 0.01),
    "claude-opus-5": (5.00, 25.00, 6.25, 0.50),
    "claude-opus-4-8": (5.00, 25.00, 6.25, 0.50),
}
# A model missing from the table (a fallback, a new release) is priced at the
# most expensive known rate, so the budget errs toward stopping early.
_UNKNOWN_PRICE = max(PRICES.values())


def cost_usd(model: str, usage) -> float:
    """Estimated cost of one response from its reported usage."""
    price = PRICES.get(model, _UNKNOWN_PRICE)
    tokens = (
        getattr(usage, "input_tokens", 0) or 0,
        getattr(usage, "output_tokens", 0) or 0,
        getattr(usage, "cache_creation_input_tokens", 0) or 0,
        getattr(usage, "cache_read_input_tokens", 0) or 0,
    )
    return sum(n * p for n, p in zip(tokens, price)) / 1_000_000


class RateLimiter:
    """Questions per client over the last hour and the last day."""

    def __init__(self, per_hour: int, per_day: int, clock=time.monotonic):
        self.limits = ((3600, per_hour), (86400, per_day))
        self.clock = clock
        self.hits: dict[str, deque] = defaultdict(deque)
        self.lock = threading.Lock()

    def check(self, key: str) -> str | None:
        """Record a question from `key`; return why it is refused, or None if allowed."""
        now = self.clock()
        with self.lock:
            hits = self.hits[key]
            while hits and now - hits[0] > self.limits[-1][0]:
                hits.popleft()
            for window, limit in self.limits:
                if sum(1 for t in hits if now - t <= window) >= limit:
                    span = "hour" if window == 3600 else "day"
                    return f"You have asked {limit} questions in the last {span}, the limit for this chat."
            hits.append(now)
            return None


class DailyBudget:
    """Estimated spend since UTC midnight, against a cap."""

    def __init__(self, usd_per_day: float, today=lambda: _dt.datetime.now(_dt.timezone.utc).date()):
        self.cap = usd_per_day
        self.today = today
        self.day = today()
        self.spent = 0.0
        self.lock = threading.Lock()

    def _roll(self):
        if self.today() != self.day:
            self.day, self.spent = self.today(), 0.0

    def exhausted(self) -> bool:
        with self.lock:
            self._roll()
            return self.spent >= self.cap

    def add(self, usd: float) -> None:
        with self.lock:
            self._roll()
            self.spent += usd
