"""The GitHub side of the autopilot: the status issue (with the heartbeat) and alert issues.

Plain HTTPS through urllib, not the gh CLI, which Klone compute nodes may not
have. urllib honours https_proxy, which Klone jobs need to reach the internet.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://api.github.com"
HEARTBEAT_RE = re.compile(r"<!-- aree-autopilot last_tick=(\S+) -->")


class GitHubError(RuntimeError):
    pass


def read_token(token_file: str | None = None) -> str | None:
    for var in ("GH_TOKEN", "GITHUB_TOKEN"):
        if os.environ.get(var):
            return os.environ[var]
    if token_file and Path(token_file).exists():
        return Path(token_file).read_text().strip() or None
    return None


class GitHub:
    def __init__(self, repo: str, token: str):
        self.repo = repo
        self.token = token

    def request(self, method: str, path: str, body: dict | None = None):
        url = path if path.startswith("http") else f"{API}{path}"
        data = json.dumps(body).encode() if body is not None else None
        req = urllib.request.Request(url, data=data, method=method, headers={
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "aree-autopilot",
        })
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                raw = resp.read()
        except urllib.error.HTTPError as exc:
            raise GitHubError(f"{method} {path}: HTTP {exc.code} {exc.read()[:300]!r}") from exc
        return json.loads(raw) if raw else None

    def ensure_label(self, name: str, color: str) -> None:
        try:
            self.request("POST", f"/repos/{self.repo}/labels", {"name": name, "color": color})
        except GitHubError as exc:
            if "HTTP 422" not in str(exc):  # 422: the label already exists
                raise

    def open_issues(self, label: str) -> list[dict]:
        encoded = urllib.parse.quote(label, safe="")
        issues = self.request("GET", f"/repos/{self.repo}/issues?state=open&per_page=100&labels={encoded}") or []
        return [i for i in issues if "pull_request" not in i]

    def create_issue(self, title: str, body: str, labels: list[str]) -> dict:
        return self.request("POST", f"/repos/{self.repo}/issues",
                            {"title": title, "body": body, "labels": labels})

    def update_issue_body(self, number: int, body: str) -> None:
        self.request("PATCH", f"/repos/{self.repo}/issues/{number}", {"body": body})


def heartbeat_marker(now: dt.datetime) -> str:
    return f"<!-- aree-autopilot last_tick={now.strftime('%Y-%m-%dT%H:%M:%SZ')} -->"


def parse_heartbeat(body: str) -> dt.datetime | None:
    m = HEARTBEAT_RE.search(body or "")
    if not m:
        return None
    return dt.datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)


def publish_status(gh, cfg: dict, body: str) -> int:
    """Edit the status issue's body (created on first use). Body edits notify nobody."""
    issues = gh.open_issues(cfg["status_label"])
    if issues:
        number = issues[0]["number"]
        gh.update_issue_body(number, body)
        return number
    gh.ensure_label(cfg["status_label"], "0e8a16")
    return gh.create_issue(cfg["status_issue_title"], body, [cfg["status_label"]])["number"]


def status_paused(gh, cfg: dict) -> bool:
    """True when the open status issue carries the paused label."""
    for issue in gh.open_issues(cfg["status_label"])[:1]:
        names = {lab["name"] if isinstance(lab, dict) else lab for lab in issue.get("labels", [])}
        return cfg["paused_label"] in names
    return False


def open_alert(gh, cfg: dict, key: str, body: str, existing_titles: set[str] | None = None) -> bool:
    """Open `[autopilot-alert] <key>` unless an open alert with that title exists.

    Returns True when a new issue was opened. A person closes the issue once the
    problem is handled; the next occurrence then opens a new one. Pass
    `existing_titles` (a mutable set) to avoid re-fetching open alert issues on
    each call in a batch.
    """
    title = f"[autopilot-alert] {key}"
    if existing_titles is None:
        existing_titles = {issue["title"] for issue in gh.open_issues(cfg["alert_label"])}
    if title in existing_titles:
        return False
    gh.ensure_label(cfg["alert_label"], "d93f0b")
    gh.create_issue(title, body, [cfg["alert_label"]])
    existing_titles.add(title)
    return True


def check_heartbeat(gh, cfg: dict, max_age_hours: float, now: dt.datetime | None = None) -> tuple[str, str]:
    """Dead-man switch, run off-cluster (GitHub Actions): is the status issue's
    last_tick recent? Returns (outcome, message); outcome is "ok", "absent"
    (no open status issue: the autopilot is not deployed, nothing to check) or
    "stale" (an alert was opened, or one is already open)."""
    now = now or dt.datetime.now(dt.timezone.utc)
    issues = gh.open_issues(cfg["status_label"])
    if not issues:
        return "absent", "no open status issue; the autopilot is not deployed"
    last = parse_heartbeat(issues[0].get("body") or "")
    age = (now - last).total_seconds() / 3600 if last else None
    if age is not None and age <= max_age_hours:
        return "ok", f"last tick {age:.1f} h ago"
    seen = f"{age:.1f} h ago" if age is not None else "never"
    open_alert(gh, cfg, "heartbeat missing",
               f"The Klone autopilot last reported {seen} (limit {max_age_hours:g} h): the watchdog "
               f"job on Klone is not running. Check `squeue -u $USER -n aree-autopilot` and the "
               f"scrontab / watchdog log on Klone, then restart it "
               f"(docs/running_on_klone.md#autopilot). Status: #{issues[0]['number']}.\n")
    return "stale", f"last tick {seen}"
