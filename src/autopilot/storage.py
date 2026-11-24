"""Free space on scrubbed.

statvfs reports the filesystem as a whole. Whether that tracks the scrubbed
fileset closely enough, or `hyakstorage` must be parsed instead, is an open
Phase 0 question (docs/klone_autopilot_plan.md).
"""
from __future__ import annotations

import os

TB = 1e12


def free_tb(path: str) -> float | None:
    try:
        st = os.statvfs(path)
    except (OSError, AttributeError):
        return None
    return st.f_bavail * st.f_frsize / TB
