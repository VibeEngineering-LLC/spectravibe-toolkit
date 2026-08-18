"""
Anchor patterns for bootstrap energy calibration.

Loads `data/anchor_patterns.json` once on first access. Each pattern is a
group of library lines whose relative spacing is invariant under any
monotonic calibration. Used by the energy-calibration bootstrap module
(Phase 1.2) when the file's stored calibration fails the residual test.

Key API:
    list_patterns() -> list[dict]
    get_pattern(name) -> dict | None
    patterns_by_priority(max_priority=None) -> list[dict]
    multi_line_patterns() -> list[dict]
    single_line_anchors() -> list[dict]
"""

from __future__ import annotations

import json
from typing import Optional

from gamma.data import DATA_DIR

_ANCHORS_PATH = DATA_DIR / "anchor_patterns.json"

_CACHE: Optional[list] = None


def _load() -> list:
    """Return the list of anchor pattern dicts."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if not _ANCHORS_PATH.is_file():
        raise FileNotFoundError(
            f"anchor_patterns.json not found at {_ANCHORS_PATH}. "
            f"Expected: <skill_root>/data/anchor_patterns.json"
        )
    with _ANCHORS_PATH.open("r", encoding="utf-8") as f:
        data = json.load(f)
    _CACHE = list(data.get("patterns", []))
    return _CACHE


def reset_cache() -> None:
    global _CACHE
    _CACHE = None


def list_patterns() -> list:
    """All patterns."""
    return list(_load())


def get_pattern(name: str) -> Optional[dict]:
    """One pattern by name."""
    for p in _load():
        if p.get("name") == name:
            return p
    return None


def patterns_by_priority(max_priority: Optional[int] = None) -> list:
    """
    Return patterns sorted by priority (1=best, 5=worst).
    If max_priority given, exclude weaker ones.
    """
    patterns = sorted(_load(), key=lambda p: p.get("priority", 99))
    if max_priority is not None:
        patterns = [p for p in patterns if p.get("priority", 99) <= max_priority]
    return patterns


def multi_line_patterns() -> list:
    """Patterns with 2 or more lines (usable for shape calibration)."""
    return [p for p in _load() if len(p.get("lines", [])) >= 2]


def single_line_anchors() -> list:
    """Single-line patterns (offset/gain correction only, not curvature)."""
    return [p for p in _load() if p.get("single_line") is True]
