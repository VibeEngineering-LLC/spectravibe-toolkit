"""
Q3 defensive doc-test: assert active-state carve-out language is absent from
living docs (KNOWN_AND_FIXED_ISSUES.md, handoff_ru.md, handoff.md).

Historical mentions in 1_Version/**/RELEASE_NOTES.md and README.md are
EXPLICITLY ALLOWED — those are factual records, not active-state language.

Background:
  DEEP-06 Step B-mini closed with v1.27.0 release (2026-06-06).
  The carve-out on compute.py / staged_pipeline.py is HISTORICAL as of v1.27.0+.
  Active-state phrases in living docs would be misleading and confusing to future
  sessions — this test guards against accidental re-introduction.

Note: No pre-fix red state was possible because all living docs were already clean
before this test was written (carve-out language was never committed to these docs
in active-state form). This is a preventive/defensive test only.
"""

import re
from pathlib import Path

# Project root is two levels up from tests/docs/
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Living docs to check (MUST NOT contain active-state phrases)
LIVING_DOCS = [
    PROJECT_ROOT / "KNOWN_AND_FIXED_ISSUES.md",
    PROJECT_ROOT / "handoff_ru.md",
    PROJECT_ROOT / "handoff.md",
]

# Active-state phrases that must NOT appear in living docs
# These indicate the carve-out is still in effect — which is false post-v1.27.0
ACTIVE_STATE_PHRASES = [
    "do not touch compute.py",
    "carve-out active",
    "carve-out in effect",
    "HARD-LOCK orchestrator on compute.py",
    "HARD-LOCK on staged_pipeline.py",
    "do not edit compute.py",
    "do not edit staged_pipeline.py",
]

# Paths that are explicitly ALLOWED to contain carve-out mentions
# (historical/archival records — do not touch these)
ALLOWED_PATHS_PATTERNS = [
    r"1_Version[/\\].+RELEASE_NOTES\.md",
    r"README\.md",
]


def _is_allowed_path(path: Path) -> bool:
    """Return True if path is in the allowed-historical-mentions list."""
    path_str = str(path)
    return any(re.search(pattern, path_str) for pattern in ALLOWED_PATHS_PATTERNS)


def test_no_active_carve_out_in_living_docs():
    """
    Assert that active-state carve-out language is absent from living docs.

    This test MUST:
    - Pass on post-fix state (v1.27.0+)
    - Fail if any active-state phrase is introduced into the listed living docs

    Historical mentions in 1_Version/ and README.md are out of scope for this test
    (see ALLOWED_PATHS_PATTERNS).
    """
    violations = []

    for doc_path in LIVING_DOCS:
        if not doc_path.exists():
            # If a doc doesn't exist yet, it can't contain the phrase — skip
            continue

        content = doc_path.read_text(encoding="utf-8")
        content_lower = content.lower()

        for phrase in ACTIVE_STATE_PHRASES:
            phrase_lower = phrase.lower()
            if phrase_lower in content_lower:
                # Find line numbers for the violation
                lines = content.splitlines()
                matched_lines = [
                    (i + 1, line.strip())
                    for i, line in enumerate(lines)
                    if phrase_lower in line.lower()
                ]
                for lineno, line_text in matched_lines:
                    violations.append(
                        f"{doc_path.name}:{lineno}: found active-state phrase "
                        f"{phrase!r} → {line_text!r}"
                    )

    assert not violations, (
        "Active-state carve-out language found in living docs. "
        "These docs should use past-tense framing (carve-out lifted by v1.27.0). "
        "Violations:\n" + "\n".join(violations)
    )


def test_allowed_paths_excluded_from_check():
    """
    Sanity: confirm that RELEASE_NOTES.md and README.md are not in LIVING_DOCS
    (they are allowed to have historical carve-out mentions and must not be checked).
    """
    living_doc_strs = [str(p) for p in LIVING_DOCS]
    for allowed_pattern in ALLOWED_PATHS_PATTERNS:
        for doc_str in living_doc_strs:
            assert not re.search(allowed_pattern, doc_str), (
                f"Allowed-historical path {doc_str!r} must not appear in LIVING_DOCS. "
                "It would cause false positives on legitimate historical records."
            )
