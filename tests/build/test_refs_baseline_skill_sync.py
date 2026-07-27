"""
tests/build/test_refs_baseline_skill_sync.py
DEEP-04 + DEEP-05 — SKILL.md doc/code drift gates.

test_refs_baseline_skill_matches_build   (DEEP-04)
    Cross-checks SKILL.md refs baseline claim against build.py baseline_refs
    set parsed via AST.  RED before SKILL.md fix, GREEN after.

test_f326_skill_documents_bg_only_exception  (DEEP-05)
    Verifies SKILL.md F-326 paragraph does not claim "always rendered" without
    documenting the background_only suppression (F-UX-04 / 2026-06-04).
    RED before SKILL.md fix, GREEN after.

Red-without-fix evidence (recorded at time of writing):
  DEEP-04: SKILL.md line 193 claims {1, 2, 7, 12, 19};
           build.py:364 has baseline_refs = {2, 7, 12, 19, 24}
           (ref 1 removed per F-337.4/v1.18.19.1; ref 24 added per BUG-16).
  DEEP-05: SKILL.md line 194 says "Section is always rendered";
           build.py:1219 + 1256 skip section when
           measurement_environment == "background_only" (F-UX-04/2026-06-04).
"""

import ast
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _extract_baseline_refs_from_build() -> frozenset:
    """Parse baseline_refs = {...} from build.py using AST."""
    src = (REPO_ROOT / "scripts" / "gamma" / "reporting" / "build.py").read_text(
        encoding="utf-8"
    )
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(t, ast.Name) and t.id == "baseline_refs"
                for t in node.targets
            )
            and isinstance(node.value, ast.Set)
        ):
            return frozenset(
                elt.value  # ast.Constant.value (replaces deprecated .n in Python 3.8+)
                for elt in node.value.elts
                if isinstance(elt, ast.Constant)
            )
    raise AssertionError("baseline_refs assignment not found in build.py")


def _extract_skill_baseline_refs() -> frozenset:
    """Parse 'refs {N, ...} always included' from SKILL.md."""
    text = (REPO_ROOT / "SKILL.md").read_text(encoding="utf-8")
    m = re.search(r"refs\s*\{([0-9,\s]+)\}\s*always included", text)
    assert m, "SKILL.md: 'refs {...} always included' pattern not found"
    return frozenset(int(x.strip()) for x in m.group(1).split(",") if x.strip())


def test_refs_baseline_skill_matches_build():
    """SKILL.md refs baseline claim must match build.py baseline_refs set.

    DEEP-04: Verifies that the set of references documented in SKILL.md as
    'always included' matches the actual baseline_refs = {...} assignment in
    build.py.  The code is authoritative (per F-337.4 + BUG-16 inline comments).
    """
    build_refs = _extract_baseline_refs_from_build()
    skill_refs = _extract_skill_baseline_refs()
    assert build_refs == skill_refs, (
        f"SKILL.md claims baseline refs {sorted(skill_refs)} "
        f"but build.py uses {sorted(build_refs)}. "
        "Update SKILL.md to match the code. "
        "(ref 1 removed per F-337.4/v1.18.19.1; ref 24 added per BUG-16/ГОСТ 26874-86)"
    )


def test_f326_skill_documents_bg_only_exception():
    """SKILL.md must document that F-326 is NOT rendered for background_only spectra.

    DEEP-05: Verifies two conditions after F-UX-04 correction (2026-06-04):
      1. SKILL.md mentions the background_only suppression condition somewhere.
      2. The F-326 paragraph specifically does NOT claim 'always rendered'
         without also documenting the background_only exception.
    """
    text = (REPO_ROOT / "SKILL.md").read_text(encoding="utf-8")

    # Condition 1: SKILL.md must mention background_only suppression for F-326
    assert "background_only" in text, (
        "SKILL.md does not document F-326 background_only suppression (F-UX-04). "
        "Update SKILL.md: F-326 section is skipped for background_only spectra "
        "(build.py:1219 + 1256 conditional, F-UX-04 / 2026-06-04)."
    )

    # Condition 2: F-326 paragraph must not say 'always rendered' without qualification
    f326_match = re.search(r"F-326.*?(?=\n-\s+\*\*F-|\Z)", text, re.DOTALL)
    if f326_match:
        f326_text = f326_match.group(0)
        assert "always rendered" not in f326_text or "background_only" in f326_text, (
            "F-326 SKILL.md paragraph says 'always rendered' without mentioning "
            "the background_only exception. Update to reflect F-UX-04 / 2026-06-04: "
            "section is suppressed for measurement_environment == 'background_only'."
        )
