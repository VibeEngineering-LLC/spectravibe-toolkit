"""
tests/build/test_requirements_lock_consistency.py
DEEP-02 · REL-3 — Lockfile pin consistency gate.

Asserts:
  test_lock_file_exists_and_nonempty       — lock file present, non-empty.
  test_lock_has_header_comment             — first non-blank line is a comment.
  test_lock_uses_strict_equals             — every non-comment, non-blank line
                                             is a strict == pin (no >=, ~=, <).
  test_lock_is_superset_of_manifest        — every top-level package in
                                             requirements.txt appears in lock.
  test_lock_versions_satisfy_manifest_bounds — every lock pin SATISFIES the
                                             manifest's own version constraints
                                             (no <upper-bound> that the validated
                                             lock violates). DEEP-02 follow-up:
                                             a lock that contradicts the manifest
                                             defeats reproducibility — `pip install
                                             -r requirements.txt` would resolve a
                                             DIFFERENT environment than the one the
                                             snapshots were validated on.

Red-without-fix: delete or empty requirements.lock → first four tests fail.
                 manifest upper bound that the lock violates → fifth test fails.
"""

import pathlib
import re

import pytest
from packaging.specifiers import SpecifierSet
from packaging.version import Version

PROJECT_ROOT = pathlib.Path(__file__).parent.parent.parent
LOCK_FILE = PROJECT_ROOT / "scripts" / "requirements.lock"
MANIFEST_FILE = PROJECT_ROOT / "scripts" / "requirements.txt"

# Regex: name==version, allowing post/local segments like 2.9.0.post0
_STRICT_PIN_RE = re.compile(
    r"^[A-Za-z0-9_\-\.]+==[0-9]+\.[0-9]+(\.[0-9A-Za-z\.\-]+)*$"
)

# Regex to strip version specifiers from a requirement line
_PKG_NAME_RE = re.compile(r"^([A-Za-z0-9_\-\.]+)")

# Regex: split a manifest line into (name, specifier-tail), e.g.
#   "scipy>=1.10,<1.18"  -> ("scipy", ">=1.10,<1.18")
#   "defusedxml>=0.7.1"  -> ("defusedxml", ">=0.7.1")
_MANIFEST_SPEC_RE = re.compile(
    r"^([A-Za-z0-9_\-\.]+)\s*([<>=!~].*)?$"
)


def _lock_lines():
    """Return all non-comment, non-blank lines from the lock file."""
    return [
        line.strip()
        for line in LOCK_FILE.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def _manifest_top_level_names():
    """
    Return the set of normalised package names declared in requirements.txt.
    Skips comment lines, blank lines, and commented-out packages (# foo>=...).
    """
    names = set()
    for line in MANIFEST_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = _PKG_NAME_RE.match(line)
        if m:
            # Normalise: lower-case, hyphens → underscores
            names.add(m.group(1).lower().replace("-", "_"))
    return names


def _lock_names():
    """Return the set of normalised package names present in the lock file."""
    names = set()
    for line in _lock_lines():
        m = _PKG_NAME_RE.match(line)
        if m:
            names.add(m.group(1).lower().replace("-", "_"))
    return names


def _norm(name: str) -> str:
    """Normalise a package name: lower-case, hyphens → underscores."""
    return name.lower().replace("-", "_")


def _manifest_specifiers():
    """
    Return {normalised_name: specifier_string} for every ACTIVE (non-commented)
    requirement line in requirements.txt that carries a version specifier.
    Lines without a specifier (bare package name) are skipped — there is no
    bound to check.
    """
    specs = {}
    for line in MANIFEST_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Drop any trailing inline comment.
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        m = _MANIFEST_SPEC_RE.match(line)
        if m and m.group(2):
            specs[_norm(m.group(1))] = m.group(2).strip()
    return specs


def _lock_versions():
    """Return {normalised_name: version_string} from the lock file."""
    versions = {}
    for line in _lock_lines():
        if "==" not in line:
            continue
        name, _, ver = line.partition("==")
        versions[_norm(name.strip())] = ver.strip()
    return versions


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_lock_file_exists_and_nonempty():
    """Lock file must exist and contain at least one pin entry."""
    assert LOCK_FILE.exists(), (
        f"requirements.lock not found at {LOCK_FILE}. "
        "Generate it via: pip install -r scripts/requirements.txt && pip freeze > ... "
        "(see DEEP-02 notes for the BFS filter procedure)."
    )
    content = LOCK_FILE.read_text(encoding="utf-8")
    assert content.strip(), "requirements.lock is empty."
    # Must have at least one pin line
    pin_lines = _lock_lines()
    assert len(pin_lines) > 0, "requirements.lock has no non-comment entries."


def test_lock_has_header_comment():
    """Lock file must start with a header comment (first non-blank line is '#...')."""
    assert LOCK_FILE.exists(), "requirements.lock not found."
    for line in LOCK_FILE.read_text(encoding="utf-8").splitlines():
        if line.strip():
            assert line.strip().startswith("#"), (
                f"First non-blank line of requirements.lock is not a comment: {line!r}. "
                "Lock file must start with a header comment explaining its purpose."
            )
            break
    else:
        pytest.fail("requirements.lock has no non-blank lines at all.")


def test_lock_uses_strict_equals():
    """
    Every non-comment, non-blank line in requirements.lock must be a strict
    == pin (e.g. 'scipy==1.17.1').  No >=, ~=, <, >, or bare package names.
    Violations cause snapshot tests to drift across environments.
    """
    assert LOCK_FILE.exists(), "requirements.lock not found."
    bad_lines = []
    for line in _lock_lines():
        if not _STRICT_PIN_RE.match(line):
            bad_lines.append(line)
    assert not bad_lines, (
        f"requirements.lock contains {len(bad_lines)} non-strict-pin line(s):\n"
        + "\n".join(f"  {l}" for l in bad_lines)
        + "\nAll entries must be name==version (exact == only)."
    )


def test_lock_is_superset_of_manifest():
    """
    Every top-level package declared in requirements.txt must appear in
    requirements.lock.  The lock is the resolved superset of the manifest.
    """
    assert LOCK_FILE.exists(), "requirements.lock not found."
    manifest_names = _manifest_top_level_names()
    lock_names = _lock_names()
    missing = manifest_names - lock_names
    assert not missing, (
        f"The following top-level packages from requirements.txt are missing "
        f"from requirements.lock: {sorted(missing)}. "
        "Regenerate the lock file to include all transitive deps."
    )


def test_lock_versions_satisfy_manifest_bounds():
    """
    DEEP-02 follow-up — every pinned version in requirements.lock MUST satisfy
    the version constraints declared for that package in requirements.txt.

    Why this matters: the lock and the manifest are two views of the SAME
    dependency contract. If the manifest says ``scipy>=1.10,<1.15`` but the
    lock pins ``scipy==1.17.1``, then `pip install -r requirements.txt`
    (the developer / `pip install -e .` path) resolves a DIFFERENT environment
    than `pip install -r requirements.lock` (the CI path the snapshots are
    validated against). That is the exact reproducibility hole DEEP-02 set out
    to close — an upper bound the validated lock violates is a trap, not a guard.

    The correct invariant: upper bounds must sit JUST ABOVE the validated minor
    (guarding against an unvalidated future bump) and therefore the lock pin
    must always satisfy them.
    """
    specs = _manifest_specifiers()
    versions = _lock_versions()

    violations = []
    for name, spec_str in specs.items():
        lock_ver = versions.get(name)
        if lock_ver is None:
            # superset test already covers presence; skip here.
            continue
        spec = SpecifierSet(spec_str)
        # prereleases=True so e.g. post-releases are not silently excluded.
        if not spec.contains(Version(lock_ver), prereleases=True):
            violations.append(
                f"  {name}: lock pins =={lock_ver} but manifest requires "
                f"'{spec_str}' (lock VIOLATES manifest)"
            )

    assert not violations, (
        "requirements.lock pins versions that violate requirements.txt bounds:\n"
        + "\n".join(violations)
        + "\n\nFix: widen the manifest upper bound to sit just above the validated "
        "minor (so the lock pin satisfies it), OR regenerate the lock inside the "
        "manifest's declared range. A manifest bound the validated lock violates "
        "defeats DEEP-02 reproducibility."
    )
