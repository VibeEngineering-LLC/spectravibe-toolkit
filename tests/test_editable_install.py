"""
tests/test_editable_install.py — Regression for task #99 (editable install bug).

Guards against:
1. pyproject.toml parse failure (structurally invalid TOML).
2. Missing [build-system] table — pip install -e . would produce
   "Multiple top-level packages discovered" or BackendUnavailable errors.
3. Missing [project] table — no installable metadata.
4. Smoke import of the gamma package itself (verifies [tool.setuptools.packages.find]
   is configured to include scripts/gamma/).

Run in any env that has gamma installed (editable or otherwise); CI enforces
cold-venv editable install before executing the test suite.
"""
from __future__ import annotations

import importlib
import sys
import tomllib
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"


# ──────────────────────────────────────────────────────────────────────────────
# 1. TOML parse (catches broken syntax before pip even sees the file)
# ──────────────────────────────────────────────────────────────────────────────


def test_pyproject_toml_is_valid_toml():
    """pyproject.toml must parse without error (tomllib stdlib, Python 3.11+)."""
    assert PYPROJECT.exists(), f"pyproject.toml not found at {PYPROJECT}"
    with PYPROJECT.open("rb") as fh:
        data = tomllib.load(fh)
    assert isinstance(data, dict), "Expected a TOML table at top level"


# ──────────────────────────────────────────────────────────────────────────────
# 2. Required PEP 517 tables
# ──────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def pyproject_data() -> dict:
    with PYPROJECT.open("rb") as fh:
        return tomllib.load(fh)


def test_build_system_table_present(pyproject_data):
    """[build-system] must exist — absence causes pip install -e . to fail
    with 'Multiple top-level packages discovered' in flat-layout repos."""
    assert "build-system" in pyproject_data, (
        "pyproject.toml is missing [build-system] — "
        "pip install -e . will fail on cold venv"
    )


def test_build_system_backend(pyproject_data):
    """build-backend must be 'setuptools.build_meta' (the canonical, portable path).
    Using 'setuptools.backends.legacy:build' fails on older bootstrapped setuptools
    in cold venvs (BackendUnavailable error observed 2026-06-06)."""
    bs = pyproject_data.get("build-system", {})
    backend = bs.get("build-backend", "")
    assert backend == "setuptools.build_meta", (
        f"Expected build-backend='setuptools.build_meta', got {backend!r}"
    )


def test_project_table_present(pyproject_data):
    """[project] table must exist with at minimum 'name' and 'version'."""
    assert "project" in pyproject_data, "pyproject.toml is missing [project] table"
    proj = pyproject_data["project"]
    assert "name" in proj, "[project] must contain 'name'"
    assert "version" in proj, "[project] must contain 'version'"


def test_setuptools_packages_find_configured(pyproject_data):
    """[tool.setuptools.packages.find] must set where=['scripts'] and
    include=['gamma*'] to avoid the flat-layout auto-discovery error."""
    find = (
        pyproject_data
        .get("tool", {})
        .get("setuptools", {})
        .get("packages", {})
        .get("find", {})
    )
    assert find.get("where") == ["scripts"], (
        "Expected tool.setuptools.packages.find.where=['scripts'], "
        f"got {find.get('where')!r}"
    )
    include = find.get("include", [])
    assert any(p.startswith("gamma") for p in include), (
        f"Expected 'gamma*' in packages.find.include, got {include!r}"
    )


# ──────────────────────────────────────────────────────────────────────────────
# 3. Smoke import — proves editable install wired scripts/gamma/ correctly
# ──────────────────────────────────────────────────────────────────────────────


def test_gamma_importable():
    """gamma package must be importable (editable install regression)."""
    gamma = importlib.import_module("gamma")
    assert hasattr(gamma, "__version__"), "gamma.__version__ not found"


def test_gamma_spectrum_importable():
    """gamma.spectrum.Spectrum must be importable (primary public surface)."""
    mod = importlib.import_module("gamma.spectrum")
    assert hasattr(mod, "Spectrum"), "gamma.spectrum.Spectrum not found"


def test_gamma_version_matches_pyproject(pyproject_data):
    """gamma.__version__ must match [project].version in pyproject.toml."""
    gamma = importlib.import_module("gamma")
    toml_version = pyproject_data["project"]["version"]
    assert gamma.__version__ == toml_version, (
        f"gamma.__version__={gamma.__version__!r} != "
        f"pyproject.toml version={toml_version!r}"
    )
