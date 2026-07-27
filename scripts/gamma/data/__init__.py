"""
gamma.data — JSON data loaders for nuclides, XRF lines, anchor patterns.

These loaders implement the token-economy principle: the AI does not see
the JSON files; it sees only the result of lookups by name or energy.

The JSON files live in `data/` at the skill root, alongside `references/`.
This module resolves that path regardless of where the package is imported
from.
"""

from __future__ import annotations

from pathlib import Path

# Path to data/ at skill root:
#   <skill_root>/data/         <- where JSON files live
#   <skill_root>/scripts/gamma/data/__init__.py  <- this file
# So data/ is at ../../../data relative to this file.
DATA_DIR = Path(__file__).resolve().parents[3] / "data"

__all__ = ["DATA_DIR"]
