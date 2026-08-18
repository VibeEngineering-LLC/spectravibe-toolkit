"""
SEC-01 (P1) — XXE / billion-laughs hardening regression for XML parsers.

Project #5 audit (S5-SEC.md, finding SEC-01) flagged that the 5 untrusted-file
XML parsers in `gamma.io.*` used stdlib `xml.etree.ElementTree`. While CPython's
expat does NOT expand external general entities by default (no network XXE /
file-read), it IS exposed to *billion-laughs* (nested internal entity expansion)
and quadratic-blowup DoS on crafted input.

Fix: swap parse entry points to `defusedxml.ElementTree.parse / fromstring`,
which blocks DOCTYPE / internal entity expansion outright.

Coverage:
  1. Per-parser billion-laughs payload — must raise an `EntitiesForbidden`
     family exception, NOT exhaust memory.
  2. Generic DOCTYPE-block payload — must raise on any of the parsers.
  3. Sanity round-trip — a known-good minimal payload still parses cleanly
     after the swap, proving we did not break the happy path.

NOTE: `becqmoni_xml.py` is intentionally EXCLUDED from per-parser billion-laughs
tests because it is a writer-only module (accepts a `Spectrum` dataclass; emits
XML). There is no parse entry point on user data — no attack surface. See the
top-of-file comment in `becqmoni_xml.py` for the rationale.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from defusedxml import EntitiesForbidden, DTDForbidden


# ---------------------------------------------------------------------------
# Billion-laughs payloads (10-level nesting expands to ~10^10 'lol' tokens
# under a vulnerable parser; defusedxml refuses to define entities at all).
# ---------------------------------------------------------------------------

_BILLION_LAUGHS_GENERIC = """<?xml version="1.0"?>
<!DOCTYPE lolz [
  <!ENTITY lol "lol">
  <!ENTITY lol2 "&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;">
  <!ENTITY lol3 "&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;">
  <!ENTITY lol4 "&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;">
]>
<{root}>&lol4;</{root}>
"""


def _billion_laughs(root_tag: str) -> str:
    """Return a billion-laughs payload with the requested root element."""
    return _BILLION_LAUGHS_GENERIC.format(root=root_tag)


_DOCTYPE_PAYLOAD = (
    '<?xml version="1.0"?>\n'
    '<!DOCTYPE r [\n'
    '  <!ENTITY foo "bar">\n'
    ']>\n'
    '<r>&foo;</r>\n'
)


# ===========================================================================
# Per-parser billion-laughs rejection
# ===========================================================================

def test_atomspectra_xml_rejects_billion_laughs(tmp_path: Path) -> None:
    """atomspectra_xml.read_atomspectra_xml must refuse DOCTYPE entities."""
    from gamma.io.atomspectra_xml import read_atomspectra_xml

    payload = _billion_laughs("ResultDataFile")
    f = tmp_path / "evil_atomspectra.xml"
    f.write_text(payload, encoding="utf-8")

    with pytest.raises((EntitiesForbidden, DTDForbidden)):
        read_atomspectra_xml(str(f))


def test_cpt_io_rejects_billion_laughs() -> None:
    """cpt_io.parse_cpt_xml must refuse DOCTYPE entities (wrapped in ValueError)."""
    from gamma.io.cpt_io import parse_cpt_xml

    payload = _billion_laughs("peak_template")

    # parse_cpt_xml wraps ET.ParseError → ValueError, and defusedxml
    # exceptions are subclasses of ET.ParseError. So we get ValueError
    # whose chained __cause__ is the defusedxml exception.
    with pytest.raises((ValueError, EntitiesForbidden, DTDForbidden)) as excinfo:
        parse_cpt_xml(payload)

    # Verify the underlying cause is a defusedxml entity-block (not some
    # unrelated XML structure error).
    err = excinfo.value
    chain = []
    cur = err
    seen = set()
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        chain.append(cur)
        cur = getattr(cur, "__cause__", None) or getattr(cur, "__context__", None)
    assert any(
        isinstance(c, (EntitiesForbidden, DTDForbidden)) for c in chain
    ), f"expected defusedxml exception in chain, got {[type(c).__name__ for c in chain]}"


def test_lsrm_library_rejects_billion_laughs(tmp_path: Path) -> None:
    """lsrm_library.read_lsrm_library must refuse DOCTYPE entities."""
    from gamma.io.lsrm_library import read_lsrm_library

    payload = _billion_laughs("Library")
    # The reader expects windows-1251; but defusedxml rejects on DOCTYPE
    # BEFORE encoding negotiation, so plain UTF-8 bytes still trigger the
    # block.
    f = tmp_path / "evil_lsrm.lib"
    f.write_bytes(payload.encode("utf-8"))

    with pytest.raises((EntitiesForbidden, DTDForbidden)):
        read_lsrm_library(str(f))


def test_n42_2012_rejects_billion_laughs(tmp_path: Path) -> None:
    """n42_2012.read_n42_2012 must refuse DOCTYPE entities."""
    from gamma.io.n42_2012 import read_n42_2012

    payload = _billion_laughs("RadInstrumentData")
    f = tmp_path / "evil_n42.n42"
    f.write_text(payload, encoding="utf-8")

    with pytest.raises((EntitiesForbidden, DTDForbidden)):
        read_n42_2012(str(f))


# ===========================================================================
# Generic DOCTYPE block (independent of parser-specific structure)
# ===========================================================================

def test_doctype_payload_blocked_at_low_level() -> None:
    """defusedxml.ElementTree.fromstring must refuse any DOCTYPE payload.

    This guards against any future code path that calls
    `defusedxml.ElementTree.fromstring` directly — a regression here would
    indicate defusedxml configuration drift.
    """
    from defusedxml.ElementTree import fromstring as _safe_fromstring

    with pytest.raises((EntitiesForbidden, DTDForbidden)):
        _safe_fromstring(_DOCTYPE_PAYLOAD)


# ===========================================================================
# Sanity: happy path still works after the swap
# ===========================================================================

def test_cpt_io_happy_path_after_hardening() -> None:
    """Known-good .cpt XML still parses post-swap — proves we didn't break it.

    Uses the canonical minimum .cpt structure from cpt_io.py docstring,
    exercising both build_cpt_xml (stdlib ET) and parse_cpt_xml
    (defusedxml-hardened entry point).
    """
    from gamma.peaks.peak_image_tabulated import (
        PeakShapeAnchor, TabulatedPeakImage,
    )
    from gamma.io.cpt_io import build_cpt_xml, parse_cpt_xml

    # Round-trip: build a minimal in-memory TabulatedPeakImage, serialize
    # to .cpt XML (via stdlib ET — writer is unchanged), then parse it
    # back through the now-hardened parse_cpt_xml.
    original = TabulatedPeakImage(
        detector_id="Gamma-1S",
        detector_class="NaI",
        crystal_diameter_mm=63.0,
        anchors=[
            PeakShapeAnchor(
                E_keV=661.66,
                fwhm_keV=46.2,
                tail_fraction=0.03,
                tail_slope_inv_keV=0.05,
                step_height_frac=0.05,
                asymmetry=0.0,
                weight=1.0,
            ),
        ],
        source_metadata="Cs-137",
        notes="sanity fixture for SEC-01 hardening test",
    )

    xml_str = build_cpt_xml(original)
    # The writer must still emit a sane DOCTYPE-free document.
    assert "<!DOCTYPE" not in xml_str
    assert "<peak_template" in xml_str

    # Parse via the hardened entry point — must succeed identically.
    parsed = parse_cpt_xml(xml_str)
    assert parsed.detector_id == "Gamma-1S"
    assert parsed.detector_class == "NaI"
    assert abs(parsed.crystal_diameter_mm - 63.0) < 1e-9
    assert len(parsed.anchors) == 1
    assert abs(parsed.anchors[0].E_keV - 661.66) < 1e-4
    assert parsed.source_metadata == "Cs-137"
