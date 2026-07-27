"""
ANSI/IEEE N42.42-2012 reader and writer.

Reference: https://www.nist.gov/programs-projects/ansiieee-n4242-standard
SpecUtils canonical implementation:
    https://github.com/sandialabs/SpecUtils/blob/master/src/SpecFile_n42.cpp

Format synopsis (the minimal valid document, attributes shown inline):

    <?xml version="1.0" encoding="UTF-8"?>
    <RadInstrumentData xmlns="http://physics.nist.gov/N42/2011/N42"
                       n42DocUUID="..." ... >
      <RadInstrumentInformation id="instr-1">
        <RadInstrumentManufacturerName>...</RadInstrumentManufacturerName>
        <RadInstrumentIdentifier>...</RadInstrumentIdentifier>
        <RadInstrumentModelName>...</RadInstrumentModelName>
        <RadInstrumentClassCode>Spectroscopic Personal Radiation Detector
        </RadInstrumentClassCode>
      </RadInstrumentInformation>
      <RadDetectorInformation id="det-1">
        <RadDetectorCategoryCode>Gamma</RadDetectorCategoryCode>
        <RadDetectorKindCode>NaI</RadDetectorKindCode>
      </RadDetectorInformation>
      <EnergyCalibration id="ec-1">
        <CoefficientValues>0.0 2.5 0.0</CoefficientValues>
      </EnergyCalibration>
      <RadMeasurement id="rm-1">
        <MeasurementClassCode>Foreground</MeasurementClassCode>
        <StartDateTime>2026-01-01T12:00:00Z</StartDateTime>
        <RealTimeDuration>PT300.00S</RealTimeDuration>
        <Spectrum id="rm-1-sp"
                  radDetectorInformationReference="det-1"
                  energyCalibrationReference="ec-1">
          <LiveTimeDuration>PT295.50S</LiveTimeDuration>
          <ChannelData>0 0 5 12 34 ...</ChannelData>
        </Spectrum>
      </RadMeasurement>
    </RadInstrumentData>

Key facts:

  1. Root element is `RadInstrumentData` in the default namespace
     `http://physics.nist.gov/N42/2011/N42`. We accept files with or
     without the namespace on read (tolerant), but we always emit it.
  2. `CoefficientValues` is whitespace-separated low-to-high.
  3. `ChannelData` is whitespace-separated integers. The optional
     `compressionCode="CountedZeroes"` attribute means: when a `0`
     value appears, the next value is the count of consecutive zeros.
  4. `LiveTimeDuration` and `RealTimeDuration` are ISO-8601 durations
     `PT<seconds>S` (we accept and emit fractional seconds).
  5. `MeasurementClassCode` is one of {Background, Foreground,
     Calibration, IntrinsicActivity, NotSpecified}. We map our
     `is_background` flag to Background ↔ Foreground.
  6. ID-based cross-references: `Spectrum/@energyCalibrationReference`
     points at `EnergyCalibration/@id`, and analogously for detectors.
"""

from __future__ import annotations

import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from xml.etree import ElementTree as ET

# SEC-01 (P1) hardening: parse untrusted .n42 content via defusedxml to block
# billion-laughs / DOCTYPE entity DoS. Stdlib ET is kept for namespace registry,
# Element/SubElement construction in write_n42_2012 (defusedxml is parse-only).
from defusedxml.ElementTree import parse as _safe_parse

import numpy as np

from gamma.spectrum import (
    ENERGY_CEILING_KEV,
    Spectrum,
)


N42_NAMESPACE = "http://physics.nist.gov/N42/2011/N42"
_NS_MAP = {"n42": N42_NAMESPACE}


# ============================================================================
# Reader
# ============================================================================

def read_n42_2012(
    path: str,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """
    Read an N42-42-2012 XML file. Returns the primary (Foreground)
    Spectrum; if the file contains a Background measurement it is
    attached via `spec.background_embedded`.

    Args mirror the other readers in `gamma.io` — see SKILL.md for the
    energy-ceiling contract. `apply_energy_ceiling` defaults to False
    since BUG-9 (v1.18.32, 2026-06-03): the reader keeps every decoded
    channel by default; the 3 MeV trim is opt-in.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"N42-2012 file not found: {path}")

    # SEC-01: _safe_parse blocks DOCTYPE / billion-laughs entity expansion.
    tree = _safe_parse(str(p))
    root = tree.getroot()

    # Detect namespace prefix used in this file (may be the default ns
    # `{http://physics.nist.gov/N42/2011/N42}RadInstrumentData` or no
    # namespace at all, depending on writer).
    ns = _detect_namespace(root)
    if _strip_ns(root.tag) != "RadInstrumentData":
        raise ValueError(
            f"Not an N42-2012 RadInstrumentData (root is {root.tag!r}): {path}"
        )

    # ----- collect cross-referenced blocks -----
    instr_info = _first_local(root, "RadInstrumentInformation", ns)
    detectors = _all_local(root, "RadDetectorInformation", ns)
    calibrations = {
        c.get("id"): c for c in _all_local(root, "EnergyCalibration", ns)
    }
    measurements = _all_local(root, "RadMeasurement", ns)
    if not measurements:
        raise ValueError(f"No RadMeasurement entries in {path}")

    # ----- pick foreground and (optional) background measurement -----
    fg_meas = None
    bg_meas = None
    for m in measurements:
        cls = _text(_first_local(m, "MeasurementClassCode", ns), "")
        if cls.lower() == "background":
            bg_meas = bg_meas or m
        else:
            fg_meas = fg_meas or m
    if fg_meas is None:
        fg_meas = measurements[0]

    fg = _measurement_to_spectrum(
        fg_meas, calibrations, ns, str(p),
        is_background=False,
        apply_energy_ceiling=apply_energy_ceiling,
        ceiling_keV=ceiling_keV,
    )
    if bg_meas is not None:
        bg = _measurement_to_spectrum(
            bg_meas, calibrations, ns, str(p),
            is_background=True,
            apply_energy_ceiling=apply_energy_ceiling,
            ceiling_keV=ceiling_keV,
        )
        fg.background_embedded = bg

    # ----- merge instrument metadata into Spectrum + extras -----
    if instr_info is not None:
        manu = _text(_first_local(instr_info, "RadInstrumentManufacturerName", ns))
        model = _text(_first_local(instr_info, "RadInstrumentModelName", ns))
        ident = _text(_first_local(instr_info, "RadInstrumentIdentifier", ns))
        if manu or model:
            fg.detector_id = (f"{manu} {model}").strip() or fg.detector_id
        if ident:
            fg.device_guid = ident

    if detectors:
        det_kinds = [
            _text(_first_local(d, "RadDetectorKindCode", ns)) for d in detectors
        ]
        det_kinds = [k for k in det_kinds if k]
        if det_kinds:
            fg.extras["n42_detector_kinds"] = det_kinds

    fg.source_format = "n42_2012"
    return fg


def _measurement_to_spectrum(
    meas: ET.Element,
    calibrations: dict,
    ns: str,
    source_path: str,
    *,
    is_background: bool,
    apply_energy_ceiling: bool,
    ceiling_keV: Optional[float],
) -> Spectrum:
    """Convert a single <RadMeasurement> into a Spectrum dataclass."""
    real_time = _parse_duration(_text(_first_local(meas, "RealTimeDuration", ns)))
    start_dt = _parse_iso(_text(_first_local(meas, "StartDateTime", ns)))

    spec_elem = _first_local(meas, "Spectrum", ns)
    if spec_elem is None:
        raise ValueError("RadMeasurement without <Spectrum>")

    live_time = _parse_duration(
        _text(_first_local(spec_elem, "LiveTimeDuration", ns))
    )
    if live_time <= 0 and real_time > 0:
        # Some emitters write LiveTimeDuration only on the parent
        # RadMeasurement; fall back to real time.
        live_time = real_time

    # Channel data
    cd = _first_local(spec_elem, "ChannelData", ns)
    compression = (cd.get("compressionCode") if cd is not None else "") or ""
    raw_text = (cd.text or "") if cd is not None else ""
    counts = _decode_channel_data(raw_text, compression)

    # Energy calibration (cross-referenced)
    cal_ref = spec_elem.get("energyCalibrationReference") or ""
    coeffs = ()
    deg = None
    cal_node = calibrations.get(cal_ref)
    if cal_node is None and calibrations:
        # Single calibration in the file? Use it implicitly.
        if len(calibrations) == 1:
            cal_node = next(iter(calibrations.values()))
    if cal_node is not None:
        cv = _first_local(cal_node, "CoefficientValues", ns)
        if cv is not None and cv.text:
            coeffs = tuple(
                float(x) for x in cv.text.split() if x
            )
            deg = max(0, len(coeffs) - 1)

    # Apply energy ceiling
    n_full = int(len(counts))
    if apply_energy_ceiling and coeffs:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        keep_n, e_max_kept = _apply_ceiling(counts, coeffs, ceiling)
        if keep_n < n_full:
            counts = counts[:keep_n].copy()
        else:
            e_max_kept = _energy_at(n_full - 1, coeffs) if n_full else None
    else:
        e_max_kept = (
            _energy_at(n_full - 1, coeffs) if (coeffs and n_full) else None
        )

    spec = Spectrum(
        counts=counts,
        live_time=float(live_time),
        real_time=float(real_time),
        source_path=source_path,
        source_format="n42_2012",
        is_background=is_background,
        start_datetime=start_dt,
        n_channels_raw=n_full,
        n_channels=int(len(counts)),
        channel_pitch=1,
        energy_cal=coeffs if coeffs else None,
        energy_cal_degree=deg,
        energy_cal_source="stored" if coeffs else "",
        energy_max_keV_kept=e_max_kept,
    )
    if compression:
        spec.extras["n42_channel_compression"] = compression
    return spec


# ============================================================================
# Writer
# ============================================================================

def write_n42_2012(
    spec: Spectrum,
    path: str,
    *,
    pretty: bool = True,
) -> None:
    """
    Write a Spectrum to N42-42-2012 XML.

    Emits a single-detector, single-measurement document. If the spectrum
    carries an embedded background (`spec.background_embedded`), it is
    added as a second <RadMeasurement> with MeasurementClassCode=Background.
    """
    ET.register_namespace("", N42_NAMESPACE)

    root = ET.Element(_qname("RadInstrumentData"))
    root.set("n42DocUUID", str(uuid.uuid4()))
    root.set("n42DocDateTime", _iso_now())

    # --- Instrument information ---
    instr_id = "instr-1"
    instr = ET.SubElement(root, _qname("RadInstrumentInformation"))
    instr.set("id", instr_id)
    manu = ET.SubElement(instr, _qname("RadInstrumentManufacturerName"))
    manu.text = "Unknown"
    model = ET.SubElement(instr, _qname("RadInstrumentModelName"))
    model.text = spec.detector_id or "Unknown"
    ident = ET.SubElement(instr, _qname("RadInstrumentIdentifier"))
    ident.text = spec.device_guid or spec.sample_id or "unknown"
    cls = ET.SubElement(instr, _qname("RadInstrumentClassCode"))
    cls.text = "Spectroscopic Personal Radiation Detector"

    # --- Detector information ---
    det_id = "det-1"
    det = ET.SubElement(root, _qname("RadDetectorInformation"))
    det.set("id", det_id)
    det_cat = ET.SubElement(det, _qname("RadDetectorCategoryCode"))
    det_cat.text = "Gamma"
    det_kind = ET.SubElement(det, _qname("RadDetectorKindCode"))
    det_kind.text = "NaI"

    # --- Energy calibration (single shared) ---
    cal_id = "ec-1"
    if spec.energy_cal:
        cal = ET.SubElement(root, _qname("EnergyCalibration"))
        cal.set("id", cal_id)
        cv = ET.SubElement(cal, _qname("CoefficientValues"))
        cv.text = " ".join(f"{c:.10g}" for c in spec.energy_cal)
    else:
        cal_id = ""

    # --- Foreground measurement ---
    _append_measurement(
        root,
        spec=spec,
        meas_id="rm-1",
        det_ref=det_id,
        cal_ref=cal_id,
        cls_code="Background" if spec.is_background else "Foreground",
    )

    # --- Background measurement (optional) ---
    if spec.background_embedded is not None:
        bg = spec.background_embedded
        bg_cal_id = ""
        if bg.energy_cal and bg.energy_cal != spec.energy_cal:
            bg_cal_id = "ec-bg"
            cal_bg = ET.SubElement(root, _qname("EnergyCalibration"))
            cal_bg.set("id", bg_cal_id)
            cv_bg = ET.SubElement(cal_bg, _qname("CoefficientValues"))
            cv_bg.text = " ".join(f"{c:.10g}" for c in bg.energy_cal)
        elif bg.energy_cal:
            bg_cal_id = cal_id
        _append_measurement(
            root,
            spec=bg,
            meas_id="rm-bg",
            det_ref=det_id,
            cal_ref=bg_cal_id,
            cls_code="Background",
        )

    # --- Serialize ---
    if pretty:
        _indent(root)
    tree = ET.ElementTree(root)
    tree.write(str(path), encoding="utf-8", xml_declaration=True)


def _append_measurement(
    root: ET.Element,
    *,
    spec: Spectrum,
    meas_id: str,
    det_ref: str,
    cal_ref: str,
    cls_code: str,
) -> None:
    """Add one <RadMeasurement> block."""
    m = ET.SubElement(root, _qname("RadMeasurement"))
    m.set("id", meas_id)

    cls = ET.SubElement(m, _qname("MeasurementClassCode"))
    cls.text = cls_code

    if spec.start_datetime is not None:
        sdt = ET.SubElement(m, _qname("StartDateTime"))
        sdt.text = _format_iso(spec.start_datetime)

    rt = ET.SubElement(m, _qname("RealTimeDuration"))
    rt.text = _format_duration(spec.real_time)

    sp = ET.SubElement(m, _qname("Spectrum"))
    sp.set("id", f"{meas_id}-sp")
    sp.set("radDetectorInformationReference", det_ref)
    if cal_ref:
        sp.set("energyCalibrationReference", cal_ref)

    lt = ET.SubElement(sp, _qname("LiveTimeDuration"))
    lt.text = _format_duration(spec.live_time)

    cd = ET.SubElement(sp, _qname("ChannelData"))
    counts = np.asarray(spec.counts, dtype=np.int64)
    # Standard recommends CountedZeroes for sparse spectra. We always emit
    # uncompressed for clarity — readers per-spec must accept this.
    cd.text = " ".join(str(int(x)) for x in counts)


# ============================================================================
# Helpers — namespaces / element access
# ============================================================================

def _detect_namespace(root: ET.Element) -> str:
    """Return the namespace prefix used in the document, '' if none."""
    if root.tag.startswith("{"):
        return root.tag[1:].split("}", 1)[0]
    return ""


def _strip_ns(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def _qname(local: str) -> str:
    return f"{{{N42_NAMESPACE}}}{local}"


def _first_local(parent: ET.Element, local: str, ns: str) -> Optional[ET.Element]:
    """Find first descendant whose local name is `local`, ignoring namespace."""
    if parent is None:
        return None
    # Try fully-qualified search first (faster)
    if ns:
        node = parent.find(f"{{{ns}}}{local}")
        if node is not None:
            return node
    # Fallback: namespace-agnostic walk
    for child in parent.iter():
        if _strip_ns(child.tag) == local:
            return child
    return None


def _all_local(parent: ET.Element, local: str, ns: str) -> list:
    if parent is None:
        return []
    out = []
    for child in parent.iter():
        if _strip_ns(child.tag) == local:
            out.append(child)
    return out


def _text(elem: Optional[ET.Element], default: str = "") -> str:
    if elem is None or elem.text is None:
        return default
    return elem.text.strip()


# ============================================================================
# Helpers — duration / datetime parsing
# ============================================================================

_DURATION_RE = re.compile(r"PT([0-9.]+)S", re.IGNORECASE)


def _parse_duration(s: str) -> float:
    """Parse ISO-8601 duration `PT<n>S` → seconds (float). 0.0 on failure."""
    if not s:
        return 0.0
    m = _DURATION_RE.match(s.strip())
    if not m:
        return 0.0
    try:
        return float(m.group(1))
    except ValueError:
        return 0.0


def _format_duration(seconds: float) -> str:
    return f"PT{float(seconds):.3f}S"


def _parse_iso(text: str) -> Optional[datetime]:
    if not text:
        return None
    t = text.strip()
    # Replace trailing Z with +00:00 for fromisoformat (Py 3.11+ tolerates Z)
    if t.endswith("Z"):
        t = t[:-1] + "+00:00"
    # Trim fractional seconds beyond microseconds
    m = re.match(
        r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(\.\d+)?([+-]\d{2}:\d{2})?$",
        t,
    )
    if m:
        base = m.group(1)
        frac = (m.group(2) or "")[:7]  # ".XXXXXX"
        tz = m.group(3) or ""
        t = base + frac + tz
    try:
        return datetime.fromisoformat(t)
    except ValueError:
        return None


def _format_iso(dt: datetime) -> str:
    if dt.tzinfo is None:
        # Treat as UTC if no tz attached
        return dt.replace(tzinfo=timezone.utc).isoformat()
    return dt.isoformat()


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ============================================================================
# Helpers — channel data encoding
# ============================================================================

#: SEC-01 hardening cap. CountedZeroes expansion of an adversary-supplied
#: <100-byte N42 XML can otherwise force a ~16 GB allocation
#: (`0 2000000000` token pair). 1e7 channels is three orders of magnitude
#: above realistic NaI / HPGe upper bounds (NaI 1024-16384, HPGe 8192-65536)
#: so legitimate spectra are unaffected. Breach raises `ValueError` so the
#: upstream reader fails loud rather than silently OOM-ing the workstation.
_COUNTED_ZEROES_MAX_CHANNELS = 10_000_000  # 1e7, fail-loud on breach


def _decode_channel_data(text: str, compression: str) -> np.ndarray:
    """Decode whitespace-separated integers; expand CountedZeroes if used."""
    if not text:
        return np.array([], dtype=np.int64)
    tokens = text.split()
    if compression.lower() == "countedzeroes":
        out = []
        it = iter(tokens)
        for tok in it:
            try:
                v = int(float(tok))
            except ValueError:
                continue
            if v == 0:
                try:
                    n = int(float(next(it)))
                except (StopIteration, ValueError):
                    n = 1
                run_length = max(1, n)
                # SEC-01: bound cumulative expansion. Fail-loud on breach so
                # the upstream N42 reader cannot silently OOM the host on
                # adversary-supplied input.
                if len(out) + run_length > _COUNTED_ZEROES_MAX_CHANNELS:
                    raise ValueError(
                        f"CountedZeroes expansion exceeds bound "
                        f"{_COUNTED_ZEROES_MAX_CHANNELS} channels "
                        f"(would reach {len(out) + run_length})"
                    )
                out.extend([0] * run_length)
            else:
                out.append(v)
        return np.array(out, dtype=np.int64)
    # Uncompressed: fast path
    # AUDIT-F7 (2026-06-25): np.fromstring → np.array(text.split(), dtype=np.int64).
    # np.fromstring deprecated since numpy 1.14; bit-identical для int64-cast того же
    # whitespace-split (split на whitespace тождествен sep=" " для контента N42).
    try:
        return np.array(text.split(), dtype=np.int64)
    except Exception:
        return np.array([int(float(x)) for x in tokens], dtype=np.int64)


# ============================================================================
# Helpers — energy ceiling (mirrors lsrm_spe.py)
# ============================================================================

def _energy_at(ch: int, coefs) -> float:
    return sum(a * (ch ** i) for i, a in enumerate(coefs))


def _apply_ceiling(counts, coefs, ceiling_keV: float):
    """Return (keep_n, e_max_kept). keep_n is the channel count to retain."""
    n = int(len(counts))
    if n == 0 or not coefs:
        return n, None
    energies = np.zeros(n, dtype=np.float64)
    channels = np.arange(n, dtype=np.float64)
    for c in reversed(coefs):
        energies = energies * channels + c
    mask = energies <= ceiling_keV
    if mask.all():
        return n, float(energies[-1])
    first_above = int(np.argmax(~mask))
    if first_above == 0:
        return 0, None
    return first_above, float(energies[first_above - 1])


# ============================================================================
# Helpers — pretty-print (ElementTree.indent only on 3.9+)
# ============================================================================

def _indent(elem: ET.Element, level: int = 0, step: str = "  ") -> None:
    """Backport of ElementTree.indent for portable pretty-printing."""
    pad = "\n" + step * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = pad + step
        for i, child in enumerate(elem):
            _indent(child, level + 1, step)
            if i == len(elem) - 1:
                if not child.tail or not child.tail.strip():
                    child.tail = pad
            else:
                if not child.tail or not child.tail.strip():
                    child.tail = pad + step
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = pad
