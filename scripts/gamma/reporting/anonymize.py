"""
F-115 — Anonymization (D-10, confidentiality).

Strips personal / confidential identifiers from a JSON report dict
before any artefact (JSON / Markdown / HTML / PDF) is written.

Rules:
  * ``operator`` → ``None``
  * ``device_guid`` → ``None``
  * ``sample_id`` — if it matches a known certified-source S/N pattern
    (e.g. ``420-7-17``), the field is reset to ``None``.
    Otherwise we keep only the symbolic prefix.
  * ``detector_id`` — strip a trailing ``№NNNN-NN`` suffix; keep only
    the type name (``УДС-ГЦ-63х63-USB``).
  * ``source_path`` → basename only.
  * Any ``*_path`` / ``*_source`` value that looks like an absolute
    filesystem path → basename only.
  * ``efficiency_source`` inside ``mda`` entries / ``calibration`` →
    basename.

Apply by calling :func:`anonymize_report_inplace` at the very end of
:func:`gamma.reporting.json_report.build_report` / :func:`build_json_report`.
"""
from __future__ import annotations

import os
import re
from typing import Any, Dict


# Certified-source S/N (e.g. "420-7-17", "SN-01", "086-16-1").
_CERT_SN_RE = re.compile(r"^\d{2,4}[-_]\d+(?:[-_]\d+)?$")

# Trailing "№NNNN-NN" or "No NNNN-NN" detector serial token.
# Matches a "№" or "No" followed by digits/dashes; we strip the
# whole token including the preceding whitespace.
_DETECTOR_SN_RE = re.compile(r"\s*(?:№|No\.?|S/?N\.?|s/?n\.?)\s*[\w\-]+\s*$")

# Filesystem path detector — Windows drive letter or POSIX absolute /
# UNC path / contains a path separator deeper than one level.
_FS_PATH_RE = re.compile(
    r"""^(?:
        [A-Za-z]:[\\/]      # C:\ or C:/
        | \\\\               # UNC \\
        | /                  # POSIX absolute
    )""",
    re.VERBOSE,
)


def _basename(p: str) -> str:
    """Return the file basename across forward and back slashes."""
    if not p:
        return p
    # Normalise both separators
    s = p.replace("\\", "/")
    return s.rsplit("/", 1)[-1]


# S/N tokens embedded in basenames — e.g. "..._SN-01_-_..." or
# "..._420-7-17_..." — we strip those before returning.
_EMBEDDED_SN_RE = re.compile(r"[_\- ]+\d{3,4}-\d+(?:-\d+)?")


def _scrub_sn_in_basename(name: str) -> str:
    """Strip embedded S/N tokens from a filename basename."""
    if not name:
        return name
    return _EMBEDDED_SN_RE.sub("", name)


def _is_path(value: str) -> bool:
    """True if value looks like an absolute or deeply-nested path."""
    if not isinstance(value, str):
        return False
    if _FS_PATH_RE.match(value):
        return True
    # Heuristic: contains a backslash or has 2+ forward slashes
    if "\\" in value:
        return True
    if value.count("/") >= 2:
        return True
    return False


def _looks_like_cert_sn(value: str) -> bool:
    if not isinstance(value, str):
        return False
    return bool(_CERT_SN_RE.match(value.strip()))


def _anonymize_sample_id(value: Any) -> Any:
    """Drop S/N from sample_id but keep the symbolic prefix.

    "Th232_420-7-17" → "Th232"  (symbolic)
    "420-7-17"        → None    (pure S/N)
    "Cs137_Marinelli" → "Cs137_Marinelli" (no S/N detected)
    """
    if not isinstance(value, str) or not value:
        return None if value == "" else value
    s = value.strip()
    if _looks_like_cert_sn(s):
        return None
    # Strip embedded S/N tokens like "_420-7-17" / "-420-7-17"
    cleaned = re.sub(r"[_\- ]+\d{2,4}-\d+(?:-\d+)?", "", s)
    if not cleaned or cleaned == s and _CERT_SN_RE.search(s):
        return None
    return cleaned or None


def _anonymize_detector_id(value: Any) -> Any:
    """Strip trailing №NNNN-NN serial token from detector_id."""
    if not isinstance(value, str) or not value:
        return None if value == "" else value
    stripped = _DETECTOR_SN_RE.sub("", value).strip()
    return stripped or None


def _anonymize_path_value(value: Any) -> Any:
    """Convert an absolute path to its basename, then strip any S/N tokens
    embedded inside that basename.  Leave plain strings unchanged unless
    they themselves carry an embedded S/N (e.g. ``foo_SN-01.efr``).
    """
    if not isinstance(value, str) or not value:
        return value
    if _is_path(value):
        return _scrub_sn_in_basename(_basename(value))
    # Even plain filenames may carry S/N tokens.
    if _EMBEDDED_SN_RE.search(value):
        return _scrub_sn_in_basename(value)
    return value


def _scrub_dict_paths(d: Dict[str, Any]) -> None:
    """Recursively basename-ize any *_path / *_source / *_file string fields."""
    if not isinstance(d, dict):
        return
    for k, v in list(d.items()):
        if isinstance(v, dict):
            _scrub_dict_paths(v)
        elif isinstance(v, list):
            for item in v:
                if isinstance(item, dict):
                    _scrub_dict_paths(item)
        elif isinstance(v, str):
            kl = k.lower()
            if kl.endswith("_path") or kl.endswith("_source") or kl.endswith("_file"):
                d[k] = _anonymize_path_value(v)


def anonymize_report_inplace(report: Dict[str, Any]) -> None:
    """Strip personal / confidential identifiers from ``report`` in place.

    F-115. Removes operator names, detector S/N, source IDs, absolute
    paths with embedded S/N, device GUIDs. Replaces with type-only
    strings or None.
    """
    if not isinstance(report, dict):
        return

    h = report.get("header")
    if isinstance(h, dict):
        # Operator and device GUID — always None
        if "operator" in h:
            h["operator"] = None
        if "device_guid" in h:
            h["device_guid"] = None
        # Sample ID — drop S/N
        if "sample_id" in h:
            h["sample_id"] = _anonymize_sample_id(h.get("sample_id"))
        # Detector ID — strip trailing serial
        if "detector_id" in h:
            h["detector_id"] = _anonymize_detector_id(h.get("detector_id"))
        # source_path → basename + S/N scrubbing
        if "source_path" in h:
            sp = h.get("source_path")
            if isinstance(sp, str) and sp:
                h["source_path"] = _scrub_sn_in_basename(_basename(sp))
        # filename also carries S/N occasionally
        if "filename" in h:
            fn = h.get("filename")
            if isinstance(fn, str) and _EMBEDDED_SN_RE.search(fn):
                h["filename"] = _scrub_sn_in_basename(fn)
        # F-144 / v1.17.7 — sample_filename / background_filename тоже
        # могут содержать S/N токены; чистим их симметрично.
        for fkey in ("sample_filename", "background_filename",
                     "background_path"):
            v = h.get(fkey)
            if isinstance(v, str) and v:
                base = _basename(v)
                if _EMBEDDED_SN_RE.search(base):
                    base = _scrub_sn_in_basename(base)
                h[fkey] = base

    # Calibration block — efficiency_source / source paths
    calib = report.get("calibration")
    if isinstance(calib, dict):
        _scrub_dict_paths(calib)

    # Diagnostics — fwhm_model_source / efficiency_source
    diag = report.get("diagnostics")
    if isinstance(diag, dict):
        # Replace by basename if these look like paths
        for key in ("efficiency_source", "fwhm_model_source"):
            if key in diag:
                diag[key] = _anonymize_path_value(diag.get(key))

    # MDA entries — each entry may have a path-bearing notes field
    mda = report.get("mda")
    if isinstance(mda, list):
        for entry in mda:
            if isinstance(entry, dict):
                _scrub_dict_paths(entry)

    # Any other *_path / *_source / *_file in top-level
    _scrub_dict_paths(report)

    # plot_files block — every leaf path → basename + S/N scrubbed
    plots = report.get("plot_files")
    if isinstance(plots, dict):
        if isinstance(plots.get("spectrum"), str):
            plots["spectrum"] = _scrub_sn_in_basename(_basename(plots["spectrum"]))
        if isinstance(plots.get("multiplets"), list):
            plots["multiplets"] = [
                _scrub_sn_in_basename(_basename(p))
                if isinstance(p, str) else p
                for p in plots["multiplets"]
            ]


__all__ = ["anonymize_report_inplace"]
