"""
Lsrm `.src` certificate file reader.

The `.src` file format stores the **certified activity** of reference
γ-ray sources used for spectrometer calibration. It is the source of
ground truth for activity validation (Phase 2.1d, F-29 — the
`compute_activity` result must agree with the certificate after
decay correction to the spectrum's measurement date).

═══════════════════════════════════════════════════════════════════
Format reference
═══════════════════════════════════════════════════════════════════

The format is **INI-style**, CP-1251 encoded, with the following
peculiarities:

  • Decimal separator is **comma** (Russian convention):
      Thick,mm=10,1   ⟶  Thick,mm = 10.1
      Cs-137=106000,3 ⟶  A = 106 000 Bq, σ = 3% (at the confidence
                          declared in [General] Sigma=...)
    These differ contextually: outside Act blocks, the comma is a
    decimal point; inside Act blocks (for nuclide-keyed lines), the
    comma separates value from σ%.

  • Section names may contain commas, so a hand-rolled parser is
    used (configparser is awkward with such headers).

  • Each main source is a hierarchy of three levels:
      [<source>]                 — metadata (geometry, date, ...)
      [<source>,structure]       — list of sub-sources contributing
      [<source>,<sub>,Act]       — activity values per sub-source

  • `[General]` carries `Sigma=N` where N is the multiplier on σ:
    `Sigma=2` means the uncertainty values are stated at 2σ (95.4%);
    `Sigma=1` would mean 1σ (68%). Callers wanting 1σ values must
    divide by `confidence_sigma`.

═══════════════════════════════════════════════════════════════════
Schema (observed in 7 reference .src files)
═══════════════════════════════════════════════════════════════════

[General]
Sigma=2

[Sets]
<source_1>=
<source_2>=
...

[<source>]
Geometry=Точечная | Маринелли | Петри-60мл | Дента-100мл | ...
Mass,g=<float-with-comma>     (may be empty)
Volume,ml=<float-with-comma>  (may be empty)
Material=<str>
Date=DD.MM.YYYY
Time=HH:MM:SS                 (often "0:00:00", treat as 00:00:00)
Units=Bq | Bq/kg              (older files: "Activity unit=Bq/kg")
Thick,mm=<float-with-comma>
Comment=<str>
UseShield=0 | 1
ShieldMaterial=<str>
ShieldThickness,mm=<float-with-comma>

[<source>,structure]
<sub_source_1>=
<sub_source_2>=
...

[<source>,<sub_source>,Act]
Mass,g=<float-with-comma>     (may be empty)
<Nuclide>=<int_Bq>,<int_sigma_pct>
<Nuclide>=...

═══════════════════════════════════════════════════════════════════
Citation
═══════════════════════════════════════════════════════════════════

Lsrm spec for `.src` is referenced in *Описание формата файла Лсрм*
§7.5.9 (the chapter heading is the only fragment supplied; the
schema above was reverse-engineered from 7 real-world fixtures).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional


# ═════════════════════════════════════════════════════════════════════
# Data classes
# ═════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class CertificateActivity:
    """One nuclide activity entry from a certificate."""

    nuclide: str               # e.g. "Cs-137"
    A_Bq: float                # activity in Bq (or Bq/kg — see `unit`)
    sigma_pct: float           # uncertainty in % (at confidence_sigma σ)
    confidence_sigma: int      # multiplier on σ: 1 for 1σ, 2 for 2σ, etc.
    unit: str                  # "Bq" or "Bq/kg"

    def sigma_1pct(self) -> float:
        """Uncertainty converted to 1σ percent."""
        if self.confidence_sigma <= 0:
            return self.sigma_pct
        return self.sigma_pct / self.confidence_sigma

    def sigma_1_Bq(self) -> float:
        """Absolute 1σ uncertainty in Bq."""
        return self.A_Bq * self.sigma_1pct() / 100.0

    def __repr__(self) -> str:
        return (f"CertificateActivity({self.nuclide}: "
                f"A={self.A_Bq:.4g} {self.unit}, "
                f"σ={self.sigma_pct:.3g}% @ {self.confidence_sigma}σ)")


@dataclass(frozen=True)
class CertificateSubSource:
    """A sub-source contributing to a main certificate source."""

    name: str
    mass_g: Optional[float] = None
    activities: tuple = ()     # tuple[CertificateActivity, ...]


@dataclass(frozen=True)
class CertificateSource:
    """One main certified source from a .src file."""

    name: str
    geometry: str
    mass_g: Optional[float] = None
    volume_ml: Optional[float] = None
    material: str = ""
    reference_datetime: Optional[datetime] = None
    activity_unit: str = "Bq"
    thickness_mm: Optional[float] = None
    comment: str = ""
    use_shield: bool = False
    shield_material: str = ""
    shield_thickness_mm: Optional[float] = None
    sub_sources: tuple = ()    # tuple[CertificateSubSource, ...]

    def all_activities(self) -> tuple:
        """Flatten all per-sub-source activities into one tuple."""
        out = []
        for sub in self.sub_sources:
            out.extend(sub.activities)
        return tuple(out)

    def get_activity(self, nuclide: str) -> Optional[CertificateActivity]:
        """Return the (first) activity for the given nuclide, or None.

        A main source may have several sub-sources, and the same
        nuclide could in principle appear in more than one of them
        — this returns the first match scanning sub-sources in
        order. For multi-sub-source aggregation, use `all_activities`
        and sum manually.
        """
        nuc_norm = nuclide.strip().lower()
        for sub in self.sub_sources:
            for act in sub.activities:
                if act.nuclide.strip().lower() == nuc_norm:
                    return act
        return None

    def __repr__(self) -> str:
        date_str = (self.reference_datetime.strftime("%Y-%m-%d")
                    if self.reference_datetime else "no date")
        return (f"CertificateSource({self.name!r}, geom={self.geometry!r}, "
                f"date={date_str}, {len(self.sub_sources)} sub-sources)")


@dataclass(frozen=True)
class Certificate:
    """A parsed .src file."""

    file_path: str
    confidence_sigma: int      # from [General] Sigma=N
    sources: dict              # dict[str, CertificateSource]

    def source_names(self) -> list:
        return list(self.sources.keys())

    def find_source(self, name_pattern: str) -> Optional[CertificateSource]:
        """Find a source by exact name (case-insensitive, whitespace-stripped).

        For fuzzy matching across multiple plausible names, use
        `find_source_fuzzy`. This method enforces strict equality
        after normalisation — useful when the caller knows the exact
        name.
        """
        key = name_pattern.strip().lower()
        for name, src in self.sources.items():
            if name.strip().lower() == key:
                return src
        return None

    def find_source_fuzzy(self, name_pattern: str) -> Optional[CertificateSource]:
        """Find a source by substring match (case-insensitive).

        Useful when the spectrum's SHIFR field carries an abbreviated
        or punctuation-variant of the certificate's source name —
        e.g. "Cs-137 SRC-02" matching certificate
        "Cs-137 №SRC-02". Returns the first hit; for ambiguous
        queries use `find_source_candidates`.
        """
        candidates = self.find_source_candidates(name_pattern)
        return candidates[0] if candidates else None

    def find_source_candidates(self, name_pattern: str) -> list:
        """All sources whose name contains the pattern (case-insensitive).

        Tokenises both pattern and candidate names on whitespace and
        common punctuation, then requires every pattern token to
        appear (as a substring) in the candidate. Robust to "№"
        variations, underscore vs space, and casing.
        """
        norm_pat = _normalise_token(name_pattern)
        pat_tokens = [t for t in norm_pat.split() if t]
        out = []
        for name, src in self.sources.items():
            norm = _normalise_token(name)
            if all(tok in norm for tok in pat_tokens):
                out.append(src)
        return out

    def __repr__(self) -> str:
        return (f"Certificate({Path(self.file_path).name!r}, "
                f"Sigma={self.confidence_sigma}σ, "
                f"{len(self.sources)} sources)")


def _normalise_token(s: str) -> str:
    """Lowercase + strip + collapse punctuation into spaces."""
    s = s.lower().strip()
    # treat №, #, /, _, , (comma in name like Маринелли-ППД ОМАНС-1750)
    # as separators by replacing with space, then collapse whitespace
    for ch in "№#/_-.,":
        s = s.replace(ch, " ")
    s = re.sub(r"\s+", " ", s)
    return s


# ═════════════════════════════════════════════════════════════════════
# Parsing
# ═════════════════════════════════════════════════════════════════════

_NUCLIDE_KEY_RE = re.compile(r"^[A-Z][a-z]?-\d{1,3}m?$")


def _decimal_comma_to_float(s: str) -> Optional[float]:
    """Parse a comma-decimal-separated number to float, or None if empty/invalid."""
    s = (s or "").strip()
    if not s:
        return None
    try:
        return float(s.replace(",", "."))
    except ValueError:
        return None


def _parse_date(date_str: str, time_str: str = "") -> Optional[datetime]:
    """Parse certificate date `DD.MM.YYYY` + optional `HH:MM:SS`."""
    date_str = (date_str or "").strip()
    time_str = (time_str or "").strip()
    if not date_str:
        return None
    for date_fmt in ("%d.%m.%Y", "%d.%m.%y", "%d-%m-%Y", "%d-%m-%y"):
        try:
            d = datetime.strptime(date_str, date_fmt)
            break
        except ValueError:
            d = None
    if d is None:
        return None
    if not time_str:
        return d
    for time_fmt in ("%H:%M:%S", "%H:%M"):
        try:
            t = datetime.strptime(time_str, time_fmt)
            return d.replace(hour=t.hour, minute=t.minute, second=t.second)
        except ValueError:
            continue
    return d  # time unparseable → date-only


def _tokenise_sections(text: str) -> list:
    """Yield (section_name, dict[key, value]) pairs in file order."""
    current_name: Optional[str] = None
    current_kv: dict = {}
    out: list = []
    for raw in text.splitlines():
        line = raw.rstrip("\r").strip()
        if not line or line.startswith(";") or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            if current_name is not None:
                out.append((current_name, current_kv))
            current_name = line[1:-1]
            current_kv = {}
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        # Note: a key may legitimately contain a comma (e.g. "Mass,g",
        # "Thick,mm"). We keep the key verbatim.
        key = key.strip()
        value = value.strip()
        # Duplicate keys within a section: keep the LAST occurrence
        # (this matches typical INI behaviour and the structure of
        # the observed .src files where duplicates are not expected).
        current_kv[key] = value
    if current_name is not None:
        out.append((current_name, current_kv))
    return out


def _parse_activity_block(kv: dict, confidence_sigma: int,
                          unit_hint: str) -> tuple:
    """Parse one [...,Act] block into (mass_g, tuple[CertificateActivity, ...])."""
    mass_g = _decimal_comma_to_float(kv.get("Mass,g", ""))
    activities = []
    for key, val in kv.items():
        if key in ("Mass,g",):
            continue
        if not _NUCLIDE_KEY_RE.match(key):
            continue
        # Activity entries: <int_value>,<int_sigma_pct>
        # Special handling: split on FIRST comma only (the value is
        # an integer here, but be defensive and accept decimal values
        # too — re-construct the "value" half by joining all but the
        # last comma-separated part).
        parts = val.split(",")
        if len(parts) < 2:
            # Lonely value with no σ — skip (or treat σ=0)
            try:
                A = float(val)
                sigma = 0.0
            except ValueError:
                continue
        else:
            # Value = everything but the last component (rejoined with .)
            # σ = last component
            val_part = ".".join(parts[:-1])  # treat first comma as decimal
            sigma_part = parts[-1]
            try:
                A = float(val_part)
                sigma = float(sigma_part)
            except ValueError:
                continue
            # For our observed files, A is always integer (no internal
            # decimal point), so the above join is overly defensive but
            # safe. The most common case is "106000,3" → A=106000, σ=3.
        activities.append(CertificateActivity(
            nuclide=key,
            A_Bq=A,
            sigma_pct=sigma,
            confidence_sigma=confidence_sigma,
            unit=unit_hint or "Bq",
        ))
    return mass_g, tuple(activities)


def _assemble_certificate(file_path: str,
                          sections: list) -> Certificate:
    """Build the Certificate object from the tokenised sections.

    Two-pass over the section list:
      1) Collect [General], [Sets], and all main source sections.
      2) Wire sub-sources and Act blocks into the right CertificateSource.
    """
    general_kv = {}
    set_names: list = []
    # main section name → metadata kv dict
    main_kv: dict = {}
    # main name → list[sub_name]
    structure: dict = {}
    # (main, sub) → (mass_g, activities tuple)
    act_blocks: dict = {}

    for name, kv in sections:
        if name == "General":
            general_kv = kv
        elif name == "Sets":
            # [Sets] lists source headers as keys with empty values
            set_names = [k for k in kv.keys()]
        elif name.endswith(",structure"):
            main_name = name[:-len(",structure")]
            structure[main_name] = [k for k in kv.keys()]
        elif name.endswith(",Act"):
            # name format: "<main>,<sub>,Act"
            inner = name[:-len(",Act")]
            # Split off the last comma: that's the sub-source. The
            # main name itself may legally not contain a comma in our
            # observed files, but be permissive.
            if "," not in inner:
                # malformed Act block — skip silently rather than crash
                continue
            main_name, sub_name = inner.rsplit(",", 1)
            act_blocks[(main_name, sub_name)] = kv
        else:
            main_kv[name] = kv

    # confidence_sigma from [General] (default 1 — assume single σ if missing)
    try:
        confidence_sigma = int(general_kv.get("Sigma", "1"))
    except ValueError:
        confidence_sigma = 1
    if confidence_sigma <= 0:
        confidence_sigma = 1

    # Now build CertificateSource for each name in [Sets]
    # (we honour the [Sets] order rather than dict insertion).
    sources: dict = {}
    for src_name in set_names:
        kv = main_kv.get(src_name, {})
        unit = (kv.get("Units")
                or kv.get("Activity unit")
                or kv.get("Activity Unit")
                or "Bq").strip()
        ref_dt = _parse_date(kv.get("Date", ""), kv.get("Time", ""))
        sub_names = structure.get(src_name, [])
        sub_objs = []
        for sub_name in sub_names:
            mass_g, activities = _parse_activity_block(
                act_blocks.get((src_name, sub_name), {}),
                confidence_sigma,
                unit,
            )
            sub_objs.append(CertificateSubSource(
                name=sub_name, mass_g=mass_g, activities=activities,
            ))

        sources[src_name] = CertificateSource(
            name=src_name,
            geometry=kv.get("Geometry", "").strip(),
            mass_g=_decimal_comma_to_float(kv.get("Mass,g", "")),
            volume_ml=_decimal_comma_to_float(kv.get("Volume,ml", "")),
            material=kv.get("Material", "").strip(),
            reference_datetime=ref_dt,
            activity_unit=unit,
            thickness_mm=_decimal_comma_to_float(kv.get("Thick,mm", "")),
            comment=kv.get("Comment", "").strip(),
            use_shield=(kv.get("UseShield", "0").strip() == "1"),
            shield_material=kv.get("ShieldMaterial", "").strip(),
            shield_thickness_mm=_decimal_comma_to_float(
                kv.get("ShieldThickness,mm", "")
            ),
            sub_sources=tuple(sub_objs),
        )

    return Certificate(
        file_path=str(file_path),
        confidence_sigma=confidence_sigma,
        sources=sources,
    )


def read_certificate_file(path) -> Certificate:
    """
    Parse one Lsrm `.src` certificate file.

    Args:
        path: path to a `.src` file (CP-1251 encoded).

    Returns:
        Certificate with all sources, sub-sources, and activities
        populated. Empty / unparseable sources still appear (with
        defaults) so caller code can detect them.

    Raises:
        FileNotFoundError: if the file does not exist.
        UnicodeDecodeError: if the file is not CP-1251 encoded.
    """
    path = Path(path)
    raw = path.read_bytes()
    try:
        text = raw.decode("cp1251")
    except UnicodeDecodeError:
        # Some files (rare) may be UTF-8 by accident; try that fallback.
        text = raw.decode("utf-8")
    sections = _tokenise_sections(text)
    return _assemble_certificate(str(path), sections)


def read_certificate_files(paths: Iterable) -> dict:
    """
    Parse multiple .src files and return a dict {filepath_str: Certificate}.

    Caller-friendly bulk-load helper for sessions with several
    certificate files (e.g. a directory containing Marinelli + Petri
    + point-source certificates).
    """
    return {str(p): read_certificate_file(p) for p in paths}


def find_certificate_for_nuclide(certs: Iterable,
                                  nuclide: str,
                                  source_hint: str = "") -> Optional[tuple]:
    """
    Convenience lookup across multiple Certificate objects.

    Args:
        certs: iterable of Certificate (e.g. .values() of the dict
            from `read_certificate_files`).
        nuclide: nuclide name to search for, e.g. "Cs-137".
        source_hint: optional substring to bias source-name matching
            (e.g. "SRC-02" to disambiguate among many Cs-137 sources).

    Returns:
        (CertificateSource, CertificateActivity) — first match — or
        None if no certificate has this nuclide.
    """
    for cert in certs:
        for src in cert.sources.values():
            if source_hint:
                # Require source_hint to match (fuzzy substring)
                if _normalise_token(source_hint) not in _normalise_token(src.name):
                    continue
            act = src.get_activity(nuclide)
            if act is not None:
                return (src, act)
    return None


__all__ = [
    "CertificateActivity",
    "CertificateSubSource",
    "CertificateSource",
    "Certificate",
    "read_certificate_file",
    "read_certificate_files",
    "find_certificate_for_nuclide",
]
