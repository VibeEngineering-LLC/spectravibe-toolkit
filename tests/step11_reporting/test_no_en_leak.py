"""F-108 (D-04, D-05, D-06) — no English leak in HTML/MD body.

Generates the Th-232 interactive HTML + Markdown reports. After
stripping CSS / JS / inline JSON data, every ASCII word of ≥4 letters
in the remaining body text must be in a permissive whitelist
(units, F-IDs, nuclide stems, file extensions).

Also asserts the absence of glossary tokens that should always be
translated:
  - ``gain drift``
  - ``trump card``
  - ``WARNING``
"""
from __future__ import annotations

import os
import re
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


# Permissive whitelist for ASCII words ≥4 letters seen in body text.
_WHITELIST = {
    # units
    "keV", "cps", "Bq", "MDA", "FWHM", "ROI", "TCS", "XRF", "HPGe", "NaI",
    "LaBr", "CeBr", "CdZnTe", "LaBr3", "ISO", "IAEA", "LSRM", "ERN",
    "Marinelli",
    # chart.js / plugin (these appear ONLY in <script>, but if they
    # leak as text we still keep them whitelisted)
    "Chart", "Annotation",
    # filename extensions
    "html", "json", "xml", "txt", "csv",
    # math symbols
    "alpha", "beta", "gamma", "sigma", "chi",
    # tier
    "fp-long", "fp-nat", "fp-phys", "fp-mp", "fp-tbl",
    # detector-name family
    "Gamma",
    # mode names / English labels we keep as-is for now
    "sample", "anchor", "rank",
    # F-160 (2026-06-20) — source-label идентификатор для FWHM-модели,
    # загруженной из references/lsrm_ground_truth/. Это технический
    # идентификатор такого же класса как `lsrm_peaks_table_quadratic`
    # и `default_NaI_63x63`, не narrative-English.
    "ground", "truth", "reference",
    # F-117 / F-118 (v1.17.5): scientific solver acronyms in multiplet
    # methodology card («NNLS / lsq_linear»). Both are method names, not
    # narrative — keep verbatim.
    "NNLS", "linear",
    # F-452 (2026-06-21) — source-label fragments for LSRM poly-4 sqrt(E)
    # FWHM model (`lsrm_ground_truth_reference_poly4_sqrtE`). Same class
    # of technical identifier as `lsrm_peaks_table_quadratic` /
    # `lsrm_ground_truth_reference` (parts already whitelisted via
    # `LSRM` + `ground/truth/reference`); these two complete the label.
    # Not narrative English — method-shape identifier.
    "poly", "sqrtE",
    # K-18 / F-117 — area-method labels surfaced in cells (имя метода
    # интегрирования — собственное имя автора алгоритма)
    "Cowell",
    # F-260 / F-321 (v1.18.15) — bilingual narrator EN equivalents
    # (русский_термин «(english_equivalent)»). Это legitimate enrichment
    # из data/glossary_budyka_2021.json — не leak.
    "Activity", "Detector", "Multiplet", "Spectrometry", "Spectrometer",
    "Full", "half", "maximum", "width",           # из FWHM expansion
    "Background", "Calibration", "Energy", "Channel", "Peak",
    "Resolution", "Efficiency", "Identification", "Deconvolution",
    "Cascade", "Coincidence", "Summing", "Annihilation",
    # F-319 (v1.18.15) — EN→ru замены в (en) скобках
    "double", "single", "escape", "plateau", "backscatter",
    "pile", "pileup",
    # F-070 W3 / v1.24.0 — visual similarity card: template ID fragments +
    # methodology reference (technical identifiers, not narrative English).
    # Template IDs follow pattern VT-NUCLIDE-GEOMCODE-YEAR; geometry codes:
    # POINT5CM, PETRI60ML (→ PETRI, POINT seen as fragments), DENTA100ML,
    # DENTA120ML (→ DENTA), MARI0CM (→ MARI), MERGED.
    "MARI", "POINT", "PETRI", "DENTA", "MERGED",
}


_F_RULE_RE = re.compile(r"\bF-\d+[a-z]?\b", re.IGNORECASE)
_NUCLIDE_RE = re.compile(
    r"\b(?:Cs|Co|Am|Eu|Pb|Bi|Tl|Ac|Ra|Rn|Th|Po|Ba|Na|K|U|I|Be|Mn|Mo|Sn|"
    r"Sr|Y|Zr|Nb|Tc|Ru|Rh|Pd|Ag|Cd|In|Sb|Te|Xe|Cs|Hf|Ta|W|Re|Os|Ir|Pt|Au|"
    r"Hg|Pa|Np|Pu)-\d+m?\b"
)
_NUMBER_RE = re.compile(r"\b\d+(?:[\.,]\d+)?\b")
_ASCII_WORD_RE = re.compile(r"[A-Za-z]{4,}")


def _strip_html_body(html: str) -> str:
    """Strip <head>, <style>, <script> and inline `const xxx = [...]` arrays."""
    # Strip <head>...</head>
    html = re.sub(r"<head\b.*?</head>", "", html, flags=re.IGNORECASE | re.DOTALL)
    # Strip <style>
    html = re.sub(r"<style\b.*?</style>", "", html, flags=re.IGNORECASE | re.DOTALL)
    # Strip <script>
    html = re.sub(r"<script\b.*?</script>", "", html, flags=re.IGNORECASE | re.DOTALL)
    # F-320 / v1.18.15 — Strip ГОСТ-references section (foreign-language
    # publication titles like "Gilmore G. R., Hemingway J. D. Practical
    # Gamma-ray Spectrometry. — Chichester: Wiley, 2024." ДОЛЖНЫ быть в
    # оригинале по ГОСТ Р 7.0.5–2008 §5.2). EN authors/publishers тут
    # — НЕ leak, а соответствие нормативу.
    html = re.sub(
        r'<section\s+class="gost-references".*?</section>',
        "", html, flags=re.IGNORECASE | re.DOTALL,
    )
    # F-326 / v1.18.18.1 — Strip passport-comparison section: contains
    # legitimate Python code block (`analyze_and_report`, `passport_activity_Bq`,
    # `sample_mass_kg`) и инструкцию для пользователя. Это deliberate
    # code example, не EN-leak в narrative.
    html = re.sub(
        r'<section\s+class="passport-comparison".*?</section>',
        "", html, flags=re.IGNORECASE | re.DOTALL,
    )
    # F-325 / v1.18.18.1 — Strip "файл фона: <filename>" pattern, потому что
    # filenames в legacy datasets могут содержать ASCII токены
    # (water, background, marinelli), которые НЕ являются narrative leak.
    html = re.sub(
        r'файл фона:\s*[\w\-.]+',
        "файл фона: «удалено»", html, flags=re.IGNORECASE,
    )
    html = re.sub(
        r'файл образца:\s*[\w\-.]+',
        "файл образца: «удалено»", html, flags=re.IGNORECASE,
    )
    # Strip HTML tags but keep text
    text = re.sub(r"<[^>]+>", " ", html)
    return text


def _strip_md_codefences(md: str) -> str:
    """Strip ```fenced code``` blocks from Markdown."""
    md = re.sub(r"```.*?```", " ", md, flags=re.DOTALL)
    md = re.sub(r"`[^`]+`", " ", md)
    # F-320 / v1.18.15 — Strip ГОСТ-references section + F-327 relabel
    # (F-327 переименовал «Список использованных источников» →
    # «Список использованной литературы»).
    md = re.sub(
        r"## Список использованной литературы.*$",
        "", md, flags=re.DOTALL,
    )
    md = re.sub(
        r"## Список использованных источников.*$",
        "", md, flags=re.DOTALL,
    )
    # F-326 / v1.18.18.1 — Strip passport comparison section (содержит
    # Python code-block с EN keywords).
    md = re.sub(
        r"## Сравнение с паспортной удельной активностью.*?(?=\n## |\n---|\Z)",
        "", md, flags=re.DOTALL,
    )
    # F-325 / v1.18.18.1 — Strip filename mentions в любом формате
    # (inline "файл фона: ..." или MD-table row "| Имя файла фона | ... |").
    # Filenames в legacy datasets могут содержать ASCII токены: water,
    # background, marinelli — не narrative leak.
    md = re.sub(
        r"(файл фона|файл образца):\s*[\w\-./\\]+",
        r"\1: «удалено»", md, flags=re.IGNORECASE,
    )
    md = re.sub(
        r"^\|.*(?:Имя файла|файл образца|файл фона|sample_filename|background_filename).*\|.*$",
        "", md, flags=re.IGNORECASE | re.MULTILINE,
    )
    return md


def _check_no_en_leak(text: str, where: str) -> None:
    probe = _F_RULE_RE.sub(" ", text)
    probe = _NUCLIDE_RE.sub(" ", probe)
    probe = _NUMBER_RE.sub(" ", probe)
    bad = []
    for m in _ASCII_WORD_RE.finditer(probe):
        word = m.group(0)
        if word in _WHITELIST:
            continue
        # Allow whitelist match by stem (case-insensitive)
        if word.lower() in {w.lower() for w in _WHITELIST}:
            continue
        bad.append(word)
    # Strict glossary tokens that must NEVER appear (post-translation).
    forbidden_glossary = ["gain drift", "trump card", "WARNING"]
    failing = []
    for tok in forbidden_glossary:
        if tok in text:
            failing.append(tok)
    if failing:
        raise AssertionError(
            f"glossary tokens leaked in {where}: {failing}"
        )
    if bad:
        # Trim duplicates for readable error
        uniq = sorted(set(bad))[:30]
        raise AssertionError(
            f"English ASCII words leaked in {where} (first 30 unique): {uniq}"
        )


def test_no_en_leak(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_no_en").
    out = str(tmp_path)
    sp = (
        "detectors/Gamma-1S/reference_spectra/"
        "archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    res = analyze_and_report(
        sp,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )

    html = open(res["html"], encoding="utf-8").read()
    body = _strip_html_body(html)
    _check_no_en_leak(body, "HTML body")

    md = open(res["markdown"], encoding="utf-8").read()
    md_body = _strip_md_codefences(md)
    _check_no_en_leak(md_body, "Markdown body")


if __name__ == "__main__":
    import tempfile, pathlib
    test_no_en_leak(pathlib.Path(tempfile.mkdtemp(prefix="_test_no_en_")))
    print("OK")
