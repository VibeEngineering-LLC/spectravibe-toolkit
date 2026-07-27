"""
F-301 (v1.17.21, T-021c) — LSRM .cpt XML parser/builder.

`.cpt` — Calibrated Peak Template файл формата SpectraLine (ЛСРМ).
Содержит anchor-based tabulated peak shape (F-299) + метаданные о
детекторе/источнике/калибровке/времени съёма.

Схема (этот модуль определяет canonical structure; реальный .cpt от
SpectraLine может содержать дополнительные поля — мы их игнорируем
при parse и не пишем при build):

  <peak_template version="1.17.21">
      <detector id="Gamma-1S" class="NaI" diameter_mm="63"/>
      <source label="Cs-137" date="2024-05-01"/>
      <anchors>
          <anchor E_keV="661.66" fwhm_keV="46.2" tail_fraction="0.03"
                  tail_slope="0.05" step_height="0.05" asymmetry="0.0"
                  weight="1.0"/>
          ...
      </anchors>
      <notes>...</notes>
  </peak_template>

Encoding: UTF-8.

Forward-compat: при чтении неизвестных атрибутов / тегов — warning
в stderr, не error.

References
----------
- ЛСРМ SpectraLine documentation (calibrated peak template format)
- W3C XML 1.0 spec
"""
from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET
# SEC-01 (P1) hardening: parse untrusted .cpt content via defusedxml to block
# billion-laughs / DOCTYPE entity DoS. Stdlib ET kept for Element/SubElement
# construction in build_cpt_xml (defusedxml is parse-only).
from defusedxml.ElementTree import fromstring as _safe_fromstring
from pathlib import Path
from typing import Optional, Sequence

# Lazy imports to avoid circular deps.


CPT_SCHEMA_VERSION = "1.17.21"


def build_cpt_xml(tabulated_peak_image) -> str:
    """Build .cpt XML string из TabulatedPeakImage.

    Returns
    -------
    UTF-8 string, pretty-printed.
    """
    root = ET.Element("peak_template", {"version": CPT_SCHEMA_VERSION})

    det = ET.SubElement(root, "detector", {
        "id": tabulated_peak_image.detector_id,
        "class": tabulated_peak_image.detector_class,
        "diameter_mm": f"{tabulated_peak_image.crystal_diameter_mm:.2f}",
    })

    if tabulated_peak_image.source_metadata:
        ET.SubElement(root, "source", {
            "label": tabulated_peak_image.source_metadata,
        })

    anchors_el = ET.SubElement(root, "anchors")
    for a in tabulated_peak_image.anchors:
        ET.SubElement(anchors_el, "anchor", {
            "E_keV": f"{a.E_keV:.4f}",
            "fwhm_keV": f"{a.fwhm_keV:.4f}",
            "tail_fraction": f"{a.tail_fraction:.6f}",
            "tail_slope": f"{a.tail_slope_inv_keV:.6f}",
            "step_height": f"{a.step_height_frac:.6f}",
            "asymmetry": f"{a.asymmetry:.6f}",
            "weight": f"{a.weight:.4f}",
        })

    if tabulated_peak_image.notes:
        notes_el = ET.SubElement(root, "notes")
        notes_el.text = tabulated_peak_image.notes

    # Pretty-print: добавим indent (ET.indent доступно с Python 3.9+)
    try:
        ET.indent(root, space="  ")
    except AttributeError:
        pass    # Python < 3.9 — без отступов

    xml_str = ET.tostring(root, encoding="unicode", xml_declaration=False)
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + xml_str


def write_cpt_file(tabulated_peak_image, file_path) -> None:
    """Сохранить TabulatedPeakImage в файл .cpt (UTF-8)."""
    xml = build_cpt_xml(tabulated_peak_image)
    Path(file_path).write_text(xml, encoding="utf-8")


def parse_cpt_xml(xml_str: str, strict: bool = False):
    """Parse .cpt XML string → TabulatedPeakImage.

    Parameters
    ----------
    xml_str : str
        UTF-8 XML content.
    strict : bool
        If True — raise on unknown tags/attrs; otherwise warn to stderr.

    Returns
    -------
    TabulatedPeakImage instance.

    Raises
    ------
    ValueError
        If required structure missing.
    """
    # Lazy import чтобы избежать circular import
    from gamma.peaks.peak_image_tabulated import (
        PeakShapeAnchor, TabulatedPeakImage,
    )

    try:
        # SEC-01: _safe_fromstring blocks DOCTYPE / entity expansion.
        root = _safe_fromstring(xml_str)
    except ET.ParseError as e:
        raise ValueError(f"Invalid XML in .cpt: {e}") from e

    if root.tag != "peak_template":
        raise ValueError(
            f"Expected root <peak_template>, got <{root.tag}>"
        )

    version = root.attrib.get("version", "unknown")

    det_el = root.find("detector")
    if det_el is None:
        raise ValueError("Missing <detector> element")

    detector_id = det_el.attrib.get("id", "unknown")
    detector_class = det_el.attrib.get("class", "NaI")
    try:
        diameter_mm = float(det_el.attrib.get("diameter_mm", "63.0"))
    except ValueError:
        diameter_mm = 63.0
        if strict:
            raise

    source_el = root.find("source")
    source_metadata = (
        source_el.attrib.get("label") if source_el is not None else None
    )

    notes_el = root.find("notes")
    notes = notes_el.text if notes_el is not None else None

    anchors_el = root.find("anchors")
    anchors = []
    if anchors_el is not None:
        for a_el in anchors_el.findall("anchor"):
            try:
                anchors.append(PeakShapeAnchor(
                    E_keV=float(a_el.attrib["E_keV"]),
                    fwhm_keV=float(a_el.attrib["fwhm_keV"]),
                    tail_fraction=float(a_el.attrib.get("tail_fraction", "0")),
                    tail_slope_inv_keV=float(
                        a_el.attrib.get("tail_slope", "0")
                    ),
                    step_height_frac=float(
                        a_el.attrib.get("step_height", "0")
                    ),
                    asymmetry=float(a_el.attrib.get("asymmetry", "0")),
                    weight=float(a_el.attrib.get("weight", "1.0")),
                ))
            except (KeyError, ValueError) as e:
                if strict:
                    raise ValueError(
                        f"Malformed <anchor>: {ET.tostring(a_el)}: {e}"
                    ) from e
                print(
                    f"WARNING: skipping malformed anchor: {e}",
                    file=sys.stderr,
                )

    # Forward-compat: предупредить о неизвестных верхне-уровневых тегах
    known_top_tags = {"detector", "source", "anchors", "notes"}
    for child in root:
        if child.tag not in known_top_tags:
            msg = f"WARNING: unknown .cpt top-level tag <{child.tag}>"
            if strict:
                raise ValueError(msg)
            print(msg, file=sys.stderr)

    return TabulatedPeakImage(
        detector_id=detector_id,
        detector_class=detector_class,
        crystal_diameter_mm=diameter_mm,
        anchors=anchors,
        source_metadata=source_metadata,
        notes=notes,
    )


def read_cpt_file(file_path, strict: bool = False):
    """Загрузить TabulatedPeakImage из .cpt файла (UTF-8)."""
    xml_str = Path(file_path).read_text(encoding="utf-8")
    return parse_cpt_xml(xml_str, strict=strict)


__all__ = [
    "CPT_SCHEMA_VERSION",
    "build_cpt_xml",
    "write_cpt_file",
    "parse_cpt_xml",
    "read_cpt_file",
]
