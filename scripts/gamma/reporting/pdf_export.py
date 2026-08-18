"""
PDF export via Microsoft Edge headless (F-114 / D-12, v1.17.4).

Generates ``{stem}_report.pdf`` next to the source HTML.  Returns
the PDF path on success, ``None`` if Edge wasn't found.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional


_EDGE_CANDIDATES = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    "msedge",
    "microsoft-edge",
]


def _find_edge() -> Optional[str]:
    for candidate in _EDGE_CANDIDATES:
        if os.path.isabs(candidate):
            if os.path.exists(candidate):
                return candidate
        else:
            p = shutil.which(candidate)
            if p:
                return p
    return None


def html_to_pdf(html_path: str, pdf_path: Optional[str] = None,
                virtual_time_budget_ms: int = 5000) -> Optional[str]:
    """Render ``html_path`` to PDF using Edge headless.

    Returns the PDF path on success, ``None`` if Edge is unavailable.
    """
    edge = _find_edge()
    if not edge:
        return None

    html_path = str(Path(html_path).resolve())
    if pdf_path is None:
        pdf_path = str(Path(html_path).with_suffix(".pdf"))
    pdf_path = str(Path(pdf_path).resolve())

    url = "file:///" + html_path.replace("\\", "/")
    cmd = [
        edge,
        "--headless=new",
        "--disable-gpu",
        f"--virtual-time-budget={int(virtual_time_budget_ms)}",
        f"--print-to-pdf={pdf_path}",
        "--no-pdf-header-footer",
        url,
    ]
    try:
        subprocess.run(
            cmd, check=False, timeout=60,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except Exception:
        return None
    if os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0:
        return pdf_path
    return None


__all__ = ["html_to_pdf"]
