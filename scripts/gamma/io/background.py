"""
Background-spectrum resolution.

The AtomSpectra primary spectrum carries an embedded BackgroundEnergySpectrum,
and that's the canonical source of the background spectrum for analysis.
The <BackgroundSpectrumFile> string in the same file is a filename hint
pointing to the source of that embedded snapshot — useful for traceability
and for the (rare) case when the embedded copy is missing or stale and we
need to re-read the original.

This module locates the external file given its hint and a set of search
directories. Empty link → returns None silently.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from gamma.spectrum import Spectrum


def resolve_external_background(
    primary: Spectrum,
    search_dirs: Optional[list] = None,
) -> Optional[Spectrum]:
    """
    Try to read the external file referenced by primary.background_link.

    Search order:
      1. Directory of the primary spectrum file
      2. Each path in search_dirs (in order)
      3. Current working directory

    Returns the parsed Spectrum (with is_background=True) on success, None
    if the file cannot be located or read. Does NOT raise on failure —
    we want analysis to continue if the embedded background is present.
    """
    if not primary.background_link:
        return None

    # Late import to avoid circular dependency (atomspectra_xml imports spectrum,
    # spectrum has no imports from io).
    from gamma.io.atomspectra_xml import read_atomspectra_xml

    candidates = []
    if primary.source_path:
        candidates.append(Path(primary.source_path).parent)
    if search_dirs:
        for d in search_dirs:
            candidates.append(Path(d))
    candidates.append(Path.cwd())

    link_name = primary.background_link.strip()

    for d in candidates:
        candidate_path = d / link_name
        if candidate_path.is_file():
            try:
                bg = read_atomspectra_xml(str(candidate_path),
                                          parse_background=False)
                bg.is_background = True
                bg.extras["resolved_from_link"] = primary.background_link
                return bg
            except (ValueError, OSError):
                continue

    return None
