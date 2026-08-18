"""Peak search and (Phase 2.1) deconvolution."""
from gamma.peaks.search import FoundPeak, mariscotti_search, estimate_fwhm_at_peak
from gamma.peaks.deconvolve import (
    MultipletComponent,
    DeconvolutionResult,
    deconvolve_multiplet,
    find_multiplet_regions,
    deconvolve_identified_multiplets,
    apply_multiplet_deconvolution,
)

__all__ = [
    "FoundPeak", "mariscotti_search", "estimate_fwhm_at_peak",
    "MultipletComponent", "DeconvolutionResult",
    "deconvolve_multiplet", "find_multiplet_regions",
    "deconvolve_identified_multiplets",
    "apply_multiplet_deconvolution",
]
