"""
gamma.reporting — Step 11 (Report with diagnostics) assembly.

Four primary deliverables per `references/06_report_format.md`:

* `build_json_report(result)` — machine-readable primary report (dict
  ready for json.dump). Schema version 0.1.
* `build_chat_summary(result, *, json_dict=None, report_path=None)` —
  3–8 line in-chat summary string.
* `build_markdown_report(result, *, json_dict=None, plots=None,
  md_dir=None)` — full Markdown report (only emitted on explicit
  request per the token-economy rule). When `plots` is provided,
  sections 9 & 10 embed the rendered PNGs.
* `build_html_report(result, *, json_dict=None, plots=None)` — same
  content as Markdown but as a self-contained HTML document with
  base64-embedded PNGs (F-86c / v1.15.0).

Plot rendering (F-86a / v1.15.0):

* `build_spectrum_plot(result, output_path)` — cps-log spectrum
  overlay with labelled FEPs, secondary peaks, and optional
  background overlay.
* `build_multiplet_plots(result, output_dir)` — one PNG per resolved
  multiplet cluster.
* `build_all_plots(result, output_dir)` — convenience bundle.

Dispatchers:

* `build_report(result, *, output_dir=None, write_json=True,
  write_markdown=False, write_plots=False, write_html=False,
  return_summary=True)` — writes the requested artefacts to disk
  and returns the paths + summary.
* `analyze_and_report(path, *, output_dir=None, sample_mass_kg=None,
  ...)` — one-call wrapper that runs the orchestrator (Round 5 on by
  default) and assembles the full report bundle. F-86e / v1.15.0.

Environment classifier:

* `classify_environment(result)` — "natural" / "low_background" /
  "unknown" based on Pb K-XRF presence, K-40 / Tl-208 rate, and
  filename hints.

Per the v1.12.0 isolation policy, this package is detector-agnostic
in its API surface but is currently tuned to Gamma-1S NaI 63×63.
"""
from __future__ import annotations

from gamma.reporting.environment import (
    classify_environment, ENV_NATURAL, ENV_LOW_BG, ENV_UNKNOWN,
)
from gamma.reporting.json_report import build_json_report
from gamma.reporting.chat_summary import build_chat_summary
from gamma.reporting.markdown_report import build_markdown_report
from gamma.reporting.build import build_report
from gamma.reporting.wrapper import analyze_and_report

# Optional (matplotlib-dependent) exports — degrade gracefully when
# matplotlib is unavailable.
try:
    from gamma.reporting.plots import (
        build_spectrum_plot,
        build_multiplet_plots,
        build_all_plots,
    )
    _PLOTS_AVAILABLE = True
except ImportError:
    _PLOTS_AVAILABLE = False
    build_spectrum_plot = None    # type: ignore
    build_multiplet_plots = None  # type: ignore
    build_all_plots = None        # type: ignore

# HTML rendering does not require matplotlib (only the data URIs do,
# at which point matplotlib was already used upstream).
from gamma.reporting.html_report import build_html_report


__all__ = [
    "classify_environment",
    "ENV_NATURAL", "ENV_LOW_BG", "ENV_UNKNOWN",
    "build_json_report",
    "build_chat_summary",
    "build_markdown_report",
    "build_html_report",
    "build_report",
    "analyze_and_report",
    "build_spectrum_plot",
    "build_multiplet_plots",
    "build_all_plots",
]
