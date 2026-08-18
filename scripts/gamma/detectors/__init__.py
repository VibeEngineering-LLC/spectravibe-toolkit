"""gamma.detectors — detector-specific path resolvers and constants.

Each spectrometric complex has an isolated subtree under ``detectors/<canonical>/``
(see project README). Algorithms in ``gamma.peaks``, ``gamma.identification``,
``gamma.calibration``, ``gamma.activity`` and ``gamma.physics`` are shared across
detectors and acquire detector-specific paths via the corresponding resolver
module exposed here.

Currently supported (v1.12.0):

- ``gamma.detectors.gamma1s`` — Gamma-1S complex (УДС-ГЦ-63×63 NaI(Tl) by Aspect +
  Lsrm SpectraLine). The only pipeline implemented end-to-end.

Deferred (F-78a roadmap):

- ``gamma.detectors.atomspectra`` (AtomSpectra GS5050)
- ``gamma.detectors.atomnano``
- ``gamma.detectors.radiacode``
"""
