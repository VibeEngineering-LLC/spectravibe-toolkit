"""
Verify the v1.6 stored_check upgrade on real fixtures.

Key claims to verify:
  - SimpleSqrtFwhm now evaluates correctly (predicts the measured
    cal-peak FWHMs to <1%)
  - Adaptive matching window catches Tl-208 2614 keV on natural-
    background NaI 50×50 (which v1.5 missed because its 10-channel
    window translated to ~6.5 keV, far narrower than the actual
    ~90 keV peak)
  - Backwards-compatible API (no kwargs changes except `match_window_fwhm`)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import math
from gamma.io.readers import read_spectrum
from gamma.peaks.search import mariscotti_search
from gamma.calibration import (
    check_stored_calibration,
    make_fwhm_at_channel_provider,
)


def fmt_match(m):
    return (f"  E_anchor={m['anchor_keV']:7.2f}  found ch={m['matched_peak_channel']:4d} "
            f"E={m['matched_peak_energy_keV']:7.2f}  Δ={m['residual_keV']:+5.2f} keV "
            f"({m['residual_over_fwhm']:.2f}·FWHM, FWHM={m['fwhm_keV']:.1f}keV)")


print("=" * 78)
print("TEST A: Кабинетный фон, 8192 ch NaI, 22.5-day экспозиция")
print("=" * 78)

spec = read_spectrum("detectors/AtomSpectra/data/fixtures/Фон_кабинет_8192к_01-01-2025.xml")
print(f"n_channels={spec.n_channels}, detector={spec.detector_id}")
print(f"Stored FWHM cal peaks: {len(spec.stored_fwhm_calibration.calibration_peaks)}")
for cp in spec.stored_fwhm_calibration.calibration_peaks:
    print(f"  ch={cp.channel:5d}  E={cp.energy_keV:7.2f} keV  FWHM={cp.fwhm_channels:.1f} ch")

# Verify the SimpleSqrtFwhm fix
sf = spec.stored_fwhm_calibration
c0, c1 = sf.coefficients[0], sf.coefficients[1]
print(f"\nSimpleSqrtFwhm coefs: c0={c0:.3f}, c1={c1:.4f}")
print("Verifying fixed interpretation FWHM(N) = √(c0 + c1·N):")
for cp in sf.calibration_peaks:
    arg = c0 + c1 * cp.channel
    predicted = math.sqrt(arg) if arg > 0 else float('nan')
    err_pct = 100.0 * abs(predicted - cp.fwhm_channels) / cp.fwhm_channels if cp.fwhm_channels else 0
    print(f"  ch={cp.channel}: predicted FWHM={predicted:.2f} vs measured {cp.fwhm_channels:.1f}  ({err_pct:.2f}% err)")

# Also verify the wrong v1.5 interpretation would fail
print("\nv1.5 (wrong) interpretation FWHM(N) = c0 + c1·√N for comparison:")
for cp in sf.calibration_peaks:
    wrong = c0 + c1 * math.sqrt(cp.channel)
    print(f"  ch={cp.channel}: v1.5 predicted FWHM={wrong:.2f} (vs actual {cp.fwhm_channels:.1f}) — clearly wrong")

# Build the provider and probe at a few energies
fwhm_at_ch = make_fwhm_at_channel_provider(spec)
print("\nAdaptive FWHM(channel) at representative channels:")
for ch in [100, 500, 1500, 3700, 6073, 7000]:
    E = spec.channel_to_energy(ch)
    fw_ch = fwhm_at_ch(ch)
    # Convert to keV via local slope
    dE_dN = sum(i*a*ch**(i-1) for i,a in enumerate(spec.energy_cal) if i > 0)
    fw_kev = fw_ch * dE_dN
    print(f"  ch={ch:5d}  E={E:7.1f} keV  FWHM={fw_ch:5.1f} ch = {fw_kev:5.1f} keV")

# Run adaptive Mariscotti with stored FWHM model
print("\nMariscotti with adaptive FWHM (from stored model):")
peaks = mariscotti_search(
    spec.counts,
    fwhm_channels=fwhm_at_ch,
    sigma_threshold=3.0,
)
print(f"Found {len(peaks)} peaks")
peaks_by_sig = sorted(peaks, key=lambda p: -p.significance)[:15]
print("Top 15 by significance:")
for p in peaks_by_sig:
    E = spec.channel_to_energy(p.channel)
    print(f"  ch={p.channel:5d}  E={E:7.2f} keV  σ={p.significance:6.1f}  "
          f"local FWHM={p.fwhm_channels:5.1f} ch")

# stored_check with adaptive window
res = check_stored_calibration(spec, peaks)
print(f"\nStored-check result:")
print(f"  Passed: {res.passed}")
print(f"  FWHM source: {res.fwhm_source}")
print(f"  Matched: {res.n_anchors_matched}/{res.n_anchors_tested}")
print(f"  Max residual: {res.max_residual_keV:.2f} keV ({res.max_residual_over_fwhm:.3f}·FWHM)")
print(f"  Reason: {res.reason}")
print("\n  Matched anchors:")
for m in res.matches:
    print(fmt_match(m))

print()
print("=" * 78)
print("TEST B: Алтайское_Зло (известный +48 keV сдвиг stored cal на 2614 keV)")
print("=" * 78)

spec2 = read_spectrum("detectors/AtomSpectra/data/fixtures/Алтайское_Зло_в_домике_маринелли_294_6г.xml")
print(f"n_channels={spec2.n_channels}, detector={spec2.detector_id}")

fwhm_at_ch2 = make_fwhm_at_channel_provider(spec2)
peaks2 = mariscotti_search(
    spec2.counts,
    fwhm_channels=fwhm_at_ch2,
    sigma_threshold=3.0,
)
print(f"Mariscotti adaptive found {len(peaks2)} peaks")

res2 = check_stored_calibration(spec2, peaks2)
print(f"\nStored-check (Алтайское):")
print(f"  Passed: {res2.passed}")
print(f"  FWHM source: {res2.fwhm_source}")
print(f"  Matched: {res2.n_anchors_matched}/{res2.n_anchors_tested}")
print(f"  Max residual: {res2.max_residual_keV:.2f} keV ({res2.max_residual_over_fwhm:.3f}·FWHM)")
print(f"  Reason: {res2.reason}")
print("\n  Matched anchors:")
for m in res2.matches:
    print(fmt_match(m))

# Specifically check Tl-208 2614 keV anchor — this is the one v1.5 misses
print("\n  Looking for Tl-208 2614 keV specifically:")
tl_match = [m for m in res2.matches if abs(m['anchor_keV'] - 2614.51) < 1.0]
if tl_match:
    m = tl_match[0]
    print(f"    ✓ Tl-208 matched: {fmt_match(m).strip()}")
    print(f"    → Δ = {m['residual_keV']:.1f} keV (TODO predicted ~48 keV here)")
else:
    print(f"    ⚠ Tl-208 2614 NOT matched at adaptive 1·FWHM window")
    # check unmatched
    if any(abs(e - 2614.51) < 1.0 for e in res2.unmatched):
        print(f"    Tl-208 IS in unmatched list — anchor was tested but no peak within window")
