# -*- coding: utf-8 -*-
"""Воспроизвести Cowell-площадь Tl-208 2614 + увидеть реальную базу под пиком."""
import sys, numpy as np
sys.path.insert(0, 'scripts')
from gamma.io import format_registry as fr
from gamma.peaks.area import cowell_area

SP = r'C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe'
BG = r'C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Фон закр кр\Фон закр кр вода_13.spe'
s = fr.read(SP); b = fr.read(BG)
cs = np.asarray(s.counts, float); cb = np.asarray(b.counts, float)
Ls = float(s.live_time); Lb = float(b.live_time)
n = len(cs)
Es = np.array([s.channel_to_energy(i) for i in range(n)], float)
net = np.zeros(n)
for i in range(n):
    bch = b.energy_to_channel(Es[i]); j = int(round(bch)) if bch is not None else -1
    net[i] = cs[i] - (cb[j] if 0 <= j < len(cb) else 0.0) * (Ls / Lb)

pk = 876
dE = Es[pk+1] - Es[pk]
fwhm_ch = 116.63958416192543 / dE
print('peak_ch=%d  dE/ch=%.3f keV  fwhm_ch=%.2f' % (pk, dE, fwhm_ch))

for label, arr in (('NET', net), ('SAMPLE', cs)):
    r = cowell_area(arr, peak_channel=pk, fwhm_channels=fwhm_ch)
    roi_x = np.arange(r.roi_low_ch, r.roi_high_ch)
    w = int(round(0.3 * (r.roi_high_ch - r.roi_low_ch)))
    wing_idx = np.concatenate([np.arange(0, w), np.arange(len(roi_x)-w, len(roi_x))])
    coefs = np.polyfit(roi_x[wing_idx], arr[r.roi_low_ch:r.roi_high_ch][wing_idx], 1)
    base_center = float(np.polyval(coefs, pk))
    base_left = float(np.polyval(coefs, r.roi_low_ch))
    base_right = float(np.polyval(coefs, r.roi_high_ch-1))
    print('\n[%s] ROI ch[%d..%d] = E[%.0f..%.0f] keV' % (
        label, r.roi_low_ch, r.roi_high_ch, Es[r.roi_low_ch], Es[min(r.roi_high_ch-1,n-1)]))
    print('   gross=%.0f  baseline_total=%.0f  NET=%.0f  conv=%s' % (
        r.gross_counts, r.baseline_counts, r.net_area_counts, r.converged))
    print('   baseline line: left=%.0f  center@pk=%.0f  right=%.0f' % (base_left, base_center, base_right))
    print('   measured net@pk=%.0f  => FEP height above Cowell base = %.0f' % (net[pk], net[pk]-base_center))
