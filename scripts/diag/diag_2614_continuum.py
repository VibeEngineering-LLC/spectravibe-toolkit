# -*- coding: utf-8 -*-
"""Диагностика fit-view континуума Tl-208 2614 (p2611)."""
import sys, numpy as np
sys.path.insert(0, 'scripts')
from gamma.io import format_registry as fr

SP = r'C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe'
BG = r'C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Фон закр кр\Фон закр кр вода_13.spe'

s = fr.read(SP); b = fr.read(BG)
cs = np.asarray(s.counts, float); cb = np.asarray(b.counts, float)
Ls = float(s.live_time); Lb = float(b.live_time)
n = len(cs)
Es = np.array([s.channel_to_energy(i) for i in range(n)], float)

net = np.zeros(n)
for i in range(n):
    bch = b.energy_to_channel(Es[i])
    j = int(round(bch)) if bch is not None else -1
    bc = cb[j] if 0 <= j < len(cb) else 0.0
    net[i] = cs[i] - bc * (Ls / Lb)

def win(elo, ehi):
    m = (Es >= elo) & (Es <= ehi)
    return np.where(m)[0]

E0 = 2610.72
print('=== SHOULDER CONTINUUM (current method, mean of raw shoulders) ===')
for sig in (49.528, 45.83, 40.0):
    Li = win(E0 - 3*sig, E0 - 2*sig); Ri = win(E0 + 2*sig, E0 + 3*sig)
    L = net[Li].mean() if len(Li) else float('nan')
    R = net[Ri].mean() if len(Ri) else float('nan')
    center = (L + R) / 2.0
    print('  sig=%6.2f  left[%.0f,%.0f]=%7.1f  right[%.0f,%.0f]=%7.1f  lin@center=%7.1f'
          % (sig, E0-3*sig, E0-2*sig, L, E0+2*sig, E0+3*sig, R, center))

sig = 49.528
pk_idx = win(E0 - sig, E0 + sig)
pk_max = net[pk_idx].max(); pk_ch = pk_idx[net[pk_idx].argmax()]
print('\nPEAK net max in +-1sig = %.1f at E=%.1f' % (pk_max, Es[pk_ch]))
area = 421.99 * sig * np.sqrt(2*np.pi)
print('amp(report)=421.99  sigma=%.3f => implied net area=%.0f' % (sig, area))
print('blue top = lin_center + amp = 262.6 + 421.99 = %.1f' % (262.6+421.99))
print('right-floor model = 194 + amp = %.1f' % (194+421.99))

fl_idx = np.concatenate([win(E0-3.2*sig, E0-2.2*sig), win(E0+2.0*sig, E0+3.2*sig)])
print('\nrobust floor (min both shoulders) = %.1f' % net[fl_idx].min())
print('robust floor (25th pct shoulders) = %.1f' % np.percentile(net[fl_idx],25))

print('\n=== NET SPECTRUM SHAPE 2360..2820 keV (per channel) ===')
print('   E_keV     net   sample')
for i in range(n):
    if 2360 <= Es[i] <= 2820:
        print('%8.1f %8.1f %8.0f' % (Es[i], net[i], cs[i]))
