# -*- coding: utf-8 -*-
"""Ground-truth: Tl-208 2614 singlet continuum_grid (rendered podstilauschaya)
vs GROSS sample spectrum. Shows where the magenta continuum line floats ABOVE
the black gross spectrum (operator 2026-06-20: 'right wing above spectrum')."""
import sys, json, glob, os
import numpy as np
sys.path.insert(0, 'scripts')
from gamma.io import format_registry as fr

SP = r'C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe'

# newest th232 rerun report json
cands = sorted(glob.glob('demo_reports/*th232*step3*/sample_v2/*_report.json'))
RJ = cands[-1]
print('report:', RJ)
rep = json.load(open(RJ, encoding='utf-8'))

# locate fit_overlay singlet peak near 2611/2614
fo = rep.get('fit_overlay') or rep.get('charts', {}).get('fit_overlay') or {}
peaks = fo.get('peaks') or []
cand = [p for p in peaks if p.get('source') == 'singlet'
        and abs(float(p.get('energy_keV', 0)) - 2611) < 12]
print('singlet near 2611:', [(p.get('label'), p.get('energy_keV'),
      p.get('sigma_keV')) for p in cand])
pk = cand[0]
grid = pk['continuum_grid']
ge = grid['energies']; gv = grid['values']

# gross spectrum
s = fr.read(SP)
cg = np.asarray(s.counts, float)
def gross_at(e):
    ch = s.energy_to_channel(e)
    if ch is None:
        return float('nan')
    return float(cg[max(0, min(len(cg)-1, int(round(ch))))])

print('\n  E_keV   cont(grid)   gross   cont-gross  ABOVE?')
for e, v in zip(ge, gv):
    g = gross_at(e)
    flag = 'ABOVE' if v > g else ''
    print('%8.1f  %10.1f  %8.1f  %+9.1f  %s' % (e, v, g, v - g, flag))
print('\ncont_left=%.1f (E=%.1f)  cont_right=%.1f (E=%.1f)'
      % (gv[0], ge[0], gv[-1], ge[-1]))
print('gross@left_edge=%.1f  gross@right_edge=%.1f'
      % (gross_at(ge[0]), gross_at(ge[-1])))
