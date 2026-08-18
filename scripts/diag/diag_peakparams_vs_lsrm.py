# -*- coding: utf-8 -*-
"""Compare MY pipeline peak parameters (FWHM, area, energy) against the LSRM
Gamma-1S ground truth tables the operator provided (2026-06-20). Goal: quantify
the FWHM mismatch ('FWHM ne sootvetstvuet spektru') peak-by-peak."""
import sys, json, glob
import numpy as np
sys.path.insert(0, 'scripts')

cands = sorted(glob.glob('demo_reports/*th232*step3*cont_clamp*/sample_v2/*_report.json'))
RJ = cands[-1]
print('report:', RJ)
rep = json.load(open(RJ, encoding='utf-8'))

# explore: top-level keys
print('TOP KEYS:', list(rep.keys()))

# fit_overlay peaks carry sigma_keV per singlet/component
fo = rep.get('fit_overlay') or {}
print('fit_overlay keys:', list(fo.keys()) if isinstance(fo, dict) else type(fo))
for p in (fo.get('peaks') or [])[:3]:
    print('  sample peak entry keys:', list(p.keys()))
    break

# Pull every peak with an energy + sigma from fit_overlay
rows = []
for p in (fo.get('peaks') or []):
    e = p.get('energy_keV'); sg = p.get('sigma_keV')
    if e and sg:
        fwhm = float(sg) * 2.355
        rows.append((float(e), fwhm, 100.0 * fwhm / float(e),
                     p.get('source'), p.get('label')))
rows.sort()
print('\nMY peaks (fit_overlay): E_keV  FWHM_keV  FWHM%   source  label')
for e, fw, pct, src, lab in rows:
    print('  %8.2f  %7.2f  %6.2f  %-18s %s' % (e, fw, pct, src, lab))
