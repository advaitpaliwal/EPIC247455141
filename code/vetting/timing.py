# Overlap/alias/harmonic checks between the two signals, and K2 thruster-firing check (logic copied from ../../212726424/scripts/thruster.py: QUALITY bits 20/21).
import numpy as np, glob, json
from astropy.io import fits
from fractions import Fraction
D = 'k2-deep/247455141'
S = {'b': (3.075223, 2989.9491, 1.30/24), 'c': (5.455991, 2990.0929, 1.49/24)}
d = fits.getdata(glob.glob(f'{D}/lc/mastDownload/K2/*/*llc.fits')[0], 1); t, q = d['TIME'], d['SAP_QUALITY']; ok = np.isfinite(t); t, q = t[ok], q[ok]
out = {'campaign_span_BKJD': [float(t.min()), float(t.max())]}
tm = {k: P0 + 0 for k, P0 in []}
for k, (P, T0, d14) in S.items():
    n = np.arange(np.ceil((t.min()-T0)/P), np.floor((t.max()-T0)/P)+1); tm[k] = T0 + n*P
    out[f'{k}_transit_times'] = tm[k].round(4).tolist()
sep = np.abs(tm['b'][:, None] - tm['c'][None, :]); hw = 0.5*(S['b'][2] + S['c'][2])
ii = np.argwhere(sep < 0.5)
out['pairs_within_0.5d'] = [dict(tb=float(tm['b'][i]), tc=float(tm['c'][j]), sep_h=float(sep[i, j]*24), overlap=bool(sep[i, j] < hw)) for i, j in ii]
out['n_overlapping_events(sep<(T14b+T14c)/2)'] = int((sep < hw).sum())
r = S['c'][0]/S['b'][0]; out['Pc_over_Pb'] = r
out['nearest_small_ratios'] = sorted([(abs(r - n/m)/r, f'{n}/{m}') for n in range(1, 13) for m in range(1, 13) if Fraction(n, m).denominator == m][:0] or
                                     [(abs(r - n/m)/r, f'{n}/{m}') for n in range(1, 13) for m in range(1, 13)])[:5]
# fraction of each signal's transits that fall inside the other's transit window (if one signal were a subset of the other)
for k, o in [('b', 'c'), ('c', 'b')]:
    P, T0, dd = S[o]; ph = ((tm[k] - T0 + 0.5*P) % P) - 0.5*P; out[f'frac_{k}_events_in_{o}_window'] = float(np.mean(np.abs(ph) < dd))
fire = ((q & (1 << 19)) | (q & (1 << 20))) > 0; tf = t[fire]; dt = np.diff(tf); dt = dt[(dt > 0.1) & (dt < 0.4)]; Pthr = float(np.median(dt))
out['thruster'] = dict(n_cad=int(len(t)), n_fire=int(fire.sum()), P_thr_d=Pthr, P_thr_h=Pthr*24)
for k, (P, T0, dd) in S.items():
    ph = ((t - T0 + 0.5*P) % P) - 0.5*P; intr = np.abs(ph) < 0.5*dd; phf = ((tf - T0 + 0.5*P) % P) - 0.5*P
    out['thruster'][k] = dict(P_over_Pthr=P/Pthr, frac_fire_in_transit_window=float(np.mean(np.abs(phf) < 0.5*dd)), expected=float(dd/P),
        flagged_frac_in=float(np.mean(fire[intr])), flagged_frac_out=float(np.mean(fire[~intr])), firing_phase_coherence_on_P=float(np.abs(np.mean(np.exp(2j*np.pi*phf/P)))),
        coherence_expected_random=float(1/np.sqrt(len(tf))), **{f'{m}P_mod_Pthr_frac': float(((m*P/Pthr) + 0.5) % 1 - 0.5) for m in [1, 2, 3]})
json.dump(out, open(f'{D}/fit/timing.json', 'w'), indent=1, default=str); print(json.dumps(out, indent=1, default=str))
