"""tls3b.py: scripts/tls3.py repeated on every reduction (raw series from lc/*_m3.npz; b and c in-transit cadences masked; k2tls_hunt.clean; same TLS settings), reporting the
top peak, the top 5 distinct peaks and the SDE at the TLS grid period nearest 9.356 d; plus an inverted-flux run (f -> 2 - f) on EVEREST as a systematics false-alarm check."""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
from transitleastsquares import transitleastsquares
sys.path.insert(0, 'k2-hunt/scripts'); from k2tls_hunt import clean
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); out = {}
for s, inv in [('k2sff_m3', 0), ('fresh_m3', 0), ('fresh40_m3', 0), ('everest_m3', 1)]:
    z = np.load(f'{D}/lc/{s}.npz'); t, f = z['t_raw'], z['f_raw']; m = np.zeros(len(t), bool)
    for p in 'bc':
        P, T0, d = RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24; m |= np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) < 0.75*d
    tt, ff = t[~m], f[~m]
    if inv: ff = 2*np.median(ff) - ff
    tc, fc = clean(tt, ff)
    r = transitleastsquares(tc, fc).power(period_min=0.5, period_max=30.0, R_star=0.39, M_star=0.382, R_star_min=0.08, R_star_max=1.0, M_star_min=0.08, M_star_max=1.0,
        oversampling_factor=2, duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
    pk = []
    for j in np.argsort(r.power)[::-1]:
        if all(abs(r.periods[j] - r.periods[k])/r.periods[k] > 0.02 for k in pk): pk.append(j)
        if len(pk) == 5: break
    near = np.abs(r.periods - 9.35604) < 0.01; jd = np.where(near)[0][np.argmax(r.power[near])]
    key = s + ('_inverted' if inv else '')
    out[key] = dict(top_P=float(r.period), top_SDE=float(r.SDE), top_T0=float(r.T0), top5=[dict(P=float(r.periods[j]), SDE=float(r.power[j])) for j in pk],
                    SDE_at_9p356=float(r.power[jd]), P_at=float(r.periods[jd]), rank_of_9p356=int(1 + sum(1 for j in pk if r.power[j] > r.power[jd])))
    print(key, json.dumps(out[key]), flush=True)
json.dump(out, open(f'{D}/fit3/tls3b.json', 'w'), indent=1)
