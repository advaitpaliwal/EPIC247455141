# Third-signal check: TLS (same clean() and settings as scripts/inj2.py null mode) on EVEREST and K2SFF raw light curves with signals b and c masked (|phase| < 0.75 T14).
# Reports the top peak's TLS ephemeris, per-transit depths/SNR, odd/even, and the top 5 distinct peaks of the SDE spectrum.
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
from transitleastsquares import transitleastsquares
sys.path.insert(0, 'k2-hunt/scripts'); from k2tls_hunt import clean
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); out = {}
for s in ['everest', 'k2sff']:
    z = np.load(f'{D}/lc/{s}.npz'); t, f = z['t_raw'], z['f_raw']; m = np.zeros(len(t), bool)
    for p in 'bc':
        P, T0, d = RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24; m |= np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) < 0.75*d
    tc, fc = clean(t[~m], f[~m])
    r = transitleastsquares(tc, fc).power(period_min=0.5, period_max=30.0, R_star=0.39, M_star=0.382, R_star_min=0.08, R_star_max=1.0, M_star_min=0.08, M_star_max=1.0,
        oversampling_factor=2, duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
    pk = []; o = np.argsort(r.power)[::-1]
    for j in o:
        if all(abs(r.periods[j] - r.periods[k])/r.periods[k] > 0.02 for k in pk): pk.append(j)
        if len(pk) == 5: break
    out[s] = dict(P=float(r.period), T0=float(r.T0), dur_h=float(r.duration*24), depth_ppm=float((1-r.depth)*1e6), SDE=float(r.SDE), snr=float(r.snr), odd_even_mismatch=float(r.odd_even_mismatch),
        transit_count=int(r.transit_count), distinct_transit_count=int(r.distinct_transit_count), transit_times=[float(x) for x in r.transit_times],
        per_transit_count=[int(x) for x in r.per_transit_count], snr_per_transit=[float(x) for x in r.snr_per_transit], transit_depths_ppm=[float((1-x)*1e6) if np.isfinite(x) else None for x in r.transit_depths],
        top5=[dict(P=float(r.periods[j]), SDE=float(r.power[j])) for j in pk], P_over_Pb=float(r.period/RJ['b']['P_d'][0]), P_over_Pc=float(r.period/RJ['c']['P_d'][0]))
    print(s, json.dumps(out[s]), flush=True)
json.dump(out, open(f'{D}/fit/tls3.json', 'w'), indent=1)
