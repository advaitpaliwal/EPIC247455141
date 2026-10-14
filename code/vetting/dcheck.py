# Third TLS peak after masking b and c (fit/tls3.json, EVEREST: P 9.35604 d, T0 2992.0428, 1.46 h): box depth, odd/even and per-transit chi2 in each reduction at that fixed ephemeris
# (b and c in-transit cadences removed), same estimators as scripts/vet2.py; plus the 300-trial transit-free null of scripts/ptnull.py for its S/N.
import sys, json, numpy as np; sys.path.insert(0, 'k2-deep/247455141/scripts'); import ptnull_stats as S
D = 'k2-deep/247455141'; lcs = S.lcs
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); T3 = json.load(open(f'{D}/fit/tls3.json'))['everest']
P, T0, T14 = T3['P'], T3['T0'], T3['dur_h']/24; out = dict(ephem=[P, T0, T14*24]); rng = np.random.default_rng(5)
for s in ['everest', 'k2sff', 'fresh', 'fresh40']:
    z = np.load(f'{D}/lc/{s}.npz'); t, f = z['t'], z['f']; k = np.ones(len(t), bool)
    for p in 'bc': k &= np.abs(lcs.phase(t, RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0])) > 0.75*RJ[p]['T14_h'][0]/24
    t, f = t[k], f[k]; ph = lcs.phase(t, P, T0); inn = np.abs(ph) < 0.35*T14; oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5)
    d = np.median(f[oot]) - np.mean(f[inn]); e = np.std(f[oot])/np.sqrt(inn.sum()); chi, oe, snr = S.stats(t, f, P, T0, T14)
    kk = np.abs(ph) > 0.75*T14; nul = np.array([S.stats(t[kk], f[kk], P, T0 + rng.uniform(0.1, 0.9)*P, T14)[2] for _ in range(300)])
    out[s] = dict(box_depth_ppm=float(d*1e6), err_ppm=float(e*1e6), n_in=int(inn.sum()), chi2_red=chi, oe_sigma=oe, wsnr=snr, null_snr_p99=float(np.percentile(nul, 99)), null_snr_max=float(nul.max()))
    print(s, out[s], flush=True)
json.dump(out, open(f'{D}/fit/dcheck.json', 'w'), indent=1)
