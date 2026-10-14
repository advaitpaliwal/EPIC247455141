"""oe_null_b.py: empirical calibration of signal b's odd/even sigma. For each reduction (lc/*_m3), all three signals' in-transit cadences removed (|phase|<0.75 T14), then
ptnull_stats.stats (verbatim estimator used for the reported odd/even sigma) at 1000 random fake epochs T0 + U(0.1,0.9) P with b's P and T14 -> distribution of |odd-even|/sigma
for a 'transit' of zero depth. Reports the fraction of fakes with sigma >= the real value (real value from fit3/vet3.json), and the same for d and c."""
import sys, json, numpy as np; sys.path.insert(0, 'k2-deep/247455141/scripts'); import ptnull_stats as S; lcs = S.lcs
D = 'k2-deep/247455141'; V = json.load(open(f'{D}/fit3/vet3.json')); E = {k: (v[0], v[1], v[2]/24) for k, v in V['ephem'].items()}
out = {}; rng = np.random.default_rng(17)
for s in ['everest_m3', 'k2sff_m3', 'fresh_m3', 'fresh40_m3c']:
    z = np.load(f'{D}/lc/{s}.npz'); t, f = z['t'], z['f']; k = np.ones(len(t), bool)
    for o in 'bcd': k &= np.abs(lcs.phase(t, E[o][0], E[o][1])) > 0.75*E[o][2]
    t, f = t[k], f[k]; out[s] = {}
    for w in 'bcd':
        P, T0, T14 = E[w]; oe = np.array([S.stats(t, f, P, T0 + rng.uniform(0.1, 0.9)*P, T14)[1] for _ in range(1000)])
        real = V[s][w]['odd_even_sigma']; out[s][w] = dict(real_oe_sigma=real, frac_fake_ge_real=float(np.mean(oe >= real)), fake_p50=float(np.median(oe)), fake_p95=float(np.percentile(oe, 95)), fake_p99=float(np.percentile(oe, 99)))
    print(s, out[s], flush=True)
json.dump(out, open(f'{D}/fit3/oe_null.json', 'w'), indent=1)
