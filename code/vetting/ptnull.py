# Noise calibration of the vet statistics: the per-transit depth chi2 and odd/even sigma of vet2.py, recomputed at 300 transit-free trial epochs per signal and reduction
# (trial T0 = T0 + U(0.1,0.9) P, computed on the light curve with both real signals' in-transit cadences (|phase| < 0.75 T14) removed); the real value's percentile in that distribution.
# Also: Gaia DR3 neighbours (catalogue/gaia_dr3_60as.csv) projected on the TPF grid with the per-pixel difference-image S/N at their positions (logic of ../../247192465/scripts/extra.py).
import sys, json, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, 'k2-deep/scripts'); import lcs
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); E = {p: (RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24) for p in 'bc'}
def stats(t, f, P, T0, T14):
    ph = lcs.phase(t, P, T0); ep = np.round((t - T0)/P).astype(int); oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5); inn = np.abs(ph) < 0.35*T14
    def dep(mi, mo): return np.median(f[mo]) - np.mean(f[mi]), np.std(f[mo])/np.sqrt(mi.sum())
    o = dep(inn & (ep % 2 == 1), oot); e = dep(inn & (ep % 2 == 0), oot); oe = abs(o[0]-e[0])/np.hypot(o[1], e[1])
    ch, ce = [], []
    for k in np.unique(ep[inn]):
        m = inn & (ep == k); mo = oot & (ep == k); d, s = dep(m, mo if mo.sum() > 5 else oot); ch.append(d); ce.append(s)
    ch, ce = np.array(ch), np.array(ce); mu = np.average(ch, weights=1/ce**2)
    return float(np.sum(((ch-mu)/ce)**2)/(len(ch)-1)), float(oe), float(mu/np.sqrt(1/np.sum(1/ce**2)))
out = {}; rng = np.random.default_rng(11)
for s in ['everest', 'k2sff', 'fresh40']:
    z = np.load(f'{D}/lc/{s}.npz'); t0, f0 = z['t'], z['f']
    for p in 'bc':
        o = 'c' if p == 'b' else 'b'; Po, To, Do = E[o]; keep = np.abs(lcs.phase(t0, Po, To)) > 0.75*Do; t, f = t0[keep], f0[keep]
        P, T0, T14 = E[p]; real = stats(t, f, P, T0, T14); null = []
        k2 = np.abs(lcs.phase(t, P, T0)) > 0.75*T14; tn, fn = t[k2], f[k2]
        while len(null) < 300:
            Tx = T0 + rng.uniform(0.1, 0.9)*P
            null.append(stats(tn, fn, P, Tx, T14))
        null = np.array(null)
        out[f'{s}_{p}'] = dict(real_chi2_red=real[0], real_oe_sigma=real[1], real_snr=real[2], null_chi2_red_pct=np.percentile(null[:, 0], [16, 50, 84, 95]).tolist(),
            frac_null_chi2_ge_real=float(np.mean(null[:, 0] >= real[0])), null_oe_pct=np.percentile(null[:, 1], [50, 84, 95]).tolist(), frac_null_oe_ge_real=float(np.mean(null[:, 1] >= real[1])),
            null_snr_pct=np.percentile(null[:, 2], [50, 99, 100]).tolist(), frac_null_snr_ge_real=float(np.mean(null[:, 2] >= real[2])))
        print(s, p, out[f'{s}_{p}'], flush=True)
import lightkurve as lk
tp = lk.read(f'{D}/lc/mastDownload/K2/ktwo247455141-c13_lc/ktwo247455141-c13_lpd-targ.fits.gz'); g = pd.read_csv(f'{D}/catalogue/gaia_dr3_60as.csv').sort_values('sep')
for p in 'bc':
    zz = np.load(f'{D}/pix/diff_{p}.npz'); snr, dm = zz['snr'], zz['depth_map']; L = []
    for _, r in g.iterrows():
        c, rr = tp.wcs.all_world2pix([[r.ra, r.dec]], 0)[0]; ic, ir = int(round(c)), int(round(rr)); ins = 0 <= ir < snr.shape[0] and 0 <= ic < snr.shape[1]
        L.append(dict(source_id=int(r.source_id), G=float(r.phot_g_mean_mag), sep=float(r.sep), col=float(c), row=float(rr), pix_snr=float(snr[ir, ic]) if ins else None, pix_depth_ppm=float(dm[ir, ic]*1e6) if ins else None))
    out[f'neighbours_{p}'] = L
out['tpf_shape'] = list(tp.shape)
json.dump(out, open(f'{D}/fit/ptnull.json', 'w'), indent=1); print(json.dumps({k: v for k, v in out.items() if k.startswith('neigh')}, indent=0)[:4000])
