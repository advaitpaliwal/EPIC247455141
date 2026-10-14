"""Three-signal vet, estimators copied from scripts/vet2.py (box depth, odd/even, per-transit chi2, secondary at phase 0.5, phase scan) and ptnull_stats.stats (weighted S/N) with a
300-trial transit-free null; every OTHER signal's in-transit cadences (|phase| < 0.75 T14) removed first. Light curves: lc/*_m3 (b, c, d all masked in detrending).
Ephemerides: b, c from fit/result_everest_rho.json; d from fit/tls3.json (EVEREST TLS). Also: thruster-firing coincidence of d's in-transit cadences (TPF QUALITY flags) and period-alias table.
usage: vet3.py   -> fit3/vet3.json"""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, 'k2-deep/scripts'); import lcs
sys.path.insert(0, 'k2-deep/247455141/scripts'); import ptnull_stats as S
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); T3 = json.load(open(f'{D}/fit/tls3.json'))['everest']
E = {p: (RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24) for p in 'bc'}; E['d'] = (T3['P'], T3['T0'], T3['dur_h']/24)
def dep(f, mi, mo): return float((np.median(f[mo]) - np.mean(f[mi]))*1e6), float(np.std(f[mo])/np.sqrt(mi.sum())*1e6), int(mi.sum())
def vet(t, f, P, T0, T14):
    ph = lcs.phase(t, P, T0); ep = np.round((t - T0)/P).astype(int); oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5); inn = np.abs(ph) < 0.35*T14
    r = dict(all=dep(f, inn, oot), odd=dep(f, inn & (ep % 2 == 1), oot), even=dep(f, inn & (ep % 2 == 0), oot))
    r['odd_even_sigma'] = abs(r['odd'][0]-r['even'][0])/np.hypot(r['odd'][1], r['even'][1]); r['odd_minus_even_ppm'] = r['odd'][0]-r['even'][0]
    per = []
    for e in np.unique(ep[inn]):
        m = inn & (ep == e); mo = oot & (ep == e); dd = dep(f, m, mo if mo.sum() > 5 else oot); per.append(dict(epoch=int(e), tmid=float(T0+e*P), depth_ppm=dd[0], err_ppm=dd[1], n=dd[2]))
    ch = np.array([p['depth_ppm'] for p in per]); ce = np.array([p['err_ppm'] for p in per]); mu = np.average(ch, weights=1/ce**2)
    r['per_transit'] = per; r['per_transit_chi2'] = float(np.sum(((ch-mu)/ce)**2)); r['per_transit_dof'] = len(ch)-1
    from scipy.stats import chi2; r['per_transit_p'] = float(chi2.sf(r['per_transit_chi2'], len(ch)-1))
    ph2 = lcs.phase(t, P, T0 + P/2); r['secondary_phase0.5'] = dep(f, np.abs(ph2) < 0.5*T14, (np.abs(ph2) > T14) & (np.abs(ph2) < 0.5) & (np.abs(ph) > T14))
    sc = []
    for g in np.arange(0.05, 0.95, 0.005):
        p3 = lcs.phase(t, P, T0 + g*P); m = (np.abs(p3) < 0.5*T14) & (np.abs(ph) > T14)
        if m.sum() >= 4: sc.append((g, *dep(f, m, oot)))
    sc = np.array(sc); snr = sc[:, 1]/sc[:, 2]
    r['phase_scan_max_snr'] = float(snr.max()); r['phase_scan_argmax'] = float(sc[np.argmax(snr), 0]); r['phase_scan_std_snr'] = float(np.std(snr)); r['primary_snr'] = r['all'][0]/r['all'][1]
    chi, oe, ws = S.stats(t, f, P, T0, T14); rng = np.random.default_rng(5); kk = np.abs(ph) > 0.75*T14
    nul = np.array([S.stats(t[kk], f[kk], P, T0 + rng.uniform(0.1, 0.9)*P, T14)[2] for _ in range(300)])
    r['wsnr'] = ws; r['null_wsnr_p99'] = float(np.percentile(nul, 99)); r['null_wsnr_max'] = float(nul.max()); r['null_frac_ge'] = float(np.mean(nul >= ws))
    return r
out = dict(ephem={k: [v[0], v[1], v[2]*24] for k, v in E.items()})
for s in ['everest_m3', 'k2sff_m3', 'fresh_m3', 'fresh40_m3c']:
    z = np.load(f'{D}/lc/{s}.npz'); out[s] = {}
    for w in 'dbc':
        t, f = z['t'], z['f']; k = np.ones(len(t), bool)
        for o in 'bcd':
            if o != w: k &= np.abs(lcs.phase(t, E[o][0], E[o][1])) > 0.75*E[o][2]
        out[s][w] = vet(t[k], f[k], *E[w])
        print(s, w, {x: out[s][w][x] for x in ['all', 'odd', 'even', 'odd_even_sigma', 'per_transit_chi2', 'per_transit_dof', 'per_transit_p', 'secondary_phase0.5', 'phase_scan_max_snr', 'wsnr', 'null_wsnr_max']}, flush=True)
# thruster firings: TPF QUALITY flags (lightkurve KeplerQualityFlags names)
from astropy.io import fits; import glob
from lightkurve.utils import KeplerQualityFlags as Q
tp = glob.glob(f'{D}/lc/mastDownload/K2/ktwo247455141-c13_lc/*lpd-targ.fits.gz')[0]; h = fits.getdata(tp, 1); tt, qq = h['TIME'], h['QUALITY']; ok = np.isfinite(tt); tt, qq = tt[ok], qq[ok]
thr = (qq & (Q.ThrusterFiring | Q.PossibleThrusterFiring)) > 0; tf = tt[thr]
dt = np.diff(tf); dtm = float(np.median(dt[dt > 0.1])) if np.any(dt > 0.1) else None
P, T0, T14 = E['d']; ph = lcs.phase(tt, P, T0); inn = np.abs(ph) < 0.5*T14
near = np.array([np.min(np.abs(tf - x)) < 0.021 for x in tt])   # within ~1 long cadence of a flagged firing
tr = []
for e in range(int(np.floor((tt.min()-T0)/P)), int(np.ceil((tt.max()-T0)/P))+1):
    tc = T0 + e*P
    if tt.min() < tc < tt.max(): tr.append(dict(epoch=e, tmid=float(tc), dt_nearest_firing_h=float(np.min(np.abs(tf - tc))*24)))
out['thruster'] = dict(n_flagged=int(thr.sum()), median_spacing_d=dtm, frac_all_cadences_near_firing=float(near.mean()), frac_d_intransit_near_firing=float(near[inn].mean()), n_d_intransit=int(inn.sum()),
    P_over_spacing=float(P/dtm) if dtm else None, per_transit=tr)
Pb, Pc, Pd = E['b'][0], E['c'][0], E['d'][0]
out['alias'] = dict(Pd_Pc=Pd/Pc, Pd_Pb=Pd/Pb, Pc_Pb=Pc/Pb, Pd_2Pb=Pd/(2*Pb), Pd_3Pb=Pd/(3*Pb), Pd_2Pc=Pd/(2*Pc), Pd_half=Pd/2, near_rational={f'{a}/{b}': Pd/Pc - a/b for a, b in [(5, 3), (7, 4), (12, 7), (19, 11)]},
    d_epochs_overlap_b=[float(np.min(np.abs(lcs.phase(np.array([T0 + e*Pd for e in range(9)]), Pb, E['b'][1]))*24))], d_epochs_overlap_c=[float(np.min(np.abs(lcs.phase(np.array([T0 + e*Pd for e in range(9)]), Pc, E['c'][1]))*24))])
print(json.dumps({k: out[k] for k in ['thruster', 'alias']}, indent=1))
json.dump(out, open(f'{D}/fit3/vet3.json', 'w'), indent=1)
