# TESS QLP light curves (MAST HLSP, sectors 43, 44, 71; tess/mast_obs.json): KSPSAP flux, quality==0, wotan biweight 0.75 d with b and c transits masked at the propagated ephemeris;
# per signal: propagated mid-time uncertainty, nominal box depth (in |phase|<T14/2, out 1-12 h), expected S/N for the K2 depth, and the max S/N of a box scan over the +-3 sigma timing window.
import json, numpy as np, warnings, lightkurve as lk; warnings.filterwarnings('ignore')
from wotan import flatten
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); E = {p: (RJ[p]['P_d'][0], max(RJ[p]['P_d'][1:]), RJ[p]['T0_BKJD'][0], max(RJ[p]['T0_BKJD'][1:]), RJ[p]['T14_h'][0]/24, RJ[p]['depth_ppm'][0]*1e-6) for p in 'bc'}
sr = lk.search_lightcurve('TIC 118598568', author='QLP'); print(sr)
open(f'{D}/tess/qlp.url', 'w').write("lightkurve.search_lightcurve('TIC 118598568', author='QLP') -> MAST HLSP QLP\n")
out = {}
for r in sr:
    lc = r.download(download_dir=f'{D}/tess'); s = int(lc.meta['SECTOR'])
    q = lc.quality.value == 0 if hasattr(lc.quality, 'value') else np.asarray(lc.quality) == 0
    t = lc.time.value + 2167.0; f = np.asarray(lc.flux.value, float); m = q & np.isfinite(f) & np.isfinite(t); t, f = t[m], f[m]/np.nanmedian(f[m])
    msk = np.zeros(len(t), bool)
    for p, (P, eP, T0, eT0, d, dep) in E.items(): msk |= np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) < 1.5*d
    tr = flatten(t, np.where(msk, np.nan, f), method='biweight', window_length=0.75, return_trend=True)[1]; ok = np.isfinite(tr); tr = np.interp(t, t[ok], tr[ok]); fl = f/tr
    cad = float(np.median(np.diff(t))*1440); sd = 1.4826*np.median(np.abs(fl - np.median(fl))); out[s] = dict(npts=int(len(t)), cadence_min=cad, mad_rms_ppm=float(sd*1e6), rms_1h_ppm=float(sd/np.sqrt(60/cad)*1e6))
    for p, (P, eP, T0, eT0, d, dep) in E.items():
        n = np.round((np.median(t) - T0)/P); eTp = float(np.hypot(eT0, n*eP)); ph = ((t - T0 + 0.5*P) % P) - 0.5*P; i = np.abs(ph) < 0.5*d; o = (np.abs(ph) > 1/24) & (np.abs(ph) < 0.5)
        e0 = sd/np.sqrt(max(i.sum(), 1)); d0 = np.median(fl[o]) - np.mean(fl[i]) if i.sum() else np.nan
        best = []
        for sh in np.arange(-3*eTp, 3*eTp, d/4):
            ii = np.abs(((t - (T0 + sh) + 0.5*P) % P) - 0.5*P) < 0.5*d
            if ii.sum() >= 3: best.append((np.median(fl) - np.mean(fl[ii]))/(sd/np.sqrt(ii.sum())))
        out[s][p] = dict(n_epochs_from_K2=int(n), sigma_Tmid_h=eTp*24, n_in=int(i.sum()), expected_snr=float(dep/e0), box_depth_ppm=float(d0*1e6), box_err_ppm=float(e0*1e6), scan_max_snr=float(max(best)) if best else None, scan_trials=len(best))
    print(s, out[s], flush=True)
for p in 'bc': out[f'combined_expected_snr_{p}'] = float(np.sqrt(sum(out[s][p]['expected_snr']**2 for s in out if isinstance(s, int))))
json.dump(out, open(f'{D}/tess/qlp.json', 'w'), indent=1, default=str); print(out)
