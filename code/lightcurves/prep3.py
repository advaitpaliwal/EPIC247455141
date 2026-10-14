# Three-signal variant of scripts/prep2.py: rebuild EVEREST, K2SFF, fresh SFF (20 and 40 windows) light curves with b, c AND the 9.356-d signal d masked
# during detrending (wotan biweight 0.75 d; SFF cadence_mask), so the third signal is not partly absorbed by the trend. Writes lc/{name}_m3.npz (new names; prep2 outputs untouched).
import sys, json, numpy as np; sys.path.insert(0, 'k2-deep/247455141/scripts'); import lcs2 as L
D = 'k2-deep/247455141'; epic = 247455141
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); T3 = json.load(open(f'{D}/fit/tls3.json'))['everest']
EPH3 = [(RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24) for p in 'bc'] + [(T3['P'], T3['T0'], T3['dur_h']/24)]
out = dict(eph=EPH3)
for name in ['everest', 'k2sff', 'fresh', 'fresh40']:
    r = L.fresh(epic, eph=EPH3, windows=40 if name == 'fresh40' else 20)[:2] if name.startswith('fresh') else getattr(L, name)(epic)
    t, f = r; tf, ff = L.flat(t, f, eph=EPH3)
    np.savez(f'{D}/lc/{name}_m3.npz', t_raw=t, f_raw=f, t=tf, f=ff)
    out[name] = dict(n=len(tf), rms_ppm=float(np.std(ff)*1e6))
    for i, nm in enumerate('bcd'):
        oth = L.intr(tf, [EPH3[j] for j in range(3) if j != i]); dd, e, n = L.boxdepth(tf[~oth], ff[~oth], *EPH3[i])
        out[name][nm] = dict(box_depth_ppm=float(dd*1e6), box_err_ppm=float(e*1e6), n_in=n, snr=float(dd/e))
    print(name, out[name], flush=True)
json.dump(out, open(f'{D}/fit3/depths_box_m3.json', 'w'), indent=1)
