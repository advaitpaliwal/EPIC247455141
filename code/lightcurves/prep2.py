# Build detrended light curves (EVEREST, K2SFF, fresh SFF from TPF) with both signals masked; box depth per signal per reduction with the other signal's in-transit points removed.
import sys, json, numpy as np; sys.path.insert(0, 'k2-deep/247455141/scripts'); import lcs2 as L
D = 'k2-deep/247455141'; epic = 247455141
out = {}
for name in ['everest', 'k2sff', 'fresh', 'fresh40']:
    r = L.fresh(epic, windows=40 if name == 'fresh40' else 20)[:2] if name.startswith('fresh') else getattr(L, name)(epic)
    if r is None: continue
    t, f = r; tf, ff = L.flat(t, f)
    np.savez(f'{D}/lc/{name}.npz', t_raw=t, f_raw=f, t=tf, f=ff)
    out[name] = dict(n=len(tf), rms_ppm=float(np.std(ff)*1e6), cdpp_like_1h_ppm=float(np.std(ff[~L.intr(tf)])*1e6/np.sqrt(2)))
    for i, (P, T0, d) in enumerate(L.EPH):
        oth = L.intr(tf, [L.EPH[1-i]]); dd, e, n = L.boxdepth(tf[~oth], ff[~oth], P, T0, d)
        out[name]['bc'[i]] = dict(box_depth_ppm=float(dd*1e6), box_err_ppm=float(e*1e6), n_in=n, snr=float(dd/e))
    print(name, out[name], flush=True)
json.dump(out, open(f'{D}/lc/depths_box.json', 'w'), indent=1)
