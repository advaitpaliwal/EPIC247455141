"""Copy of ../../scripts/pix.py with (a) the other signal's in-transit cadences (|phase|<0.75 T14) excluded from in-transit, reference and fake-epoch sets (optional args OP OT0 OT14_H; OP=0 disables), and (b) the null test of ../../220179502/scripts/pix_null.py appended (centroid offset of each transit-free fake-epoch difference image, same estimator).
Matched-pointing difference image + per-pixel depth map for a K2 TPF.
usage: pix.py TPF P T0_BKJD T14_H RA DEC OUTPREFIX
For each in-transit cadence, average the 8 out-of-transit cadences within +-1.5 d whose (POS_CORR1,POS_CORR2) is closest; diff = OOT - IN.
Noise: same estimator at 60 fake epochs (phase offsets avoiding transit), per-pixel std."""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
tp_path, P, T0, T14, RA, DEC, OUT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])/24, float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
OP, OT0, OT14 = (float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])/24) if len(sys.argv) > 10 else (0, 0, 0)
tp = lk.read(tp_path, quality_bitmask='hard')
t = tp.time.value; F = tp.flux.value; pc1 = tp.pos_corr1; pc2 = tp.pos_corr2
ok = np.isfinite(t) & np.isfinite(pc1) & np.isfinite(pc2) & np.all(np.isfinite(F.reshape(len(t), -1)) | True, axis=1)
t, F, pc1, pc2 = t[ok], F[ok], pc1[ok], pc2[ok]
oth = (np.abs(((t - OT0 + 0.5*OP) % OP) - 0.5*OP) < 0.75*OT14) if OP > 0 else np.zeros(len(t), bool)
# remove slow trend per pixel? local matching within 1.5 d keeps it small
def diffimg(T0x):
    ph = ((t - T0x + 0.5*P) % P) - 0.5*P
    phr = ((t - T0 + 0.5*P) % P) - 0.5*P
    intr = np.where((np.abs(ph) < 0.35*T14) & ~oth)[0]
    oot_ok = (np.abs(phr) > 0.75*T14) & ~oth
    D = []; O = []
    for i in intr:
        cand = np.where(oot_ok & (np.abs(t - t[i]) < 1.5) & (np.abs(ph) > 0.75*T14))[0]
        if len(cand) < 8: continue
        d = np.hypot(pc1[cand]-pc1[i], pc2[cand]-pc2[i]); sel = cand[np.argsort(d)[:8]]
        o = np.nanmean(F[sel], axis=0); D.append(o - F[i]); O.append(o)
    if len(D) < 3: return None, None, 0
    return np.nanmean(D, axis=0), np.nanmean(O, axis=0), len(D)
diff, oot, n = diffimg(T0)
rng = np.random.default_rng(1); fakes = []
for k in range(200):
    off = rng.uniform(0.15, 0.85)*P
    if abs(((off + 0.5*P) % P) - 0.5*P) < 3*T14: continue
    d, _, nn = diffimg(T0 + off)
    if d is not None and nn >= 0.6*n: fakes.append(d)
    if len(fakes) >= 60: break
fakes = np.array(fakes); sig = np.nanstd(fakes, axis=0)
snr = diff/sig; depth_map = diff/oot
# aperture (threshold) and target pixel from WCS
ap = tp.create_threshold_mask(threshold=3, reference_pixel='center')
from astropy.wcs import WCS
w = tp.wcs; cx, cy = w.all_world2pix([[RA, DEC]], 0)[0]   # column, row in cutout coordinates
yy, xx = np.mgrid[:diff.shape[0], :diff.shape[1]]
def cen(img, m):
    im = np.where(m, img, 0); im = np.clip(im, 0, None); s = im.sum(); return (xx*im).sum()/s, (yy*im).sum()/s
m3 = (np.hypot(xx-cx, yy-cy) < 3.5)
dcx, dcy = cen(diff, m3); ocx, ocy = cen(oot, m3)
# bootstrap centroid error from fakes added to diff
bc = []
for fk in fakes: bc.append(cen(diff + (fk - np.nanmean(fakes, axis=0)), m3))
bc = np.array(bc); e_dc = np.nanstd(bc, axis=0)
# depth vs aperture
def apdepth(mask):
    return float(np.nansum(diff[mask])/np.nansum(oot[mask]))
core = np.zeros_like(ap); idx = np.argsort(np.where(ap, oot, -np.inf).ravel())[::-1][:4]; core.flat[idx] = True
big = tp.create_threshold_mask(threshold=1, reference_pixel='center')
ring = big & ~core
res = dict(tpf=tp_path, n_in_cadences=int(n), n_fake=int(len(fakes)), target_pix_colrow=[float(cx), float(cy)],
           oot_centroid=[float(ocx), float(ocy)], diff_centroid=[float(dcx), float(dcy)], diff_centroid_err=e_dc.tolist(),
           offset_pix=[float(dcx-ocx), float(dcy-ocy)], offset_arcsec=float(np.hypot(dcx-ocx, dcy-ocy)*3.98),
           offset_sigma=float(np.hypot((dcx-ocx)/e_dc[0], (dcy-ocy)/e_dc[1])),
           offset_from_gaia_pix=[float(dcx-cx), float(dcy-cy)],
           peak_diff_pixel=[int(np.unravel_index(np.nanargmax(np.where(np.isfinite(snr), snr, -1e9)), snr.shape)[1]), int(np.unravel_index(np.nanargmax(np.where(np.isfinite(snr), snr, -1e9)), snr.shape)[0])],
           brightest_oot_pixel=[int(np.unravel_index(np.nanargmax(oot), oot.shape)[1]), int(np.unravel_index(np.nanargmax(oot), oot.shape)[0])],
           depth_ap_ppm=apdepth(ap)*1e6, depth_core4_ppm=apdepth(core)*1e6, depth_ring_ppm=apdepth(ring)*1e6, depth_big_ppm=apdepth(big)*1e6,
           depth_ap_err_ppm=float(np.nanstd([np.nansum(fk[ap])/np.nansum(oot[ap]) for fk in fakes])*1e6),
           depth_core4_err_ppm=float(np.nanstd([np.nansum(fk[core])/np.nansum(oot[core]) for fk in fakes])*1e6),
           depth_ring_err_ppm=float(np.nanstd([np.nansum(fk[ring])/np.nansum(oot[ring]) for fk in fakes])*1e6),
           max_pixel_snr=float(np.nanmax(snr)), snr_summed_ap=float(np.nansum(diff[ap])/np.sqrt(np.nansum(sig[ap]**2))))

real = np.hypot(dcx-ocx, dcy-ocy)*3.98
fo = np.array([np.hypot(*(np.array(cen(fk, m3)) - [ocx, ocy]))*3.98 for fk in fakes])
fs = np.array([np.nansum(fk[ap])/np.nansum(oot[ap])*1e6 for fk in fakes])
res['null'] = dict(real_offset_arcsec=float(real), n_fake=int(len(fakes)), fake_offset_arcsec_pct=dict(zip(['p16','p50','p84','p95'], np.percentile(fo, [16,50,84,95]).tolist())),
         frac_fake_offset_ge_real=float(np.mean(fo >= real)), fake_depth_ap_std_ppm=float(np.std(fs)))
json.dump(res, open(OUT + '.json', 'w'), indent=1)
np.savez(OUT + '.npz', diff=diff, oot=oot, sig=sig, snr=snr, depth_map=depth_map, ap=ap)
fig, ax = plt.subplots(1, 4, figsize=(16, 4))
for a, im, ti in zip(ax, [oot, diff, snr, depth_map*1e6], ['out-of-transit', 'difference (OOT-IN)', 'difference S/N', 'per-pixel depth ppm']):
    v = np.nanpercentile(im, [2, 98]) if ti != 'per-pixel depth ppm' else (-3000, 3000)
    c = a.imshow(im, origin='lower', cmap='viridis', vmin=v[0], vmax=v[1]); plt.colorbar(c, ax=a, fraction=0.046)
    a.plot(cx, cy, 'r*', ms=12); a.plot(dcx, dcy, 'wx', ms=10); a.contour(ap, levels=[0.5], colors='w', linewidths=0.8); a.set_title(ti)
plt.suptitle(OUT.split('/')[-1] + f"  offset {res['offset_arcsec']:.2f}\" ({res['offset_sigma']:.1f} sigma)")
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=90)
for k, v in res.items(): print(k, v)
