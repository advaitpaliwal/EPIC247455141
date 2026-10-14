# Joint TRICERATOPS(K2)+MOLUSC FPP for EPIC 247455141 (C13), one signal at a time. Copied from ../../247192465/fpp_joint/scripts/run_joint.py; ephemeris, depth, T14 from the joint two-planet fit ../fit/result_everest_rho.json; the other signal's in-transit cadences (|phase|<0.75 T14) are removed before folding. No neighbours removed (all stars TRICERATOPS finds within 10 pixels kept).
# (same read-only shim, 15-min binning of the EVEREST LC within +-0.3 d, threshold-3 aperture, filt Kepler, N=1e6, seed 2000+k, molusc_file, parallel=False).
# usage: run_joint.py VARIANT K   VARIANT in {b_molusc, c_molusc (MOLUSC prior, 10 runs each), b_none, c_none (no MOLUSC, 1 reference run each)}
import sys, os, json, pickle, time, glob, types, numpy as np, warnings; warnings.filterwarnings('ignore')
D = 'k2-deep/247455141/fpp_joint'; os.chdir(D)
import triceratops.triceratops as tr, triceratops.marginal_likelihoods as ml, lightkurve as lk
class _F:
    def __init__(s, d): s.d = d
    def __getitem__(s, k):
        if isinstance(k, str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
        return _F(s.d[np.asarray(k)])
_rc = ml.read_csv
ml.read_csv = lambda f, *a, **kw: _F(_rc(f, *a, **kw)) if 'molusc' in str(f) else _rc(f, *a, **kw)
var, k = sys.argv[1], int(sys.argv[2]); exptime = 29.4/1440
epic = 247455141
RJ = json.load(open('../fit/result_everest_rho.json')); w_ = var[0]; o_ = 'c' if w_ == 'b' else 'b'; R = RJ[w_]
P, T0, dur, depth = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24, R['depth_ppm'][0]*1e-6
Po, T0o, duro = RJ[o_]['P_d'][0], RJ[o_]['T0_BKJD'][0], RJ[o_]['T14_h'][0]/24
pk = f'out/target_{epic}.pkl'
if not os.path.exists(pk):
    if os.environ.get('GAIAPATCH'): sys.path.insert(0, D + '/scripts'); import gaia_patch
    tg = tr.target(ID=epic, sectors=np.array([13]), mission='K2', lightkurve_cache_dir=D + '/lkcache', search_radius=10)
    pickle.dump(tg, open(pk, 'wb')); tg.stars.to_csv(f'out/stars_{epic}.csv', index=False)
    if var.endswith('init'): print(tg.stars.to_string()); sys.exit()
tg = pickle.load(open(pk, 'rb'))
tpp = glob.glob(f'k2-deep/{epic}/lc/mastDownload/K2/ktwo{epic}-c13_lc/*lpd-targ.fits.gz')[0]
tp = lk.read(tpp); ap = tp.create_threshold_mask(threshold=3, reference_pixel='center')
rows, cols = np.nonzero(ap); apx = np.column_stack([cols + tp.column, rows + tp.row])
z = np.load('../lc/everest.npz'); t, f = z['t'], z['f']
ko = np.abs(t - T0o - np.round((t - T0o)/Po)*Po) > 0.75*duro; t, f = t[ko], f[ko]
ph = t - T0 - np.round((t - T0)/P)*P; w = np.abs(ph) < 0.3
edges = np.arange(-0.3, 0.3 + 1e-9, 15/1440); idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)]); fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
err = float(np.std(fb[np.abs(tb) > dur]))
tg.calc_depths(tdepth=depth, all_ap_pixels=[apx])
mf = None if var.endswith('none') else D + '/molusc_out/molusc_for_trc_N1e6.csv'
np.random.seed(2000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=None, filt='Kepler', N=1000000, parallel=False,
              exptime=exptime, nsamples=20, verbose=0, molusc_file=mf)
tag = f'{var}_{k}'
tg.probs.to_csv(f'out/probs_{tag}.csv', index=False)
json.dump(dict(epic=epic, variant=var, run=k, seed=2000+k, N=1000000, exptime_min=exptime*1440, P=P, T0=T0, depth_ppm=depth*1e6, err_ppm=err*1e6, nbins=len(tb),
               molusc=mf, filt='Kepler', FPP=float(tg.FPP), NFPP=float(tg.NFPP), n_stars=len(tg.stars), sec=time.time() - t0, aperture=apx.tolist()), open(f'out/res_{tag}.json', 'w'), indent=1)
print(tag, tg.FPP, tg.NFPP, len(tg.stars), round(time.time() - t0), flush=True)
