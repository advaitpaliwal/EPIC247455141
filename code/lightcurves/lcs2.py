"""Copy of ../../scripts/lcs.py adapted for two signals: every transit mask is the union over both ephemerides (EPH list of (P, T0, T14_d))."""
import numpy as np, glob, os, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, 'k2-deep/scripts'); import lcs
from wotan import flatten
EPH = [(3.075223, 2989.9491, 1.30/24), (5.455991, 2990.0929, 1.49/24)]
everest, k2sff, k2sc, phase, boxdepth, tpf = lcs.everest, lcs.k2sff, lcs.k2sc, lcs.phase, lcs.boxdepth, lcs.tpf
def intr(t, eph=EPH, fac=0.75):
    m = np.zeros(len(t), bool)
    for P, T0, d in eph: m |= np.abs(phase(t, P, T0)) < fac*d
    return m
def flat(t, f, eph=EPH, window=0.75):
    """wotan biweight (0.75 d) on the light curve with all transits of both signals masked; trend interpolated across masks; 4-sigma upper / 8-sigma lower clip of out-of-transit points (MAD)."""
    it = intr(t, eph)
    trend = flatten(t, np.where(it, np.nan, f), method='biweight', window_length=window, return_trend=True)[1]
    ok = np.isfinite(trend); trend = np.interp(t, t[ok], trend[ok])
    fl = f/trend; m = np.isfinite(fl); t, fl = t[m], fl[m]
    oot = ~intr(t, eph)
    med = np.median(fl[oot]); sd = 1.4826*np.median(np.abs(fl[oot]-med))
    keep = ((fl < med + 4*sd) & (fl > med - 8*sd)) | ~oot
    keep &= (fl > med - 20*sd)
    return t[keep], fl[keep]
def fresh(epic, eph=EPH, seg='c13', thresh=3, windows=20):
    import lightkurve as lk
    tp = tpf(epic, seg)
    ap = tp.create_threshold_mask(threshold=thresh, reference_pixel='center')
    lc = tp.to_lightcurve(aperture_mask=ap).remove_nans()
    lc = lc[np.isfinite(lc.centroid_col.value) & np.isfinite(lc.centroid_row.value)]
    t = lc.time.value; it = intr(t, eph)
    corr = lk.SFFCorrector(lc).correct(windows=windows, bins=15, cadence_mask=~it, restore_trend=True)
    f = corr.flux.value; m = np.isfinite(f)
    return t[m], f[m]/np.nanmedian(f[m]), ap
