"""Copy of ../../scripts/tess.py looping over both signals (b, c) of the joint fit fit/result_everest_rho.json; one TESScut download per sector.
TESScut aperture photometry + ephemeris propagation + expected-SNR / constrained search. usage: tess.py EPIC RA DEC SECTOR(s comma)"""
import sys, json, numpy as np, warnings, os; warnings.filterwarnings('ignore')
import lightkurve as lk
from wotan import flatten
epic, ra, dec, secs = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), [int(s) for s in sys.argv[4].split(',')]
D = f'k2-deep/{epic}'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); OUT = {}; _cache = {}
for which in 'bc':
  R = RJ[which]
  T0, eT0 = R['T0_BKJD'][0], max(R['T0_BKJD'][1:]); P, eP = R['P_d'][0], max(R['P_d'][1:]); T14 = R['T14_h'][0]/24; dep = R['depth_ppm'][0]*1e-6
  out = dict(T0_BKJD=T0, eT0=eT0, P=P, eP=eP, sectors={}); OUT[which] = out
  for s in secs:
      sr = lk.search_tesscut(f'{ra} {dec}', sector=s)
      tpf = sr.download(cutout_size=15, download_dir=f'{D}/tess') if s not in _cache else _cache[s]; _cache[s] = tpf
      open(f'{D}/tess/tesscut_s{s}.url', 'w').write(f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra}&dec={dec}&y=15&x=15&sector={s}\n')
      ap = tpf.create_threshold_mask(threshold=3, reference_pixel='center')
      if ap.sum() == 0 or ap.sum() > 12:
          ap = np.zeros(tpf.shape[1:], bool); c = np.array(tpf.shape[1:])//2; ap[c[0]-1:c[0]+1, c[1]-1:c[1]+1] = True
      bk = tpf.create_threshold_mask(threshold=0.001, reference_pixel=None) & False
      med = np.nanmedian(tpf.flux.value, 0); bkm = med < np.nanpercentile(med, 30)
      bkg = np.nanmedian(tpf.flux.value[:, bkm], axis=1)
      raw = np.nansum(tpf.flux.value[:, ap], axis=1) - bkg*ap.sum()
      t = tpf.time.value; q = tpf.quality == 0; m = q & np.isfinite(raw) & (raw > 0)
      t, fl = t[m], raw[m]/np.nanmedian(raw[m])
      fl = flatten(t, fl, method='biweight', window_length=0.75)
      ok = np.isfinite(fl); t, fl = t[ok], fl[ok]
      dt = np.median(np.diff(t)); nb = max(1, int(round((30/1440)/dt)))
      rms30 = np.std(fl)/np.sqrt(nb)  # white-noise approx
      tb = t + 2167.0  # BTJD -> BKJD
      n = np.round((np.median(tb) - T0)/P); Tpred = T0 + n*P; eTpred = np.hypot(eT0, n*eP)
      # transits in sector at nominal ephemeris
      epochs = np.arange(np.ceil((tb.min()-T0)/P), np.floor((tb.max()-T0)/P)+1)
      tmids = T0 + epochs*P; cov = [float(tm) for tm in tmids if np.any(np.abs(tb - tm) < T14/2)]
      nin = np.sum(np.abs(((tb - T0 + 0.5*P) % P) - 0.5*P) < 0.5*T14)
      exp_snr = dep/(np.std(fl)/np.sqrt(max(nin, 1)))
      # nominal-ephemeris box depth
      ph = ((tb - T0 + 0.5*P) % P) - 0.5*P; i = np.abs(ph) < 0.5*T14; o = (np.abs(ph) > T14) & (np.abs(ph) < 0.5)
      d0 = np.median(fl[o]) - np.mean(fl[i]) if i.sum() else np.nan; e0 = np.std(fl[o])/np.sqrt(max(i.sum(), 1))
      # constrained search: T0 shift within +-3 sigma propagated, P within +-3 sigma
      best = []
      for dP in np.linspace(-3*eP, 3*eP, 13):
          for sh in np.arange(-3*eTpred, 3*eTpred, T14/4):
              ph2 = ((tb - (Tpred + sh) + 0.5*(P+dP)) % (P+dP)) - 0.5*(P+dP); ii = np.abs(ph2) < 0.5*T14
              if ii.sum() < 3: continue
              best.append((dP, sh, (np.median(fl) - np.mean(fl[ii]))/(np.std(fl)/np.sqrt(ii.sum()))))
      best = np.array(best)
      out['sectors'][s] = dict(npts=int(len(t)), cadence_min=float(dt*1440), aperture_npix=int(ap.sum()), rms_per_cadence_ppm=float(np.std(fl)*1e6),
          rms_30min_ppm=float(rms30*1e6), n_epochs_from_K2=int(n), Tpred_mid_BKJD=float(Tpred), Tpred_mid_BTJD=float(Tpred-2167), sigma_Tpred_h=float(eTpred*24),
          transits_covered_nominal=len(cov), n_in_nominal=int(nin), expected_snr=float(exp_snr), nominal_box_depth_ppm=float(d0*1e6), nominal_box_err_ppm=float(e0*1e6),
          search_max_snr=float(best[:, 2].max()) if len(best) else None, search_trials=int(len(best)),
          tic_contamination_note='TIC contratio applies to TESS pixels (21 arcsec)')
      np.savez(f'{D}/tess/s{s}_lc_{which}.npz', t=t, f=fl)
      print(which, s, out['sectors'][s], flush=True)
json.dump(OUT, open(f'{D}/tess/tess_check.json', 'w'), indent=1)
