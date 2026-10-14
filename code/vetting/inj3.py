"""inj3.py: scripts/inj2.py adapted to the third signal d (P 9.35604 d, fit/tls3.json EVEREST). Same detrend (k2tls_hunt.clean) and TLS settings as inj2.py.
mode=real: b and c in-transit cadences masked (|phase|<0.75 T14), d left in -> TLS SDE of d itself.
mode=null: b, c and d all masked, no injection -> TLS top-peak SDE with nothing in the data.
mode=inj:  b, c, d masked; inject a batman transit at random P in [0.7, 1.3] x 9.356 d, random T0, Rp/R* = sqrt(1100 ppm) = 0.0332, b = 0.4,
           a/R* from the TIC 8.2 density (6.385 rho_sun) at each P, quadratic LD u = [0.3, 0.3], 15x supersampled over 29.42 min. Recovered = TLS top peak within 1% of P (or 2% of 2P, P/2).
usage: inj3.py N SEED OUT MODE"""
import sys, json, numpy as np, batman, warnings; warnings.filterwarnings('ignore')
from transitleastsquares import transitleastsquares
sys.path.insert(0, 'k2-hunt/scripts'); from k2tls_hunt import clean
N, seed, OUT, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
D = 'k2-deep/247455141'
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); T3 = json.load(open(f'{D}/fit/tls3.json'))['everest']
z = np.load(f'{D}/lc/everest.npz'); t0, f0 = z['t_raw'], z['f_raw']
def m(P, T0, d): return np.abs(((t0 - T0 + 0.5*P) % P) - 0.5*P) < 0.75*d
mbc = m(RJ['b']['P_d'][0], RJ['b']['T0_BKJD'][0], RJ['b']['T14_h'][0]/24) | m(RJ['c']['P_d'][0], RJ['c']['T0_BKJD'][0], RJ['c']['T14_h'][0]/24)
P0, Td = T3['P'], T3['T0']; md = m(P0, Td, T3['dur_h']/24)
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
def a_rs(P): return ((G*6.38504*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
rng = np.random.default_rng(seed); out = []
for i in range(N):
    if mode == 'real': t, f = t0[~mbc], f0[~mbc].copy(); Pi, Ti = P0, Td
    else:
        t, f = t0[~(mbc | md)], f0[~(mbc | md)].copy(); Pi, Ti = np.nan, np.nan
        if mode == 'inj':
            Pi = rng.uniform(0.7*P0, 1.3*P0); Ti = t.min() + rng.uniform(0, Pi)
            pm = batman.TransitParams(); pm.t0=Ti; pm.per=Pi; pm.rp=np.sqrt(1100e-6); pm.a=a_rs(Pi); pm.inc=np.degrees(np.arccos(0.4/pm.a)); pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
            f = f*batman.TransitModel(pm, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(pm)
    tc, fc = clean(t, f)
    r = transitleastsquares(tc, fc).power(period_min=0.5, period_max=30.0, R_star=0.39, M_star=0.382, R_star_min=0.08, R_star_max=1.0,
        M_star_min=0.08, M_star_max=1.0, oversampling_factor=2, duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
    rec = bool(np.isfinite(Pi) and (abs(r.period - Pi)/Pi < 0.01 or abs(r.period - 2*Pi)/Pi < 0.02 or abs(2*r.period - Pi)/Pi < 0.02))
    out.append(dict(i=i, mode=mode, P_inj=float(Pi), T0_inj=float(Ti), P_rec=float(r.period), SDE=float(r.SDE), snr=float(r.snr), depth_ppm=float((1-r.depth)*1e6), recovered_period=rec, n_transits_inj=int(r.distinct_transit_count) if rec else None))
    print(out[-1], flush=True)
    if mode != 'inj': break
json.dump(out, open(OUT, 'w'), indent=1)
