"""Copy of ../../scripts/inj.py for two signals (EPIC 247455141): WHICH (b or c) is the signal under test; ephemerides from fit/result_everest_rho.json.
mode=real masks only the OTHER signal; null and inj mask both signals (|phase| < 0.75 T14) before injecting.
usage: inj2.py WHICH N SEED OUT MODE
Env MASKD=1: also mask the third TLS peak of fit/tls3.json (EVEREST P 9.35604 d) in every mode.
Injection-recovery using the hunt's own detrend+TLS (k2tls_hunt.py clean/power settings).
usage: inj.py EPIC SRC RPRS ARS_AT_P0 P0 N SEED OUT [mode]
mode=inj (default): inject batman transits at random P in [P0*0.7, P0*1.3] and random T0 into the raw signal-masked light curve.
mode=null: no injection, real signal masked -> TLS false-alarm level.
mode=real: real data unmasked -> SDE of the candidate itself."""
import sys, json, numpy as np, batman, warnings; warnings.filterwarnings('ignore')
from transitleastsquares import transitleastsquares
sys.path.insert(0, 'k2-hunt/scripts'); from k2tls_hunt import clean
which, N, seed, OUT, mode = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]; epic, src = '247455141', 'everest'
RJ = json.load(open('k2-deep/247455141/fit/result_everest_rho.json')); oth = 'c' if which == 'b' else 'b'
k, aR0, P0 = RJ[which]['RpRs'][0], RJ[which]['a_Rs'][0], RJ[which]['P_d'][0]
T0r, Pr, T14r = float(sys.argv[10]) if len(sys.argv) > 10 else None, P0, None
D = f'k2-deep/{epic}'
z = np.load(f'{D}/lc/{src}.npz'); t0, f0 = z['t_raw'], z['f_raw']
Tc, Pc, T14 = RJ[which]['T0_BKJD'][0], RJ[which]['P_d'][0], RJ[which]['T14_h'][0]/24
To, Po, T14o = RJ[oth]['T0_BKJD'][0], RJ[oth]['P_d'][0], RJ[oth]['T14_h'][0]/24
intro = np.abs(((t0 - To + 0.5*Po) % Po) - 0.5*Po) < 0.75*T14o
import os
if os.environ.get('MASKD'):
    T3 = json.load(open('k2-deep/247455141/fit/tls3.json'))['everest']; intro = intro | (np.abs(((t0 - T3['T0'] + 0.5*T3['P']) % T3['P']) - 0.5*T3['P']) < 0.75*T3['dur_h']/24)
intr = (np.abs(((t0 - Tc + 0.5*Pc) % Pc) - 0.5*Pc) < 0.75*T14) | intro
rng = np.random.default_rng(seed); out = []
targ = [l for l in open('k2-hunt/targets_all.csv') if l.startswith(epic)][0].split(',')
Rs, Ms = float(targ[7]), float(targ[8])
for i in range(N):
    if mode == 'real': t, f = t0[~intro], f0[~intro].copy(); Pi = Pc; Ti = Tc
    else:
        t, f = t0[~intr], f0[~intr].copy()
        if mode == 'inj':
            Pi = rng.uniform(0.7*P0, 1.3*P0); Ti = t.min() + rng.uniform(0, Pi)
            pm = batman.TransitParams(); pm.t0=Ti; pm.per=Pi; pm.rp=k; pm.a=aR0*(Pi/P0)**(2/3); pm.inc=np.degrees(np.arccos(0.4/pm.a)); pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
            f = f*batman.TransitModel(pm, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(pm)
        else: Pi, Ti = np.nan, np.nan
    tc, fc = clean(t, f)
    r = transitleastsquares(tc, fc).power(period_min=0.5, period_max=min(30.0, (tc.max()-tc.min())/2), R_star=Rs, M_star=Ms, R_star_min=0.08, R_star_max=1.0,
        M_star_min=0.08, M_star_max=1.0, oversampling_factor=2, duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
    rec = bool(np.isfinite(Pi) and (abs(r.period - Pi)/Pi < 0.01 or abs(r.period - 2*Pi)/Pi < 0.02 or abs(2*r.period - Pi)/Pi < 0.02))
    out.append(dict(i=i, mode=mode, P_inj=float(Pi), T0_inj=float(Ti), P_rec=float(r.period), SDE=float(r.SDE), snr=float(r.snr), depth=float(1-r.depth), recovered_period=rec))
    print(out[-1], flush=True)
    if mode == 'real': break
json.dump(out, open(OUT, 'w'), indent=1)
