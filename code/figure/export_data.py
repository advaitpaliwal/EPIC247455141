# Export the light curves and posterior summaries of EPIC 247455141 to CSV for the public repository. usage: export_data.py K2DEEP_TARGET_DIR OUT_DATA_DIR
import sys, json, numpy as np, pandas as pd, batman
D, O = sys.argv[1], sys.argv[2]
R = json.load(open(f'{D}/fit3/result_everest_m3_rho.json')); z = np.load(f'{D}/lc/everest_m3.npz')
c = np.load(f'{D}/fit3/chain_everest_m3_rho.npz'); ch, lp, Tref = c['ch'], c['lp'], c['Tref']; med = np.median(ch[lp > lp.max() - 30], axis=0)
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
def model(t, i):
    t0, P, lrp, b = med[4*i:4*i+4]; rho, q1, q2 = med[12], med[13], med[14]; a = ((G*rho*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
    pm = batman.TransitParams(); pm.t0 = t0; pm.per = P; pm.rp = np.exp(lrp); pm.a = a; pm.inc = np.degrees(np.arccos(b/a)); pm.ecc = 0; pm.w = 90; pm.limb_dark = 'quadratic'; pm.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
    return batman.TransitModel(pm, t - Tref[i], supersample_factor=15, exp_time=1765.5/86400).light_curve(pm)
t, f = z['t'], z['f']; df = pd.DataFrame(dict(time_bkjd=t, flux_detrended=f))
for i, p in enumerate('bcd'):
    P, T0 = R[p]['P_d'][0], R[p]['T0_BKJD'][0]; df[f'phase_{p}_d'] = ((t - T0 + 0.5*P) % P) - 0.5*P; df[f'model_{p}'] = model(t, i)
df['model_total'] = med[16]*df.model_b*df.model_c*df.model_d
df.to_csv(f'{O}/lc/everest_detrended_m3_c13.csv', index=False, float_format='%.8f')
pd.DataFrame(dict(time_bkjd=z['t_raw'], flux_everest_normalised=z['f_raw'])).to_csv(f'{O}/lc/everest_undetrended_c13.csv', index=False, float_format='%.8f')
rows = []
for p in 'bcd':
    for k, v in R[p].items(): rows.append(dict(planet=p, parameter=k, median=v[0], minus_1sigma=v[1], plus_1sigma=v[2]))
for k in ['q1', 'q2', 'jitter_ppm']: rows.append(dict(planet='shared', parameter=k, median=R[k][0], minus_1sigma=R[k][1], plus_1sigma=R[k][2]))
pd.DataFrame(rows).to_csv(f'{O}/fit/posterior_summary_three_planet_everest.csv', index=False)
R2 = json.load(open(f'{D}/fit/result_everest_rho.json')); rows = []
for p in 'bc':
    for k, v in R2[p].items(): rows.append(dict(planet=p, parameter=k, median=v[0], minus_1sigma=v[1], plus_1sigma=v[2]))
pd.DataFrame(rows).to_csv(f'{O}/fit/posterior_summary_two_planet_everest.csv', index=False)
print(len(df), 'rows')
