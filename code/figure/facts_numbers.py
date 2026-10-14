# Recompute the comparison numbers quoted in FACTS.md from the fit outputs. usage: facts_numbers.py K2DEEP_TARGET_DIR
import sys, json, numpy as np
D = sys.argv[1]
R2 = json.load(open(f'{D}/fit/result_everest_rho.json')); R3 = json.load(open(f'{D}/fit3/result_everest_m3_rho.json'))
for p in 'bc':
    print(f'== {p}: two-planet -> three-planet')
    for k in ['P_d', 'T0_BKJD', 'RpRs', 'b', 'depth_ppm', 'T14_h', 'Rp_earth', 'a_Rs', 'inc_deg', 'Teq_K', 'S_earth']:
        a, b = R2[p][k], R3[p][k]; s = 0.5*(b[1] + b[2]); d = b[0] - a[0]
        print(f'  {k}: {a[0]:.6f} -> {b[0]:.6f} (-{b[1]:.6f}/+{b[2]:.6f}); diff {d:+.6f} = {d/s:+.3f} sigma3, {100*d/a[0]:+.3f}%')
for k in ['jitter_ppm', 'q1', 'q2']: print(k, R2[k], R3[k])
print('rho2', R2['b']['rho_sun'], 'rho3', R3['b']['rho_sun'])
for nm, R in [('2pl', R2), ('3pl', R3)]:
    t = R['tau']; post = R['nstep'] - R['nstep']//3; print(nm, 'nwalk', R['nwalk'], 'nstep', R['nstep'], 'npts', R['npts'], 'acc', round(R['acc'], 3), 'tau max', round(max(t), 1), 'post/tau', round(post/max(t), 1), 'maxlnp', round(R['maxlnp'], 1))
S = json.load(open(f'{D}/fpp_joint/summary.json'))
for p in 'bc':
    dp = S[f'{p}_molusc']['depth_ppm']; f2 = R2[p]['depth_ppm'][0]; f3 = R3[p]['depth_ppm'][0]
    print(f'FPP depth {p}: passed {dp:.4f}; 2pl fit {f2:.4f} (diff {100*(dp-f2)/f2:+.4f}%); adopted 3pl fit {f3:.4f} (diff {100*(dp-f3)/f3:+.3f}%, {(dp-f3)/(0.5*(R3[p]["depth_ppm"][1]+R3[p]["depth_ppm"][2])):+.3f} sigma)')
Sd = json.load(open(f'{D}/fpp_joint_d/summary.json')); print('FPP depth d passed', Sd['d_molusc']['depth_ppm'], 'fit', R3['d']['depth_ppm'][0])
for f in ['k2sff_rho', 'fresh40_rho', 'everest_free']:
    R = json.load(open(f'{D}/fit/result_{f}.json')); print(f, {p: {k: [round(x, 5) for x in R[p][k]] for k in ['P_d', 'depth_ppm', 'Rp_earth', 'b', 'rho_sun']} for p in 'bc'}, 'jit', [round(x, 1) for x in R['jitter_ppm']], R['priors'].get('rho'))
RF = json.load(open(f'{D}/fit3/result_everest_m3_dfree.json')); print('dfree b,c', {p: {k: [round(x, 5) for x in RF[p][k]] for k in ['P_d', 'depth_ppm', 'Rp_earth']} for p in 'bcd'})
