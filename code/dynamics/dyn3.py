"""dyn3.py: masses from the Chen & Kipping (2017) probabilistic mass-radius relation (forecaster, github.com/chenjj2/forecaster @90f36e5 via _sources/forecaster_repo/py3 print-only port, Rpost2M on the radius posterior of fit3),
then period ratios, a from Kepler's third law with M* = 0.381 +- 0.021 Msun (TIC 8.2), and mutual Hill separations of adjacent pairs
Delta = (a_out - a_in) / R_H,  R_H = [(m_in + m_out)/(3 M*)]^(1/3) (a_in + a_out)/2  (Chambers et al. 1996 convention), circular orbits.
usage (needs h5py): dyn3.py RESULT_JSON CHAIN_NPZ -> fit3/dyn3.json"""
import sys, os, json, numpy as np
RES, CH = sys.argv[1], sys.argv[2]; D = 'k2-deep/247455141'
R = json.load(open(RES)); rng = np.random.default_rng(9)
FD = f'{D}/_sources/forecaster_repo/py3'; np.random.seed(9); cwd = os.getcwd(); os.chdir(FD); sys.path.insert(0, FD); import mr_forecast as mr; os.chdir(cwd)
MS, EMS = 0.381037, 0.0205361; n = 2000
out = {}; M = {}; A = {}
Ms = MS + EMS*rng.standard_normal(n)
for p in 'bcd':
    r = R[p]['Rp_earth']; rs = np.clip(r[0] + 0.5*(r[1]+r[2])*rng.standard_normal(n), 0.3, None)
    m = np.array(mr.Rpost2M(rs, unit='Earth', grid_size=1000)); M[p] = m
    P = R[p]['P_d'][0]; A[p] = (Ms*(P/365.25)**2)**(1/3)    # AU
    q = np.percentile(m, [16, 50, 84]); out[p] = dict(Rp_earth=r, M_earth_CK17=[float(q[1]), float(q[1]-q[0]), float(q[2]-q[1])], a_AU=float(np.median(A[p])), P_d=P)
MSE = 332946.0
for i, o in [('b', 'c'), ('c', 'd'), ('b', 'd')]:
    rh = ((M[i] + M[o])/MSE/(3*Ms))**(1/3)*(A[i] + A[o])/2; dl = (A[o] - A[i])/rh
    q = np.percentile(dl, [2.5, 16, 50, 84, 97.5])
    out[f'{i}{o}'] = dict(period_ratio=R[o]['P_d'][0]/R[i]['P_d'][0], Delta_mutual_Hill=dict(zip(['p2.5', 'p16', 'p50', 'p84', 'p97.5'], q.tolist())),
                         frac_Delta_gt_10=float(np.mean(dl > 10)), frac_Delta_gt_2sqrt3=float(np.mean(dl > 2*np.sqrt(3))))
# 1-sigma-heavy case: masses at the CK17 84th percentile x 2 as a stress test
for i, o in [('b', 'c'), ('c', 'd')]:
    mm = 2*np.percentile(M[i], 84) + 2*np.percentile(M[o], 84); rh = (mm/MSE/(3*MS))**(1/3)*(np.median(A[i]) + np.median(A[o]))/2
    out[f'{i}{o}']['Delta_if_masses_2x_p84'] = float((np.median(A[o]) - np.median(A[i]))/rh)
json.dump(out, open(f'{D}/fit3/dyn3.json', 'w'), indent=1); print(json.dumps(out, indent=1))
