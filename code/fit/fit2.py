"""Joint two-planet batman+emcee fit, adapted from ../../scripts/fit.py.
usage: fit2.py SRC MODE NSTEP   MODE=rho: one shared stellar density with Gaussian prior (TIC 8.2 rho); MODE=free: each planet its own density, uniform 0.1-60 rho_sun, no prior (duration-consistency test).
Shared: quadratic LD (Kipping q1,q2 uniform), white-noise jitter sigma (log-uniform), flux offset f0. Per planet: T0, P, ln(Rp/R*), b. Circular orbits.
Likelihood: ln L = -0.5 * sum_i [ (f_i - m_i)^2 / sigma^2 + ln(2 pi sigma^2) ], over all cadences within max(3 T14, 0.25 d) of either signal's transit.
Model: batman, supersample 15 over the 29.42-min long cadence; product of the two planets' transit models."""
import sys, json, numpy as np, batman, emcee, time
from scipy.optimize import minimize
src, MODE, NSTEP = sys.argv[1], sys.argv[2], int(sys.argv[3])
D = 'k2-deep/247455141'
RHO, ERHO, RS, ERS, TEFF, ETEFF = 6.38504, 0.244728, 0.390782, 0.0120151, 3550., 157.   # TIC 8.2 (TIC 118598568), catalogue/mast_tic.json
EPH = [(3.075223, 2989.9491, 1.30/24), (5.455991, 2990.0929, 1.49/24)]
z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
w = np.zeros(len(t), bool)
for P0, T00, d in EPH: w |= np.abs(((t - T00 + 0.5*P0) % P0) - 0.5*P0) < max(3*d, 0.25)
t, f = t[w], f[w]
Tref = []; 
for P0, T00, d in EPH: n = np.round((np.median(t) - T00)/P0); Tref.append(T00 + n*P0)
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
def a_rs(rho, P): return ((G*rho*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
pms, ms = [], []
for i in range(2):
    pm = batman.TransitParams(); pm.t0=0; pm.per=EPH[i][0]; pm.rp=0.03; pm.a=20; pm.inc=89.5; pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
    pms.append(pm); ms.append(batman.TransitModel(pm, t - Tref[i], supersample_factor=15, exp_time=1765.5/86400))
# theta: [t0b,Pb,lrpb,bb, t0c,Pc,lrpc,bc, rho_b, rho_c, q1, q2, lsig, f0]  (rho_c tied to rho_b in MODE=rho)
def unpack(th):
    th = list(th)
    if MODE == 'rho': th = th[:9] + [th[8]] + th[9:]
    return th
def model(th):
    th = unpack(th); q1, q2, lsig, f0 = th[10:]; out = f0*np.ones(len(t))
    for i in range(2):
        t0, P, lrp, b = th[4*i:4*i+4]; rho = th[8+i]; a = a_rs(rho, P); pm = pms[i]
        pm.t0 = t0; pm.per = P; pm.rp = np.exp(lrp); pm.a = a; pm.inc = np.degrees(np.arccos(b/a)); pm.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
        out *= ms[i].light_curve(pm)
    return out
def lnp(th):
    u = unpack(th); q1, q2, lsig, f0 = u[10:]
    if not (0 < q1 < 1 and 0 < q2 < 1 and -12 < lsig < -5): return -np.inf
    for i in range(2):
        t0, P, lrp, b = u[4*i:4*i+4]; rho = u[8+i]; rp = np.exp(lrp)
        if not (-0.1 < t0 < 0.1 and abs(P - EPH[i][0]) < 0.05 and 0.001 < rp < 0.3 and 0 <= b < 1 + rp and 0.1 < rho < 60): return -np.inf
        if b/a_rs(rho, P) >= 1: return -np.inf
    lp = -0.5*((u[8]-RHO)/ERHO)**2 if MODE == 'rho' else 0.0
    s2 = np.exp(2*lsig); r = f - model(th)
    return lp - 0.5*np.sum(r**2/s2 + np.log(2*np.pi*s2))
x0 = [0., EPH[0][0], np.log(0.029), 0.3, 0., EPH[1][0], np.log(0.028), 0.3, RHO] + ([] if MODE == 'rho' else [RHO]) + [0.4, 0.3, np.log(np.std(f)), 1.0]
x0 = np.array(x0); best = minimize(lambda x: -lnp(x), x0, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-6))
x1 = best.x if np.isfinite(best.fun) else x0
ndim = len(x1); nw = 4*ndim
scl = np.array([1e-3, 1e-5, 0.02, 0.05, 1e-3, 1e-5, 0.02, 0.05, 0.3] + ([] if MODE == 'rho' else [0.3]) + [0.05, 0.05, 0.02, 1e-5])
p0 = []
while len(p0) < nw:
    c = x1 + scl*np.random.randn(ndim); c[3] = abs(c[3]); c[7] = abs(c[7])
    if np.isfinite(lnp(c)): p0.append(c)
s = emcee.EnsembleSampler(nw, ndim, lnp); t1 = time.time(); s.run_mcmc(np.array(p0), NSTEP, progress=False)
ch = s.get_chain(discard=NSTEP//3, thin=5, flat=True); lpch = s.get_log_prob(discard=NSTEP//3, thin=5, flat=True)
try: tau = s.get_autocorr_time(discard=NSTEP//3, quiet=True).tolist()
except Exception as e: tau = str(e)
U = np.array([unpack(c) for c in ch])
def q(x): p = np.percentile(x, [16, 50, 84]); return [float(p[1]), float(p[1]-p[0]), float(p[2]-p[1])]
res = dict(source=src, mode=MODE, nstep=NSTEP, nwalk=nw, npts=int(len(t)), acc=float(np.mean(s.acceptance_fraction)), tau=tau, maxlnp=float(lpch.max()),
           priors=dict(rho=[RHO, ERHO] if MODE == 'rho' else 'uniform 0.1-60 per planet', Rstar=[RS, ERS], Teff=[TEFF, ETEFF], ld='uniform q1,q2 (Kipping 2013), shared', ecc='fixed 0'),
           q1=q(U[:, 10]), q2=q(U[:, 11]), jitter_ppm=q(np.exp(U[:, 12])*1e6), runtime_s=time.time()-t1)
rng = np.random.default_rng(3); Rs = RS + ERS*rng.standard_normal(len(U)); Te = TEFF + ETEFF*rng.standard_normal(len(U))
for i, nm in enumerate('bc'):
    t0, P, lrp, b = U[:, 4*i:4*i+4].T; rho = U[:, 8+i]; rp = np.exp(lrp); a = a_rs(rho, P); inc = np.arccos(b/a)
    T14 = P/np.pi*np.arcsin(np.sqrt(np.clip((1+rp)**2 - b**2, 0, None))/(a*np.sin(inc)))*24
    dep = []
    for j in rng.choice(len(U), 2000, replace=False):
        p2 = batman.TransitParams(); p2.t0=0; p2.per=P[j]; p2.rp=rp[j]; p2.a=a[j]; p2.inc=np.degrees(inc[j]); p2.ecc=0; p2.w=90; p2.limb_dark='quadratic'
        p2.u=[2*np.sqrt(U[j,10])*U[j,11], np.sqrt(U[j,10])*(1-2*U[j,11])]; dep.append(1 - batman.TransitModel(p2, np.array([0.0])).light_curve(p2)[0])
    Tc = Tref[i] + t0
    res[nm] = dict(T0_BKJD=q(Tc), T0_BJD_TDB=q(Tc + 2454833.0), P_d=q(P), RpRs=q(rp), b=q(b), rho_sun=q(rho), a_Rs=q(a), inc_deg=q(np.degrees(inc)),
                   depth_ppm=q(np.array(dep)*1e6), T14_h=q(T14), Rp_earth=q(rp*Rs*109.076), Teq_K=q(Te*np.sqrt(1/(2*a))), S_earth=q(Rs**2*(Te/5772.)**4*(a*Rs*0.00465047)**-2))
if MODE == 'free': r = U[:, 8]/U[:, 9]; res['rho_ratio_b_over_c'] = q(r); res['frac_rho_c_gt_rho_b'] = float(np.mean(r < 1))
json.dump(res, open(f'{D}/fit/result_{src}_{MODE}.json', 'w'), indent=1)
np.savez(f'{D}/fit/chain_{src}_{MODE}.npz', ch=ch, lp=lpch, Tref=Tref)
bi = ch[np.argmax(lpch)]; np.savez(f'{D}/fit/bestmodel_{src}_{MODE}.npz', t=t, f=f, mod=model(bi), Tref=Tref, best=bi)
print(json.dumps(res, indent=1))
