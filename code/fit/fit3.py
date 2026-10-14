"""Joint three-planet batman+emcee fit; scripts/fit2.py extended to the 9.356-d signal d. usage: fit3.py SRC MODE NSTEP  (SRC e.g. everest_m3)
MODE=rho:   one stellar density shared by b, c, d with Gaussian prior (TIC 8.2: 6.385 +- 0.245 rho_sun).
MODE=dfree: b and c share the TIC-prior density; d gets its own density, uniform 0.1-60 rho_sun, no prior -> density implied by d's transit shape alone.
Shared: quadratic LD (Kipping q1,q2 uniform), white-noise jitter sigma (log-uniform), flux offset f0. Per planet: T0 (offset from reference epoch, +-0.1 d), P (+-0.05 d of start), ln(Rp/R*), b (uniform 0..1+k). Circular orbits.
Likelihood: Gaussian white noise, all cadences within max(3 T14, 0.25 d) of any of the three signals' transits. Model: batman, supersample 15 over the 29.42-min long cadence; product of the three planets' models."""
import sys, json, numpy as np, batman, emcee, time
from scipy.optimize import minimize
src, MODE, NSTEP = sys.argv[1], sys.argv[2], int(sys.argv[3])
D = 'k2-deep/247455141'
RHO, ERHO, RS, ERS, TEFF, ETEFF = 6.38504, 0.244728, 0.390782, 0.0120151, 3550., 157.   # TIC 8.2, catalogue/mast_tic.json (as fit2.py)
RJ = json.load(open(f'{D}/fit/result_everest_rho.json')); T3 = json.load(open(f'{D}/fit/tls3.json'))['everest']
EPH = [(RJ[p]['P_d'][0], RJ[p]['T0_BKJD'][0], RJ[p]['T14_h'][0]/24) for p in 'bc'] + [(T3['P'], T3['T0'], T3['dur_h']/24)]
NP = 3
z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
w = np.zeros(len(t), bool)
for P0, T00, d in EPH: w |= np.abs(((t - T00 + 0.5*P0) % P0) - 0.5*P0) < max(3*d, 0.25)
t, f = t[w], f[w]
Tref = [T00 + np.round((np.median(t) - T00)/P0)*P0 for P0, T00, d in EPH]
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
def a_rs(rho, P): return ((G*rho*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
pms, ms = [], []
for i in range(NP):
    pm = batman.TransitParams(); pm.t0=0; pm.per=EPH[i][0]; pm.rp=0.03; pm.a=20; pm.inc=89.5; pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
    pms.append(pm); ms.append(batman.TransitModel(pm, t - Tref[i], supersample_factor=15, exp_time=1765.5/86400))
# theta: [t0,P,lrp,b]x3, rho_bc, (rho_d if dfree), q1, q2, lsig, f0 ; unpack -> rho list of 3
def unpack(th):
    th = list(th); base = th[:12]
    if MODE == 'rho': rho = [th[12]]*3; rest = th[13:]
    else: rho = [th[12], th[12], th[13]]; rest = th[14:]
    return base, rho, rest
def model(th):
    base, rho, (q1, q2, lsig, f0) = unpack(th); out = f0*np.ones(len(t))
    for i in range(NP):
        t0, P, lrp, b = base[4*i:4*i+4]; a = a_rs(rho[i], P); pm = pms[i]
        pm.t0 = t0; pm.per = P; pm.rp = np.exp(lrp); pm.a = a; pm.inc = np.degrees(np.arccos(b/a)); pm.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
        out *= ms[i].light_curve(pm)
    return out
def lnp(th):
    base, rho, (q1, q2, lsig, f0) = unpack(th)
    if not (0 < q1 < 1 and 0 < q2 < 1 and -12 < lsig < -5): return -np.inf
    for i in range(NP):
        t0, P, lrp, b = base[4*i:4*i+4]; rp = np.exp(lrp)
        if not (-0.1 < t0 < 0.1 and abs(P - EPH[i][0]) < 0.05 and 0.001 < rp < 0.3 and 0 <= b < 1 + rp and 0.1 < rho[i] < 60): return -np.inf
        if b/a_rs(rho[i], P) >= 1: return -np.inf
    lp = -0.5*((rho[0]-RHO)/ERHO)**2
    s2 = np.exp(2*lsig); r = f - model(th)
    return lp - 0.5*np.sum(r**2/s2 + np.log(2*np.pi*s2))
x0 = []
for i, k0 in enumerate([0.0305, 0.0273, 0.032]): x0 += [0., EPH[i][0], np.log(k0), 0.3]
x0 += [RHO] + ([] if MODE == 'rho' else [RHO]) + [0.4, 0.3, np.log(np.std(f)), 1.0]
x0 = np.array(x0); best = minimize(lambda x: -lnp(x), x0, method='Nelder-Mead', options=dict(maxiter=60000, maxfev=60000, xatol=1e-8, fatol=1e-6))
x1 = best.x if np.isfinite(best.fun) else x0
ndim = len(x1); nw = 4*ndim
scl = np.array([1e-3, 1e-5, 0.02, 0.05]*3 + [0.3] + ([] if MODE == 'rho' else [0.3]) + [0.05, 0.05, 0.02, 1e-5])
rng0 = np.random.default_rng(11); p0 = []
while len(p0) < nw:
    c = x1 + scl*rng0.standard_normal(ndim)
    for j in (3, 7, 11): c[j] = abs(c[j])
    if np.isfinite(lnp(c)): p0.append(c)
s = emcee.EnsembleSampler(nw, ndim, lnp); t1 = time.time(); s.run_mcmc(np.array(p0), NSTEP, progress=False)
ch = s.get_chain(discard=NSTEP//3, thin=5, flat=True); lpch = s.get_log_prob(discard=NSTEP//3, thin=5, flat=True)
try: tau = s.get_autocorr_time(discard=NSTEP//3, quiet=True).tolist()
except Exception as e: tau = str(e)
U = [unpack(c) for c in ch]; B = np.array([u[0] for u in U]); RH = np.array([u[1] for u in U]); RE = np.array([u[2] for u in U])
def q(x): p = np.percentile(x, [16, 50, 84]); return [float(p[1]), float(p[1]-p[0]), float(p[2]-p[1])]
res = dict(source=src, mode=MODE, nstep=NSTEP, nwalk=nw, npts=int(len(t)), acc=float(np.mean(s.acceptance_fraction)), tau=tau, maxlnp=float(lpch.max()),
           priors=dict(rho_bc=[RHO, ERHO], rho_d='shared with b,c' if MODE == 'rho' else 'uniform 0.1-60, no prior', Rstar=[RS, ERS], Teff=[TEFF, ETEFF], ld='uniform q1,q2 (Kipping 2013), shared', ecc='fixed 0'),
           q1=q(RE[:, 0]), q2=q(RE[:, 1]), jitter_ppm=q(np.exp(RE[:, 2])*1e6), runtime_s=time.time()-t1)
rng = np.random.default_rng(3); Rs = RS + ERS*rng.standard_normal(len(U)); Te = TEFF + ETEFF*rng.standard_normal(len(U))
for i, nm in enumerate('bcd'):
    t0, P, lrp, b = B[:, 4*i:4*i+4].T; rho = RH[:, i]; rp = np.exp(lrp); a = a_rs(rho, P); inc = np.arccos(b/a)
    T14 = P/np.pi*np.arcsin(np.sqrt(np.clip((1+rp)**2 - b**2, 0, None))/(a*np.sin(inc)))*24
    dep = []
    for j in rng.choice(len(U), 2000, replace=False):
        p2 = batman.TransitParams(); p2.t0=0; p2.per=P[j]; p2.rp=rp[j]; p2.a=a[j]; p2.inc=np.degrees(inc[j]); p2.ecc=0; p2.w=90; p2.limb_dark='quadratic'
        p2.u=[2*np.sqrt(RE[j,0])*RE[j,1], np.sqrt(RE[j,0])*(1-2*RE[j,1])]; dep.append(1 - batman.TransitModel(p2, np.array([0.0])).light_curve(p2)[0])
    Tc = Tref[i] + t0
    res[nm] = dict(T0_BKJD=q(Tc), T0_BJD_TDB=q(Tc + 2454833.0), P_d=q(P), RpRs=q(rp), b=q(b), rho_sun=q(rho), a_Rs=q(a), inc_deg=q(np.degrees(inc)),
                   depth_ppm=q(np.array(dep)*1e6), T14_h=q(T14), Rp_earth=q(rp*Rs*109.076), Teq_K=q(Te*np.sqrt(1/(2*a))), S_earth=q(Rs**2*(Te/5772.)**4*(a*Rs*0.00465047)**-2))
if MODE == 'dfree':
    r = RH[:, 2]; res['rho_d_alone'] = q(r); res['rho_d_pct'] = dict(zip(['p2.5', 'p97.5'], np.percentile(r, [2.5, 97.5]).tolist()))
    res['rho_d_over_TIC'] = q(r/RHO); res['P(rho_d within TIC 1sigma-convolved)'] = float(np.mean(np.abs(r - RHO) < ERHO))
json.dump(res, open(f'{D}/fit3/result_{src}_{MODE}.json', 'w'), indent=1)
np.savez(f'{D}/fit3/chain_{src}_{MODE}.npz', ch=ch, lp=lpch, Tref=Tref)
bi = ch[np.argmax(lpch)]; np.savez(f'{D}/fit3/bestmodel_{src}_{MODE}.npz', t=t, f=f, mod=model(bi), Tref=Tref, best=bi)
print(json.dumps(res, indent=1))
