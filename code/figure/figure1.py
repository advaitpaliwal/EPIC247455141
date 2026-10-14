"""Figure 1 for EPIC 247455141 (K2 C13): (a) detrended EVEREST light curve over the campaign with transit ticks for b, c and d;
(b)-(d) folded b, c, d, each with the other two signals' median posterior models divided out, 20-min bins and the 16-84% band of 300 posterior models.
Model, posterior and light curve: the joint three-planet fit (fit3.py, shared stellar density with TIC 8.2 prior) on the EVEREST light curve detrended with all three signals masked.
Main mode only: draws with ln p > max ln p - 30 (the chain has no second mode; this cut removes 3 of 181,288 draws).
usage: figure1.py LIGHTCURVE_NPZ CHAIN_NPZ OUT_PREFIX"""
import sys, json, numpy as np, batman, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
LC, CH, OUT = sys.argv[1:4]
z = np.load(LC); t, f, traw, fraw = z['t'], z['f'], z['t_raw'], z['f_raw']
c = np.load(CH); ch, lp, Tref = c['ch'], c['lp'], c['Tref']
keep = lp > lp.max() - 30; ch = ch[keep]
rng = np.random.default_rng(42); draws = ch[rng.choice(len(ch), 300, replace=False)]
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
def a_rs(rho, P): return ((G*rho*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
def planet(th, i, tt):
    t0, P, lrp, b = th[4*i:4*i+4]; rho, q1, q2 = th[12], th[13], th[14]; a = a_rs(rho, P)
    pm = batman.TransitParams(); pm.t0 = t0; pm.per = P; pm.rp = np.exp(lrp); pm.a = a; pm.inc = np.degrees(np.arccos(b/a)); pm.ecc = 0; pm.w = 90
    pm.limb_dark = 'quadratic'; pm.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
    return batman.TransitModel(pm, tt - Tref[i], supersample_factor=15, exp_time=1765.5/86400).light_curve(pm)
med = np.median(ch, axis=0); f0 = med[16]
Mods = [np.median([planet(d, i, t) for d in draws], axis=0) for i in range(3)]
Tc = [Tref[i] + med[4*i] for i in range(3)]; Ps = [med[4*i+1] for i in range(3)]
names = ['b', 'c', 'd']; cols = ['#D55E00', '#009E73', '#6E6E6E']
plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in', 'font.family': 'DejaVu Sans'})
fig = plt.figure(figsize=(3.45, 5.3))
gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 1, 1, 1], hspace=0.42, left=0.16, right=0.965, top=0.98, bottom=0.07)
ax = fig.add_subplot(gs[0])
y = (f/f0 - 1)*1e3
ax.plot(t, y, '.', ms=1.0, color='0.65', rasterized=True)
bw = 2/24; e = np.arange(t.min(), t.max() + bw, bw); idx = np.digitize(t, e)
bt = np.array([t[idx == k].mean() for k in np.unique(idx)]); by = np.array([y[idx == k].mean() for k in np.unique(idx)])
ax.plot(bt, by, '-', lw=0.5, color='k')
ylo = -3.4
for i in range(3):
    tt = Tc[i] + Ps[i]*np.arange(-40, 41); tt = tt[(tt > traw.min()) & (tt < traw.max())]
    ax.plot(tt, np.full(len(tt), ylo + 0.45*(2 - i) + 0.2), '|', ms=4, mew=0.9, color=cols[i])
    ax.text(traw.max() + 0.6, ylo + 0.45*(2 - i) + 0.2, names[i], color=cols[i], va='center', fontsize=6.5)
ax.set_xlim(traw.min() - 0.5, traw.max() + 2.2); ax.set_ylim(ylo - 0.3, 3.3)
ax.set_xlabel('Time (BJD$_{\\rm TDB}$ $-$ 2454833)', labelpad=1.5); ax.set_ylabel('Rel. flux (ppt)')
ax.text(0.015, 0.93, '(a) K2 C13, EVEREST, detrended; line: 2-h bins', transform=ax.transAxes, va='top', fontsize=6.5, bbox=dict(fc='w', ec='none', alpha=0.85, pad=0.6))
lab = ['(b) b: $P$ = 3.0748 d, 1.30 R$_\\oplus$', '(c) c: $P$ = 5.4557 d, 1.16 R$_\\oplus$', '(d) d: $P$ = 9.3585 d, candidate signal, not vetted (FPP 0.07)']
W = 4.0; info = {}
for i in range(3):
    a = fig.add_subplot(gs[i+1])
    others = np.prod([Mods[j] for j in range(3) if j != i], axis=0)
    r = f/f0/others
    ph = ((t - Tc[i] + 0.5*Ps[i]) % Ps[i] - 0.5*Ps[i])*24; m = np.abs(ph) < W
    a.plot(ph[m], (r[m] - 1)*1e3, '.', ms=1.4, color='0.7', rasterized=True)
    be = np.arange(-W, W + 1e-9, 1/3); k = np.digitize(ph[m], be)
    xb, yb, eb = [], [], []
    for kk in np.unique(k):
        s = (r[m] - 1)[k == kk]*1e3
        if len(s) >= 3: xb.append(ph[m][k == kk].mean()); yb.append(s.mean()); eb.append(s.std(ddof=1)/np.sqrt(len(s)))
    a.errorbar(xb, yb, eb, fmt='o', ms=2.2, color='k', elinewidth=0.6, capsize=0, zorder=4)
    tm = np.linspace(-W, W, 600)
    mods = np.array([planet(d, i, Tc[i] + tm/24) for d in draws])
    lo, mid, hi = np.percentile((mods - 1)*1e3, [16, 50, 84], axis=0)
    a.fill_between(tm, lo, hi, color=cols[i], alpha=0.35, lw=0, zorder=2)
    a.plot(tm, mid, color=cols[i], lw=1.0, ls='--' if i == 2 else '-', zorder=3)
    a.set_xlim(-W, W); a.set_ylim(-2.6, 2.0)
    a.set_ylabel('Rel. flux (ppt)'); a.text(0.015, 0.95, lab[i], transform=a.transAxes, va='top', fontsize=6.5, color='k' if i < 2 else '0.25', bbox=dict(fc='w', ec='none', alpha=0.85, pad=0.6), zorder=6)
    if i == 2: a.set_xlabel('Hours from mid-transit', labelpad=1.5)
    info[names[i]] = dict(n_points_window=int(m.sum()), n_bins=len(xb), median_bin_err_ppt=float(np.median(eb)), model_min_depth_ppm=float(-mid.min()*1e3), Tc_BKJD=float(Tc[i]), P_d=float(Ps[i]))
fig.savefig(OUT + '.png', dpi=300); fig.savefig(OUT + '.pdf')
info.update(n_draws_total=int(len(keep)), n_draws_main_mode=int(keep.sum()), n_draws_used=300, f0=float(f0), size_in=[3.45, 5.3], dpi=300)
json.dump(info, open(OUT + '.json', 'w'), indent=1); print(json.dumps(info, indent=1))
