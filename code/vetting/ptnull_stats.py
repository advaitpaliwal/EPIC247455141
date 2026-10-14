# stats() copied verbatim from ptnull.py
import sys, json, numpy as np
sys.path.insert(0, 'k2-deep/scripts'); import lcs
def stats(t, f, P, T0, T14):
    ph = lcs.phase(t, P, T0); ep = np.round((t - T0)/P).astype(int); oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5); inn = np.abs(ph) < 0.35*T14
    def dep(mi, mo): return np.median(f[mo]) - np.mean(f[mi]), np.std(f[mo])/np.sqrt(mi.sum())
    o = dep(inn & (ep % 2 == 1), oot); e = dep(inn & (ep % 2 == 0), oot); oe = abs(o[0]-e[0])/np.hypot(o[1], e[1])
    ch, ce = [], []
    for k in np.unique(ep[inn]):
        m = inn & (ep == k); mo = oot & (ep == k); d, s = dep(m, mo if mo.sum() > 5 else oot); ch.append(d); ce.append(s)
    ch, ce = np.array(ch), np.array(ce); mu = np.average(ch, weights=1/ce**2)
    return float(np.sum(((ch-mu)/ce)**2)/(len(ch)-1)), float(oe), float(mu/np.sqrt(1/np.sum(1/ce**2)))
