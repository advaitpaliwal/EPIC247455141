# Copied from ../../247192465/fpp_joint/scripts/aggregate.py: b/c with MOLUSC (10 runs each) and without (1 run each).
import json, hashlib, numpy as np, pandas as pd
from pathlib import Path
D = Path('k2-deep/247455141/fpp_joint')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
out = {'input_sha256': {f: sha(D / f) for f in ['out/target_247455141.pkl', 'molusc_out/molusc_for_trc_N1e6.csv']}}
for v, nr in [('b_molusc', 10), ('b_none', 1), ('c_molusc', 10), ('c_none', 1)]:
    rows, sc = [], []
    for k in range(nr):
        r = json.loads((D / f'out/res_{v}_{k}.json').read_text()); p = pd.read_csv(D / f'out/probs_{v}_{k}.csv')
        assert np.isclose(p.prob.sum(), 1, atol=1e-9) and (p.prob >= 0).all()
        fpp = 1 - p.loc[p.scenario.isin(['TP', 'PTP', 'DTP']), 'prob'].sum(); nfpp = p.loc[p.scenario.isin(['NTP', 'NEB', 'NEBx2P']), 'prob'].sum()
        assert np.isclose(fpp, r['FPP'], atol=1e-9) and np.isclose(nfpp, r['NFPP'], atol=1e-9)
        rows.append({x: r[x] for x in r if x != 'aperture'}); sc.append(p.groupby('scenario').prob.sum().rename(k))
    d = pd.DataFrame(rows); s = pd.concat(sc, axis=1).fillna(0)
    st = pd.DataFrame({'mean': s.mean(axis=1), 'sd': s.std(axis=1, ddof=1), 'min': s.min(axis=1), 'max': s.max(axis=1)})
    st.to_csv(D / f'scenario_summary_{v}.csv'); d.to_csv(D / f'run_results_{v}.csv', index=False)
    out[v] = {'seeds': d.seed.tolist(), 'err_ppm': float(d.err_ppm.iloc[0]), 'nbins': int(d.nbins.iloc[0]), 'depth_ppm': float(d.depth_ppm.iloc[0]),
              'FPP': {'mean': d.FPP.mean(), 'sd': d.FPP.std(ddof=1), 'sem': d.FPP.std(ddof=1)/np.sqrt(len(d)), 'min': d.FPP.min(), 'max': d.FPP.max(), 'runs_le_0.035': int((d.FPP <= 0.035).sum()), 'all': d.FPP.tolist()},
              'NFPP': {'mean': d.NFPP.mean(), 'max': d.NFPP.max()},
              'scenarios': {i: {c: float(st.loc[i, c]) for c in st.columns} for i in st.index if st.loc[i, 'mean'] > 1e-6}}
(D / 'summary.json').write_text(json.dumps(out, indent=1, default=float) + '\n'); print(json.dumps(out, indent=1, default=float))
