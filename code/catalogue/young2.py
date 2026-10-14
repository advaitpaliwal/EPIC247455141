# Supplement to young.py: columns young.py did not auto-detect. Hyades parallax ranges (Gaia Collab. 2018 'plx'; Lodieu+2019 'Plx', pmRA*, pmDE, RV) and Taurus UVW (Krolikowski+2021 U,V,W) vs target; target source_id membership in each list.
import json, numpy as np, warnings; warnings.filterwarnings('ignore'); import astropy.units as u
from astropy.coordinates import SkyCoord; from astroquery.vizier import Vizier
c0 = SkyCoord(69.53047*u.deg, 21.84265*u.deg); out = {}
def q(cat, r=15): return Vizier(columns=['**'], row_limit=-1, timeout=300).query_region(c0, radius=r*u.deg, catalog=cat)[0]
h = q('J/A+A/616/A10/tablea1a'); out['Hyades_GC18'] = dict(n=len(h), plx_pct=np.percentile(np.asarray(h['plx'], float), [0, 1, 50, 99, 100]).tolist(), target_listed=bool(np.any(np.asarray(h['Source']).astype(str) == '144850112816196992')))
l = q('J/A+A/623/A35/tablec1'); out['Hyades_L19'] = dict(n=len(l), plx_pct=np.nanpercentile(np.asarray(l['Plx'], float), [0, 1, 50, 99, 100]).tolist(), pmra_pct=np.nanpercentile(np.asarray(l['pmRA_'], float) if 'pmRA_' in l.colnames else np.asarray(l[l.colnames[7]], float), [1, 50, 99]).tolist(), rv_pct=np.nanpercentile(np.asarray(l['RV'], float), [1, 50, 99]).tolist())
k = q('J/AJ/162/110/table1'); U, V, W = [np.asarray(k[x], float) for x in 'UVW']; ok = np.isfinite(U)
out['Taurus_K21'] = dict(n=len(k), UVW_median=[float(np.median(U[ok])), float(np.median(V[ok])), float(np.median(W[ok]))], UVW_p1_p99=[np.percentile(x[ok], [1, 99]).tolist() for x in (U, V, W)], target_listed=bool(np.any(np.asarray(k['Gaia']).astype(str) == '144850112816196992')))
json.dump(out, open('k2-deep/247455141/catalogue/young2.json', 'w'), indent=1); print(json.dumps(out, indent=1))
open('k2-deep/247455141/catalogue/young2.json.url', 'w').write('VizieR query_region (astroquery) r=15deg around 69.53047 +21.84265: J/A+A/616/A10/tablea1a, J/A+A/623/A35/tablec1, J/AJ/162/110/table1\n')
