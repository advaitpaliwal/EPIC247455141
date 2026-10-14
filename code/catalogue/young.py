# Young-group check: Taurus (Luhman 2018, J/AJ/156/271; Krolikowski+2021, J/AJ/162/110) and Hyades (Gaia Collab. 2018, J/A+A/616/A10 tablea1a; Lodieu+2019, J/A+A/623/A35)
# members within 15 deg of the target via VizieR; compare target parallax, proper motion and LAMOST RV with member distributions; also UVW (astropy) of target and members with RV.
import json, numpy as np, warnings; warnings.filterwarnings('ignore')
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
from astroquery.vizier import Vizier
D = 'k2-deep/247455141'
c0 = SkyCoord(69.53047*u.deg, 21.84265*u.deg)
T = dict(plx=6.7936, pmra=57.653, pmdec=2.141, rv=np.mean([59.7, 54.7]))   # Gaia DR3 (catalogue/gaia_astrom.csv); LAMOST LRS RVs 59.7+-2.1, 54.7+-2.2 (J/ApJS/253/19, vizier_all_5as.tsv)
out = dict(target=T)
cats = {'Taurus_Luhman2018': 'J/AJ/156/271', 'Taurus_Krolikowski2021': 'J/AJ/162/110', 'Hyades_GaiaCollab2018': 'J/A+A/616/A10/tablea1a', 'Hyades_Lodieu2019': 'J/A+A/623/A35'}
urls = []
for k, cat in cats.items():
    try:
        v = Vizier(columns=['**', '_RAJ2000', '_DEJ2000'], row_limit=-1, timeout=300); r = v.query_region(c0, radius=15*u.deg, catalog=cat)
    except Exception as e: out[k] = str(e); continue
    urls.append(f'VizieR query_region {cat} r=15deg around 69.53047 +21.84265 (astroquery)')
    res = {}
    for t in r:
        cols = t.colnames
        pk = next((c for c in cols if c.lower() in ('plx', 'plxg', 'plx_g', 'plx3')), None); ra = next((c for c in cols if c.lower() in ('pmra', 'pmra3', 'pmrag')), None); de = next((c for c in cols if c.lower() in ('pmde', 'pmdec', 'pmde3', 'pmdeg')), None)
        sep = SkyCoord(np.asarray(t['_RAJ2000'], float)*u.deg, np.asarray(t['_DEJ2000'], float)*u.deg).separation(c0).deg
        d = dict(n=len(t), cols=cols[:40], min_sep_deg=float(np.nanmin(sep)))
        if pk and ra and de:
            p, a, b = [np.asarray(t[x], float) for x in (pk, ra, de)]; ok = np.isfinite(p) & np.isfinite(a) & np.isfinite(b)
            d.update(plx_pct=np.percentile(p[ok], [1, 50, 99]).tolist(), pmra_pct=np.percentile(a[ok], [1, 50, 99]).tolist(), pmdec_pct=np.percentile(b[ok], [1, 50, 99]).tolist(),
                     n_within_target_box=int(np.sum(ok & (np.abs(p - T['plx']) < 1.5) & (np.hypot(a - T['pmra'], b - T['pmdec']) < 10))))
        res[t.meta.get('name', str(len(res)))] = d
    out[k] = res
g = SkyCoord(ra=c0.ra, dec=c0.dec, distance=(1000/T['plx'])*u.pc, pm_ra_cosdec=T['pmra']*u.mas/u.yr, pm_dec=T['pmdec']*u.mas/u.yr, radial_velocity=T['rv']*u.km/u.s).galactic
v = g.velocity; out['target_UVW_kms_heliocentric'] = [float(v.d_x.to(u.km/u.s).value), float(v.d_y.to(u.km/u.s).value), float(v.d_z.to(u.km/u.s).value)]
out['target_vtan_kms'] = float(4.74047*np.hypot(T['pmra'], T['pmdec'])/T['plx'])
json.dump(out, open(f'{D}/catalogue/young.json', 'w'), indent=1, default=str); open(f'{D}/catalogue/young.json.url', 'w').write('\n'.join(urls) + '\n')
print(json.dumps(out, indent=1, default=str)[:6000])
