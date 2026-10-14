# EPIC 247455141: two small planet candidates at 3.07 and 5.46 days around an M dwarf from K2 Campaign 13

Advait Paliwal (independent researcher, advaitpaliwal@gmail.com), 2026-10-13.

Two dips in brightness repeat in K2 Campaign 13 (2017) photometry of EPIC 247455141, an M dwarf 147 pc away: one of about 0.11% every 3.07 days (26 transits) and one of about 0.09% every 5.46 days (15 transits). If they are planets, they are about 1.3 and 1.2 times the size of Earth. A third signal of about 0.13% every 9.36 days (9 transits) is an unvetted candidate signal, listed for completeness only.

| | b | c | d (unvetted signal) |
|---|---|---|---|
| Period (d) | 3.074847 (+0.000242/−0.000245) | 5.455730 (+0.000549/−0.000580) | 9.358462 (+0.001234/−0.001063) |
| Mid-transit (BJD_TDB) | 2457859.8542 ± 0.0018 | 2457861.2849 ± 0.0025 | 2457862.4613 ± 0.0026 |
| Depth (ppm) / duration (h) | 1091 ± 97 / 1.34 (+0.10/−0.12) | 892 (+97/−99) / 1.70 (+0.06/−0.12) | 1339 (+161/−156) / 1.51 (+0.17/−0.14) |
| Radius (R⊕) | 1.30 ± 0.09 | 1.16 ± 0.08 | 1.54 (+0.14/−0.13) |
| Equilibrium temperature (K) / insolation (S⊕) | 617 ± 27 / 24 (+5/−4) | 510 ± 23 / 11.3 (+2.1/−1.9) | 426 ± 19 / 5.5 (+1.0/−0.9) |
| False-positive probability (TRICERATOPS + MOLUSC, 10 runs) | 0.014 ± 0.002 | 0.014 ± 0.002 | 0.072 ± 0.006 |

Host: EPIC 247455141 = TIC 118598568 = Gaia DR3 144850112816196992 = 2MASS J04380724+2150334; TIC Teff 3550 ± 157 K, R 0.391 ± 0.012 R☉, M 0.381 ± 0.021 M☉; Gaia DR3 parallax 6.794 ± 0.040 mas, G = 15.56, RUWE 1.03; Kp = 15.29.

Parameters come from one joint fit of all three signals to the EVEREST light curve, with a shared stellar density (prior from TIC v8.2). Fitting b and c alone changes no parameter of either by more than 0.08σ.

Status: **planet candidates** (b, c); d is an **unvetted candidate signal**, not a planet candidate. None is a TOI, CTOI or K2 candidate in the NASA Exoplanet Archive or ExoFOP, and the star is not in SIMBAD, as of 2026-10-13. For b and c the transit depths agree across the EVEREST, K2SFF and our own SFF light curves, there is no secondary eclipse, and difference-image centroids are no further from the star than in transit-free difference images (the same test on the known planet K2-133 b, also in C13, gives the same result). Signal b shows a 1.7σ odd/even difference in EVEREST; 14% of transit-free trial epochs give at least as much. Signal d is the top periodogram peak once b and c are removed and appears in our own SFF light curves, but it is not found in the K2SFF light curve, its false-positive probability stays above 0.035 in all 10 runs, and its duration fits the catalogue stellar density only with a high impact parameter (0.73), and a fit with its own density prefers a star 3.3 times denser. With the Chen & Kipping (2017) mass–radius relation, adjacent pairs are separated by 17 (b–c) and 15 (c–d) mutual Hill radii. No multiplicity factor is applied to any false-positive probability. TESS (Sectors 43, 44, 71) is not sensitive enough to recover these transits at T = 14.46.

## Contents

- `note/rnaas_EPIC247455141.md` — the Research Note draft
- `note/figure1_EPIC247455141.png` (and `.pdf`) — (a) the detrended C13 light curve with transit times of b, c and d; (b)–(d) b, c and d folded, each with the other two signals removed, 20-minute bins and the posterior model band; d is marked as an unvetted candidate signal
- `data/lc/` — the detrended light curve with phases and models for b, c, d (`everest_detrended_m3_c13.csv`), the undetrended EVEREST light curve, box depths per reduction, and the note on removing a correction spike from one SFF reduction
- `data/fit/` — posterior summaries of the three-planet fit (adopted) and the two-planet fit, all fit outputs (including free-density and other-reduction fits), transit-timing overlaps and Hill-stability results
- `data/vetting/`, `data/injection/` — odd/even, secondary-eclipse, phase-scan and per-transit checks with their transit-free calibrations; periodogram searches for d; injection-recovery and null searches
- `data/pixels/` — difference images and their null tests for b, c and d, with the same tests on K2-133 b for reference
- `data/fpp/` — false-positive probability runs (per-run results, scenario tables, MOLUSC summary, TRICERATOPS star list) for b, c (`b_c/`) and d (`d/`)
- `data/tess/`, `data/catalogue/`, `data/search/` — TESS coverage and sensitivity, catalogue queries with their URLs (including the 2026-10-13 recheck), and the search detections
- `code/` — search, light-curve, fit, vetting, pixel, FPP, TESS, catalogue, dynamics and figure scripts

Data: K2 (NASA, via MAST), EVEREST and K2SFF high-level science products, TESS (TESScut and QLP), ESA Gaia DR3, TIC v8.2, LAMOST. Software: lightkurve, batman, emcee, transitleastsquares, wotan, TRICERATOPS, MOLUSC, forecaster, tess-point, astropy, matplotlib. Licences: code MIT; text, figures and data CC BY 4.0.

Scripts are kept as run. File locations inside them are given relative to a single working directory (for example `k2-deep/247455141/`); set that directory up, or edit the paths, to rerun them.
