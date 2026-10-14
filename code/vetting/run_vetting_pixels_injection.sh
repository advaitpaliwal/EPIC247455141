#!/bin/sh
# vets (3 reductions x 2 signals), difference images + centroid null tests (b, c, and K2-133 b same-campaign comparison), injection-recovery.
B=k2-deep; D=$B/247455141; PY=python; cd $D
E=$($PY -I -c "import json;r=json.load(open('fit/result_everest_rho.json'));print(*[r[p][k][0] for p in 'bc' for k in ('P_d','T0_BKJD','T14_h')])")
set -- $E; PB=$1; TB=$2; DB=$3; PC=$4; TC=$5; DC=$6
for s in everest k2sff fresh40; do
  $PY -I scripts/vet2.py $s b $PB $TB $DB $PC $TC $DC > fit/vet_${s}_b.log 2>&1
  $PY -I scripts/vet2.py $s c $PC $TC $DC $PB $TB $DB > fit/vet_${s}_c.log 2>&1
done
T=$D/lc/mastDownload/K2/ktwo247455141-c13_lc/ktwo247455141-c13_lpd-targ.fits.gz
$PY -I scripts/pix2.py $T $PB $TB $DB 69.53046951614664 21.842649811547997 pix/diff_b $PC $TC $DC > pix/diff_b.log 2>&1
$PY -I scripts/pix2.py $T $PC $TC $DC 69.53046951614664 21.842649811547997 pix/diff_c $PB $TB $DB > pix/diff_c.log 2>&1
C=$(ls pix/cmp_K2-133/mastDownload/K2/*/*lpd-targ.fits.gz)
RD=$($PY -I -c "import lightkurve as lk;t=lk.read('$C');print(t.ra,t.dec)")
# K2-133 b: Wells+2018 via NEA ps (pix/cmp_K2-133/nea_K2-133.csv): P 3.0712 d, T0 BKJD 2988.3168, T14 1.22 h
$PY -I scripts/pix2.py $C 3.0712 2988.3168 1.22 $RD pix/cmp_K2-133/diff_K2-133b > pix/cmp_K2-133/run_b.log 2>&1
for w in b c; do
  $PY -I scripts/inj2.py $w 1 1 inj/real_$w.json real > inj/real_$w.log 2>&1
  $PY -I scripts/inj2.py $w 1 1 inj/null_$w.json null > inj/null_$w.log 2>&1
  $PY -I scripts/inj2.py $w 40 7 inj/inj_$w.json inj > inj/inj_$w.log 2>&1
done
echo LANEA_DONE > scripts/laneA.done
