#!/bin/zsh
# MOLUSC 2 x 100k (RUWE + Gaia contrast; target has a 5-parameter Gaia DR3 solution, astrometric_params_solved=31, RUWE 1.030) = 200,000 companions,
# M = 0.381 Msun (TIC 118598568, TIC 8.2 = TRICERATOPS M_s), default age. Copied from ../../211038138/fpp_joint/scripts/run_molusc.sh; run sequentially.
F=k2-deep/247455141/fpp_joint
PYM=python
cd $F
for s in 21 22; do
  SEED=$s $PYM -I -W ignore molusc_run.py -v --ruwe --gaia -a $F/molusc_out/m100k_s$s 04h38m07.2466s -- +21d50m33.504s 100000 0.381 > molusc_out/m100k_s$s.log 2>&1
done
echo MOLUSC_DONE >> molusc_out/pids
