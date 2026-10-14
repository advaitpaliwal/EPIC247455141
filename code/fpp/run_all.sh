#!/bin/sh
# Per signal: 1 reference run without MOLUSC, then 10 joint runs with MOLUSC (seeds 2000-2009); 1 worker.
D=k2-deep/247455141/fpp_joint
(echo b_none 0; echo c_none 0; seq 0 9 | sed "s/^/b_molusc /"; seq 0 9 | sed "s/^/c_molusc /") | xargs -P ${NP:-1} -n 2 sh -c 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -I -u '"$D"'/scripts/run_joint.py "$1" "$2" > '"$D"'/out/run_$1_$2.log 2>&1' sh
echo ALL_DONE > $D/out/done
