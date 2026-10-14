#!/bin/sh
# Signal d: 1 reference run without MOLUSC, then 10 joint runs with MOLUSC (seeds 2000-2009); NP workers (default 2).
(echo d_none 0; seq 0 9 | sed "s/^/d_molusc /") | xargs -P ${NP:-2} -n 2 sh -c 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -I -u k2-deep/247455141/fpp_joint_d/scripts/run_joint_d.py "$1" "$2" > k2-deep/247455141/fpp_joint_d/out/run_$1_$2.log 2>&1' sh
echo ALL_DONE > k2-deep/247455141/fpp_joint_d/out/done
