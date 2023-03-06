#!/bin/bash

if [ $# -lt 1 ]; then
      NBTHREADS=32
else
      NBTHREADS=$1
fi

SECU=2048
ITER=11

SQFIL=bench_sq_vespo_P254.txt
for DEGREE in 5816 18390 58154 186093
do
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${SQFIL}
done

REFIL=bench_re_vespo_P254.txt
for DEGREE in 5125 46551 426519 4026778
do
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${REFIL}
done
