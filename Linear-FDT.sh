#!/bin/bash

if [ $# -lt 1 ]; then
      NBTHREADS=32
else
      NBTHREADS=$1
fi

SECU=2048
ITER=11

LINFIL=bench_lin_vespo_P254.txt


for LogDeg in {8..17}
do
  OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench $((2**LogDeg)) ${SECU} ${ITER} &>> ${LINFIL}
done
