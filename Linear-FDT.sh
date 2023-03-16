#!/bin/bash

#####################################################
### VESPo benchmarks with increasing degrees
### Copyright(c) 2023 Jean-Guillaume Dumas
#####################################################


if [ $# -lt 1 ]; then
      NBTHREADS=32
else
      NBTHREADS=$1
fi

# Number of iterations
ITER=3

# Paillier security
SECU=2048

# Polynomial degrees from 2^8 to 2^16
BEGLD=8
ENDLD=16
SIZLD=$((ENDLD - BEGLD + 1))

LINFIL=bench_lin_vespo_P254.txt

echo "##### VESPO linear benchmarks"
echo "##### Doubling degrees from $((2**BEGLD)) to $((2**ENDLD))"
PASSED=0
PERCEN=0
TOTALT=$((SIZLD * ITER))

for (( LogDeg=$BEGLD; LogDeg<=$ENDLD; LogDeg++ ))
do
  DEGREE=$((2**LogDeg))
  echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, current degree: ${DEGREE} ...\r"

  # Actual benchmark:
  OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${LINFIL}

  PASSED=`grep VAUDIT ${LINFIL} | grep OK | wc -l`
  PERCEN=$(( (100 * PASSED) / TOTALT ))
done
echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m                                 \n"
