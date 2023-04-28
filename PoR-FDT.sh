#!/bin/bash

#####################################################
### VESPo benchmarks for PoR on 1GB, 10GB, 100GB, 1TB
### Copyright(c) 2023 Jean-Guillaume Dumas
### usage: PoR-FDT.sh [threads] [iterations]
#####################################################


# Number of threads
if [ $# -lt 1 ]; then
      NBTHREADS=32
else
      NBTHREADS=$1
fi

# Number of iterations
if [ $# -lt 3 ]; then
      ITER=3
else
      ITER=$3
fi

# Paillier security
SECU=2048

SQFIL=bench_sq_vespo_P254.txt

echo "##### VESPO PoR benchmarks"
echo "##### 8 degrees from 5816 to 4026778"

PSINIT=0
if [ -e ${SQFIL} ]
then
  PSINIT=`grep VAUDIT ${SQFIL} | grep OK | wc -l`
fi
PASSED=0
PERCEN=0
TOTALT=$(( (4+4) * ITER))

for DEGREE in 5816 18390 58154 186093
do
    echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, current degree: ${DEGREE} ...\r"

    # Actual benchmark:
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${SQFIL}
    PASSED=`grep VAUDIT ${SQFIL} | grep OK | wc -l`
    PERCEN=$(( (100 * (PASSED - PSINIT)) / TOTALT ))
done

REFIL=bench_re_vespo_P254.txt
PRINIT=0
if [ -e ${REFIL} ]
then
  PRINIT=`grep VAUDIT ${REFIL} | grep OK | wc -l`
fi

for DEGREE in 5125 46551 426519 4026778
do
    echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, current degree: ${DEGREE} ...\r"

    # Actual benchmark:
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${REFIL}

    PASSED=`grep VAUDIT ${SQFIL} ${REFIL} | grep OK | wc -l`
    PERCEN=$(( (100 * (PASSED - PSINIT - PRINIT)) / TOTALT ))
done

echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m                                 \n"
