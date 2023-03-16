#!/bin/bash

if [ $# -lt 1 ]; then
      NBTHREADS=32
else
      NBTHREADS=$1
fi

SECU=2048
ITER=3

SQFIL=bench_sq_vespo_P254.txt

echo "##### VESPO PoR benchmarks"
echo "##### 8 degrees from 5816 to 4026778"

PASSED=0
PERCEN=0
TOTALT=$(( (4+4) * ITER))

for DEGREE in 5816 18390 58154 186093
do
    echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, current degree: ${DEGREE} ...\r"

    # Actual benchmark:
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${SQFIL}
    PASSED=`grep VAUDIT ${SQFIL} | grep OK | wc -l`
    PERCEN=$(( (100 * PASSED) / TOTALT ))
done

REFIL=bench_re_vespo_P254.txt

for DEGREE in 5125 46551 426519 4026778
do
    echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, current degree: ${DEGREE} ...\r"

    # Actual benchmark:
    OMP_NUM_THREADS=${NBTHREADS} ./vespo_bench ${DEGREE} ${SECU} ${ITER} &>> ${REFIL}

    PASSED=`grep VAUDIT ${SQFIL} ${REFIL} | grep OK | wc -l`
    PERCEN=$(( (100 * PASSED) / TOTALT ))
done

echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m                                 \n"
