#!/bin/bash

#####################################################
### VESPo testing range of degrees with given tasks
### Copyright(c) 2023 Jean-Guillaume Dumas
#####################################################

START=$1
ENDLO=$2

if [ $# -lt 3 ]; then
      NBTASKS=4
else
      NBTASKS=$3
fi

if [ $# -lt 4 ]; then
      THREADS=${NBTASKS}
else
      THREADS=$4
fi

SECU=2048
ITER=2

LINFIL=/tmp/bench_range_vespo_P254.txt

SIZLD=$((ENDLO - START + 1))
PINIT=0
if [ -e $file ]
then
  PINIT=`grep VAUDIT ${LINFIL} | grep OK | wc -l`
fi
PCURR=${PINIT}
PASSED=$((PCURR-PINIT))
PERCEN=0
TOTALT=$((SIZLD * ITER))

i=${START}
while [ $i -le $ENDLO ]
do
  echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m, OMP_NUM_THREADS=${THREADS} ./vespo_bench $i ${SECU} ${ITER} ${NBTASKS}\r"

  # Benchmaring
  OMP_NUM_THREADS=${THREADS} ./vespo_bench $i ${SECU} ${ITER} ${NBTASKS} &>> ${LINFIL}

  i=$(( $i + 1 ))
  PCURR=`grep VAUDIT ${LINFIL} | grep OK | wc -l`
  PASSED=$((PCURR-PINIT))
  PERCEN=$(( (100 * PASSED) / TOTALT ))
done
echo -ne "##### Passed: \e[32m${PERCEN}%\e[0m                                                   \n"
