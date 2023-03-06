#!/bin/bash

for FILE in $*
do
	echo "############### VESPO library benchmarks ###############"
	echo -en "$FILE, \e[32mPASSED: "`fgrep VAUDIT ${FILE} | wc -l`
	echo -e "\e[0m"

	egrep '(STATE|TIMING)' ${FILE} | xargs -n2 -d'\n' | awk 'BEGIN { print "Degree Setup CStore CTime STime Utime Htime" } {printf("%d %.2fs %.2fkB %.2fms %.2fs %.2fms %.2fms\n",$8,$14,$3/8192.,$18*1000.,$22,$30*1000.,$26*1000.)}'

	fgrep DEVIATIO ${FILE} | sed 's/\%//g' | awk 'BEGIN {c=0;s=0} { if (c<$9) { c=$9 } ; if (s<$13) { s=$13 } } END {m=c; if (m<s) m=s; printf("%.2f%% (%f%% %f%%)\n",m,c,s)}'
done
