echo "PASSED: "`fgrep VAUDIT bench_* | wc -l`


egrep '(STATE|TIMING)' bench_sq_vespo_P254.txt bench_re_vespo_P254.txt | xargs -n2 -d'\n' | awk 'BEGIN { print "Degree Setup CStore CTime STime Utime Htime" } {printf("%d %.2fs %.2fkB %.2fms %.2fs %.2fms %.2fms\n",$8,$14,$3/8192.,$18*1000.,$22,$30*1000.,$26*1000.)}' | tee res_vespo_P254.txt

fgrep DEVIATIO bench_sq_vespo_P254.txt bench_re_vespo_P254.txt | sed 's/\%//g' | awk 'BEGIN {c=0;s=0} { if (c<$9) { c=$9 } ; if (s<$13) { s=$13 } } END {m=c; if (m<s) m=s; printf("%.2f%% (%f%% %f%%)\n",m,c,s)}'
