#fgrep SERVER bench_sq_i7_2048_P254.txt bench_re_i7_2048_P254.txt | sed 's/(//g' | xargs -n2 -d'\n' | awk '{print $6+1,$5,$14,$5+$14}' | sort -k4n
fgrep VAUDIT bench_sq_i7_2048_P254.txt bench_re_i7_2048_P254.txt | wc


egrep '(STATE|TIMING)' bench_sq_i7_2048_P254.txt | xargs -n2 -d'\n' | awk 'BEGIN { print "\\multirow{1}{*}{Matrix view} &\\multirow{1}{*}{Client Audit (dotproduct part)} &\\multirow{1}{*}{Client Audit (polynomial part)} &\\multirow{1}{*}{Server Audit (matrix-vector part)} &\\multirow{1}{*}{Server Audit (polynomial part)} &\\multirow{1}{*}{Client storage} &\\multirow{1}{*}{Communications} \\\\" } {printf("$ {\\times}%d$ & & $%.1f$ms & & $%.1f$s & $%.1f$KB & \\\\\n",$8,$18*1000.,$22,($3+2*$9)/1024./8.);}' | Transpose_Latex_Array.csh | tee res_sq_i7_2048_P254.txt

egrep '(STATE|TIMING)' bench_re_i7_2048_P254.txt | xargs -n2 -d'\n' | awk 'BEGIN { print "\\multirow{1}{*}{Matrix view} &\\multirow{1}{*}{Client Audit (dotproduct part)} &\\multirow{1}{*}{Client Audit (polynomial part)} &\\multirow{1}{*}{Server Audit (matrix-vector part)} &\\multirow{1}{*}{Server Audit (polynomial part)} &\\multirow{1}{*}{Client storage} &\\multirow{1}{*}{Communications} \\\\" } {printf("$ {\\times}%d$ & & $%.1f$ms & & $%.1f$s & $%.1f$KB & \\\\\n",$8,$18*1000.,$22,($3+2*$9)/1024./8.);}' | Transpose_Latex_Array.csh | tee res_re_i7_2048_P254.txt

fgrep DEVIATIO bench_sq_i7_2048_P254.txt bench_re_i7_2048_P254.txt | sed 's/\%//g' | awk 'BEGIN {c=0;s=0} { if (c<$9) { c=$9 } ; if (s<$13) { s=$13 } } END {m=c; if (m<s) m=s; printf("%.1f%% (%f%% %f%%)\n",m,c,s)}'
