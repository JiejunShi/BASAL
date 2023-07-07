# This script is called by 'epiallele' function. DONOT run it alone!

sort -k1b,1 -k2n,2 $1".unsrt" > $1".srt";

rm $1".unsrt";

uniq -c $1".srt" | awk '{for(i=2;i<=7;i++){printf $i"\t"};if(NF==7){print $1}else{print $1"\t"$8"\t"$9}}' > $1;

rm $1".srt";

bgzip $1;

tabix -b 2 -e 3 $1".gz";
