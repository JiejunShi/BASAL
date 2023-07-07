# This script is called by 'merge' function. DONOT run it alone!

cat $1"_t2g.tsv" $2 | sort -k1b,1 -k2n,2 -k3b,3 | awk '{if($1~/^chr/){print}}' > $1".mergetmp";

awk 'BEGIN{FS="\t";OFS="\t";chr="";pos=""}{if(NR==1){print}else{if($4!=""){if($1==chr && $2==pos){eff_coverage+=$6;N_mod+=$7;N_total+=$8;n+=1}else{if(chr!="" && N_total>0){ratio=N_mod/N_total;print chr,pos,strand,context,ratio,eff_coverage,N_mod,N_total};chr=$1;pos=$2;strand=$3;context=$4;eff_coverage=$6;N_mod=$7;N_total=$8}}}}END{if(chr!="" && N_total>0){ratio=N_mod/N_total;print chr,pos,strand,context,ratio,eff_coverage,N_mod,N_total}}' $1".mergetmp" > $1"_AvgMod_merged.tsv";

gzip $1"_AvgMod_merged.tsv";

rm $1".mergetmp";

rm $1"_t2g.tsv";
