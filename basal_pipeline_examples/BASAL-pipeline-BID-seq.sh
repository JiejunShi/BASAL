#!/bin/bash
# Below is the BASAL pipeline for analyzing BID-seq

# 1 Prepare your reads after trimming

# 2 Mapping reads
## 2.1 Mapping reads to genome
basal -p {CORES} \
-a {input.fq.gz} -d {genome.fa} \
-o {map2genome.bam} \
-M T:- -n 1 -g 3 -R -u
# {CORES}: core number used
# {input.fq.gz}: pre-processed read file (.fq.gz)
# {genome.fa}: genome fasta file
# {map2genome.bam}: output bam file (.bam)

# extract reads aligned to genome
samtools view -b -F 3588 -@ {CORES} -o {tmp.bam} {map2genome.bam}
sambamba sort -m 8GB -t {CORES} -o {genomeAlign.bam} {tmp.bam}
# {CORES}: core number
# {map2genome.bam}: input reads mapped to genome in last step (.bam)
# {tmp.bam}: temperory file
# {genomeAlign.bam}: output reads mapped to genome (.bam)

# extract reads not aligned to genome
samtools view -b --include-flags 4 -@{CORES} -o {unmap2genome.bam} {map2genome.bam}
samtools fastq {unmap2genome.bam} > {unmap2genome.fq}
gzip {unmap2genome.fq}
# {CORES}: core number used
# {map2genome.bam}: input reads mapped to genome in last step (.bam)
# {unmap2genome.bam}: output reads unmapped to genome (.bam)
# {unmap2genome.fq}: output reads unmapped to genome (.fq)

## 2.2 Mapping reads to transcriptome
basal -p {CORES} \
-a {unmap2genome.fq.gz} -d {transcriptome.fa} \
-o {tmp.bam} \
-M T:- -n 1 -g 3 -R;
# {CORES}: core number used
# {unmap2genome.fq.gz}: input reads unmapped to genome (.fq.gz)
# {transcriptome.fa}: transcriptome fasta file
# {tmp.bam}: temperory file

# sort reads mapping to transcriptome
sambamba sort -m 8GB -t {CORES} -o {trxptomeAlign.bam} {tmp.bam};
# {CORES}: core number used
# {trxptomeAlign.bam}: output sorted reads map to transcriptome (.bam)
# {tmp.bam}: temperory file from last step

## 2.3 CIGAR correction for consecutive pU ('-R' must be used in basal)
python basalkit.py shiftD {genomeAlign.bam} -o {tmp}
sambamba sort -m 8GB -t {CORES} -o {genomeAlign.corrected.bam} {tmp}
python basalkit.py shiftD {trxptomeAlign.bam} -o {tmp}
sambamba sort -m 8GB -t {CORES} -o {trxptomeAlign.corrected.bam} {tmp}
# {CORES}: core number
# {output.tmp}: temperory file
# {genomeAlign.corrected.bam}: output corrected genome reads (.bam)
# {trxptomeAlign.corrected.bam}: output corrected transcriptome reads (.bam)

## 2.4 merge genome and transcriptome bam file
python basalkit.py mergeBAM {trxptomeAlign.corrected.bam} {genomeAlign.corrected.bam} \
{gtf} \
-o {output_prefix}
# {trxptomeAlign.corrected.bam}: transcriptome bam file(.bam)
# {genomeAlign.corrected.bam}: genome bam file(.bam)
# {merged}: output merged bam file prefix (.bam)
# {gtf}: annotation gtf file

# 3 measure sites ratio
## 3.1 measure average modification with basalkit "avgmod"
python basalkit.py avgmod {merge.bam} {genome.fa} \ 
-o {output_prefix} \ 
-M T:- -D M -T RNA -y 7
# {merge.bam}: input merged bam file (.bam)
# {genome.fa}: genome fasta file
# {treat|ctrl}: output avgmod file prefix (_AvgMod.tsv)

## 3.2 perform significant test of Treat vs Ctrl and calculate FDR
python basalkit.py fdr {treat_AvgMod.tsv.gz} -c {ctrl_AvgMod.tsv.gz} \
-o {output_FDR}
# {treat_AvgMod.tsv.gz}: input AvgMod file in treat
# {ctrl_AvgMod.tsv.gz}: input AvgMod file in ctrl
# {output_FDR}: output file prefix (_FDR.tsv.gz)
