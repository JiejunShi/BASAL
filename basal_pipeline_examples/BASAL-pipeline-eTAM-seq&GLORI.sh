#!/bin/bash
# Below is the BASAL pipeline for analyzing eTAM-seq and GLORI

# 1 Prepare your reads after trimming

# 2 Mapping reads
## 2.1 Mapping reads to genome
basal -p {CORES} \
-a {input.fq.gz} -d {genome.fa} \
-o {map2genome.bam} \
-M A:G -u
# {CORES}: core number
# {input.fq.gz}: pre-processed read file (.fq.gz)
# {genome.fa}: genome fasta file
# {map2genome.bam}: output bam file (.bam)

# reads successfully aligned to genome were extracted
samtools view -b -F 3588 -@ {CORES} -o {genomeAlign.bam} {map2genome.bam}
sambamba sort -m 8GB -t {CORES} -o {tmp.bam} {genomeAlign.bam};mv {tmp.bam} {genomeAlign.bam};
samtools view -b --include-flags 4 -@ {CORES} -o {unmap2genome.bam} {map2genome.bam}
# {CORES}: core number used
# {map2genome.bam}: input reads mapped to genome in last step (.bam)
# {tmp.bam}: temperory file
# {genomeAlign.bam}: output reads mapped to genome (.bam)
# {unmap2genome.bam}: output reads unmapped to genome (.bam)

## 2.2 Mapping reads to transcriptome
basal -p {CORES} \
-a {unmap2genome.bam} -d {transcriptome.fa} \
-o {trxptomeAlign.bam} \
-M A:G
# {CORES}: core number
# {unmap2genome.bam}: reads unmapped to genome (.bam)
# {transcriptome.fa}: transcriptome fasta file
# {trxptomeAlign.bam}: output reads map to transcriptome (.bam)

## 2.3 merge genome and transcriptome alignments
python basalkit.py mergeBAM {trxptomeAlign.bam} {genomeAlign.bam} \
{gtf} \
-o {merged}
# {trxptomeAlign.bam}: input transcriptome bam file(.bam)
# {genomeAlign.bam}: input genome bam file(.bam)
# {merged}: output merged bam file prefix (.bam)
# {gtf}: annotation gtf file

# 3 modification sites detection
## 3.1 measure average modification level of each site
python basalkit.py avgmod {merged.sorted.bam} {genome.fa} \ 
-o {treat|ctrl} \ 
-M A:G
# {merged.bam}: input merge bam file (.bam)
# {genome.fa}: genome fasta file
# {treat|ctrl}: output avgmod file prefix (_AvgMod.tsv)

## 3.2 perform significant test of Treat vs Ctrl and calculate FDR
python basalkit.py fdr {treat_AvgMod.tsv.gz} -c {ctrl_AvgMod.tsv.gz} \
-o {output_FDR}
# {treat_AvgMod.tsv.gz}: input AvgMod file in treat
# {ctrl_AvgMod.tsv.gz}: input AvgMod file in ctrl
# {output_FDR}: output file prefix (_FDR.tsv.gz)
