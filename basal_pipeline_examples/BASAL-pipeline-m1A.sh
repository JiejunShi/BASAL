#!/bin/bash
#===for m1A=====

# 1 Prepare your reads after trimming
# 2 Mapping reads
## 2.1 Mapping reads to genome
basal -p {CORES} \
-a {input.processed_Read} -d {genome.fa} \
-o {output.map2genome.bam} \
-M A:CGT -n 1 -S 2024 -u
# {CORES}: core number used
# {input.processed_Read}: pre-processed read file (.fq.gz)
# {genome.fa}: genome fasta file
# {output.map2genome.bam}: output bam file (.bam)


# reads that were successfully matched to the genome were extracted
samtools view -b -F 3588 -@ {CORES} -o {output.genomeAlign.bam} {input.map2genome.bam}
sambamba sort -m 8GB -t {CORES} -o {output.srt.bam} {input.genomeAlign.bam};mv {input.srt.bam} {output.genomeAlign.bam};
samtools view -b --include-flags 4 -@ {CORES} -o {output.unmap2genome.bam} {input.map2genome.bam}
# {CORES}: core number used
# {input.map2genome.bam}: input map to genome bam file with BASAL (.bam)
# {output.tmp.bam}: temperory file
# {output.genomeAlign.bam}: output map to genome sorted BAM file (.bam)

## 2.2 Mapping reads to transcriptome
basal -p {CORES} \
-a {input.unmap2genome.bam} -d {transcriptome.fa} \
-M A:CGT -S 2024 \
-o {ouput.trxptomeAlign.bam}
# {CORES}: core number used
# {input.unmap2genome.bam}: input unmap to genome bam file (.bam)
# {transcriptome.fa}: transcriptome fasta file
# {output.trxptomeAlign.bam}: output map to transcriptome sorted BAM file (.bam)


## 2.3 merge genome and transcriptome bam file
basalkit mergeBAM {input.trxptomeAlign.bam} {input.genomeAlign.bam} \
-o {output_prefix} \
{gtf} \
# {input.trxptomeAlign.bam}: input transcriptome bam file(.bam)
# {input.genomeAlign.bam}: input genome bam file(.bam)
# {output_prefix}: output merge bam file prefix (.bam)
# {gtf}: annotation gtf file


# 3 measure sites ratio
## 3.1 measure average modification with basalkit "avgmod"
basalkit avgmod {input.merge.sorted.bam} {genome.fa} \ 
-o {output_prefix} \
-M A:CGT -D M
# {input.merge.sorted.bam}: input merge bam file (.bam)
# {genome.fa}: genome fasta file
# {output_prefix}: output avgmod file prefix (_AvgMod.tsv)

## 3.2 measure p-value and FDR with basalkit "fdr"
basalkit fdr {input.treat_AvgMod.tsv.gz} -c {input.ctrl_AvgMod.tsv.gz} \
-o {output_prefix}
# {input.treat_AvgMod.tsv}: input treat AvgMod file
# {input.ctrl_AvgMod.tsv}: input ctrl AvgMod file
# {output_prefix}: output file prefix (_FDR.tsv.gz)
