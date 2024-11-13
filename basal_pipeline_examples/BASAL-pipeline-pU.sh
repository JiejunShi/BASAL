#!/bin/bash
#===for Ψ=====

# 1 Prepare your reads after trimming
# 2 Mapping reads
## 2.1 Mapping reads to genome
basal -p {CORES} \
-a {input.processed_Read} -d {genome.fa} \
-o {output.map2genome.bam} \
-M T:- -n 1 -g 3 -R -S 2024 -u
# {CORES}: core number used
# {input.processed_Read}: pre-processed read file (.fq.gz)
# {genome.fa}: genome fasta file
# {output.map2genome.bam}: output bam file (.bam)

# reads that were successfully matched to the genome were extracted
samtools view -b -F 3588 -@ {CORES} -o {output.tmp.bam} {input.map2genome.bam}
sambamba sort -m 8GB -t {CORES} -o {output.genomeAlign.bam} {output.tmp.bam}
# {CORES}: core number used
# {input.map2genome.bam}: input map to genome bam file with BASAL (.bam)
# {output.tmp.bam}: temperory file
# {output.genomeAlign.bam}: output map to genome sorted BAM file (.bam)

# reads that were not successfully matched to the genome were extracted
samtools view -b --include-flags 4 -@{CORES} -o {output.unmap2genome.bam} {input.map2genome.bam}
samtools fastq {output.unmap2genome.bam} > {output.unmap2genome.fq}
gzip {output.unmap2genome.fq}
# {CORES}: core number used
# {input.map2genome.bam}: input map to genome bam file with BASAL (.bam)
# {output.unmap2genome.bam}: output unmap to genome bam file (.bam)
# {output.unmap2genome.fq}: output unmap to genome fastq file (.fq)

## 2.2 Mapping reads to transcriptome
basal -p {CORES} \
-a {input.unmap2genome.fq.gz} -d {transcriptome.fa} \
-o {ouput.tmp.bam} \
-M T:- -n 1 -g 3 -R -S 2024;
# {CORES}: core number used
# {input.unmap2genome.fq.gz}: input unmap to genome fastq file (.fq.gz)
# {transcriptome.fa}: transcriptome fasta file
# {ouput.tmp.bam}: temperory file

# sort reads mapping to transcriptome
sambamba sort -m 8GB -t {CORES} -o {output.trxptomeAlign.bam} {input.tmp.bam};
# {CORES}: core number used
# {output.trxptomeAlign.bam}: output map to transcriptome sorted BAM file (.bam)
# {input.tmp.bam}: input temperory file from BASAL map to transcriptome result

## 2.3 correct CIGAR with basalkit "shiftD": for this module parameters '-R' is necessary in basal
basalkit shiftD {input.genomeAlign.bam} -o {output.tmp}
sambamba sort -m 8GB -t {CORES} -o {output.genomeAlign.corrected.bam} {input.tmp}
basalkit shiftD input.trxptomeAlign.bam -o {output.tmp}
sambamba sort -m 8GB -t {CORES} -o {output.trxptomeAlign.corrected.bam} {input.tmp}
# {CORES}: core number used
# {output.tmp}: temperory file
# {output.genomeAlign.corrected.bam}: output corrected genome bam file(.bam)
# {output.trxptomeAlign.corrected.bam}: output corrected transcriptome bam file(.bam)

## 2.4 merge genome and transcriptome bam file
basalkit mergeBAM {input.trxptomeAlign.corrected.bam} {input.genomeAlign.corrected.bam} \
{gtf} \
-o {output_prefix}
# {input.trxptomeAlign.corrected.bam}: input transcriptome bam file(.bam)
# {input.genomeAlign.corrected.bam}: input genome bam file(.bam)
# {output_prefix}: output merge bam file prefix (.bam)
# {gtf}: annotation gtf file


# 3 measure sites ratio
## 3.1 measure average modification with basalkit "avgmod"
basalkit avgmod {input.merge.sorted.bam} {genome.fa} \ 
-o {output_prefix} \ 
-M T:- -D M -T RNA -y 7
# {input.merge.sorted.bam}: input merge bam file (.bam)
# {genome.fa}: genome fasta file
# {output_prefix}: output avgmod file prefix (_AvgMod.tsv)

## 3.2 measure p-value and FDR with basalkit "fdr"
basalkit fdr {input.treat_AvgMod.tsv.gz} -c {input.ctrl_AvgMod.tsv.gz} \
-o {output_prefix};
# {input.treat_AvgMod.tsv.gz}: input treat AvgMod file
# {input.ctrl_AvgMod.tsv.gz}: input ctrl AvgMod file
# {output_prefix}: output file prefix (_FDR.tsv.gz)