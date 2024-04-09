# BAse-conversion Sequencing ALigner (BASAL)
Nucleotide modifications, encompassing both RNA and DNA modifications, play a pivotal role in gene transcription regulation. To pinpoint the modified sites, a variety of Base Conversion (BC) methods have been introduced. These BC methods can be categorized into three types: (1) one-way conversion, such as C-to-T for 5mC detection, or A-to-G for m6A detection; (2) multi-way conversion, such as converting A to C/G/T for m1A detection; and (3) deletion-induced conversion, such as U-to-deletion conversion for pseudouridine detection. These methods, particularly effective when coupled with sequencing, offer single-base resolution, surpassing immunoprecipitation-based techniques.

However, these evolving BC methods present significant data analysis challenges, with no single bioinformatic tool capable of handling the diversity of data produced. The primary hurdle is reads mapping, with two main strategies: the "mutation-rate approach" and the "conversion-sensitive approach."(**Fig 1**) The former often results in misalignment or the erroneous discard of reads due to treating converted bases as mismatches, while the latter, though more logical, lacks tools support the wide array of BC methods, especially those involving multi-way or deletion-induced conversions.
<div align=center><img src="https://github.com/JiejunShi/BASAL/blob/main/media/mapping_strategies.png" /></div>  

**Fig 1. Illustration of the differences between mutation-rate strategy and conversion-sensitive strategy in mapping BC sequencing reads using m6A sequencing (A is converted to G, while m6A is not converted) as an example.**

To address these challenges, we have introduced BASAL (BAse-conversion Sequencing ALigner), leveraging bitwise masking technology to support the analysis of diverse BC methods (**Fig 2**). BASAL has demonstrated superior performance in mapping accuracy and efficiency over existing tools, excelling at identifying reliable modification sites, and uncovering cell clusters and trajectories in single-cell epitranscriptomic data that align with biological functions. This breakthrough positions BASAL as a universal tool for analyzing various RNA and DNA modification detection technologies, facilitating groundbreaking discoveries in epigenomics and epitranscriptomics.
<div align=center><img src="https://github.com/JiejunShi/BASAL/blob/main/media/BASAL_algorithm.png" /></div>  

**Fig 2. The principle and compatibility of BASAL.**  
**(a)** Illustrations of different conversion modes used by base conversion (BC) sequencing techniques. Since RNA is reverse transcribed into DNA during sequencing, uridines are represented by thymidines (T). **(b)** The principle of the BASAL algorithm. **(c)** An example demonstrating how BASAL calculates the penalty score for one-way conversion data. **(d)** Using the same example as in (c) to explain how BASAL maps the multi-way conversion data. **(e)** Comparison of compatibility among aligners currently used in BC sequencing data mapping.

## Authors
- Jiejun Shi
- Moping Xu
- Miao Wang
## Dependencies
- Python3 with following packages
  - numpy; pandas; copy; collections; multiprocessing; pysam
- samtools (lastest version)
## Installation
BASAL is designed for linux64. In linux, simply `cd` to the `BASAL` directory and type `make` to make the executable binary.
No installation needed for BASALkit.
## Usage
### BASAL
	   ___    __    __    __    _
	  | |_)  / /\  ( (`  / /\  | |
	  |_|_) /_/--\ _)_) /_/--\ |_|__

	Usage:  basal [options]
	  Options for input/output files:
	       -a  <str>    input reads in FASTA/FASTQ/BAM format [Required option]
	       -b  <str>    input reads which is paired with -a, (default: none, single-end)
	       -d  <str>    reference sequences in FASTA format [Required option]
	       -o  <str>    output alignment in SAM/BAM format, if omitted, the output will be written to STDOUT in SAM format.

	  Options for base-conversion:
	       -M  <str>    the convert-from and convert-to base(s) seperated by ':' [Required option]
	                    the convert-from base must be single letter from [A,T,C,G],
	                    the convert-to base(s) can be single or multiple letters from [A,T,C,G,-], '-' represents deletion.
	                   examples:
	                    -M C:T, can detect C>T conversion(e.g. DNA bisulfite seq)
	                    -M A:G, can detect A>G conversion in RNA m6A seq(e.g. GLORI) or DNA 6mA seq(e.g. NT-seq)
	                    -M A:CGT, can detect RNA m6A in m6A-SAC-seq, which convert A to C/G/T
	                    -M T:-, can detect pseudouridine in BID-seq, which convert pseudouridine to deletion

	  Options for alignment:
	       -v  <float>  maximum percentage/number of mismatch bases in each read. (default: 0.1)
	                    The float value(between 0 and 1) is interpreted as the percentage of read length.
	                    The integer value is interpreted as absolute number of mismatches.
	                    The maximum mismatches will be reduced to 15 if it exceed 15.
	       -g  <int>    maximum size of gap (deletion/insertion), <=3 bp. default: 0
	       -w  <int>    maximum number of equal best hits to count, <=1000
	       -B  <int>    start from the Nth read or read pair, default: 1
	       -E  <int>    end at the Nth read or read pair, default: 4,294,967,295
	       -I  <int>    index interval (1~16), the reference genome will be indexed every Nbp, default: 4. Larger -I uses less memory.
	       -k  <float>  the cut-off ratio for over-represented kmers, default: 5e-07
	                    example: -k 1e-6 means the top 0.0001% over-represented kmer will be skipped in alignment
	       -s  <int>    seed size (8~16), default: 16.
	       -S  <int>    seed for random number generation used in selecting multiple hits
	                    set identical values to allow reproducible mapping results.
	                    (default: 0, get seed from system clock, mapping results not resproducible)
	       -p  <int>    number of processors to use, default: 1

	  Options for pair-end alignment:
	       -m  <int>    minimal insert size allowed, default: 28
	       -x  <int>    maximal insert size allowed, default: 1000

	  Options for reads trimming:
	       -q  <int>    quality threshold in trimming, 0-40, default: 0 (no trim)
	       -z  <int>    base quality, default: 33 [set 64 for Illumina, 33 for Sanger]
	       -f  <int>    reads containing more than this number of Ns will be skipped, default=5
	       -A  <str>    3' end adapter sequence to be trimmed, default: none (no trim)
	       -L  <int>    map the first N bases of the read, the max is 480 (default).

	  Options for mapping strand:
	       -n  [0,1,2]  -n 0: directional protocol, map single-end(SE) reads to forward strands, i.e. ++(same as OT in bismark) and -+(same as OB in bismark). For pair-end(PE), map read#1 to ++ and -+, map read#2 to +-(same as CTOT in bismark) and --(same as CTOB in bismark).
	                    -n 1: non-directional protocol, map reads to all 4 strands.
	                    -n 2: PBAT protocol, map SE reads to reverse strands, i.e. +- and --. For PE, map read#1 to +- and --, read#2 to ++ and -+.

	  Options for reporting:
	       -r  [0,1,2]  how to report repeat hits, 0=none(unique hit/pair); 1=random one; 2=all(slow), default:1.
	       -R           print corresponding reference sequences in SAM output, default: off
	       -u           report unmapped reads, default: off
	       -H           do not print header information in SAM format output
	       -V  [0,1,2]  verbose level: 0=no message displayed (quiet mode); 1=major message (default); 2=detailed message.
	       -h           help

### BASALkit
This is the downstream analysis toolkit for BASAL, designed to summarize sequence alignment results and identifying significantly converted sites. The executable script is `basalkit.py`. The other one, `basalkit_functions.py`, is not executable and contains the functions required by `basalkit.py`.

	$ python basalkit.py -h

	For help information of each function, try:

	  python basalkit.py <Function> -h

	Availible Functions:
		mergeBAM	Transfer the transcriptome BAM file to genome positions, and then merge it with the genome BAM file. This function is designed for RNA modification sequencing.
		avgmod	Calculate average modification level(AvgMod) of tested nucleotide(e.g. 5mC/6mA)
		fdr	Perform significance test between treatment and control/background, and report FDR for each sites
		regmod	Summarise the modification level of given regions
		shiftD	Shift the position of D in CIGAR in bam/sam. For deletion-induced techniques(e.g. BID-seq), if a deletion is detected in a polymer of convert-from bases, it is re-assigned to the rightmost base of the polymer.

#### BASALkit - mergeBAM
The `mergeBAM` module, specialized for RNA modification data, merges two-step alignment results and converts transcriptome alignment coordinates to genomic coordinates using the transcriptome annotation file. Reads spanning introns have their CIGAR values in the BAM file adjusted to reflect RNA splicing. CIGAR, an acronym for Concise Idiosyncratic Gapped Alignment Report, documents alignment discrepancies. The converted transcriptome alignments are then merged with genomic alignments into a single BAM file.

#### BASALkit - avgmod
This is the core module of BASALkit and calculates the conversion rate at potential modification sites, i.e., the frequency of the ‘convert-to’ base relative to sequencing depth. It's important to note that the conversion rate may not directly reflect the modification rate, depending on whether the converted base is modified or not (take note of the `-D` option): for methods like TAPS, the modification rate equals the conversion rate (use `-D M`); for BS-seq, it is 1 minus the conversion rate (use `-D U`).

#### BASALkit - fdr
Identifying credible modification sites often involves statistical comparisons between treatment and control groups, where the control could be untreated or unmodified samples. For RNA m6A detection, controls could be regular RNA-seq or FTO-treated samples. The `fdr` module, designed for this purpose, employs unilateral statistical tests (e.g., Fisher Test, Binomial Test) to assess modification significance, aiding in the selection of credible sites.

#### BASALkit - regmod
Beyond individual site analysis, DNA/RNA modification studies often assess the average modification level across genes or regions. The `regmod` module facilitates this, extracting convert-to base frequency and sequencing depth for all modification sites within a specified region, determining the region's average modification level as the ratio of their respective sums.
