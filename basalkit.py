#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, sys, time, argparse, array
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import OrderedDict
from basalkit_functions import *
import pysam

version="1.8"

ScriptName=re.split('/', sys.argv[0])[-1]
cmds=OrderedDict([("avgmod","Calculate average modification level(AvgMod) of tested nucleotide(e.g. 5mC/6mA)"),\
                  ("shiftD","Shift the position of D in CIGAR in bam/sam. For deletion-induced techniques(e.g. BID-seq), if a deletion is detected in a polymer of convert-from bases, it is re-assigned to the rightmost base of the polymer."),\
                  ("mergeBAM","Transfer the transcriptome BAM file to genome positions, and then merge it with the genome BAM file. This function is designed for RNA modification sequencing."),\
                  ("fdr","Perform significance test between treatment and control/background, and report FDR for each sites"),\
                  ("regmod","Summarise the modification level of given regions")
                  ])

def printHelp():
    print("BASAL Toolkit v"+version+"\n")
    print("For help information of each function, try:\n")
    print("  python "+ScriptName+" <Function> -h\n")
    print("Availible Functions:\n")
    for i in cmds.keys():
        print('  '+i+'\t'+cmds[i]+'\n')

def main(cmd=''):
    # check whether the requested function is included
    if cmd not in cmds.keys():
        print("Function ",cmd," not found, please see the help below:\n")
        printHelp()
        #sys.exit()
        return False

    # parse parameters
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage='python '+ScriptName+' '+cmd+' <positional arguments> [optional arguments]\n',\
                                     description=cmds[cmd]+"\n",\
                                     epilog='')
    # designate the first positional argument as ' '
    parser.add_argument(' ',default=None,help='')

    # parse input parameters of each function
    if cmd=="avgmod":
        parser.add_argument('Alignments',default=None,\
                            help="Aligned reads stored in SAM/BAM files. Multiple files generated by the same aligner are seperated by ','.")
        parser.add_argument('Reference',default=None,\
                            help="Reference genome sequence in fasta format.")
        parser.add_argument("-M", "--converted_base",dest="converted_base",default="C:T",\
                            help='The convert-from and convert-to base(s) seperated by ":". Use the same -M option for BASAL mapping. The convert-from base must be single letter from [A,T,C,G], the convert-to base(s) can be single/multiple letters from [A,T,C,G,-], "-" represents deletion. For example, -M C:T can detect C>T conversion in DNA bisulfite sequencing; -M A:G can detect A>G conversion in RNA m6A GLORI seq; -M T:- can detect pseudouridine in BID-seq.')
        parser.add_argument("-D", "--conversion_mode",dest="conversion_mode",default="U",\
                            help='Base-conversion mode, whether the modified(-D M) or unmodified(-D U, the default) nucleobase is converted. For example, "-M C:T -D U" can detect 5mC in bisulfite-seq or EM-seq; "-M C:T -D M" can detect 5mC in TAPS.')
        parser.add_argument("-T", "--molecule_type",dest="molecule_type",default="DNA",choices=['DNA', 'RNA'],\
                            help='The type of molecule which is tested.')
        parser.add_argument("-a", "--aligner",dest="aligner",default="BASAL",\
                            help='Because different aligners store mapping strand information in different ways, currently only the alignments of following aligners are supported: BASAL, Bismark and gemBS.')
        parser.add_argument('-o', '--tsv_prefix',dest="tsv_prefix",default='output',\
                            help="Prefix of output tsv files saving AvgMod of each tested nucleobase.")
        parser.add_argument('-w', '--wig_prefix',dest="wig_prefix",default=None,\
                            help="Prefix of output wiggle files saving AvgMod. As default, NO wiggle file is saved.")
        parser.add_argument('-b', '--wig_bin',dest="wig_bin",metavar='BIN',type=int,default=25,\
                            help="Bin size for wiggle files.")
        parser.add_argument('-c', '--chroms',dest="chroms",metavar='CHR',default=None,\
                            help="Only process given chromosomes, which separated by ','. For example, 'chr1,chr2'. As default, all chromosomes are processed.")
        parser.add_argument('-s', '--sam_path',dest="sam_path",metavar='PATH',default=None,\
                            help="Path to samtools.")
        parser.add_argument('-u', '--unique',action="store_true",dest="unique",default=False,\
                            help="Only process unique mappings/pairs.")
        parser.add_argument("-p", "--pair", action="store_true", dest="pair",default=False,\
                            help="Only process properly paired mappings.")
        parser.add_argument("-r", "--rm_dup", action="store_true", dest="rm_dup",default=False,\
                            help="Remove duplicated reads. This option may too stringent, it only keep one for reads with same left-most pos and strand info.")
        parser.add_argument("-t", "--trim_fillin", dest="trim_fillin", type=int, metavar='N', default=0,\
                            help="Trim N end-repairing fill-in nucleotides.")
        parser.add_argument("-g", "--combine", action="store_true", dest="combine", default=False,\
                            help="If set, reads from both Waston and Crick strands covering the same CpG are combined. This option is only valid if -M C:T.")
        parser.add_argument("-m", "--min_depth", dest="min_depth", type=int, metavar='FOLD', default=4,\
                            help="Report loci with sequencing depth>=FOLD.")
        parser.add_argument('-z', '--converted_site',dest="converted_site",type=float,default=0,\
                            help="Only use reads cover more than this number of converted nucleotides(e.g. 5mC/m6A). A value >=1 means the absolute number, while a value between 0 and 1 means the ratio to all convert-from bases covered by the read.")
        parser.add_argument("-i", "--handle_SNP", dest="handle_SNP", default="no-action",\
                            help='How to handle CT(for -M C:T) or AG(for -M A:G) SNP? Must be one of ["no-action","correct","skip"].')
        parser.add_argument("-x", "--context", dest="context", default=None,\
                            help='Select the sequence context of modified bases, only availabel for -M C:T, available options are "CG","CHG","CHH", multiple options are separated by ",". If -M is not C:T or no -x specified, all context are included and the sequence around the reported sites will be reported in "context" column of output, and the reported sequence length is controled by -y option.')
        parser.add_argument("-y", "--motif_length",dest="motif_length",type=int,default=5,\
                            help='The length of reported sequence around each reported site. As the reported site is in the middle of reported sequence, this option must be an odd number.')
    if cmd=="shiftD":
        parser.add_argument('Alignments',default=None,\
                            help="Aligned reads generated by BASAL. Note: the XR tag is required, so BASAL -R option must be used for reads mapping.")
        parser.add_argument("-c", "--convert_from_base",dest="convert_from_base",default="T",\
                            help='The convert-from base. It must be single letter from [A,T,C,G].')
        parser.add_argument('-s', '--sam_path',dest="sam_path",metavar='PATH',default=None,\
                            help="Path to samtools.")
        parser.add_argument('-o', '--out',dest="out",default='corrected',\
                            help="Prefix of output sam/bam file whose CIGAR has been corrected.")
    if cmd=="mergeBAM":
        parser.add_argument('transcriptomeAlignmentBam',default=None,\
                            help="Transcriptome BAM file generated by BASAL")
        parser.add_argument('genomeAlignmentBam',default=None,\
                            help="Genome BAM file generated by BASAL")
        parser.add_argument("gtf",default=None, help="gtf file contains the gene annotation")
        parser.add_argument('-o', '--output',dest="output",metavar='',default='sample1',\
                            help="The prefix of output merged bam file.")
        parser.add_argument("-u","--unlift", dest="unlift", default=False, action="store_true",
                            help="Whether to output unlifted transcriptomeAlignment into unlifted.bam")
    if cmd=="fdr":
        parser.add_argument('treat', default=None,\
                            help='AvgMod.tsv file of Treatment group.')
        parser.add_argument('-c','--ctrl', dest="ctrl", metavar='', default=None,\
                            help='AvgMod.tsv file of Ctrl group. If omitted, the overall conversion ratio(CR) of Treatment group will be used as background.')
        parser.add_argument("-m", "--min_depth", dest="min_depth", type=int, metavar='FOLD', default=4,\
                            help="Report loci with sequencing depth>=FOLD in both Treatment and Ctrl(if valid) group. This option is used before pval calculation.")
        parser.add_argument("-d","--method", dest="method", default="binomial", choices=['binomial', 'poisson','fisher'],\
                            help="statistical method: binomial, poisson, or fisher test")
        parser.add_argument("-r","--fdr_method", dest="fdr_method", default="fdr_bh",\
                            choices=["bonferroni","sidak","holm-sidak","holm","simes-hochberg","hommel","fdr_bh","fdr_by","fdr_tsbh","fdr_tsbky"],\
                            help="choose method for FDR correction: bonferroni,sidak,holm-sidak,holm,simes-hochberg,hommel,fdr_bh(the default),fdr_by,fdr_tsbh,fdr_tsbky. Check https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for details.")
        parser.add_argument('-o', '--output_prefix',dest="output_prefix",metavar='',default='output',\
                            help="prefix of output file name")
    if cmd=="regmod":
        parser.add_argument('Bed',default=None,\
                            help="Genomic regions in bed file with at least the first 3 columns")
        parser.add_argument('AvgMod',default=None,\
                            help="AvgMod file generated by 'avgmod' command")
        parser.add_argument('-s', '--usestrand',action="store_true",dest="usestrand",default=False,\
                            help="Use the strand information(6th column) in Bed file, then reads on forward and reverse strands will not be combined.")
        parser.add_argument('-o', '--output',dest="OUT",metavar='',default='region_ratio.tsv',\
                            help="Output file name")

    # print help information for the requested function
    if '-h' in sys.argv or '--help' in sys.argv:
        parser.print_help()
        return False
    else:
        if len(sys.argv)==2:
            print("Too few arguments. Try:\npython "+ScriptName+" "+cmd+" -h")
            return False

    # save all paramter values in args
    args=parser.parse_args()

    # run function
    if cmd=='avgmod':
        disp("avgmod Started")
        # check options
        converted_bases=args.converted_base.split(':')
        if len(converted_bases)!=2:parser.error('Invalid -M value')
        else:
            convert_from_base=converted_bases[0]
            convert_to_base=list(converted_bases[1])
            if len(convert_from_base)!=1 or convert_from_base not in ["A", "C", "G", "T"]:
                parser.error('Invalid -M value, the convert-from base must be single letter from [A,T,C,G]')
            for i in convert_to_base:
                if i not in ["A", "C", "G", "T", "-"]:
                    parser.error('Invalid -M value, the convert-to base(s) should be single/multiple letters from [A,T,C,G,-]')
        if args.conversion_mode not in ["U", "M"]:
            parser.error('Invalid -D value, select "U" or "M"')
        if args.aligner not in ["BASAL", "Bismark", "gemBS"]:
            parser.error('Invalid -a value, select "BASAL" or "Bismark" or "gemBS"')
        if args.sam_path!=None:
            if args.sam_path[-1] != '/': args.sam_path += '/'
        else:
            args.sam_path=""
        if args.chroms!=None:
            args.chroms = set(args.chroms.split(','))
        else:
            args.chroms = set()
        if args.combine==True:
            if args.converted_base!="C:T":
                parser.error('-g is only valid for -M C:T')
        if args.converted_site<0:
            parser.error('-z option should >= 0')
        handle_SNP_val = {"no-action": 0, "correct": 1, "skip": 2}
        try: args.handle_SNP = handle_SNP_val[args.handle_SNP.lower()]
        except KeyError: parser.error('Invalid -i value, select "no-action", "correct" or "skip"')
        if args.min_depth <= 0: parser.error('Invalid -m value, must >= 1')
        if args.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
        if args.converted_base=="C:T":
            seq_context_str = ['CG','CHG','CHH']
        else:
            seq_context_str = []
        if args.context!=None:
            args.context = set(args.context.upper().split(','))
            try: seq_context = set([seq_context_str.index(cc)+1 for cc in args.context])
            except ValueError: parser.error('Invalid -x value. For -M C:T, -x must be one or multiple in ["CG", "CHG", "CHH"].')
        else:
            seq_context=[]
        if args.motif_length<1 or args.motif_length % 2==0:
            parser.error('-y must be an odds number')
        else:
            motif_length=args.motif_length//2
        # Load Reference Genome
        ref=Load_Reference(ifile=args.Reference,chroms=args.chroms);
        args.chroms = set(ref.keys())
        # Create arrays
        coverage, meth0 = {}, {}
        for cr in ref:
            meth0[cr] = array.array('H', [0]) * len(ref[cr])
            if args.rm_dup==True:
                coverage[cr] = array.array('B', [0]) * len(ref[cr])
        depth=deepcopy(meth0)
        if args.handle_SNP > 0:meth_ct=deepcopy(meth0);depth_ct=deepcopy(meth0)
        else:meth_ct=[];depth_ct=[]
        # Mark Reference Genome
        refmark={}
        if seq_context!=[]:
            for cr in ref:refmark[cr] = array.array('b', [0]) * len(ref[cr])
            refmark=Mark_Reference(ref=ref,refmark=refmark,convert_from_base=convert_from_base,convert_to_base=convert_to_base)
        # Load Alignment
        args.Alignments=set(args.Alignments.split(','))
        meth0,depth,meth_ct,depth_ct,nmap=Load_Alignment(\
            ifiles=args.Alignments,convert_from_base=convert_from_base,convert_to_base=convert_to_base,conversion_mode=args.conversion_mode,molecule_type=args.molecule_type,aligner=args.aligner,ref=ref,refmark=refmark,coverage=coverage,\
            meth0=meth0,depth=depth,meth_ct=meth_ct,depth_ct=depth_ct,\
            sam_path=args.sam_path,unique=args.unique,pair=args.pair,rm_dup=args.rm_dup,\
            trim_fillin=args.trim_fillin,chroms=args.chroms,\
            seq_context=seq_context,handle_SNP=args.handle_SNP,converted_site=args.converted_site)
        # Combine cytosine methylation from both strands
        if args.combine==True:
            disp('Combining cytosine methylation from both strands')
            meth0=Combine_Methylation_Both_Strands(ref=ref,uncombined=meth0)
            depth=Combine_Methylation_Both_Strands(ref=ref,uncombined=depth)
            if args.handle_SNP > 0:
                meth_ct=Combine_Methylation_Both_Strands(ref=ref,uncombined=meth_ct)
                depth_ct=Combine_Methylation_Both_Strands(ref=ref,uncombined=depth_ct)
        # Output AvgMod files and wiggle file
        Out_base_ratio(tsv_prefix=args.tsv_prefix,wig_prefix=args.wig_prefix,wig_bin=args.wig_bin,\
            min_depth=args.min_depth,ref=ref,refmark=refmark,handle_SNP=args.handle_SNP,convert_from_base=convert_from_base,\
            seq_context_str=seq_context_str,seq_context=seq_context,motif_length=motif_length,meth0=meth0,depth=depth,meth_ct=meth_ct,depth_ct=depth_ct,nmap=nmap)
        # Delete big array to release memory
        del ref, refmark, coverage, depth, meth0, depth_ct, meth_ct
        disp("avgmod Finished")

    if cmd=="shiftD":
        disp("shiftD Started")
        # check options
        if len(args.convert_from_base)!=1 or args.convert_from_base not in ["A", "C", "G", "T"]:
            parser.error('Invalid -c value, the convert-from base must be single letter from [A,T,C,G]')
        if args.sam_path!=None:
            if args.sam_path[-1] != '/': args.sam_path += '/'
        else:
            args.sam_path=""
        shiftD(Alignments=args.Alignments,output=args.out,convert_from_base=args.convert_from_base,sam_path=args.sam_path)
        tmp = os.system("{}samtools view -o {} {}".format(args.sam_path,args.out+".bam",args.out+".sam"))
        if tmp!=0:
            disp('shiftD Finished. Failed to write bam file. Corrected alignments are still in sam file: {}.'.format(args.out+".sam"))
        else:
            tmp = os.system("rm {}".format(args.out+".sam"))
            disp('shiftD Finished. Corrected alignments are saved in {}'.format(args.out+".bam"))

    if cmd=='mergeBAM':
        disp("mergeBAM Started")
        gtf = read_gtf(args.gtf)
        disp("gtf loaded")

        in_mode = 'rb'
        out_mode = 'wb'

        disp("Buidling genome header")
        header_trans = {}
        header_trans['HD'] = {'VN': '1.0'}
        header_trans['SQ'] = []

        #genomeAlignmentBam="test2_map2genome.bam"
        if args.genomeAlignmentBam:
            with pysam.AlignmentFile(args.genomeAlignmentBam, 'rb') as INPUT:
                for header in INPUT.header["SQ"]:
                    #print(header["SN"])
                    header_trans['SQ'].append({"SN": header["SN"], "LN": header["LN"]})
        else:
            raise Warning("Please provide genomeAlignmentBam")
        outputprefix=args.output
        file_outBAM=outputprefix + ".t2g.bam"
        if args.unlift == True:
            unlift_fn = args.transcriptomeAlignmentBam.replace(".bam", ".unlift.bam")
            with pysam.AlignmentFile(args.transcriptomeAlignmentBam, in_mode) as INPUT, pysam.AlignmentFile(file_outBAM, out_mode,header=header_trans) as OUTPUT,\
                    pysam.AlignmentFile(unlift_fn, 'wb', template=INPUT) as UNLIFT:
                header_dict = {}
                n_header = 0
                for headerr in OUTPUT.header["SQ"]:
                    header_dict[headerr['SN']] = n_header
                    n_header += 1
                for segment in INPUT:
                    segment_output = map_to_genome(header_dict,gtf,segment,args.unlift,UNLIFT)
                    if segment_output:
                        OUTPUT.write(segment_output)
        else:
            with pysam.AlignmentFile(args.transcriptomeAlignmentBam, in_mode) as INPUT, pysam.AlignmentFile(file_outBAM, out_mode,header=header_trans) as OUTPUT:
                header_dict = {}
                n_header = 0
                UNLIFT=[]
                for headerr in OUTPUT.header["SQ"]:
                    header_dict[headerr['SN']] = n_header
                    n_header += 1
                for segment in INPUT:
                    segment_output = map_to_genome(header_dict,gtf,segment,args.unlift,UNLIFT)
                    if segment_output:
                        OUTPUT.write(segment_output)

        disp("merging transcriptome and genome BAM")
        hid = 0
        hid_dict = {}
        lift_over = {}
        new_header = {}
        new_header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}
        new_header['SQ'] = []

        file_sortBAM=outputprefix + ".t2g.bam"
        file_genomeBAM=args.genomeAlignmentBam
        file_mergeBAM=outputprefix + ".merge.bam"
        current_path = os.getcwd()
        input2=[os.path.join(current_path, file_sortBAM),os.path.join(current_path, file_genomeBAM)]
        for fn in input2:
            hid,new_header,hid_dict,lift_over = read_headers(fn,hid,new_header,hid_dict,lift_over)

        with pysam.AlignmentFile(file_mergeBAM, 'wb', header = new_header) as OUTPUT:
            for fn in input2:
                merge_bam(lift_over,fn,OUTPUT)
        pysam.sort("-o", file_mergeBAM.replace(".bam",".sorted.bam"), file_mergeBAM)
        pysam.index(file_mergeBAM.replace(".bam",".sorted.bam"))
        disp("mergeBAM Finished")

    if cmd=='fdr':
        disp("fdr Started")
        calc_pval(args.treat,args.ctrl, output_prefix=args.output_prefix,min_depth=args.min_depth, method=args.method, fdr_method=args.fdr_method)
        disp("fdr Finished")

    if cmd=='regmod':
        disp("regmod Started")
        Ratio_df=read_methy_files(ifile=args.AvgMod, cols=[0,1,2,6,7])
        o1=open(args.OUT,'w')
        if args.usestrand==True:
            Bed=pd.read_csv(args.Bed, sep="\t",usecols=[0,1,2,5],header=None)
            Bed.columns=["chr","start","end","strand"]
            Bed.sort_values(['chr','strand','start','end'], inplace=True, ascending=True)
            disp("Generating AvgMod ratio for {} Regions ...".format(Bed.shape[0]))
            o1.write('\t'.join(["chr","start","end","strand","AvgMod","site","coverage"])+'\n')
            chr0="";strand0="";
            for row in Bed.iterrows():
                chr1=row[1]['chr'];
                start1=int(row[1]['start']);
                end1=int(row[1]['end']);
                strand1=row[1]['strand']
                if chr1 != chr0 or strand1!=strand0:
                    Ratio_sub=Ratio_df[((Ratio_df['chr']==chr1) & (Ratio_df['strand']==strand1))]
                ratio_list=Region_weighted_Ratio(ratio_sub=Ratio_sub,start=start1,end=end1)
                aline=[chr1,start1,end1,strand1];aline.extend(ratio_list);
                o1.write('\t'.join(map(str,aline))+'\n')
                chr0=chr1;strand0=strand1;
        else:
            Bed=pd.read_csv(args.Bed, sep="\t",usecols=[0,1,2],header=None)
            Bed.columns=["chr","start","end"]
            Bed.sort_values(['chr','start','end'], inplace=True, ascending=True)
            disp("Generating AvgMod ratio for {} regions ...".format(Bed.shape[0]))
            o1.write('\t'.join(["chr","start","end","AvgMod","site","coverage"])+'\n')
            chr0="";
            for row in Bed.iterrows():
                chr1=row[1]['chr'];
                start1=row[1]['start'];
                end1=row[1]['end'];
                if chr1 != chr0:
                    Ratio_sub=Ratio_df[Ratio_df['chr']==chr1]
                ratio_list=Region_weighted_Ratio(ratio_sub=Ratio_sub,start=start1,end=end1)
                aline=[chr1,start1,end1];aline.extend(ratio_list);
                o1.write('\t'.join(map(str,aline))+'\n')
                chr0=chr1;
        o1.close()
        disp("regmod Finished")

if len(sys.argv)>1:
    main(cmd=sys.argv[1])
else:
    printHelp()
