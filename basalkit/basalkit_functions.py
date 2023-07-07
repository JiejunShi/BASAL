#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import multiprocessing as mp
import os,re,sys,time
from collections import Counter,defaultdict,OrderedDict
import numpy as np
import pandas as pd
import math
import scipy.stats
from statsmodels.stats.multitest import multipletests
from copy import copy

def disp(text):
    print('[BASALKIT @{}] \t{}'.format(time.asctime(), text), file=sys.stderr)

def Load_Reference(ifile,chroms):
    ref, cr, seq = {}, '', ''
    disp('Loading Reference Genome: {}'.format(ifile))
    for line in open(ifile):
        if line[0] == '>':
            if len(cr) > 0:
                if len(chroms)==0 or cr in chroms: ref[cr] = seq.upper()
            cr, seq = line[1:-1].split()[0], ''
        else: seq += line.strip()
    if len(chroms) == 0 or cr in chroms: ref[cr] = seq.upper()
    del seq
    return ref

def Mark_Reference(ref,refmark,convert_from_base,convert_to_base):
    disp('Marking Reference Genome')
    if convert_from_base=="C" and convert_to_base==["T"]:
        CG, CHG, CHH = 1, 2, 3
        for cr in ref:
            refcr, refmarkcr = ref[cr], refmark[cr]
            index = refcr.find('C', 0, len(refcr)-2)
            while index >= 0:
                if refcr[index+1] == 'G': refmarkcr[index] = CG
                elif refcr[index+2] == 'G': refmarkcr[index] = CHG
                else: refmarkcr[index] = CHH
                index = refcr.find('C', index+1, len(refcr)-2)
            index = refcr.find('G', 2, len(refcr))
            while index >= 0:
                if refcr[index-1] == 'C': refmarkcr[index] = CG
                elif refcr[index-2] == 'C': refmarkcr[index] = CHG
                else: refmarkcr[index] = CHH
                index = refcr.find('G', index+1, len(refcr))
    return refmark

def parseCigar(seq, cigar):
    cigarRE = re.compile(r'\d+[a-zA-Z]')
    index = 0
    #originalLen = len(seq)
    cigarMatch = cigarRE.findall(cigar)
    for align in cigarMatch:
        length = int(align[:-1])
        op = align[-1]
        if op == 'M' or op == '=' or op == 'X':
            index += length
        elif op == 'I' or op == 'S':
            seq = seq[:index]+seq[index+length:]
        elif op == 'D' or op == 'N':
            seq = seq[:index]+'-'*length+seq[index:]
            index += length
        elif op == 'H' or op == 'P':
            continue
        #else:
            #raise ValueError("%c not a valid CIGAR character"%(op))
    #assert originalLen == index, "String length does not match CIGAR"
    return seq

def strand_bismark2bsmap(XR,XG):
    if XR == "CT" and XG == "CT": return "++"
    elif XR == "CT" and XG == "GA": return "-+"
    elif XR == "GA" and XG == "CT": return "+-"
    elif XR == "GA" and XG == "GA": return "--"
    else: return ""

def strand_gemBS2bsmap(XB,sam_flag):
    # read forward strand
    if (32&sam_flag==32):
        # XB:A:C means read from C2T strand, i.e. original strand
        if XB == "C": return "++"
        # XB:A:G means read from G2A strand, i.e. complement to original strand
        elif XB == "G": return "--"
    # read reverse strand
    elif (16&sam_flag==16):
        # XB:A:C means read from C2T strand, i.e. original strand
        if XB == "C": return "+-"
        # XB:A:G means read from G2A strand, i.e. complement to original strand
        elif XB == "G": return "-+"
    else: return ""

def Load_One_Read(line,ref,coverage,sam_format,molecule_type,aligner,unique,pair,rm_dup,trim_fillin,chroms,ZSselect):
    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = int(col[1])
        # unmapped 0x4=4; not uniq 0x100=256; proper paired 0x2=2;
        if (4&flag==4): return []# unmapped
        if unique and (256&flag==256): return []# not unique mapped
        if pair and (2&flag!=2): return []# not proper paired

        cr, pos, cigar, seq, strand, insert, pos_mate = col[2], int(col[3])-1, col[5], col[9], "", int(col[8]), int(col[7])-1
        # 1-based (sam) -> 0-based
        if cr not in chroms: return []
        seq = parseCigar(seq, cigar)
        if seq=="": return []

        # strand information
        if aligner == "BASAL":
            ZS_index = line.find('ZS:Z:')
            strand = line[ZS_index+5:ZS_index+7]
        elif aligner == "Bismark":
            XR_index = line.find('XR:Z:');XG_index = line.find('XG:Z:');
            XR_tag=line[XR_index+5:XR_index+7]
            XG_tag=line[XG_index+5:XG_index+7]
            strand=strand_bismark2bsmap(XR=XR_tag,XG=XG_tag)
        elif aligner == "gemBS":
            XB_index = line.find('XB:A:')
            XB_tag = line[XB_index+5:XB_index+6]
            strand=strand_gemBS2bsmap(XB=XB_tag,sam_flag=flag)
    else:
        #BSP format
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if unique and flag != 'UM': return []
        if pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference (deletion?)
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]

    if strand not in ZSselect:return []
    pos2 = pos + len(seq)
    if pos2 >= len(ref[cr]): return []
    if strand == '+-' or strand == '-+': frag_end, direction = pos2, 2
    elif strand == '++' or strand == '--': frag_end, direction = pos, 1
    if rm_dup:  # remove duplicate hits
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[trim_fillin:], pos+trim_fillin
    # remove overlapped regions in paired hits, SAM/BAM format only
    # if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # by yxi
    # pair on +
    #         |------->
    # pair on -, 6 senarios
    # 1: <--|                   # keep pair on +/- (reported as unpaired in BASAL)
    # 2: <-------|              # trim pair on + (left part), keep pair on -
    # 3: <--------------|       # delete pair on +, keep pair on -
    # 4:        <--|            # delete pair on - (not feasible now), keep pair on +
    # 5:        <-------|       # trim pair on +(right part), keep pair on -
    # 6:               <--|     # keep pair on +/-
    if sam_format:
        if strand == '++' or strand == '--':
            if insert > 0: # pair on + has smaller 5', senario 2/3/4/5/6
                if pos < pos_mate:# senario 4/5/6
                    if insert > len(seq):# senario 5/6
                        seq = seq[:pos_mate-pos]
                else:# senario 2/3
                    if insert >= len(seq):# senario 3
                        return []
                    else:# senario 2
                        seq = seq[insert:]
                        pos = pos+insert
    # unable to calculate seq_mate from seq, cannot delete pair on - for senario 4, only feasible with name-sorted bam
    #    elif strand == '+-' or strand == '-+':
    #        if insert < 0: # pair on - has bigger 5', senario 2/3/4/5/6
    #            if pos > pos_mate:# senario 4/5/6
    #                if abs(insert) > len(seq_mate):# senario 4
    #                    return []
    if molecule_type=="DNA":
        return (seq, strand[0], cr, pos)
    else:
        if direction==1:
            return (seq, "+", cr, pos)
        elif direction==2:
            return (seq, "-", cr, pos)

def Load_Alignment(ifiles,convert_from_base,convert_to_base,conversion_mode,molecule_type,aligner,camda,ref,refmark,coverage,meth0,meth1,depth,meth_ct,depth_ct,sam_path,unique,pair,rm_dup,trim_fillin,chroms,seq_context,handle_SNP,converted_site,mapping_strand):
    pipes = []
    for ifile in ifiles:
        if ifile[-4:].upper() == '.SAM': sam_format, fin = True, open(ifile)
        elif ifile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view %s' % (sam_path, ifile))
        elif ifile[-5:].upper() == '.CRAM': sam_format, fin = True, os.popen('%ssamtools view %s' % (sam_path, ifile))
        else: sam_format, fin = False, open(ifile)
        pipes.append((ifile,sam_format,fin))
    # conversion_rule: raw, modified, unmodified, raw_reverse_compl, modified_reverse_compl, unmodified_reverse_compl
    complement={"A":"T","C":"G","G":"C","T":"A","-":"-"}
    convert_to_base_cp=[complement[i] for i in convert_to_base]
    if conversion_mode == "U":# BS mode
        conversion_rule = {'+': (convert_from_base,[convert_from_base],convert_to_base,complement[convert_from_base],[complement[convert_from_base]],convert_to_base_cp), '-': (complement[convert_from_base],[complement[convert_from_base]],convert_to_base_cp,convert_from_base,[convert_from_base],convert_to_base)}
    elif conversion_mode == "M":# TAPS mode
        conversion_rule = {'+': (convert_from_base,convert_to_base,[convert_from_base],complement[convert_from_base],convert_to_base_cp,[complement[convert_from_base]]), '-': (complement[convert_from_base],convert_to_base_cp,[complement[convert_from_base]],convert_from_base,convert_to_base,[convert_from_base])}
    if mapping_strand==0:
        ZSselect=["++","-+"]
    elif mapping_strand==1:
        ZSselect=["++","-+","+-","--"]
    else:
        ZSselect=["+-","--"]

    nmap = 0
    for ifile, sam_format, fin in pipes:
        disp('Load Alignment: {}'.format(ifile))
        nline = 0
        for line in fin:
            nline += 1
            map_info = Load_One_Read(line,ref,coverage,sam_format,molecule_type,aligner,unique,pair,rm_dup,trim_fillin,chroms,ZSselect)
            if len(map_info) == 0: continue
            seq, strand, cr, pos = map_info
            pos2 = pos + len(seq)
            nmap += 1
            raw, modified, unmodified, raw_rc, modified_rc, unmodified_rc = conversion_rule[strand]
            refseq, refmarkcr = ref[cr], {}
            if refmark!={}:refmarkcr=refmark[cr]
            indexes=[];n_covered_1read=0;n_converted_1read=0;
            if camda==True:read_modified=False
            for index in re.finditer(raw,refseq[pos:pos2]):
                n_covered_1read += 1
                index0=index.span()[0]
                indexes.append(index0)
                if seq[index0] in unmodified and conversion_mode == "U":n_converted_1read += 1
                if seq[index0] in modified:
                    if conversion_mode == "M":n_converted_1read += 1
                    if camda==True and read_modified==False:
                        if refmarkcr=={} or (refmarkcr[index0+pos] in seq_context):read_modified=True
            if converted_site>=1:
                if n_converted_1read < converted_site:continue
            else:
                if n_converted_1read < converted_site*n_covered_1read:continue
            if n_covered_1read>0:
                depth_cr = depth[cr];
                meth0_cr = meth0[cr];
                if camda==True:meth1_cr = meth1[cr];
                for index in indexes:
                    if (refmarkcr=={} or (refmarkcr[index+pos] in seq_context)) and depth_cr[index+pos] < (2**16-1):
                        if seq[index] in unmodified:
                            depth_cr[index+pos] += 1
                            if camda==True:
                                if read_modified==True:
                                    meth1_cr[index+pos] += 1
                        elif seq[index] in modified:
                            depth_cr[index+pos] += 1
                            meth0_cr[index+pos] += 1
                            if camda==True:
                                meth1_cr[index+pos] += 1
            if handle_SNP == 0: continue
            # use G/GA on reverse strand to adjust eff_CT_counts = CT_counts * (rev_G_counts / rev_GA_counts)
            indexes=[]
            for index in re.finditer(raw_rc,refseq[pos:pos2]):
                index0=index.span()[0]
                indexes.append(index0)
            if len(indexes)>0:
                depth_ct_cr = depth_ct[cr];
                meth_ct_cr = meth_ct[cr];
                for index in indexes:
                    if (refmarkcr=={} or (refmarkcr[index+pos] in seq_context)) and depth_ct_cr[index+pos] < (2**16-1):
                        if seq[index] in unmodified_rc:
                            depth_ct_cr[index+pos] += 1
                        elif seq[index] in modified_rc:
                            depth_ct_cr[index+pos] += 1
                            meth_ct_cr[index+pos] += 1
        fin.close()
        disp('Read {} lines'.format(nline))
    return meth0,meth1,depth,meth_ct,depth_ct,nmap

def correctCigar(cigar,XR,convert_from_base):
    cigarRE = re.compile(r'\d+[a-zA-Z]')
    cigar_old=[]
    for align in cigarRE.findall(cigar):
        if align[-1]=="D":Dsize=int(align[:-1])
        cigar_old.append(int(align[:-1]));cigar_old.append(align[-1])

    cigar_new=copy(cigar_old)
    D1 = cigar_old[0]+cigar_old[2]-1
    if XR[D1]==convert_from_base:
        shift=0
        for j in range(1,len(XR)-D1):
            if XR[D1+j]==convert_from_base:shift=j
            else:break
        if shift>0:
            if Dsize==1 or \
            (Dsize==2 and XR[D1-1]==convert_from_base) or \
            (Dsize==3 and XR[D1-1]==convert_from_base and XR[D1-2]==convert_from_base):#deletion=T/TT/TTT
                cigar_new[0]=cigar_old[0]+shift;cigar_new[4]=cigar_old[4]-shift
            else:
                if (Dsize==2 and XR[D1-1]!=convert_from_base) or \
                (Dsize==3 and XR[D1-1]==convert_from_base and XR[D1-2]!=convert_from_base):#deletion=XT/XTT (X means not T)
                    cigar_new[2]=1
                elif Dsize==3 and XR[D1-1]!=convert_from_base:#deletion=TXT/XXT
                    cigar_new[2]=2
                cigar_new.insert(4,shift);cigar_new.insert(5,'M')
                cigar_new.insert(6,Dsize-cigar_new[2]);cigar_new.insert(7,'D')
                cigar_new[8]=cigar_old[4]-shift
    if cigar_new[-2]<=0:cigar_new.pop();cigar_new.pop()
    cigar_new="".join(map(str,cigar_new))
    return cigar_new

def shiftCIGAR(Alignments,output,convert_from_base,sam_path):
    # load Alignment
    if Alignments[-4:].upper() == '.SAM':
        fin=open(Alignments)
    elif Alignments[-4:].upper() == '.BAM':
        fin=os.popen('%ssamtools view -h %s' % (sam_path, Alignments))
    else:
        parser.error('Input Alignment must be sam or bam file')
    fout = open(output+".sam", 'w')
    disp('Load Alignments in: {}'.format(Alignments))
    for line in fin:
        if line[0] == '@':fout.write(line);continue
        col = line.split('\t')
        cigar=col[5]
        if 'D' not in cigar:fout.write(line);continue
        XR=col[12][7:-2]#ref
        ZS=col[13][5:7]
        which_base=convert_from_base
        if ZS=="+-" or ZS=="-+":which_base=reverse_complement(convert_from_base)
        cigar_new=correctCigar(cigar=cigar,XR=XR,convert_from_base=which_base)
        col[5]=cigar_new
        fout.write("\t".join(col))
    fout.close()
    fin.close()

def Combine_Methylation_Both_Strands(ref,uncombined):
    for cr in uncombined:
        uncombined_cr, ref_cr = uncombined[cr], ref[cr]
        pos = ref_cr.find('CG')
        while pos >= 0:
            try:
                uncombined_cr[pos] += uncombined_cr[pos+1]
            except OverflowError:
                uncombined_cr[pos] = int((uncombined_cr[pos] + uncombined_cr[pos+1]) / 2)
            uncombined_cr[pos+1] = 0
            pos = ref_cr.find('CG', pos+2)
    return uncombined

def reverse_complement(seq):
    complement={"A":"T","C":"G","G":"C","T":"A","N":"N"}
    seq=seq.upper()
    seq_rc=[complement[i] for i in seq]
    seq_rc.reverse()
    seq_rc="".join(seq_rc)
    return(seq_rc)

def Out_base_ratio(tsv_prefix,wig_prefix,wig_bin,camda,min_depth,ref,refmark,handle_SNP,convert_from_base,seq_context_str,seq_context,motif_length,meth0,meth1,depth,meth_ct,depth_ct,nmap):

    header=['chr','pos','strand','context','ratio','eff_coverage','N_mod','N_total']
    if handle_SNP > 0: header.extend(['N_mod_rev','N_total_rev'])

    fo_mr = open(tsv_prefix+"_AvgMod.tsv", 'w')
    fo_mr.write('\t'.join(header)+'\n')
    if camda==True:
        fo_camda = open(tsv_prefix+"_CAMDA.tsv", 'w')
        fo_camda.write('\t'.join(header)+'\n')
    if wig_prefix!=None:
        fo_mr_wig = open(wig_prefix+"_AvgMod.wig", 'w');fo_mr_wig.write('track type=wiggle_0 name='+wig_prefix+'_MethRatio\n')
        if camda==True:
            fo_camda_wig = open(wig_prefix+"_CAMDA.wig", 'w');fo_camda_wig.write('track type=wiggle_0 name='+wig_prefix+'_CAMDA\n')
        disp('Output ratios in tsv files and wiggle files')
    else:
        disp('Output ratios in tsv files')

    #z95, z95sq = 1.96, 1.96 * 1.96
    nc, nd, dep0 = 0, 0, min_depth
    for cr in sorted(depth.keys()):
        depth_cr, meth0_cr, refcr, refmarkcr = depth[cr], meth0[cr], ref[cr], {}
        if refmark!={}:refmarkcr=refmark[cr]
        if camda==True: meth1_cr=meth1[cr]
        if handle_SNP > 0: depth_ct_cr, meth_ct_cr = depth_ct[cr], meth_ct[cr]
        if wig_prefix!=None:
            fo_mr_wig.write('variableStep chrom={} span={}\n'.format(cr, wig_bin))
            if camda==True:fo_camda_wig.write('variableStep chrom={} span={}\n'.format(cr, wig_bin))
            bin = wigd = wigm0 = wigm1 = 0
        for i, dd in enumerate(depth_cr):
            if dd < dep0: continue
            if handle_SNP > 0:
                m1, d1 = meth_ct_cr[i], depth_ct_cr[i]
                if m1 != d1:
                    if handle_SNP == 2: continue
                    d = float(dd) * m1 / d1
                else: d = dd#float(dd)
            else: d = dd#float(dd)

            if seq_context != []:
                if refmarkcr[i] not in seq_context: continue
                else:seq = seq_context_str[refmarkcr[i]-1]
            else:
                if refcr[i] == convert_from_base:seq=refcr[i-motif_length:i+motif_length+1]
                else:seq=reverse_complement(refcr[i-motif_length:i+motif_length+1])
            if refcr[i] == convert_from_base:strand = '+'
            else:strand = '-'

            m_mr = meth0_cr[i];
            if camda==True:m_full=meth1_cr[i];m_camda=m_full-m_mr;
            try: ratio_mr = min(m_mr, d)*1.0 / d
            except ZeroDivisionError: continue
            if camda==True:
                try: ratio_full = min(m_full, d)*1.0 / d
                except ZeroDivisionError: continue
                ratio_camda = ratio_full - ratio_mr
            nc += 1
            nd += d
            if wig_prefix!=None:
                if i // wig_bin != bin:
                    if wigd > 0:
                        fo_mr_wig.write('{:.0f}\t{:.3f}\n'.format(bin*wig_bin+1, min(wigm0/wigd,1)))# wiggle is 1-based
                        if camda==True:fo_camda_wig.write('{:.0f}\t{:.3f}\n'.format(bin*wig_bin+1, min((wigm1-wigm0)/wigd,1)))# wiggle is 1-based
                    bin = i // wig_bin#use integer division
                    wigd = wigm0 = 0.0
                    if camda==True: wigm1 = 0.0
                wigd += d
                wigm0 += m_mr
                if camda==True: wigm1 += m_full
            #denorminator = 1 + z95sq / d
            #pmid_mr = ratio_mr + z95sq / (2 * d)
            #sd_mr = z95 * ((ratio_mr*(1-ratio_mr)/d + z95sq/(4*d*d)) ** 0.5)
            #CIl_mr, CIu_mr = (pmid_mr - sd_mr) / denorminator, (pmid_mr + sd_mr) / denorminator
            #pmid_camda = ratio_camda + z95sq / (2 * d)
            #sd_camda = z95 * ((ratio_camda*(1-ratio_camda)/d + z95sq/(4*d*d)) ** 0.5)
            #CTl_camda, CTu_camda = (pmid_camda - sd_camda) / denorminator, (pmid_camda + sd_camda) / denorminator
            if handle_SNP > 0:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd, m1, d1))# 1-based
            else:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd))# 1-based
            if camda==True:
                if handle_SNP > 0:
                    fo_camda.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_camda, d, m_camda, dd, m1, d1))# 1-based
                else:
                    fo_camda.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_camda, d, m_camda, dd))# 1-based
    fo_mr.close();
    if camda==True:fo_camda.close();
    if wig_prefix!=None:
        fo_mr_wig.close();
        if camda==True:fo_camda_wig.close();
    if nc == 0:
        fold="NA"
    else:
        fold=round(float(nd)/nc,2)
    disp('Total {} valid mappings, cover {} base {}, average depth: {} fold.'.format(nmap, nc, convert_from_base, fold))

def read_gtf(fn):
    fin=open(fn,'r')
    output = defaultdict(dict)
    for line in fin:
        if line.startswith("#"):continue
        line = line.strip().split("\t")
        #if line[2] not in ["exon", "UTR"]:continue
        if line[2] != "exon":continue
        chr=line[0];pleft=int(line[3]);pright=int(line[4]);strand=line[6]
        for field in line[8].split("; "):
            if field.startswith("transcript_id"):
                trans_id=field.replace("transcript_id ","")
                trans_id=trans_id.replace("\"","")
                break
        if trans_id not in output:
            output[trans_id]['strand'] = strand
            output[trans_id]['chr'] = chr
            output[trans_id]['exonStarts'] = []
            output[trans_id]['exonEnds'] = []
        if strand=="+":
            output[trans_id]['exonStarts'].append(pleft)
            output[trans_id]['exonEnds'].append(pright)
        else:
            output[trans_id]['exonStarts'].append(pright)
            output[trans_id]['exonEnds'].append(pleft)
    fin.close()

    # 1-based, fully-closed
    for trans_id in output.keys():
        start_t=1
        exons=OrderedDict()
        if output[trans_id]['strand'] == "+":
            output[trans_id]['exonStarts'].sort()
            output[trans_id]['exonEnds'].sort()
        else:
            output[trans_id]['exonStarts'].sort(reverse=True)
            output[trans_id]['exonEnds'].sort(reverse=True)
        bins = list(zip(output[trans_id]['exonStarts'], output[trans_id]['exonEnds']))
        for start_g, end_g in bins:
            end_t = abs(end_g - start_g) + start_t
            exons[(start_t, end_t)] = (start_g, end_g)
            start_t=end_t+1
        output[trans_id]['exons']=exons
        output[trans_id].pop('exonStarts')
        output[trans_id].pop('exonEnds')

    return output

def mergeFile(transcriptomeAlignment,output,gtf,unlift,rvstrand):
    outfile = open(output+"_t2g.tsv", 'w')
    if unlift == True:unlifted=open(output+'_unlifted.tsv','w')

    transcriptomeAlignment=open(transcriptomeAlignment)
    for line in transcriptomeAlignment:
        row = line.strip().split('\t')
        if row[0]=='':continue
        if rvstrand==False and row[2]=="-":continue
        if row[0]=="chr":
            if unlift == True:unlifted.write(line)
            continue
        if "|" in row[0]:
            row[0]=row[0].split("|")[0]
        transcript = gtf[row[0]]
        if transcript != {}:
            row[0] = transcript['chr']
            row[2] = transcript['strand']
            pos_t = int(row[1])

            pos_g = None
            genome_info_iter = list(transcript["exons"].items())
            for t, g in genome_info_iter:
                start_t, end_t = t
                start_g, end_g = g
                if start_t <= pos_t <= end_t:
                    if row[2] == "+":
                        pos_g = start_g + (pos_t - start_t)
                    elif row[2] == "-":
                        pos_g = start_g - (pos_t - start_t)
                    continue
            if pos_g is None:
                if unlift == True:unlifted.write(line)
            else:
                row[1]=pos_g
                outfile.write('\t'.join(map(str,row))+'\n')
        else:
            if unlift == True:unlifted.write(line)
    outfile.close()
    if unlift == True:unlifted.close()

def calc_pval(treat,ctrl,output_prefix,min_depth,method,fdr_method,fdr_cutoff,min_diff):
    treat_df=pd.read_csv(treat,sep='\t',compression='infer')
    treat_df=treat_df[treat_df.N_total>=min_depth]
    if ctrl is not None:
        ctrl_df=pd.read_csv(ctrl,sep='\t',compression='infer')
        ctrl_df=ctrl_df[ctrl_df.N_total>=min_depth]
    pvalue_col=[]
    outfile = open(output_prefix + "_pval.tsv", 'w')
    if ctrl is None:
        header=['chr','pos','strand','context','ratio','eff_coverage','N_mod','N_total','ratio_ctrl','pvalue']
        outfile.write('\t'.join(header)+'\n')
        N_mod_ctrl=treat_df['N_mod'].sum()
        N_total_ctrl=treat_df['N_total'].sum()
        ctrl_CR=N_mod_ctrl/N_total_ctrl
        for i in range(0,len(treat_df)):
            row=treat_df.iloc[i]
            #if row['ratio'] - ctrl_CR < min_diff:continue
            N_mod=row['N_mod']
            N_total=row['N_total']
            if N_mod>N_total:continue
            if method == "binomial":
                pvalue = scipy.stats.binom_test(x=N_mod, n=N_total, p=ctrl_CR, alternative='greater')
            elif method == "poisson":
                pvalue = scipy.stats.poisson.sf(N_mod, int(math.ceil(ctrl_CR * N_total)))
            elif method == "fisher":
                pvalue = scipy.stats.fisher_exact(table=[[N_mod,N_total-N_mod],[N_mod_ctrl,N_total_ctrl-N_mod_ctrl]],alternative='greater')
                pvalue = pvalue.pvalue
            pvalue_col.append(pvalue)
            outfile.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{:.3f}\t{:.3e}\n'.format(row['chr'],row['pos'],row['strand'],row['context'],row['ratio'],row['eff_coverage'],row['N_mod'],row['N_total'],ctrl_CR,pvalue))
    else:
        header=['chr','pos','strand','context','ratio','eff_coverage','N_mod','N_total','N_mod_ctrl','N_total_ctrl','ratio_ctrl','pvalue']
        outfile.write('\t'.join(header)+'\n')
        # common sites of treat & ctrl
        matched_rows = pd.merge(treat_df.iloc[:, :3], ctrl_df.iloc[:, :3], how='inner')
        matched_treat = pd.merge(matched_rows, treat_df, on=treat_df.columns[:3].tolist())
        matched_ctrl = pd.merge(matched_rows, ctrl_df, on=treat_df.columns[:3].tolist())
        disp("{} common sites found between treat and ctrl".format(len(matched_treat)))
        for i in range(0,len(matched_treat)):
            row_treat=matched_treat.iloc[i]
            row_ctrl=matched_ctrl.iloc[i]
            N_mod=row_treat['N_mod']
            N_total=row_treat['N_total']
            N_mod_ctrl=row_ctrl['N_mod']
            N_total_ctrl=row_ctrl['N_total']
            if N_mod>N_total or N_mod_ctrl>N_total_ctrl:continue
            ctrl_CR=N_mod_ctrl/N_total_ctrl
            #if row_treat['ratio'] - ctrl_CR < min_diff:continue
            if method == "binomial":
                pvalue = scipy.stats.binom_test(x=N_mod, n=N_total, p=ctrl_CR, alternative='greater')
            elif method == "poisson":
                pvalue = scipy.stats.poisson.sf(N_mod, int(math.ceil(ctrl_CR * N_total)))
            elif method == "fisher":
                pvalue = scipy.stats.fisher_exact(table=[[N_mod,N_total-N_mod],[N_mod_ctrl,N_total_ctrl-N_mod_ctrl]],alternative='greater')
                pvalue = pvalue.pvalue
            pvalue_col.append(pvalue)
            outfile.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3e}\n'.format(row_treat['chr'],row_treat['pos'],row_treat['strand'],row_treat['context'],row_treat['ratio'],row_treat['eff_coverage'],row_treat['N_mod'],row_treat['N_total'],N_mod_ctrl,N_total_ctrl,ctrl_CR,pvalue))
    outfile.close()

    fdr_col=multipletests(pvals=pvalue_col,method=fdr_method)[1]
    outfile = open(output_prefix + "_FDR_col.tsv", 'w')
    outfile.write('FDR\n')
    for i in fdr_col:outfile.write('{:.3e}\n'.format(i))
    outfile.close()

    tmp = os.system("paste {} {} > {}".format(output_prefix + "_pval.tsv",output_prefix + "_FDR_col.tsv",output_prefix + "_FDR.tsv"))
    if tmp!=0:disp('Failed to paste FDR column')
    else:
        os.system("rm {} {};gzip {}".format(output_prefix + "_pval.tsv",output_prefix + "_FDR_col.tsv",output_prefix + "_FDR.tsv"))
        disp('FDR values are saved in {}'.format(output_prefix + "_FDR.tsv.gz"))
        if ctrl is None:
            tmp=os.system("zcat "+output_prefix + "_FDR.tsv.gz | awk -v fdr_cutoff="+str(fdr_cutoff)+" -v min_diff="+str(min_diff)+" '{if(NR==1){print}else{d=$5-$9;if(($11<=fdr_cutoff)&&(d>=min_diff)){print}}}' > "+output_prefix + "_SelectedSites.tsv")
        else:
            tmp=os.system("zcat "+output_prefix + "_FDR.tsv.gz | awk -v fdr_cutoff="+str(fdr_cutoff)+" -v min_diff="+str(min_diff)+" '{if(NR==1){print}else{d=$5-$11;if(($13<=fdr_cutoff)&&(d>=min_diff)){print}}}' > "+output_prefix + "_SelectedSites.tsv")
        if tmp!=0:disp('Failed to Select Sites')
        else:disp('Selected Sites are saved in {}'.format(output_prefix + "_SelectedSites.tsv"))

def read_methy_files(ifile, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'modified', 'total']
    disp('Loading ratios in tsv file: {}'.format(ifile))
    meth_file = pd.read_csv(ifile, sep='\t', header=0, usecols=cols, names=names, compression='infer',low_memory=False)
    meth_file.index = meth_file['pos']
    meth_file.drop(['pos'], axis=1, inplace=True)
    return meth_file

def merge_strand_each_chr(df):
    df_p = df[df['strand']=='+']
    df_n = df[df['strand']=='-']
    df_n.index =  df_n.index.values - 1
    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['modified','total']].values
    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['modified','total']].values
    df_merge.columns = ['modified','total']
    df_merge = df_merge.loc[0:,:] # remove the minus index pos -1
    return df_merge

def merge_strand(df):
    chs = df["chr"].unique().tolist()
    df_merge = pd.DataFrame()
    for ch in chs:
        chr_sub = df[df["chr"] == ch]
        if chr_sub.shape[0] > 0:
            chr_sub = merge_strand_each_chr(chr_sub)
            chr_sub['chr']=pd.Series([ch] * chr_sub.shape[0], index=chr_sub.index)
            df_merge=pd.concat([df_merge, chr_sub])#df_merge=df_merge.append(chr_sub)
    return df_merge

def Region_weighted_Ratio(ratio_sub,start=0,end=0):
    #ratio_sub: ratio df of one chrom
    #region_meth=ratio_sub.loc[start:end,:]
    region_meth=ratio_sub[(ratio_sub.index >= start) & (ratio_sub.index <= end)]
    count_C=region_meth.shape[0]
    if count_C > 0:
        region_meth=merge_strand(df=region_meth)
        methy_C=region_meth["modified"].sum()
        total_C=region_meth["total"].sum()
        region_methratio=methy_C*1.0/total_C
    else:
        region_methratio=np.nan
        total_C=np.nan
    return [region_methratio,count_C,total_C]

def bam2epiallele(ifiles,convert_from_base,convert_to_base,conversion_mode,aligner,ref,refmark,coverage,sam_path,unique,pair,rm_dup,trim_fillin,seq_context,chroms,output,read_pos):
    pipes = []
    for ifile in ifiles:
        if ifile[-4:].upper() == '.SAM': sam_format, fin = True, open(ifile)
        elif ifile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view %s' % (sam_path, ifile))
        elif ifile[-5:].upper() == '.CRAM': sam_format, fin = True, os.popen('%ssamtools view %s' % (sam_path, ifile))
        else: sam_format, fin = False, open(ifile)
        pipes.append((ifile,sam_format,fin))
    complement={"A":"T","C":"G","G":"C","T":"A","-":"-"}
    convert_to_base_cp=[complement[i] for i in convert_to_base]
    if conversion_mode == "U":# BS mode
        conversion_rule = {'+': (convert_from_base,[convert_from_base],convert_to_base,complement[convert_from_base],[complement[convert_from_base]],convert_to_base_cp), '-': (complement[convert_from_base],[complement[convert_from_base]],convert_to_base_cp,convert_from_base,[convert_from_base],convert_to_base)}
    elif conversion_mode == "M":# TAPS mode
        conversion_rule = {'+': (convert_from_base,convert_to_base,[convert_from_base],complement[convert_from_base],convert_to_base_cp,[complement[convert_from_base]]), '-': (complement[convert_from_base],convert_to_base_cp,[complement[convert_from_base]],convert_from_base,convert_to_base,[convert_from_base])}

    nmap = 0
    fo_epa = open(output+".unsrt", 'w')
    for ifile, sam_format, fin in pipes:
        disp('Load Alignment: {}'.format(ifile))
        nline = 0; n_NotMatch = 0;
        for line in fin:
            nline += 1
            map_info = Load_One_Read(line,ref,coverage,sam_format,aligner,unique,pair,rm_dup,trim_fillin,chroms)
            if len(map_info) == 0: continue
            seq, strand, cr, pos = map_info
            pos2 = pos + len(seq)
            nmap += 1
            raw, modified, unmodified = conversion_rule[strand][0:3]
            refseq, refmarkcr = ref[cr], {}
            if refmark!={}:refmark=refmark[cr]
            indexes=[]
            for index in re.finditer(raw,refseq[pos:pos2]):
                index0=index.span()[0]
                indexes.append(index0)
            MU_pos_seq={};
            if len(indexes)>0:
                for index in indexes:
                    if refmarkcr=={} or (refmarkcr[index+pos] in seq_context):
                        # report M/U(modified/unmodified) string
                        if seq[index] in modified:
                            MU_pos_seq[index+pos]="M"
                        elif seq[index] in unmodified:
                            MU_pos_seq[index+pos]="U"
            if MU_pos_seq=={}:
                n_NotMatch += 1
            else:
                MU_pos=sorted(MU_pos_seq);MU_seq=[];
                for key in MU_pos:
                    MU_seq.append(MU_pos_seq[key])
                MU_pos=[x+1 for x in MU_pos]# 1-based
                MU_start=MU_pos[0];MU_end=MU_pos[-1];
                # report absolute position of each nucleotide
                if read_pos:
                    fo_epa.write("\t".join([cr,str(MU_start),str(MU_end),strand,"".join(MU_seq), ";".join(map(str,MU_pos)),str(pos),str(pos2)])+"\n")
                else:
                    fo_epa.write("\t".join([cr,str(MU_start),str(MU_end),strand,"".join(MU_seq), ";".join(map(str,MU_pos))])+"\n")
                # report relative position to the first nucleotide of each nucleotide, not convenient for downstream analysis
                #dist_to_CT1=[x-MU_start for x in MU_pos];fo_epa.write("\t".join([cr,str(MU_start),str(MU_end),strand,"".join(MU_seq), ";".join(map(str,dist_to_CT1))])+"\n")
        fin.close()
    fo_epa.close()
    disp('Read {} lines'.format(nline))
    disp('Total {} valid mappings, {} of them donot cover {} base in requested context'.format(nmap,n_NotMatch,convert_from_base))

def PDR(MU_str=[]):
    x=dict(Counter(MU_str))
    x1=0
    for pat in x.keys():
        if "M" in pat and "U" in pat:
            x1+=x[pat]
    return x1*1.0/len(MU_str)

def Entropy(MU_str=[]):
    x=dict(Counter(MU_str))
    x1=0
    for pat in x.keys():
        p=x[pat]*1.0/len(MU_str)
        x1+=p*np.log2(p)
    return 0-x1

def Epipolymorphism(MU_str=[]):
    x=dict(Counter(MU_str))
    x1=0
    for pat in x.keys():
        p=x[pat]*1.0/len(MU_str)
        x1+=p*p
    return 1-x1

def MHL(MU_str=[]):
    x0=0;x1=0;
    for i in range(1,len(MU_str[0])):
        x1=x1+i
        s0=list(map(lambda x: x[0:i], MU_str))
        x0=x0+i*s0.count("M" * i)*1.0/len(MU_str)
    return x0/x1

def epiallele2score(Bed,region_i,epiallele,usestrand=True,min_depth=4,window_size=4,window_depth=1,nucleotide=False):
    region0=Bed.iloc[region_i]
    chr=region0['chr'];
    start=int(region0['start']);
    end=int(region0['end']);
    epa_col = ["start","end","strand","MU_seq","MU_pos","read_count"]
    with os.popen("tabix --verbosity 1 {} {}:{}-{}".format(epiallele,chr,start,end)) as file_in:
        epa_df=pd.read_csv(file_in,sep='\t',header=None,usecols=[1,2,3,4,5,6],names=epa_col)
    if usestrand==True:
        strand=region0['strand']
        epa_df=epa_df[epa_df['strand']==strand]
    else:
        strand=""

    output_lines=[]
    line0=[chr,start,end,strand,0];
    # 8 types of quantifications
    line0.extend([np.nan] * 8);
    # For single base, the weighted and unweighted CAMDA values are the same.
    if nucleotide==True:del(line0[-1])
    if strand=="":del(line0[3])
    # for region, report NA value if no read cover
    if epa_df.shape[0] == 0:
        if nucleotide==False:output_lines.append('\t'.join(map(str,line0))+'\n')
        return output_lines
    read_count_total=epa_df['read_count'].sum()

    if nucleotide==False:
        region_MUX=[];region_PDR=[];region_Entropy=[];region_Epipolymorphism=[];region_MHL=[]

    # If -s not specified, epiallele pattern on reverse strand is converted to forward(pos-1) strand position, and combined with forward strand epialleles.
    epa_forward=epa_df[epa_df['strand']=="+"]
    if epa_forward.shape[0] > 0:
        allc_forward=list(map(int,map(float,set(";".join(map(str,list(epa_forward['MU_pos']))).split(";")))))
    else:allc_forward=[]
    epa_reverse=epa_df[epa_df['strand']=="-"]
    if epa_reverse.shape[0] > 0:
        allc_reverse=list(map(int,map(float,set(";".join(map(str,list(epa_reverse['MU_pos']))).split(";")))))
        allc_reverse= [ x-1 for x in allc_reverse ]
    else:allc_reverse=[]
    allc=list(set(allc_forward+allc_reverse));allc.sort();
    allc=[val for idx, val in enumerate(allc) if val >= start and val <= end]
    del allc_forward, allc_reverse, epa_forward, epa_reverse

    window_dict={};
    if len(allc)>=window_size:
        window_count=len(allc)-window_size+1
        for w in range(0,window_count):window_dict[w]=[]
    if nucleotide==True:
        nucleotide_MUX={}
        for c in range(0,len(allc)):nucleotide_MUX[c]=""

    for row in epa_df.iterrows():
        MU_start=int(row[1]['start']);
        MU_end=int(row[1]['end']);
        MU_strand=row[1]['strand'];
        MU_seq=row[1]['MU_seq'];
        MU_pos=str(row[1]['MU_pos']);
        MU_pos=list(map(int,map(float,MU_pos.split(";"))));MU_pos.sort();
        read_count=int(row[1]['read_count'])
        # merge +/- strand
        if MU_strand=="-":
            MU_start-=1;MU_end-=1;MU_pos=[x-1 for x in MU_pos];
        # for window-based scores: entropy, ...
        if len(allc)>=window_size:
            for w in range(0,window_count):
                start_w=allc[w]
                end_w=allc[w+window_size-1]
                pos0 = [idx for idx, val in enumerate(MU_pos) if val >= start_w and val <= end_w]
                if len(pos0)==window_size:
                    seq0=MU_seq[min(pos0):(max(pos0)+1)]
                    window_dict[w].extend([seq0] * read_count)
        # X: concurrence cytosine
        if "M" in MU_seq and "U" in MU_seq:
            seq0=MU_seq.replace("U","X")
        else:
            seq0=MU_seq;
        if MU_start < start or MU_end > end:
            pos0 = [idx for idx, val in enumerate(MU_pos) if val >= start and val <= end]
            if len(pos0)>0:
                seq0=seq0[min(pos0):(max(pos0)+1)]
                pos0=MU_pos[min(pos0):(max(pos0)+1)]
            else:seq0="";pos=[]
        else:pos0=MU_pos
        if seq0!="":
            if nucleotide==False:
                region_MUX.extend([seq0] * read_count)
            else:
                for i in range(0,len(pos0)):
                    nucleotide_MUX[allc.index(pos0[i])] += seq0[i] * read_count

    if nucleotide==False:
        if window_dict!={}:
            for w in range(0,window_count):
                if len(window_dict[w])>=window_depth:
                    region_PDR.append(PDR(window_dict[w]))
                    region_Entropy.append(Entropy(window_dict[w]))
                    region_Epipolymorphism.append(Epipolymorphism(window_dict[w]))
                    region_MHL.append(MHL(window_dict[w]))
    else:
        for c in range(0,len(allc)):
            MUX=nucleotide_MUX[c]
            if len(MUX)>=min_depth:
                AvgMod=MUX.count("M")*1.0/len(MUX)
                CAMDA_w=MUX.count("X")*1.0/len(MUX)
                CHALM=1-MUX.count("U")*1.0/len(MUX)
                nucleotide_PDR=[];nucleotide_Entropy=[];nucleotide_Epipolymorphism=[];nucleotide_MHL=[]
                c1=max(c-window_size+1,0);c2=min(c,len(allc)-window_size);
                for i in range(c1,c2+1):
                    if len(window_dict[i])>=window_depth:
                        nucleotide_PDR.append(PDR(window_dict[i]))
                        nucleotide_Entropy.append(Entropy(window_dict[i]))
                        nucleotide_Epipolymorphism.append(Epipolymorphism(window_dict[i]))
                        nucleotide_MHL.append(MHL(window_dict[i]))
                v1=v2=v3=v4=np.nan
                if nucleotide_PDR!=[]:v1=np.mean(nucleotide_PDR)
                if nucleotide_Entropy!=[]:v2=np.mean(nucleotide_Entropy)
                if nucleotide_Epipolymorphism!=[]:v3=np.mean(nucleotide_Epipolymorphism)
                if nucleotide_MHL!=[]:v4=np.mean(nucleotide_MHL)

                raw_pos=allc[c]
                if strand=="":
                    output_lines.append('{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(chr,raw_pos,len(MUX),AvgMod,CAMDA_w,CHALM,v1,v2,v3,v4))
                else:
                    if strand=="-":raw_pos=raw_pos+1
                    output_lines.append('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(chr,raw_pos,strand,len(MUX),AvgMod,CAMDA_w,CHALM,v1,v2,v3,v4))

    if nucleotide==False:
        # region AvgMod/CAMDA
        AvgMod=CAMDA_w=CAMDA_unw=CHALM=np.nan
        if region_MUX!=[]:
            M_frag=U_frag=X_frag=0
            for frag0 in region_MUX:
                M_frag+=len(re.findall("M+", frag0))
                U_frag+=len(re.findall("U+", frag0))
                X_frag+=len(re.findall("X+", frag0))
            CAMDA_unw=X_frag*1.0/(M_frag+U_frag+X_frag)
            region_MUX="".join(region_MUX)
            AvgMod=region_MUX.count("M")*1.0/len(region_MUX)
            CAMDA_w=region_MUX.count("X")*1.0/len(region_MUX)
            CHALM=1-region_MUX.count("U")*1.0/len(region_MUX)
        # region window-based scores
        v1=v2=v3=v4=np.nan
        if region_PDR!=[]:v1=np.mean(region_PDR)
        if region_Entropy!=[]:v2=np.mean(region_Entropy)
        if region_Epipolymorphism!=[]:v3=np.mean(region_Epipolymorphism)
        if region_MHL!=[]:v4=np.mean(region_MHL)

        if strand=="":
            output_lines.append('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(chr,start,end,read_count_total,AvgMod,CAMDA_w,CAMDA_unw,CHALM,v1,v2,v3,v4))
        else:
            output_lines.append('{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(chr,start,end,strand,read_count_total,AvgMod,CAMDA_w,CAMDA_unw,CHALM,v1,v2,v3,v4))

    return output_lines
