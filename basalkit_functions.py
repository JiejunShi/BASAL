#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os,re,sys,time
from collections import defaultdict,OrderedDict
import numpy as np
import pandas as pd
import math
import scipy.stats
from statsmodels.stats.multitest import multipletests
from copy import copy
import pysam

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
        elif op == 'D':
            seq = seq[:index]+'-'*length+seq[index:]
            index += length
        elif op == 'N':
            seq = seq[:index]+'+'*length+seq[index:]
            index += length
        elif op == 'H' or op == 'P':
            continue
        #else:
            #raise ValueError("%c not a valid CIGAR character"%(op))
    #assert originalLen == index, "String length does not match CIGAR"
    return seq

def strand_bismark2basal(XR,XG):
    if XR == "CT" and XG == "CT": return "++"
    elif XR == "CT" and XG == "GA": return "-+"
    elif XR == "GA" and XG == "CT": return "+-"
    elif XR == "GA" and XG == "GA": return "--"
    else: return ""

def strand_gemBS2basal(XB,sam_flag):
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

def Load_One_Read(line,ref,coverage,sam_format,molecule_type,aligner,unique,pair,rm_dup,trim_fillin,chroms):
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
            strand=strand_bismark2basal(XR=XR_tag,XG=XG_tag)
        elif aligner == "gemBS":
            XB_index = line.find('XB:A:')
            XB_tag = line[XB_index+5:XB_index+6]
            strand=strand_gemBS2basal(XB=XB_tag,sam_flag=flag)
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

    if strand not in ['++',"-+","+-","--"]:return []
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
    if molecule_type=="DNA":
        return (seq, strand[0], cr, pos)
    else:
        if (16&flag!=16):
        #if direction==1:
            return (seq, "+", cr, pos)
        else:
        #elif direction==2:
            return (seq, "-", cr, pos)

def Load_Alignment(ifiles,convert_from_base,convert_to_base,conversion_mode,molecule_type,aligner,ref,refmark,coverage,meth0,depth,meth_ct,depth_ct,sam_path,unique,pair,rm_dup,trim_fillin,chroms,seq_context,handle_SNP,converted_site):
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

    nmap = 0
    for ifile, sam_format, fin in pipes:
        disp('Load Alignment: {}'.format(ifile))
        nline = 0
        for line in fin:
            nline += 1
            map_info = Load_One_Read(line,ref,coverage,sam_format,molecule_type,aligner,unique,pair,rm_dup,trim_fillin,chroms)
            if len(map_info) == 0: continue
            seq, strand, cr, pos = map_info
            pos2 = pos + len(seq)
            nmap += 1
            raw, modified, unmodified, raw_rc, modified_rc, unmodified_rc = conversion_rule[strand]
            refseq, refmarkcr = ref[cr], {}
            if refmark!={}:refmarkcr=refmark[cr]
            indexes=[];n_covered_1read=0;n_converted_1read=0;
            for index in re.finditer(raw,refseq[pos:pos2]):
                n_covered_1read += 1
                index0=index.span()[0]
                indexes.append(index0)
                if seq[index0] in unmodified and conversion_mode == "U":n_converted_1read += 1
                if seq[index0] in modified:
                    if conversion_mode == "M":n_converted_1read += 1
            if converted_site>=1:
                if n_converted_1read < converted_site:continue
            else:
                if n_converted_1read < converted_site*n_covered_1read:continue
            if n_covered_1read>0:
                depth_cr = depth[cr];
                meth0_cr = meth0[cr];
                for index in indexes:
                    if (refmarkcr=={} or (refmarkcr[index+pos] in seq_context)) and depth_cr[index+pos] < (2**16-1):
                        if seq[index] in unmodified:
                            depth_cr[index+pos] += 1
                        elif seq[index] in modified:
                            depth_cr[index+pos] += 1
                            meth0_cr[index+pos] += 1
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
    return meth0,depth,meth_ct,depth_ct,nmap

def rightmostD(cigar,XR,convert_from_base):
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

def shiftD(Alignments,output,convert_from_base,sam_path):
    # load Alignment
    if Alignments[-4:].upper() == '.SAM':
        fin=open(Alignments)
    elif Alignments[-4:].upper() == '.BAM':
        fin=os.popen('%ssamtools view -h %s' % (sam_path, Alignments))
    else:
        disp('Input Alignment must be sam or bam file');sys.exit()
    fout = open(output+".sam", 'w')
    disp('Load Alignments in: {}'.format(Alignments))
    for line in fin:
        if line[0] == '@':fout.write(line);continue
        col = line.split('\t')
        XR="";ZS=""
        for col0 in col:
            if col0.startswith("XR:Z:"):
                XR=re.sub(r"[a-z]","",col0[5:])#ref
            if col0.startswith("ZS:Z:"):
                ZS=col0[5:7]
        if XR=="" or ZS=="":
            disp("The XR and ZS tag of input alignments are required. Please make sure it is generated by BASAL with -R option. Exit.")
            sys.exit()
        cigar=col[5]
        if 'D' not in cigar:fout.write(line);continue
        which_base=convert_from_base
        if ZS=="+-" or ZS=="-+":which_base=reverse_complement(convert_from_base)
        cigar_new=rightmostD(cigar=cigar,XR=XR,convert_from_base=which_base)
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

def Out_base_ratio(tsv_prefix,wig_prefix,wig_bin,min_depth,ref,refmark,handle_SNP,convert_from_base,seq_context_str,seq_context,motif_length,meth0,depth,meth_ct,depth_ct,nmap):

    header=['chr','pos','strand','context','ratio','eff_coverage','N_mod','N_total']
    if handle_SNP > 0: header.extend(['N_mod_rev','N_total_rev'])

    fo_mr = open(tsv_prefix+"_AvgMod.tsv", 'w')
    fo_mr.write('\t'.join(header)+'\n')
    if wig_prefix!=None:
        fo_mr_wig = open(wig_prefix+"_AvgMod.wig", 'w');fo_mr_wig.write('track type=wiggle_0 name='+wig_prefix+'_MethRatio\n')
        disp('Output ratios in tsv files and wiggle files')
    else:
        disp('Output ratios in tsv files')

    #z95, z95sq = 1.96, 1.96 * 1.96
    nc, nd, dep0 = 0, 0, min_depth
    for cr in sorted(depth.keys()):
        depth_cr, meth0_cr, refcr, refmarkcr = depth[cr], meth0[cr], ref[cr], {}
        if refmark!={}:refmarkcr=refmark[cr]
        if handle_SNP > 0: depth_ct_cr, meth_ct_cr = depth_ct[cr], meth_ct[cr]
        if wig_prefix!=None:
            fo_mr_wig.write('variableStep chrom={} span={}\n'.format(cr, wig_bin))
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

            # for C2T, candidate site not in defined context (-x) is excluded
            if seq_context != [] and refmarkcr[i] not in seq_context:continue

            # for any conversion, output motif sequence
            if refcr[i] == convert_from_base:seq=refcr[i-motif_length:i+motif_length+1]
            else:seq=reverse_complement(refcr[i-motif_length:i+motif_length+1])

            if refcr[i] == convert_from_base:strand = '+'
            else:strand = '-'

            m_mr = meth0_cr[i];
            try: ratio_mr = min(m_mr, d)*1.0 / d
            except ZeroDivisionError: continue
            nc += 1
            nd += d
            if wig_prefix!=None:
                if i // wig_bin != bin:
                    if wigd > 0:
                        fo_mr_wig.write('{:.0f}\t{:.3f}\n'.format(bin*wig_bin+1, min(wigm0/wigd,1)))# wiggle is 1-based
                    bin = i // wig_bin#use integer division
                    wigd = wigm0 = 0.0
                wigd += d
                wigm0 += m_mr
            #denorminator = 1 + z95sq / d
            #pmid_mr = ratio_mr + z95sq / (2 * d)
            #sd_mr = z95 * ((ratio_mr*(1-ratio_mr)/d + z95sq/(4*d*d)) ** 0.5)
            #CIl_mr, CIu_mr = (pmid_mr - sd_mr) / denorminator, (pmid_mr + sd_mr) / denorminator
            if handle_SNP > 0:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd, m1, d1))# 1-based
            else:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd))# 1-based
    fo_mr.close();
    if wig_prefix!=None:
        fo_mr_wig.close();
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

def generate_new_cigar(all_bins, start, end, old_cigar, trans_dir):
    ''' order: small --> big corrdinate '''
    new_cigar_tmp = []  # no del and insert, with intron
    if trans_dir == "-":
        old_cigar = old_cigar[::-1]
        all_bins = all_bins[::-1]
        start, end = end, start
    all_bins_iter = iter(all_bins)
    while (1):
        try:
            x, y = next(all_bins_iter)
            if trans_dir == "-": x, y = y, x
            if x <= start <= y < end:
                new_cigar_tmp.append([0, y - start + 1])
                exon_edge = y
            elif x <= start <= end <= y:
                new_cigar_tmp.append([0, end - start + 1])
                break
            elif start < x <= y < end:
                if x - exon_edge - 1 > 0:
                    new_cigar_tmp.append([3, x - exon_edge - 1])
                new_cigar_tmp.append([0, y - x + 1])
                exon_edge = y
            elif start < x <= end <= y:
                if x - exon_edge - 1 > 0:
                    new_cigar_tmp.append([3, x - exon_edge - 1])
                new_cigar_tmp.append([0, end - x + 1])
                break
        except StopIteration:
            sys.exit()

    new_cigar_tmp_tmp = []

    new_cigar_tmp_iter = iter(new_cigar_tmp)
    cigar_type, number = next(new_cigar_tmp_iter)
    while (1):
        try:
            cigar_type_1, number_1 = next(new_cigar_tmp_iter)
            if cigar_type == cigar_type_1:
                number = number + number_1
            else:
                new_cigar_tmp_tmp.append([cigar_type, number])
                cigar_type, number = cigar_type_1, number_1
        except StopIteration:
            new_cigar_tmp_tmp.append([cigar_type, number])
            break
    new_cigar_tmp = new_cigar_tmp_tmp
    new_cigar = []
    # debug

    # old_M, old = cal(old_cigar)
    # NT_M, NT = cal(new_cigar_tmp)

    new_cigar_tmp_iter = iter(new_cigar_tmp)
    block = next(new_cigar_tmp_iter)
    for cigar_type, num in old_cigar:
        try:
            if block[0] == 3:
                new_cigar.append((block[0], block[1]))
                block = next(new_cigar_tmp_iter)
            if cigar_type == 0:  # matched
                if num < block[1]:  # samller than the original block
                    new_cigar.append((0, num))
                    block[1] = block[1] - num
                elif num == block[1]:  # remove a block
                    new_cigar.append((0, num))
                    block = next(new_cigar_tmp_iter)
                    if block[0] == 3:  # intron
                        new_cigar.append((block[0], block[1]))
                        block = next(new_cigar_tmp_iter)
                else:
                    while num > block[1]:
                        new_cigar.append((0, block[1]))
                        num = num - block[1]
                        block = next(new_cigar_tmp_iter)
                        new_cigar.append((block[0], block[1]))  # intron
                        block = next(new_cigar_tmp_iter)  # until the last exon
                    if num == block[1]:
                        new_cigar.append((0, num))
                        block = next(new_cigar_tmp_iter)
                    elif num < block[1]:
                        block[1] = block[1] - num
                        new_cigar.append((0, num))
                    if block[1] < 0:
                        raise Warning("Block start <0")
            elif cigar_type == 1:  # insert
                new_cigar.append((1, num))
            elif cigar_type == 2:  # del
                if num < block[1]:
                    new_cigar.append((2, num))
                    block[1] = block[1] - num
                elif num == block[1]:
                    new_cigar.append((2, num))
                    block = next(new_cigar_tmp_iter)
                    if block[0] == 3:
                        new_cigar.append((block[0], block[1]))
                        block = next(new_cigar_tmp_iter)
                else:
                    while num > block[1]:
                        new_cigar.append((2, block[1]))
                        num = num - block[1]
                        block = next(new_cigar_tmp_iter)
                        new_cigar.append((block[0], block[1]))  # intron
                        block = next(new_cigar_tmp_iter)  # until the last exon
                    if num == block[1]:
                        new_cigar.append((2, num))
                        block = next(new_cigar_tmp_iter)
                    elif num < block[1]:
                        block[1] = block[1] - num
                        new_cigar.append((2, num))
                    if block[1] < 0:
                        raise Warning("Block start <0")
            elif cigar_type == 3:
                new_cigar.append((3, num))
            elif cigar_type == 4:
                new_cigar.append((4, num))
            elif cigar_type == 5:
                new_cigar.append((5, num))
            elif cigar_type == 6:
                new_cigar.append((6, num))
        except StopIteration:
            continue
    # new_M,new = cal(new_cigar)
    # if new_M != old_M:

    # debug

    return new_cigar

def map_to_genome(header_dict,gtf,segment,unlift,UNLIFT):
    try:
        if '|' in segment.reference_name:
            reference_name = segment.reference_name.split("|")[0]
        else:
            reference_name = segment.reference_name
        genome_info = gtf.get(reference_name)
        new_ref_id = header_dict.get(genome_info['chr'])
        trans_dir = genome_info['strand']
    except TypeError:
        genome_info = None
        new_ref_id = None
    if genome_info and new_ref_id is not None:
        old_start = segment.reference_start  # 0-based
        old_end = segment.reference_end - 1  # 1-based
        modified_exons = OrderedDict(((key[0] - 1, key[1] - 1), (value[0] - 1, value[1] - 1))
            for key, value in genome_info["exons"].items())
        #genome_info["exons"] = modified_exons

        new_start = None
        new_end = None
        if trans_dir == "+":
            genome_info_iter = list(modified_exons.items())
        elif trans_dir == "-":
            genome_info_iter = list(modified_exons.items())[::-1]
        list_maxend = []
        for key, values in genome_info_iter:
            list_maxend += [key[0], key[1]]
        len_transcript = max(list_maxend)
        if old_end <= len_transcript:
            while new_start is None or new_end is None:
                for key, values in genome_info_iter:
                    start, end = key
                    geno_start, geno_end = values
                    if trans_dir == "+":
                        if start <= old_start <= end:
                            new_start = geno_start + old_start - start
                        if start <= old_end <= end:
                            new_end = geno_start + old_end - start
                    elif trans_dir == "-":
                        geno_start, geno_end = geno_end, geno_start
                        if start <= old_end <= end:
                            new_end = geno_start + (end - old_end)
                        if start <= old_start <= end:
                            new_start = geno_start + (end - old_start)
                    else:
                        raise Warning("Transcription direction loss.")
                new_cigar = generate_new_cigar(list(modified_exons.values()), new_start, new_end, segment.cigar,
                                               genome_info['strand'])

                qual = segment.query_qualities
                mpq = segment.mapping_quality
                seq = segment.query_sequence
                tags = segment.get_tags()
                flag = segment.flag

                segment_output = pysam.AlignedSegment()
                segment_output.tags = segment.tags
                if trans_dir == "-":
                    new_start, new_end = new_end, new_start
                    qual = qual[::-1]
                    seq = reverse_complement(segment.query_sequence)
                    if segment.is_reverse:
                        segment.is_reverse = False
                        segment.mate_is_reverse = True
                    else:
                        segment.is_reverse = True
                        segment.mate_is_reverse = False
                    if flag == 0:
                        new_flag = 16
                    else:
                        flags = []
                        bit = 1
                        while flag >= bit:
                            if flag & bit:
                                flags.append(bit)
                            bit <<= 1
                        if 16 in flags:
                            flags.remove(16)
                        else:
                            flags.append(16)
                        new_flag = 0
                        for f in flags:
                            new_flag |= f
                else:new_flag=segment.flag
                if trans_dir == "-":
                    if "ZS" in [tag[0] for tag in tags]:
                        for tag, value in tags:
                            if tag == "ZS":
                                ts_value = value
                                if ts_value=='++':new_ts='-+'
                                elif ts_value=='+-':new_ts='--'
                                elif ts_value=='-+':new_ts='++'
                                elif ts_value=='--':new_ts='+-'
                    segment_output.set_tag("ZS", new_ts, value_type="Z")
                    if "XR" in [tag[0] for tag in tags]:
                        for tag, value in tags:
                            if tag == "XR":
                                xr_value = value
                                xr_value = xr_value.upper()
                                new_xr=reverse_complement(xr_value)
                                new_xr=new_xr[:2].lower() + new_xr[2:-2] + new_xr[-2:].lower()
                        segment_output.set_tag("XR", new_xr, value_type="Z")

                segment_output.query_name = segment.query_name
                segment_output.flag = new_flag
                segment_output.reference_id = new_ref_id
                segment_output.reference_start = new_start
                segment_output.cigar = new_cigar
                segment_output.query_sequence = seq
                segment_output.query_qualities = qual
                segment_output.mapping_quality = mpq
                segment_output.set_tag("TN", segment.reference_name)
                if segment_output:
                    return segment_output
        else:
            if unlift == True:
                UNLIFT.write(segment)
    else:
        if unlift == True:
            UNLIFT.write(segment)

def read_headers(fn,hid,new_header,hid_dict,lift_over):
    lift_over[fn] = {}
    with pysam.AlignmentFile(fn, 'rb') as INPUT:
        n = 0
        for header in INPUT.header["SQ"]:
            if header['SN'] not in hid_dict:
                # print(header['SN'])
                hid_dict[header['SN']] = hid
                new_header['SQ'].append(header)
                hid += 1
            lift_over[fn][n] = hid_dict[header['SN']]
            n += 1
    return hid,new_header,hid_dict,lift_over

def merge_bam(lift_over,fn,fout):
    with pysam.AlignmentFile(fn, 'rb') as INPUT:
        for read in INPUT:
            read.reference_id = lift_over[fn][read.reference_id]
            read.next_reference_id = -1
            read.next_reference_start = 0
            fout.write(read)

def calc_pval(treat,ctrl,output_prefix,min_depth,method,fdr_method):
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
                #pvalue = scipy.stats.fisher_exact(table=[[N_mod, N_total-N_mod], [N_mod_ctrl, N_total_ctrl-N_mod_ctrl]], alternative='greater')
                test_statistic, pvalue = scipy.stats.fisher_exact(table=[[N_mod, N_total-N_mod], [N_mod_ctrl, N_total_ctrl-N_mod_ctrl]], alternative='greater')
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
        del matched_rows,treat_df,ctrl_df
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
                #pvalue = scipy.stats.fisher_exact(table=[[N_mod, N_total-N_mod], [N_mod_ctrl, N_total_ctrl-N_mod_ctrl]], alternative='greater')
                test_statistic, pvalue = scipy.stats.fisher_exact(table=[[N_mod, N_total-N_mod], [N_mod_ctrl, N_total_ctrl-N_mod_ctrl]], alternative='greater')
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
