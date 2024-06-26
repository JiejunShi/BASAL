#ifndef _PARAM_H_
#define _PARAM_H_

#define SEGLEN 32

#include <unistd.h>
#include<string>
#include<cstring>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include<vector>
#include<algorithm>

using namespace std;

const unsigned int FIXELEMENT=16; // 480/32+1
const unsigned int MAXSNPS=15;
const unsigned int MAXGAPS=3;//max length of gap, only one gap permitted

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;

typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;
typedef ref_loc_t* NewIndex;

struct Hit {
	ref_id_t chr;       //index of chr
	ref_loc_t loc;     //location of first bp on reference seq, count from 0
};

struct gHit {
    ref_loc_t loc;
    bit32_t chr: 18;
    bit32_t strand: 2; // strand infomation 00: '++', 01: '+-', 10: '-+', 11: '--'
    //bit16_t gap_size:4;  // positive for insert on read, negative for insert on ref, 0 for no-gap
    int gap_size;//use int type, or there will be bug for insertions
    bit16_t gap_pos:9;//511 at most
};

class Param {
public:
	Param();
	void SetSeedSize(int n);
    void InitMapping();
	void SetAdaptors(int n);
    void SetDigestionSite(const char *a);
    //void SetAlign(char readnt, char refnt);
    void SetAlign(string &convert_rule0);

public:

    //hot data section
    
    //cache line 1
    bit32_t seed_size __attribute__((aligned(64))); // 4 byte
    bit32_t seed_bits, max_snp_num, max_num_hits; //16 byte
    bit32_t min_insert, max_insert, max_seedseg_num, chains; //32 byte
    bit32_t RRBS_flag, pairend, index_interval, randseed; //48 byte
    bit32_t gap, gap_edge, trim_lowQ, max_ns; //64 byte

    //cache line 2
    bit32_t min_read_size, n_adapter, read_start, read_end; // 16 byte
    bit32_t out_sam, out_ref, out_unmap, report_repeat_hits; //32 byte   
    bit32_t sam_header, max_readlen; bit16_t *_T; //48 byte
    int input_format, output_format, gz_input, gz_ref; // 64 byte    

    //cache line 3
    bit32_t max_dbseq_size, append_dbseq_size, total_ref_seq;
    vector <bit32_t> digest_pos; // 16 byte
    vector <string> digest_site; 
    string useful_nt; // 32 byte
    string nx_nt; int num_procs; //48 byte
    bit32_t stdout; bit8_t zero_qual, qual_threshold, default_qual; //51 byte
    bit8_t N_mis;//bit8_t read_nt, ref_nt, N_mis;
    int readnt_cnt;
    char refnt,readnt;
    char readnts[5];
    bit32_t max_kmer_num, verbose_level; 
    float max_kmer_ratio;
    bit32_t nt3, pipe_out, seed_bits_lz; 
    
    //cache line 4+    
    bit32_t profile[MAXSNPS+1][16] __attribute__((aligned(64)));
    
    //int chains;   //0: forward strands only ; 1: forward and reverse strands	
	string adapter[10];
    //int input_format; // 0: fasta, 1:fastq, 2: SAM, 3: BAM, -1: auto detect
    //int output_format; // 1:SAM, 2: BAM, -1: auto detect
    //int gz_input, gz_ref; // 0: no, 1: yes, -1: auto detect

    void disp_bin64(bit64_t t, int len=32) {
        for(int i=len*2-1;i>=0;i--) {
            if(t&(1ULL<<i)) cout<<'1';
            else cout<<'0';
            if(i%4==0) cout<<' ';
        }
        cout<<endl;
    }

    inline bit64_t XT64(register bit64_t tt) {tt-=(tt<<1)&tt&0xAAAAAAAAAAAAAAAAULL; return tt;}
    inline bit32_t XT32(register bit32_t tt) {tt-=(tt<<1)&tt&0xAAAAAAAAUL; return tt;}

    inline bit32_t XT(register bit32_t tt) {
        register bit32_t ss;
        tt-=(tt<<1)&tt&0xAAAAAAAAUL;
        tt-=(tt>>2)&0x33333333UL; // 4 bit transform
        ss=(tt&0xF0F0F0F0UL)>>1;
        tt-=ss-(ss>>3); //8 bit transform, ss+(ss<<3) = s*9
        ss=(tt&0xFF00FF00UL)>>2; 
        tt=(tt&0x00FF00FFUL)+ss+(ss>>2)+(ss>>6);
        return (tt&0xFFFFUL)+(tt>>16)*6561;
    }                                                                                                       

    inline bit32_t XC(bit32_t tt) {return ((~tt)<<1)|tt|0x55555555UL;}  // generate T2C mask according to C locations
    inline bit64_t XC64(register bit64_t tt) {return ((~tt)<<1)|tt|0x5555555555555555ULL;}
    inline bit32_t swap_endian(bit32_t tt) {return ((tt>>16)|(tt<<16));}
    inline bit64_t swap_endian64(bit64_t tt) {return ((tt>>32)|(tt<<32));}

    inline bit32_t XM(bit32_t tt) {
        tt=(tt|(tt>>1))&0x55555555;
        tt=(tt+(tt>>2))&0x33333333;
        return (((bit32_t)(tt*0x1111111))>>28)+(tt&0x3);
    }

    inline bit32_t XM64(register bit64_t tt) {
        tt|=tt>>1;
        tt&=0x5555555555555555ULL;
        tt+=tt>>2;
        tt&=0x3333333333333333ULL;
        tt+=tt>>4;
        tt&=0x0F0F0F0F0F0F0F0FULL;
        tt*=0x0101010101010101ULL;
        tt>>=56;
        return tt;
    }

    // 01 -> 00, 11 unchanged
    inline bit64_t M2_judge(register bit64_t tt) {tt&=((tt & 0xAAAAAAAAAAAAAAAAULL) >> 1) | ((tt & 0x5555555555555555ULL) << 1); return tt;}
    
    bit32_t map3to4(bit32_t tt){
      	int s=0, i; for(i=0;i<16;i++) {s|=(tt%3)<<i*2; tt/=3;};
	return s;
    }
};

#endif //_PARAM_H_
