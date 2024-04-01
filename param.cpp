#include "param.h"
#include "utilities.h"
#include<iostream>

using namespace std;

Param::Param()
{
	//num_procs=sysconf(_SC_NPROCESSORS_ONLN);
    num_procs=1;
	//if(num_procs>8) num_procs=8;

	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
	
	max_ns = 5;
	trim_lowQ=0;
	
	zero_qual= '!';
	qual_threshold= 0;
	default_qual=40;
	
	min_insert= 28;
	max_insert= 1000;

	SetSeedSize(16);	
	//seed_size= 16;
	//seed_bits=(1<<(seed_size*2))-1;
	
	max_snp_num = 110;
	max_num_hits = MAXHITS>100?100:MAXHITS;
	max_kmer_ratio = 5e-7;
	
	min_read_size=seed_size;
	input_format=gz_input=gz_ref=-1;	
    n_adapter=0;

	report_repeat_hits = 1;

	useful_nt="ACGTacgt";
	nx_nt="NXnx";
	
	for(bit32_t i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;

    out_sam=1;
	read_start=1;
	read_end=~0;

    out_ref=0;
    out_unmap=0;
    RRBS_flag=0;
    index_interval=4;
    randseed=0;
    chains=0;
    pairend=0;
    gap=0;
    gap_edge=6;
    max_readlen=(FIXELEMENT-1)*SEGLEN;
    refnt='C';readnt='T';
    readnt_cnt=0;
    sam_header=1;
    stdout=1;
    N_mis=0;
    nt3=0;
    pipe_out=0;
    verbose_level=1;

};

void Param::InitMapping(){
	for(bit32_t i=0; i<index_interval; i++) for(bit32_t j=0; j<=MAXSNPS; j++) {  //for 4 binary seqs
	    profile[j][i]=((j*seed_size+i+index_interval-1)/index_interval)*index_interval;
	}
}

void Param::SetDigestionSite(const char *a) {
    string ds=a, ds1; bit32_t dp;
    int IUPAC_count[256], i, j, ii[256], l;
    char IUPAC[256][4];
    for(i=0;i<256;++i) {IUPAC_count[i]=0; for(j=0;j<4;++j) IUPAC[i][j]=(char) 0;}
    IUPAC['A'][0]='A'; IUPAC['C'][0]='C'; IUPAC['G'][0]='G'; IUPAC['T'][0]='T'; 
    IUPAC['N'][0]='A'; IUPAC['N'][1]='C'; IUPAC['N'][2]='G'; IUPAC['N'][3]='T';
    strcpy(IUPAC['R'], "AG"); strcpy(IUPAC['Y'], "CT"); strcpy(IUPAC['S'], "CG");
    strcpy(IUPAC['W'], "AT"); strcpy(IUPAC['K'], "GT"); strcpy(IUPAC['M'], "AC");
    strcpy(IUPAC['B'], "CGT"); strcpy(IUPAC['D'], "AGT"); strcpy(IUPAC['H'], "ACT"); strcpy(IUPAC['V'], "ACG");
    for(i=0;i<256;++i) for(j=0;j<4;++j) IUPAC_count[i]+=(IUPAC[i][j]>0);
    //cout<<digest_site<<endl;
    if((dp=ds.find('-'))<0) {
        cout<<"Digestion position not marked, use \'-\' to mark. example: \'C-CGG\'\n";  
        exit(1);
    }
    ds.erase(dp,1); l=ds.size(); ds1=ds;
    for(i=0;i<l;++i) ii[i]=0; j=0;
    while(j<l) {
        for(i=0;i<l;++i) ds1[i]=IUPAC[ds[i]][ii[i]];
        digest_site.push_back(ds1);
        digest_pos.push_back(dp);
        j=0; ii[j]++; 
        while(ii[j]>=IUPAC_count[ds[j]]&&j<l) {
            ii[j]=0; j++; ii[j]++;
        }
    }
    RRBS_flag=1;
    index_interval=1;
    //SetSeedSize(12);
}

void Param::SetSeedSize(int n) {
    if(n>16||n<10) {cerr<<"seed size must be between 10 and 16\n"; exit(1);}
	seed_size=n;
	seed_bits_lz=(SEGLEN-seed_size)*2;
	min_read_size=seed_size+index_interval-1;	
	seed_bits=0;
	for(bit32_t i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;
}

bit8_t alphabet[256];
bit8_t rev_alphabet[256];
bit8_t alphabet0[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

bit8_t reg_alphabet[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

//mask convert_to_bases as 01, others are 11
bit8_t alphabet_Mread[256];
bit8_t rev_alphabet_Mread[256];

//bit8_t nv3='N'; //set any unknown char as 'N'
char nv3='N';
char rev_char[256] = {
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'A' */ 
'T',nv3,'G',nv3,nv3,nv3,'C',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,'A',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'a' */ 
't',nv3,'g',nv3,nv3,nv3,'c',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,'a',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3
};

bit8_t bit_nt[4];

char nt_code[5] = {'A', 'C', 'G', 'T', '-'};
char revnt_code[5] = {'T', 'G', 'C', 'A', '-'};

void Param::SetAlign(string &convert_rule0){
    if (convert_rule0[1]!=':'){
        cerr<<"invalid -M, ref base(one letter in A/C/G/T) should be assigned first before :"<<endl;exit(1);
    }
    //char readnts[5] = {' ', ' ', ' ', ' ', ' '};//convert_to_bases
    for(int i=0;i<5;i++) readnts[i]=' ';
    //char refnt,readnt;
    char readnt_valid = 0;char readnt_used = 0;
    refnt = toupper(convert_rule0[0]);
    if(!reg_alphabet[refnt]){
        cerr<<"invalid -M, ref base "<< convert_rule0[0] <<" not in A/C/G/T"<<endl;exit(1);
    }
    cout << "[BASAL @"<<Curr_Time()<<"] convert-from base: " << refnt << endl;
    string convert_rule=convert_rule0.substr(2,convert_rule0.length());
    for (size_t i = 0; i < convert_rule.size(); i++){
        readnt = toupper(convert_rule[i]);
        for(int j = 0; j<5; j++){
            if (nt_code[j] == readnt){
                readnt_valid = 1;
            };
            if (readnts[j] == readnt){
                readnt_used = 1;
            };
        };
        if(readnt==refnt){
            cerr<<"invalid -M, read base "<< convert_rule[i] <<" should not be equal to ref base "<< refnt <<endl;exit(1);
        }
        if(readnt_valid==0){
            cerr<<"invalid -M, read base "<< convert_rule[i] <<" not in A/C/G/T/-"<<endl;exit(1);
        }
        if(readnt_used==0){
            readnts[readnt_cnt] = readnt;
            readnt_cnt++;
        }
        readnt_valid = 0;
        readnt_used = 0;
    }
    cout << "[BASAL @"<<Curr_Time()<<"] convert-to base(s):"<< readnts[0] << readnts[1] << readnts[2] << readnts[3] << readnts[4] << endl;
    //mask convert_to_bases as 01, others are 11
    for(int i=0;i<256;i++){
        alphabet_Mread[i] = reg_alphabet[i];
        rev_alphabet_Mread[i] = reg_alphabet[i];
    };
    for(int i=0;i<readnt_cnt;i++){
        alphabet_Mread[readnts[i]] = 1;
        alphabet_Mread[tolower(readnts[i])] = 1;
        for(int j = 0; j<5; j++){
            if ((nt_code[j] == readnts[i]) && (readnts[i]!='-')){
                rev_alphabet_Mread[revnt_code[j]] = 1;
                rev_alphabet_Mread[tolower(revnt_code[j])] = 1;
            };
        };
    };    
    int i, j; j = 0;
    for(i=0;i<4;i++) bit_nt[i]=100;
    // convert_from_base must be 01, consistent to XC64(keep 01 as 01, mask 00/10/11 as 11)
    bit_nt[alphabet0[refnt]]=1;
    // if only one convert_to_base(not '-'), mark convert_to_base as 11, otherwise mark other bases as 0/2/3(order dowsn't matter)
    int other_bit[3] = {0, 2, 3};
    if(readnt_cnt==1){
        if(readnts[0]!='-'){
            bit_nt[alphabet0[readnts[0]]]=3;
            int other_bit[2] = {0, 2};
        };
    };
    for(i=0;i<4;i++){
        if(bit_nt[i]==100){
            bit_nt[i]=other_bit[j];
            j++;
        };
    };
    //for(i=0;i<5;i++) cout<<" "<<nt_code[i]<<":"<<(int) bit_nt[i]; cout<<endl;
    // for -M C:T, bit_nt[4] = {0,1,2,3}; for -M A:G, bit_nt[4] = {1,0,2,3}

    // alphabet indicates the bits, for BS: A = 00, C = 01, G = 10, T = 11
    for(i=0;i<256;i++) alphabet[i]=0;
    alphabet[(unsigned char)'a']=bit_nt[0]; alphabet[(unsigned char)'A']=bit_nt[0]; 
    alphabet[(unsigned char)'c']=bit_nt[1]; alphabet[(unsigned char)'C']=bit_nt[1]; 
    alphabet[(unsigned char)'g']=bit_nt[2]; alphabet[(unsigned char)'G']=bit_nt[2];
    alphabet[(unsigned char)'t']=bit_nt[3]; alphabet[(unsigned char)'T']=bit_nt[3];
    //alphabet[(unsigned char)'-']=bit_nt[4];
    //for(i=0;i<5;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<5;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(alphabet[nt_code[i]]); cout<<endl;
    // for -M C:T, alphabet = alphabet0 = {...,0(which is A),...,1(which is C),...,2(which is G),...,3(which is T),...}, others are all 0
    // for -M A:G, alphabet[256] = {...,1(which is A),...,0(which is C),...,3(which is G),...,2(which is T),...}, others are all 1

    for(i=0;i<256;i++) rev_alphabet[i]=0;
    rev_alphabet[(unsigned char)'a']=bit_nt[3]; rev_alphabet[(unsigned char)'A']=bit_nt[3]; 
    rev_alphabet[(unsigned char)'c']=bit_nt[2]; rev_alphabet[(unsigned char)'C']=bit_nt[2]; 
    rev_alphabet[(unsigned char)'g']=bit_nt[1]; rev_alphabet[(unsigned char)'G']=bit_nt[1];
    rev_alphabet[(unsigned char)'t']=bit_nt[0]; rev_alphabet[(unsigned char)'T']=bit_nt[0];
    //for(i=0;i<5;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) rev_alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<5;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(rev_alphabet[nt_code[i]]); cout<<endl;
    // for -M C:T, rev_alphabet[256] = {...,3(which is A),...,2(which is C),...,1(which is G),...,0(which is T),...}, others are all 3
    // for -M A:G, rev_alphabet[256] = {...,2(which is A),...,3(which is C),...,0(which is G),...,1(which is T),...}, others are all 2

    for(i=0;i<4;i++) useful_nt[bit_nt[i]]=nt_code[i];
    for(i=0;i<4;i++) useful_nt[bit_nt[i]+4]=tolower(nt_code[i]);
    //cout<<useful_nt<<endl;
    // for -M C:T, useful_nt = 'ACGTacgt'; for -M A:G, useful_nt = 'CAGTcagt'
}

char chain_flag[2] = {'+', '-'};
