#### libraries
library("getopt");library("plotrix");#library("dplyr");
#### arguments
#list of options:long name, short name, flag(0=no argument, 1=required argument, 2=optional argument), data type, description
OptList=c('help',       'h', 0, 'logical',   'useage',
          'epiallele',  'e', 1, 'character', 'Epiallele file. REQUIRED.',
          'Regions',    'r', 1, 'character', 'Bed file of regions in which epiallele will be plotted. REQUIRED.',
          'UseStrand',  's', 0, 'logical', 'If -s is specified, strand infomation(6th column) in Regions file will be used.',
          'eLen',       'l', 2, 'character', 'epialleles whose length(cytosine count) lower than -l will be omitted. [Default=3].',
          'eFreq',      'f', 2, 'character', 'epialleles whose frequency lower than -f will be omitted. [Default=0].',
          'DrawX',      'x', 0, 'logical', 'If -x is specified, concurrence CpGs will be plotted as red circles. By default, it only plots black(methylated) and white(unmethylated) circles.',
          'RealDist',   'd', 0, 'logical', 'If -d is specified, adjacent cytosines will be plotted in real distance. If not(default), adjacent cytosines are plotted one after another.',
          'Output',     'o', 2, 'character', 'filename prefix of output pdf files saving the epiallele plots, one region per file. [Default="Epiallele_plot"].'
);
OptList=matrix(OptList, byrow=TRUE, ncol=5)
opt=getopt(OptList);
#### check options
if(!is.null(opt$help)){
  cat(getopt(OptList, command = paste("Rscript",get_Rscript_filename()), usage=TRUE));
  q(status=1);
};
if(is.null(opt$epiallele)){
  print("[ERROR] Epiallele file (-e) is not specified.");
  q(status=1)
}else{
  epiallele=as.character(opt$epiallele)
};
if(!file.exists(epiallele)){
  print("[ERROR] Epiallele file does not exist.");
  q(status=1)
};
if(is.null(opt$Regions)){
  print("[ERROR] Regions file (-r) is not specified.");
  q(status=1)
}else{
  Regions=as.character(opt$Regions)
};
if(!file.exists(Regions)){
  print("[ERROR] Region bed file does not exist.");
  q(status=1)
};
if(!is.null(opt$UseStrand)){UseStrand=TRUE}else{UseStrand=FALSE};
if(is.null(opt$eLen)){
  eLen=3;
}else{
  eLen=as.numeric(opt$eLen)
};
if(is.null(opt$eFreq)){
  eFreq=0;
}else{
  eFreq=as.numeric(opt$eFreq)
};
if(!is.null(opt$DrawX)){DrawX=TRUE}else{DrawX=FALSE};
if(!is.null(opt$RealDist)){RealDist=TRUE}else{RealDist=FALSE};
if(is.null(opt$Output)){
  Output="Epiallele_plot";
}else{
  Output=as.character(opt$Output)
};
#### functions
Compress_CpG_Mat<-function(x,rc,gap=3){
  leftb=as.numeric(colnames(x)[1])
  rightb=as.numeric(colnames(x)[ncol(x)])
  tmp=as.numeric(colnames(x));tmp=which(tmp>=leftb & tmp<=rightb);x=x[,tmp];
  tmp=which(rowSums(x,na.rm = T)>0);x=x[tmp,];rc=rc[tmp];
  for(i0 in 1:nrow(x)){
    j0=as.numeric(which(!is.na(x[i0,])));
    k0=j0[length(j0)];
    compressed=F;
    if(i0==1){
      rc_mat=c(rc[i0],i0,k0)
    }else{
      for(i1 in 1:(i0-1)){
        j1=which(!is.na(x[i1,]));
        l=c();for(l1 in j1){l=c(l,seq(l1-gap,l1+gap))};
        j1=unique(l);j1=j1[which(j1>0)];
        if(length(intersect(j0,j1))==0){
          compressed=T;
          x[i1,j0]=x[i0,j0];x[i0,j0]=NA;
          rc_mat=rbind(rc_mat,c(rc[i0],i1,k0));
          break;
        }
      };
      if(compressed==F){
        rc_mat=rbind(rc_mat,c(rc[i0],i0,k0));
      }
    }
  };
  rc_mat=apply(rc_mat,c(1,2),as.numeric);colnames(rc_mat)=c("read_count","row_n","col_n");
  x=x[which(rowSums(x,na.rm = T)>0),];
  return(list(x=x,rc_mat=rc_mat));
};
Plot_CpG_Mat<-function(x,rc_mat,RealDist=FALSE){
  if(RealDist==FALSE){
    x_plot=1:ncol(x);
  }else{
    x_plot=as.numeric(colnames(x))-as.numeric(colnames(x)[1])+1;
  };

  plot(0,type="n",xlim=c(0,x_plot[length(x_plot)]+2),ylim=c(0,nrow(x)+1),xaxt="n",yaxt="n",bty="n",xaxs="i",yaxs="i",xlab="",ylab="")
  for(i in 1:nrow(x)){
    aY=nrow(x)+1-i;
    if(RealDist==TRUE){
      xs=which(!is.na(as.numeric(x[i,])));
      for(s in xs[1:(length(xs)-1)]){
        if((s+1)%in%xs){
          segments(x0=x_plot[s],y0=aY,x1=x_plot[s+1],y1=aY,lwd=0.5);
        }
      }
    };

    for(j in 1:ncol(x)){
      aX=x_plot[j];
      if(is.na(x[i,j])){
        next
      }else if(x[i,j]==1){
        draw.circle(x=aX,y=aY,radius=0.48,border="black",col="white",lwd=1)
      }else if(x[i,j]==2){
        draw.circle(x=aX,y=aY,radius=0.48,border="black",col="black",lwd=1)
      }else if(x[i,j]==3){
        draw.circle(x=aX,y=aY,radius=0.48,border="black",col="tomato",lwd=1)
      }
    }
  }
  if(RealDist==TRUE){
    rc_mat[,3]=as.numeric(colnames(x)[rc_mat[,3]])-as.numeric(colnames(x)[1])+1;
  };
  for(i in 1:nrow(rc_mat)){
    text(x=rc_mat[i,3]+1,y=nrow(x)+1-rc_mat[i,2],labels=rc_mat[i,1],adj=c(0,0.5),cex=0.5)
  }
  axis(3,at=x_plot,labels=colnames(x),las=2,cex.axis=0.5,tick=F,line=0);
  #mr=length(which(x==2))/length(which(!is.na(x)));
  #mtext(text=round(mr,3),side=3,at=ncol(x),adj=0.9,cex=0.1)
};
#### main process
print(paste("[",Sys.time(),"] Starting"))
if(UseStrand==FALSE){
  Regions=read.table(Regions,sep="\t",header=F)[,1:3];
  names(Regions)=c("chr","start","end");
}else{
  Regions=read.table(Regions,sep="\t",header=F)[,c(1:3,6)];
  names(Regions)=c("chr","start","end","strand");
};
#Regions=dplyr::distinct(Regions);
Regions=Regions[order(Regions$chr,Regions$start),];
print(paste("[",Sys.time(),"]",nrow(Regions),"regions imported"))

for(r in 1:nrow(Regions)){
  chr0=as.character(Regions$chr[r]);
  start0=as.character(Regions$start[r]);
  end0=as.character(Regions$end[r]);
  if(UseStrand==TRUE){strand0=as.character(Regions$strand[r])};
  c0=paste0("tabix --verbosity 1 ",epiallele," ",chr0,":",start0,"-",end0);
  c0=pipe(c0);d0=readLines(c0);close(c0);
  m0=c("strand","seq","pos","read_count");
  if(length(d0)<=1){
    print(paste0("[",Sys.time(),"] Only ",length(d0)," epiallele detected in ",chr0,ifelse(UseStrand==TRUE,strand0,""),":",start0,"-",end0,". No need to plot."));
    next;
  }else{
    e=0;
    for(m in 1:length(d0)){
      d1=unlist(strsplit(x=d0[m],split="\t",fixed=T))[4:7];
      if(UseStrand==TRUE){if(d1[1]!=strand0){next}};
      if(eLen>1){if(nchar(d1[2])<eLen){next}};
      if(eFreq>1){if(as.numeric(d1[4])<eFreq){next}};
      m0=rbind(m0,d1);e=e+1;
    };
    if(e<=1){
      print(paste0("[",Sys.time(),"] Only ",e," epiallele detected in ",chr0,ifelse(UseStrand==TRUE,strand0,""),":",start0,"-",end0,". No need to plot."));
      next;
    };
    colnames(m0)=m0[1,];m0=m0[-1,];
  };
  read_count=as.numeric(m0[,4]);

  pos_f=as.numeric(unique(unlist(strsplit(as.character(m0[m0[,1]=="+",3]),","))));
  pos_r=as.numeric(unique(unlist(strsplit(as.character(m0[m0[,1]=="-",3]),","))));
  all_pos=sort(unique(c(pos_f,pos_r-1)));
  epa_mat=matrix(NA,ncol=length(all_pos),nrow=nrow(m0));
  colnames(epa_mat)=all_pos;

  for(i in 1:nrow(m0)){
    strand_i=as.character(m0[i,1]);
    if(UseStrand==TRUE){if(strand_i!=strand0){next}};
    p0=unlist(strsplit(as.character(m0[i,3]),","));
    if(strand_i=="-"){p0=as.character(as.numeric(p0)-1)};
    # skip the read with mismatch
    q0=which(colnames(epa_mat)==p0[length(p0)])-which(colnames(epa_mat)==p0[1])+1;
    if(q0>length(p0)){next};
    c0=unlist(strsplit(as.character(m0[i,2]),""));
    if(DrawX==TRUE){
      if(("M" %in% c0) & ("U" %in% c0)){c0[which(c0=="U")]="X"};
    };
    epa_mat[i,p0]=c0;
  };
  epa_mat[which(epa_mat=="U",arr.ind = T)]=1;
  epa_mat[which(epa_mat=="M",arr.ind = T)]=2;
  if(DrawX==TRUE){epa_mat[which(epa_mat=="X",arr.ind = T)]=3};
  epa_mat=apply(epa_mat,c(1,2),as.numeric);
  
  x_comp=Compress_CpG_Mat(x=epa_mat,rc=read_count);
  x=x_comp$x;rc_mat=x_comp$rc_mat;

  if(RealDist==FALSE){
    w=ncol(x)
  }else{
    w=as.numeric(colnames(x)[ncol(x)])-as.numeric(colnames(x)[1])+1;
  }
  pdf(paste0(Output,"_",chr0,ifelse(UseStrand==TRUE,strand0,""),"_",start0,"_",end0,".pdf"),
      width=w/10,height=nrow(x)/5+1)
  par(mai=c(0,0,1,0));
  Plot_CpG_Mat(x=x,rc_mat=rc_mat,RealDist=RealDist);
  dev.off();
};
print(paste("[",Sys.time(),"] Finished"))
