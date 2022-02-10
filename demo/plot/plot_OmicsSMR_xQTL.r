# ===============================================================================
# Author: Yang Wu, Futao Zhang, Zhihong Zhu
# Date started: 21/05/2016
# Date last updated: 05/07/2018
# R script to draw regional plot and effect size plot for Omics-SMR analysis
#
#   For regional plot, users should specify the parameters such as eSMR threshold,mSMR threshold,
# m2eSMR threshold,Heidi threshold, plot window size et al. We also provide parameters esmr_thresh_plot, msmr_thresh_plot
# and probeNEARBY to draw eQTL and mQTL plots of probes of interest. We set these three parameters
# to NULL as default. Only gene expression probes passed (specified) SMR threshold and Heidi threshold would be plotted.
# Only methylation probes that associates with the target gene expression probe and also passed SMR threshold and Heidi threshold of the trait would be drawn.
#If you want to plot other probes, please use parameter interst_only.
#
#   we also can integrate Chromatin state into the regional plot. 
#   The plot file can be a little big which would be up to hundreds of megabytes.
# This causes a long time to read by this R script. In order to reduce the redundant information,
# you can use these two parameters: --psmr and --phet which are 1 and 0 as default respectively.
#
# Amendment:
#  1. In GWAS layer, we use solid and hollow rhombus to indicate gene expression probes.
#                    we use solid and hollow cycle to indicate methylation probes.
#                    we use maroon color to indicate  gene expression probes,
#                         navy color to indicate methylation probes.
#                    we use solid to indicate this probe passed HEIDI threshold,
#                        hollow to indicate this probe did not pass HEIDI threshold.
#  2. In eQTL layers, we use cross symbol with maroon color to indicate this probe passed eSMR threshold
#
#  3. In mQTL layers, we use dot symbol with navy color to indicate this probe passed mSMR threshold
#
#  3. shifted label of SMR threshold a little right in case of shading GWAS signals.
#  4. enlarged vetical axis of eQTL layer in order to prevent the label of probe name shading
#       eQTL signals.
#  5. added the funtion that users can specify probes of interest to plot.
#  6. only annotate probes passing both SMR and heidi test
#  05. 07. 2018
#  1. add plot position flag.
#  2. implement the highlight SNPs flag
#  3. annotate the significant probes only
#  4. put the trait name in the plot if given
#  example:
#        source("plot_OmicsSMR.r")
#        SMRData=ReadomicSMRData("height.ILMN_1676393.txt")
#        omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-6,msmr_thresh=1e-9,m2esmr_thresh=1e-8,window=300)
#        omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-6, esmr_thresh_plot=2e-6,msmr_thresh=1e-9,msmr_thresh_plot=1e-8,m2esmr_thresh=1e-8,window=300)
#        omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-6, esmr_thresh_plot=2e-6,msmr_thresh=1e-9,msmr_thresh_plot=1e-8,m2esmr_thresh=1e-8,window=300,max_plot_mprobe=10)
#        omicSMRLocusPlot(data=SMRData, probeNEARBY=c("ILMN_1783771", "cg11783901"),esmr_thresh=1e-6,msmr_thresh=1e-9,m2esmr_thresh=1e-8,window=300,max_plot_mprobe=10)
#        omicSMRLocusPlot(data=SMRData,probeNEARBY = c("cg27201301"),esmr_thresh=1e-6,msmr_thresh=1e-9,m2esmr_thresh=1e-8,window=300,interest_only = TRUE)
#        omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-6,msmr_thresh=1e-9,m2esmr_thresh=1e-8,window=300,anno_methyl = TRUE)
#        omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-6,msmr_thresh=1e-9,m2esmr_thresh=1e-8,window=1000,max_anno_probe = 30,anno_methyl = TRUE, anno_self = TRUE) #using library: anno_self=FALSE 
#  6.YW modified: 1. genes in isoforms 2.epigenome start
#  22.02.2019, 1. replace HEIDI NA as 1e-300; 2. add flag to show eQTL + GWAS only
# ===============================================================================

is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
}
# check if package "TeachingDemos" is installed
if (!is.installed("TeachingDemos")){
    install.packages("TeachingDemos");
}
library("TeachingDemos")

# parameters for plot
genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;


GeneRowNum = function(GENELIST) {
    BP_THRESH = 0.03; MAX_ROW = 6
    # get the start and end position
    #GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
    START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
    STRLENGTH = nchar(as.character(GENELIST$GENE))
    MIDPOINT = (START1 + END1)/2
    START2 = MIDPOINT-STRLENGTH/250; END2 = MIDPOINT+STRLENGTH/250
    START = cbind(START1, START2); END = cbind(END1, END2);
    START = apply(START, 1, min); END = apply(END, 1, max)
    GENELIST = data.frame(GENELIST, START, END)
    GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
    START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
    # get the row index for each gene
    NBUF = dim(GENELIST)[1]
    ROWINDX = rep(1, NBUF)
    ROWEND = as.numeric(rep(0, MAX_ROW))
    MOVEFLAG = as.numeric(rep(0, NBUF))
    if(NBUF>1) {
        for( k in 2 : NBUF ) {
            ITERFLAG=FALSE
            if(START[k] < END[k-1]) {
                INDXBUF=ROWINDX[k-1]+1
            } else INDXBUF = 1
            if(INDXBUF>MAX_ROW) INDXBUF=1;
            REPTIME=0
            repeat{
                if( ROWEND[INDXBUF] > START[k] ) {
                    ITERFLAG=FALSE
                    INDXBUF=INDXBUF+1
                    if(INDXBUF>MAX_ROW) INDXBUF = 1
                } else {
                    ITERFLAG=TRUE
                }
                if(ITERFLAG) break;
                REPTIME = REPTIME+1
                if(REPTIME==MAX_ROW) break;
            }
            ROWINDX[k]=INDXBUF;
            
            if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
            | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
                MOVEFLAG[k] = 1
                SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
                MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
            }
            if(ROWEND[ROWINDX[k]]<END[k]) {
                ROWEND[ROWINDX[k]] = END[k]  }
        }
    }
    GENEROW = data.frame(as.character(GENELIST$GENE),
    as.character(GENELIST$ORIENTATION),
    as.numeric(GENELIST$GENESTART),
    as.numeric(GENELIST$GENEEND),
    ROWINDX, MOVEFLAG)
    colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
    return(GENEROW)
}
GeneRowNumYW = function(GENELIST) {
    BP_THRESH = 0.03; MAX_ROW = 6
    # get the start and end position
    #GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
    START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
    STRLENGTH = nchar(as.character(GENELIST$GENE))
    MIDPOINT = (START1 + END1)/2
    START2 = MIDPOINT-STRLENGTH/250; END2 = MIDPOINT+STRLENGTH/250
    START = cbind(START1, START2); END = cbind(END1, END2);
    START = apply(START, 1, min); END = apply(END, 1, max)
    GENELIST = data.frame(GENELIST, START, END)
    GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
    START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
    # get the row index for each gene
    NBUF = dim(GENELIST)[1]
    ROWINDX = rep(1, NBUF)
    ROWEND = as.numeric(rep(0, MAX_ROW))
    MOVEFLAG = as.numeric(rep(0, NBUF))
    if(NBUF>1) {
        for( k in 2 : NBUF ) {
            ITERFLAG=FALSE
            if(START[k] < END[k-1]) {
                INDXBUF=ROWINDX[k-1]+1
            } else INDXBUF = 1
            if(INDXBUF>MAX_ROW) INDXBUF=1;
            REPTIME=0
            repeat{
                if( ROWEND[INDXBUF] > START[k] ) {
                    ITERFLAG=FALSE
                    INDXBUF=INDXBUF+1
                    if(INDXBUF>MAX_ROW) INDXBUF = 1
                } else {
                    ITERFLAG=TRUE
                }
                if(ITERFLAG) break;
                REPTIME = REPTIME+1
                if(REPTIME==MAX_ROW) break;
            }
            ROWINDX[k]=INDXBUF;
            
            if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
            | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
                MOVEFLAG[k] = 1
                SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
                MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
            }
            if(ROWEND[ROWINDX[k]]<END[k]) {
                ROWEND[ROWINDX[k]] = END[k]  }
        }
    }
    GENEROW = data.frame(as.character(GENELIST$GENE),
    as.character(GENELIST$ORIENTATION),
    as.numeric(GENELIST$GENESTART),
    as.numeric(GENELIST$GENEEND),
    ROWINDX, MOVEFLAG)
    colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
    return(GENEROW)
}

plot_probe = function(probeinfobuf, k, colplot, x.min, x.max, y.min, y.max,pchbuf,heidi) {
    #xstart = as.numeric(probeinfobuf[k,5])
    #xend = as.numeric(probeinfobuf[k,6])
    #xcenter = (xstart+xend)/2
    xcenter = as.numeric(probeinfobuf[k,3])
    pvalbuf = as.numeric(probeinfobuf[k,8])
    strbuf = probeinfobuf[k,1]    
    par(new=TRUE)
    if(heidi==TRUE) {
        plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
        xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
    } else {
        plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
        xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
    }
}
ReadomicSMRData = function(plotfile)
{
    SMRData = list();
    key=c("$eprobe","$mprobe","$m2eprobe","$SNP","$GWAS","$eQTL","$mQTL","$Gene");
    skiplines=0;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[1])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    neprobes=as.numeric(keywords[2]);
    SMRData$eprobeID=keywords[3];
  
    
    skiplines=skiplines+1;
    SMRData$eSMR=read.table(plotfile, header=F, nrows=neprobes, skip=skiplines);
    skiplines=skiplines+neprobes;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[2])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nmprobes=as.numeric(keywords[2]);
    
    skiplines=skiplines+1;
    SMRData$mSMR=read.table(plotfile, header=F, nrows=nmprobes, skip=skiplines);
    skiplines=skiplines+nmprobes;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[3])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nme2e=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$me2e=read.table(plotfile, header=F, nrows=nme2e, skip=skiplines);
    SMRData$mprobeID=unique(as.character(SMRData$me2e[,1]));
    skiplines=skiplines+nme2e;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[4])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nrs=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$SNP=read.table(plotfile, header=F, nrows=nrs, skip=skiplines);
    skiplines=skiplines+nrs;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[5])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    ngwas=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$GWAS=read.table(plotfile, header=F, nrows=ngwas, skip=skiplines);
    skiplines=skiplines+ngwas;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[6])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    neqtl=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    prbname=keywords[1];
    neqtlsnp=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$eQTL=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
    SMRData$eQTL=cbind(prbname,SMRData$eQTL)
    skiplines=skiplines+neqtlsnp;
    if(neqtl>1)
    {
        for(i in 2:neqtl)
        {
            keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
            prbname=keywords[1];
            neqtlsnp=as.numeric(keywords[2]);
            skiplines=skiplines+1;
            raweQTLtmp=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
            raweQTLtmp=cbind(prbname,raweQTLtmp);
            SMRData$eQTL=rbind(SMRData$eQTL,raweQTLtmp);
            skiplines=skiplines+neqtlsnp;
        }
    }
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[7])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nmeqtl=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    mprbname=keywords[1];
    nmeqtlsnp=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$meQTL=read.table(plotfile, header=F, nrows=nmeqtlsnp, skip=skiplines);
    SMRData$meQTL=cbind(mprbname,SMRData$meQTL)
    skiplines=skiplines+nmeqtlsnp;
    if(nmeqtl>1)
    {
        for(i in 2:nmeqtl)
        {
            keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
            mprbname=keywords[1];
            nmeqtlsnp=as.numeric(keywords[2]);
            skiplines=skiplines+1;
            rawmeQTLtmp=read.table(plotfile, header=F, nrows=nmeqtlsnp, skip=skiplines);
            rawmeQTLtmp=cbind(mprbname,rawmeQTLtmp);
            SMRData$meQTL=rbind(SMRData$meQTL,rawmeQTLtmp);
            skiplines=skiplines+nmeqtlsnp;
        }
    }
   
   keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
   if(keywords[1]!=key[8])
   {
       print("ERROR: plot file is not correct!");
       quit();
   }
   ngenes=as.numeric(keywords[2]);
   skiplines=skiplines+1;
   SMRData$Gene=read.table(plotfile, header=F, nrows=ngenes, skip=skiplines);

    return(SMRData)
}

omicSMRLocusPlot = function(data=SMRData, plotBP=NULL, eprobeNEARBY=NULL, mprobeNEARBY=NULL,esmr_thresh, esmr_thresh_plot=NULL,msmr_thresh,msmr_thresh_plot=NULL,m2esmr_thresh=1,esmr_heidi=0.01,msmr_heidi=0.01,m2esmr_heidi=0,window=500,pointsize=20,max_anno_probe=10,ASO=TRUE,max_plot_mprobe=4, epi_plot=FALSE,interest_only=FALSE,anno_methyl=FALSE,anno_self=TRUE,anno_dist=1.05,highlight=NULL,trait_name=NULL,annoSig_only=TRUE,rmDNAm=FALSE,funcAnnoFile=NULL)
{
    ######################## epigenome parts (YW) ############################
    if(epi_plot==TRUE)
    {           
        # library(data.table)
        # dir="/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/"
        # samplefile=read.table("/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/legend_YW.txt",head=F,sep="\t")
        # colfile=read.delim("/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/annotation_25_imputed12marks.txt",head=T,sep="\t")    
        # collegend=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")
        # type1=c("IMR90","ESC","iPSC","ES-deriv.","Blood & T cell","HSC & B cell","Mesench.","Myosat.","Epithelial","Neurosph.","Thymus","Brain","Adipose","Muscle","Heart","Sm. Muscle","Digestive","Other","ENCODE")
        # type2=c("","ESC","iPSC","ES-deriv","Blood & T-cell","HSC & B-cell","Mesenchymal","","Epithelial","","","Brain","","Muscle","Heart","","Digestive","Other","ENCODE")
        # colorset=c("red","brown","slateblue","royalblue","green3","green4","lightcoral","darkorange3","darkorange","gold","gold3","goldenrod4","indianred","sienna","lightsalmon3","hotpink","lightpink","grey66","black")
        # samples=samplefile[,1]
        # rectplot<-function(data,colfile,ymin,ymax) {
        #     for(i in 1:25){
        #     tmp=subset(data,data$V4==i)
        #         if(dim(tmp)[1]>0){
        #         red=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][1]);green=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][2]);blue=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][3]);
        #         rect(as.numeric(tmp$V2)/1e6,ymin,as.numeric(tmp$V3)/1e6,ymax,col=rgb(red,green,blue,max=255),border=NA)
        #         }
        #     }
        # }
        # file = list();
        # for(i in 1:length(samples)) {
        #     filename=paste(as.character(samples[i]),"_25_imputed12marks_stateno.bed.gz",sep="")
        #     filetmp=fread(paste("gunzip -c ",dir,filename,sep=""),head=F)
        #     filetmp=as.data.frame(filetmp[-dim(filetmp)[1],])
        #     file[[i]] = filetmp;
        # }
        # save.image(file = "/home/uqywu16/90days/SMRdata/funcAnno.RData")
        if(is.null(funcAnnoFile)) {
	       return("ERROR: please input .RData for epigenome annotation, which can be downloaded at.");
        } else {
            load(funcAnnoFile)
        }        
    }

    ########## checking threshold flags ############# 
    if(is.na(esmr_thresh))
    {
        return("ERROR: please input the p-value threshold of your eSMR test.");
    }
    if(is.na(msmr_thresh))
    {
        return("ERROR: please input the p-value threshold of your mSMR test.");
    }
    if(is.na(m2esmr_thresh))
    {
        return("ERROR: please input the p-value threshold of your m2eSMR test.");
    }
    cex_coeff=5/7 * pointsize/15;
    ########## checking the physical positions ##############
    if(length(which(is.na(data$eSMR[,3])))>0)
    {
        return("ERROR: Some eprobes' physical positon is missing!");
    }
    if(length(which(is.na(data$mSMR[,3])))>0)
    {
        return("ERROR: Some mprobes' physical positon is missing!");
    }
    ########## checking the plot threshold #############
    if(length(esmr_thresh_plot)==0){
        esmr_thresh_plot=esmr_thresh;
    }
    if(length(msmr_thresh_plot)==0){
        msmr_thresh_plot=msmr_thresh;
    }
    idx=match(data$eprobeID,data$eSMR[,1]);
    #data$mprobeID from m2eSMR result
    idx2=match(data$mprobeID,data$mSMR[,1]);
    if(length(idx)==0 | length(idx2)==0){
        return("ERROR: Plot file is not generated correctly, can't find target probe!");
    }
    ######### set the plot center and bp #############
    if(length(plotBP)!=0) {
        plot_start=plotBP-window*1000
        if(plot_start<0) plot_start=0
        plot_end=plotBP+window*1000;
    } else {
        plot_start=data$eSMR[idx,3]-window*1000
        if(plot_start<0) plot_start=0
        plot_end=data$eSMR[idx,3]+window*1000;
    }
    ######### set the plot window for different data  #############
    idx=which(data$eSMR[,3]>=plot_start & data$eSMR[,3]<=plot_end)
    data$eSMR=data$eSMR[idx,]
    idx=which(data$mSMR[,3]>=plot_start & data$mSMR[,3]<=plot_end)
    data$mSMR=data$mSMR[idx,]
    idx=which(data$me2e[,3]>=plot_start & data$me2e[,3]<=plot_end)
    data$me2e=data$me2e[idx,]
    idx=match(data$GWAS[,1],data$SNP[,1])
    tmpsnpbp=data$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data$GWAS=data$GWAS[idx,]
    idx=match(data$eQTL[,2],data$SNP[,1])
    tmpsnpbp=data$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data$eQTL=data$eQTL[idx,]
    idx=match(data$meQTL[,2],data$SNP[,1])
    tmpsnpbp=data$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data$meQTL=data$meQTL[idx,]
    data$mSMR=data$mSMR[order(data$mSMR[,8]),]
    idx=which(data$Gene[,2]>=plot_start & data$Gene[,3]<=plot_end )
    data$Gene=data$Gene[idx,]
    data$Gene=unique(data$Gene)

    ########## replace the missing HEIDI pvalue as 1e-300 ###########
    data$eSMR$V9[is.na(data$eSMR[,9])] = 1e-300
    data$mSMR$V9[is.na(data$mSMR[,9])] = 1e-300
    data$me2e$V9[is.na(data$me2e[,9])] = 1e-300

    ########## filtering eprobes with eSMR & HEIDI threshold ##########
    esmrindx_plot = which(data$eSMR[,8] <= esmr_thresh_plot)
    if(length(esmrindx_plot)>0) { esmrprobes_plot =  as.character(data$eSMR[esmrindx_plot,1]) }
    
    esmrindx = which(data$eSMR[,8] <= esmr_thresh);
    eheidiindx = which(data$eSMR[,9] >= esmr_heidi);
    esmrprobes = NA; eheidiprobes = NA;
    if(length(esmrindx)>0) { 
        esmrprobes =  as.character(data$eSMR[esmrindx,1]) 
    } else { 
        return("ERROR: No expression probe pass the eSMR threshold.");
    }
    if(length(eheidiindx)>0) { 
        eheidiprobes = as.character(data$eSMR[eheidiindx,1]) 
    } else { 
        return("ERROR: No expression probe further pass the eSMR HEIDI threshold.");
    } 
    esmrprobes_plot=intersect(esmrprobes_plot,eheidiprobes) # only plot the probes passed (specified) smr threshold and heidi threshold
    
    ########## filtering mprobes with mSMR & HEIDI threshold ##############
    msmrindx_plot = which(data$mSMR[,8] <= msmr_thresh_plot)
    if(length(msmrindx_plot)>0) { msmrprobes_plot =  as.character(data$mSMR[msmrindx_plot,1]) }
    
    msmrindx = which(data$mSMR[,8] <= msmr_thresh)
    mheidiindx = which(data$mSMR[,9] >= msmr_heidi)
    msmrprobes = NA; mheidiprobes = NA;
    if(length(msmrindx)>0) { 
        msmrprobes =  as.character(data$mSMR[msmrindx,1]) 
    } else { 
        return("ERROR: No methylation probe pass the mSMR threshold.");
    }
    if(length(mheidiindx)>0) { 
        mheidiprobes = as.character(data$mSMR[mheidiindx,1]) 
    } else { 
        return("ERROR: No methylation probe further pass the mSMR HEIDI threshold.");
    }
    msmrprobes_plot=intersect(msmrprobes_plot,mheidiprobes) # only plot the mprobes passed (specified) msmr threshold and mheidi threshold
    
    ########## filtering mprobes with m2e SMR & HEIDI threshold ##############
    m2esmrindx = which(data$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx = which(data$me2e[,8] <= m2esmr_thresh & data$me2e[,9] >= m2esmr_heidi)
    m2esmrprobes = NA; m2eheidiprobes = NA;
    if(length(m2esmrindx)>0) { 
        m2esmrprobes =  as.character(data$me2e[m2esmrindx,1]) 
    } else { 
        return("ERROR: No methylation and expression association pass the m2eSMR threshold.");
    }
    if(length(m2eheidiindx)>0) { 
        m2eheidiprobes = as.character(data$me2e[m2eheidiindx,1]) 
    } else { 
        return("ERROR: No methylation and expression association further pass the m2eSMR HEIDI threshold.");
    }
    msmrprobes_plot=intersect(msmrprobes_plot,m2eheidiprobes) # only plot the mprobes passed m2esmr threshold and m2eheidi threshold
    
    ########## only select the probes of interest ##########
    if(length(eprobeNEARBY)>0)
    {
        idx=match(eprobeNEARBY,data$eSMR[,1])
        idxx=which(is.na(idx))
        if(length(idxx)>0)
        {
            eprobeNEARBY=eprobeNEARBY[-idxx]
        } else {
            eprobeNEARBY=eprobeNEARBY
        }
        eprobePLOT=data$eprobeID
        eprobePLOT=unique(c(eprobePLOT,eprobeNEARBY))
    }   else {
        eprobePLOT=esmrprobes_plot
    }

    if(length(mprobeNEARBY)>0)
    {
        idx=match(mprobeNEARBY,data$mSMR[,1])
        idxx=which(is.na(idx))
        if(length(idxx)>0)
        {
            mprobeNEARBY=mprobeNEARBY[-idxx]
        } else {
            mprobeNEARBY=mprobeNEARBY
        }
        mprobePLOT=mprobeNEARBY       
    }   else {
        mprobePLOT=msmrprobes_plot
    }

    neprobePLOT = length(eprobePLOT)  
    nmprobePLOT = length(mprobePLOT)
    if (neprobePLOT==0) {
        return("Warnings: Specified expression probes were not found in the extracted data.");
        eprobePLOT=esmrprobes_plot       
    }
    if (nmprobePLOT==0) {
        return("Warnings: Specified methylation probes were not found in the extracted data.");
        mprobePLOT=msmrprobes_plot 
    }
    num_mprb_plot=min(length(mprobePLOT), max_plot_mprobe)

    ########## flag to determine whether plot the DNAm ########
    if (rmDNAm == T) num_mprb_plot = 0;
    
	idx=which(is.na(data$GWAS[,2]) | is.na(data$GWAS[,3]))
    if(length(idx)>0) data$GWAS=data$GWAS[-idx,]
    ########## highlight the specified SNPs ############
    if(length(highlight)>0) {
        highlightindx=which(data$GWAS[,1]==highlight)
        if (length(highlightindx)==0) {
            return("Warnings: Specified highlight SNP were not found in the GWAS data.");
        }
    } else {
        highlightindx=NULL
    }

    ########## replace the infinite value #############
	pZY=-log10(pchisq((data$GWAS[,2]/data$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY)))!=0) {
        maxpZY = max(pZY[-which(is.infinite(pZY))])
        if (maxpZY > 320) {pZY[which(is.infinite(pZY))] = maxpZY + 5
            } else {pZY[which(is.infinite(pZY))] = 325}
    }

    ########## find the plot target expression probe #############
    idx=match(data$eprobeID,data$eSMR[,1]);
    if(length(idx)>0){
        chrPLOT = data$eSMR[idx,2]
    } else {
        print("ERROR: Plot file is not generated correctly, please report this bug!");
        quit();
    }
    
    ########## remove the unannotated eprobes  #############
    idx=which(is.na(data$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO=data$eSMR[-idx,];
    }else{
        eprobeINFO=data$eSMR;
    }
    idx=which(is.na(data$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO=data$mSMR[-idx,];
    }else{
        mprobeINFO=data$mSMR;
    }

    idx=which(is.na(eprobeINFO[,5]) | is.na(eprobeINFO[,6]));
    idx2=which(is.na(eprobeINFO[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO[idx,5]=eprobeINFO[idx,3]-7500;
    eprobeINFO[idx,6]=eprobeINFO[idx,3]+7500;
    eprobeINFO[,8]=-log10(eprobeINFO[,8]);
    eprobeINFO[,3]=eprobeINFO[,3]/1e6;
    eprobeINFO[,5]=eprobeINFO[,5]/1e6;
    eprobeINFO[,6]=eprobeINFO[,6]/1e6;
    
    ########## remove the unannotated mprobes  #############
    idx=which(is.na(mprobeINFO[,5]) | is.na(mprobeINFO[,6]));
    idx2=which(is.na(mprobeINFO[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO[idx,5]=mprobeINFO[idx,3]-7500;
    mprobeINFO[idx,6]=mprobeINFO[idx,3]+7500;
    mprobeINFO[,8]=-log10(mprobeINFO[,8]);
    mprobeINFO[,3]=mprobeINFO[,3]/1e6;
    mprobeINFO[,5]=mprobeINFO[,5]/1e6;
    mprobeINFO[,6]=mprobeINFO[,6]/1e6;
    
    ##########  plot start here  ##########
    epXY=eprobeINFO[,8];
    mpXY=mprobeINFO[,8];
    yMAX = ceiling(max(c(pZY, epXY,mpXY), na.rm=T)) + 1;
    #glist=cbind(eprobeINFO[,2],eprobeINFO[,5:6],as.character(eprobeINFO[,4]),eprobeINFO[,7]);#old
    glist=data$Gene;
    glist[,2]=glist[,2]/1e6;
    glist[,3]=glist[,3]/1e6;
    colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
    idx=which(is.na(glist[,2]) | is.na(glist[,3]));
    if(length(idx>0)) glist=glist[-idx,];
    generow = GeneRowNum(glist);
    num_row = max(as.numeric(generow$ROW));
    offset_map = ceiling(yMAX);
    
    offset_map=max(offset_map,num_row*2.5)
    
    offset_probe = yMAX / 2.5;
   
    offset_eqtl = ceiling(yMAX / 2.5) + 0.5;
    dev_axis = 0.2*yMAX;# control distance between panles
    if(dev_axis<1.5) dev_axis = 1.5;
    if(epi_plot==TRUE) {
        yaxis.min = -offset_map - dev_axis - offset_eqtl*neprobePLOT - offset_eqtl*num_mprb_plot- dev_axis*(neprobePLOT)- dev_axis*(num_mprb_plot)-offset_eqtl*6;
    } else {
        yaxis.min = -offset_map - dev_axis - offset_eqtl*neprobePLOT - offset_eqtl*num_mprb_plot- dev_axis*(neprobePLOT)- dev_axis*(num_mprb_plot);
    }
    yaxis.max = yMAX + ceiling(offset_probe) + 1;
    
    # scales of x-axis
    idx=match(data$GWAS[,1],data$SNP[,1]);
    gwasBP = as.numeric(data$SNP[idx,3])/1e6;
    min.pos = min(gwasBP);
    max.pos = max(gwasBP);
    start = min(as.numeric(glist[,2]));
    end = max(as.numeric(glist[,3]));
    bp = c(min.pos, max.pos, start, end);
    xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
    xmax=xmax+(xmax-xmin)*0.15 #extend
  
    ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
    xlab = paste("Chromosome", chrPLOT, "Mb");
    # plot GWAS p value
    par(mar=c(5,5,3,2), xpd=TRUE)
    plot(gwasBP, pZY, yaxt="n",xaxt="n",bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    if(length(highlightindx)>0) points(gwasBP[highlightindx], pZY[highlightindx],,col="peru",pch="*",cex=2)
    # x axis
    #axis(1, at=seq(xmin,xmax,(xmax-xmin)/5), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
    axis(1, at=round(seq(xmin,xmax,(xmax-xmin)/5),2))
    # y1 axis
    devbuf1 = yMAX/4
    axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
    mtext(ylab, side=2, line=3, at=(yMAX*4/7), cex=cex_coeff);
    if (length(trait_name)>0) text(xmin, max(pZY,na.rm=T) , label=trait_name,col="blue", cex=1, adj=0)
    
    eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
    axis.start = 0; axis.down = offset_eqtl + dev_axis;
    for( k in 1 : neprobePLOT ) {
        axis.start = axis.start - axis.down
        eqtlinfobuf = data$eQTL[which(data$eQTL[,1]==eprobePLOT[k]),]
        if(dim(eqtlinfobuf)[1]==0) next;
        pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
        if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
        }

        if(length(which(esmrprobes==eprobePLOT[k]))==0) { #probes of eQTL are maroon
            col_eqtl = "maroon"
        } else col_eqtl = "maroon"
        eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
        eqtl.max =ceiling(eqtl.max *1.5) #extend
        pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
        idx=match(eqtlinfobuf[,2],data$SNP[,1]);
        eqtlbp = as.numeric(data$SNP[idx,3])/1e6;
        probegene = unique(as.character(data$eSMR[which(data$eSMR[,1]==eprobePLOT[k]),4]))
        par(new=TRUE)
        pchbuf = 4;
        #if(k%%2==0) pchbuf = 20;
        plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the eQTLs
        colme2e="black"
        #if(eprobePLOT[k]==data$eprobeID) colme2e="darkorchid2"
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=eprobePLOT[k], geneid=probegene)),col=colme2e, cex=1, adj=0)
        if(k==1) labpos=axis.start+offset_eqtl-dev_axis/2;
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,eqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    ypos = labpos-abs(labpos - axis.start)/12;
    mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff,adj=1)
    
    #plot meQTL
    if (num_mprb_plot>0 & rmDNAm == F) {
        idx=which(data$me2e[,8]<=m2esmr_thresh & data$me2e[,9]>=m2esmr_heidi)
        me2e_probes=c()
        if(length(idx>0)) me2e_probes=as.character(data$me2e[idx,1])
        
        meqtl.lab = expression(paste("-", log[10], "(", italic(P), " xQTL)", sep=""));
        for( k in 1 : num_mprb_plot) {
            axis.start = axis.start - axis.down
            meqtlinfobuf = data$meQTL[which(data$meQTL[,1]==mprobePLOT[k]),]
            if(dim(meqtlinfobuf)[1]==0) next;
            pvalbuf=-log10(pchisq((meqtlinfobuf[,3]/meqtlinfobuf[,4])^2,1,lower.tail=F));
            if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
            }

            if(length(which(msmrprobes==mprobePLOT[k]))==0) {
                col_eqtl = "navy"
            } else col_eqtl = "navy"
            meqtl.min = 0; meqtl.max = ceiling(max(pvalbuf))
            meqtl.max =ceiling(meqtl.max *1.5) #extend
            pvalbuf = pvalbuf/meqtl.max * offset_eqtl + axis.start
            idx=match(meqtlinfobuf[,2],data$SNP[,1]);
            meqtlbp = as.numeric(data$SNP[idx,3])/1e6;
            mprobegene = unique(as.character(data$mSMR[which(data$mSMR[,1]==mprobePLOT[k]),4]))
            par(new=TRUE)
            pchbuf = 20;
            #if(k%%2==0) pchbuf = 4;
            plot(meqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
            ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
            # annotate the meQTLs
            colme2e="black"
            #if(length(me2e_probes)>0)
            #{
                #if(length(which(me2e_probes==mprobePLOT[k]))) colme2e="darkorchid2"
            #}
            #text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=mprobePLOT[k], geneid=mprobegene)),col=colme2e, cex=1, adj=0)
            text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, sep=""),list(probeid=mprobePLOT[k], geneid=mprobegene)),col=colme2e, cex=1, adj=0)
            # axis
            devbuf1 = offset_eqtl/3; devbuf2 = meqtl.max/3
            axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
            labels=round(seq(0,meqtl.max,devbuf2),0),
            las=1, cex.axis=axis)
            # add separator line
            segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
            col="dim grey", lty="24", lwd=1)
        }
        segments(xmin, axis.start-dev_axis/2, xmax, axis.start-dev_axis/2,col="dim grey", lty="24", lwd=1)
        ypos=(axis.start - dev_axis - axis.down)*16/24
        mtext(meqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)
    }
    
    
    if(epi_plot==TRUE)
    {
        ########### plot the 127 chromtin states samples #################
        plot_start=xmin*1e6
        axis.start=axis.start-dev_axis
        epi.start=color.start=axis.start
        regionchr=paste("chr",chrPLOT,sep="")
        unit=offset_eqtl*6/130
        dev=0.02*unit
        axis.down=unit+dev
        for(i in 1:length(samples)){            
            regionindx=which(file[[i]]$V1==regionchr & as.numeric(file[[i]]$V2)>=plot_start & as.numeric(file[[i]]$V3)<=plot_end)
            plotfile=file[[i]][regionindx,]
            par(new=T)
            rectplot(plotfile,colfile,axis.start-unit,axis.start)
            lines(c(plot_start/1e6,plot_end/1e6),c(axis.start,axis.start),col="dim grey", lty="24", lwd=1)
            lines(c(plot_start/1e6,plot_end/1e6),c(axis.start-unit,axis.start-unit),col="dim grey", lty="24", lwd=1)
            axis.start=axis.start-axis.down
        }

        ########################## annotation of cell types ###################
        labwidth=(plot_end/1e6-plot_start/1e6)/60
        for(i in 1:length(unique(samplefile[,3]))){
            dup=length(which(samplefile[,3]==i))
            epi.down=epi.start-dup*(unit+dev)
            par(new=T)
            rect(plot_start/1e6-labwidth,epi.start,plot_start/1e6,epi.down,col=colorset[i],border=NA)
            par(new=T)
            #rect(plot_end/1e6,epi.start,plot_end/1e6+labwidth,epi.down,col=colorset[i],border=NA)
            epiylab=type2[i]
            labpos=(epi.down+epi.start)/2
            mtext(epiylab, side=2, line=0, at=labpos,cex=axis*0.8,las=1);
            epi.start=epi.down-dev
        }

        ########################## color legend ###################
        labwidth=(plot_end/1e6-plot_start/1e6)/60
        unit=(offset_eqtl*6)/length(collegend)
        legendcol=unique(colfile[,5])
        for(i in c(1:7,14,8:13)) {
            color.down=color.start-unit
            red=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][1]);green=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][2]);blue=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][3]);
            par(new=T)
            rect(plot_end/1e6+4*labwidth,color.start,plot_end/1e6+5*labwidth,color.down,col=rgb(red,green,blue,max=255),border=NA)
            colorylab=collegend[i]
            labpos=(color.down+color.start)/2
            #mtext(collegend[i], side=4, line=-5, at=labpos,cex=0.8,las=1);
            text(plot_end/1e6+6*labwidth, labpos, labels=collegend[i], adj=0, cex=axis*0.8);
            color.start=color.down-dev
            }
    }
    ########################## plot the gene labels #######################
    # plot p value of bTG
    # all the probes
    num_gene = dim(generow)[1]
    dist = offset_map/num_row
    for( k in 1 : num_row ) {
        generowbuf = generow[which(as.numeric(generow[,5])==k),]
        xstart = as.numeric(generowbuf[,3])
        xend = as.numeric(generowbuf[,4])
        snbuf = which(xend-xstart< 1e-3)
        if(length(snbuf)>0) {
            xstart[snbuf] = xstart[snbuf] - 0.0025
            xend[snbuf] = xend[snbuf] + 0.0025
        }
        xcenter = (xstart+xend)/2
        xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
        num_genebuf = dim(generowbuf)[1]
        for( l in 1 : num_genebuf ) {
            ofs=0.3
            if(l%%2==0) ofs=-0.8
            m = num_row - k
            ypos = m*dist + yaxis.min
            code = 1
            if(generowbuf[l,2]=="+") code = 2;
            arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
            col=colors()[75], lwd=1)
            movebuf = as.numeric(generowbuf[l,6])*genemove
            text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.8)
        }    
    }

    # plot the probes
    eprobeINFO=eprobeINFO[order(eprobeINFO[,8],decreasing = TRUE),];
    neprobeINFO=dim(eprobeINFO)[1];
    mprobeINFO=mprobeINFO[order(mprobeINFO[,8],decreasing = TRUE),];
    nmprobeINFO=dim(mprobeINFO)[1];
    
    ########## flag to annotate DNAm probes on the top layer ###########
    if (rmDNAm == T) anno_methyl=FALSE
    ########## flag to plot DNAm probes on the top layer ###########
    if (rmDNAm == F) {
        for( k in 1 : nmprobeINFO)  {
            hitflag=FALSE
            if(length(which(mheidiprobes==mprobeINFO[k,1]))>0 & length(which(msmrprobes==mprobeINFO[k,1]))>0) {
                hitflag=TRUE
                colplot = "navy"; colfont=2; pchbuf=21
            } else if( length(which(msmrprobes==mprobeINFO[k,1]))>0) {
                colplot = "navy"; colfont=2; pchbuf=1
            }else {
                colplot = "navy"; colfont=1; pchbuf=1
            }
            if( as.numeric(mprobeINFO[k,8]) < 0 ) {
                colplot = "black"; colfont=1;
            }
            # plot p value of bxy            
            plot_probe(mprobeINFO, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
            
        }
    }

    for( k in 1 : neprobeINFO)  {
        hitflag=FALSE
        if(length(which(eheidiprobes==eprobeINFO[k,1]))>0 & length(which(esmrprobes==eprobeINFO[k,1]))>0)
        {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else if(length(which(esmrprobes==eprobeINFO[k,1]))>0) {
            colplot = "maroon"; colfont=2; pchbuf=5
        } else {
            colplot = "maroon"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobeINFO[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        plot_probe(eprobeINFO, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
    }
    

    if (annoSig_only==T) {
    eannoidx = which(eprobeINFO$V8 >= -log10(esmr_thresh) & eprobeINFO$V9 >= esmr_heidi)
    eprobeINFO = eprobeINFO[eannoidx,]; neprobeINFO=dim(eprobeINFO)[1];
    mannoidx = which(mprobeINFO$V8 >= -log10(msmr_thresh) & mprobeINFO$V9 >= msmr_heidi)
    mprobeINFO = mprobeINFO[mannoidx,]; nmprobeINFO=dim(mprobeINFO)[1]; }
     #anno probe
    if(anno_methyl==TRUE) {tmp_anno_probe=ceiling(max_anno_probe/2)
        } else { tmp_anno_probe = max_anno_probe }
    if(neprobeINFO<=tmp_anno_probe) {
        probeINFO=eprobeINFO
    } else {
        probeINFO=eprobeINFO[c(1:tmp_anno_probe),]
        idx=match(eprobePLOT,eprobeINFO[,1])
        remainidx=setdiff(c(1:tmp_anno_probe),idx)
        #anno all the probes drawn in eQTL layer
        if(length(remainidx)>0)
        {
            hit=0
            for(k in 1:length(idx))
            {
                if(idx[k]>tmp_anno_probe) {
                    hit=hit+1
                    probeINFO[remainidx[hit],]=eprobeINFO[idx[k],]
                }
            }
        }
    }
   if(ASO)
    {
        idx=which(probeINFO[,8] >= -log10(esmr_thresh))
        if(length(idx)>0) probeINFO=probeINFO[idx,]
        else {
            print("No Significant gene expression probes found for SMR test!")
            probeINFO=c()
        }
    }
    if(anno_methyl==TRUE){ tmp_anno_probe=dim(probeINFO)[1]; tmp_anno_probe=max_anno_probe-tmp_anno_probe
        if(nmprobeINFO<=tmp_anno_probe) {
            probeINFOm=mprobeINFO
        } else {
            probeINFOm=mprobeINFO[c(1:tmp_anno_probe),]
            idx=match(mprobePLOT[1:num_mprb_plot],mprobeINFO[,1])
            remainidx=setdiff(c(1:tmp_anno_probe),idx)
            #anno all the probes drawn in eQTL layer
            if(length(remainidx)>0)
            {
                hit=0
                for(k in 1:length(idx))
                {
                    if(idx[k]>tmp_anno_probe) {
                        hit=hit+1
                        probeINFOm[remainidx[hit],]=mprobeINFO[idx[k],]
                    }
                }
            }
        }
        if(ASO)
        {
            idx=which(probeINFOm[,8] >= -log10(msmr_thresh))
            if(length(idx)>0) probeINFOm=probeINFOm[idx,]
            else {
                print("No Significant methylation probes found for SMR test!")
                probeINFOm=c()
            }
        }
        if(length(probeINFO)==0 & length(probeINFOm)==0) {
             print("No Significant methylation/gene expression probes found for SMR test! Please disbale ASO then rerun the plot.")
             quit()
        }
        probeINFO=rbind(probeINFO,probeINFOm)   
    }
    ################### 2018.04.25 #######################
    if(length(highlightindx)>0) { 
        tmpinfo = data.frame(matrix(NA,length(highlightindx),dim(probeINFO)[2]))
        tmpinfo[,1] = highlight[which(highlightindx!=0)]; tmpinfo[,2] = chrPLOT; tmpinfo[,3] = gwasBP[highlightindx]; tmpinfo[,8] = pZY[highlightindx]
        colnames(tmpinfo) = colnames(probeINFO); probeINFO = rbind(probeINFO, tmpinfo)
    }
    ################### 2018.04.25 #######################
    if(anno_self) probeINFO=probeINFO[order(probeINFO[2],probeINFO[3]),] ####20170217
    num_anno_probe=dim(probeINFO)[1]
    xcenter = as.numeric(probeINFO[,3])
    xcbuf = xcenter
    ####20170217####
    if(anno_self)
    {
        reginlength=(xmax-(xmax-xmin)*0.15)-xmin
        leftspot=xmin+reginlength/20
        rightspot=(xmax-(xmax-xmin)*0.15)-reginlength/20
        itvl=(rightspot-leftspot)/dim(probeINFO)[1]
        if(dim(probeINFO)[1]==1) {
            xcenter=as.numeric(probeINFO[,3])
        } else {
            xcenter=leftspot+itvl/2
            for( k in 2:dim(probeINFO)[1]) xcenter=c(xcenter,leftspot+k*itvl)
        }
    } else {
        xcenter = spread.labs(xcenter[1:num_anno_probe], mindiff=0.08, maxiter=1000, min = xmin, max = (xmax-(xmax-xmin)*0.15))
        # adjust the line position
    }
    ########
    
    adjflag = rep(0, num_anno_probe)
    if(num_anno_probe>1) {
        dbuf = c(0, xcbuf[1:(num_anno_probe-1)])
        mflag = as.numeric(abs(xcbuf[1:(num_anno_probe)] - dbuf) < 0.01)
        adjflag = as.numeric( mflag | c(mflag[2:num_anno_probe],0) )
    }
    
    for( k in 1 : num_anno_probe)  {
      
        # annotate the probes
            hitflag=FALSE
            if(length(which(eheidiprobes==probeINFO[k,1]))>0 & length(which(esmrprobes==probeINFO[k,1]))>0)
            {
                hitflag=TRUE
                colplot = "maroon"; colfont=2; pchbuf=23
            } else if(length(which(esmrprobes==probeINFO[k,1]))>0) {
                colplot = "maroon"; colfont=2; pchbuf=5
            } else if (length(which(eprobeINFO[,1]==as.character(probeINFO[k,1])))>0) {
                colplot = "maroon"; colfont=1; pchbuf=5
            } else if(length(which(mheidiprobes==probeINFO[k,1]))>0 & length(which(msmrprobes==probeINFO[k,1]))>0) {
                hitflag=TRUE
                colplot = "navy"; colfont=2; pchbuf=23
            } else if(length(which(msmrprobes==probeINFO[k,1]))>0) {
                colplot = "navy"; colfont=2; pchbuf=5
            } else {
                colplot = "navy"; colfont=1; pchbuf=5
            }
            if( as.numeric(probeINFO[k,8]) < 0 ) {
                colplot = "black"; colfont=1;
            }
            if( as.numeric(probeINFO[k,8]) < 0 ) {
                colplot = "black"; colfont=1;
            }
            if(length(which(highlight==probeINFO[k,1]))>0){
                hitflag=TRUE; colplot = "peru"; colfont=2; pchbuf=23
            }

            ypos = anno_dist*yMAX
            if(colplot == "maroon")strbuf=text(xcenter[k],ypos,labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),list(probeid=as.character(probeINFO[k,1]),genename=as.character(probeINFO[k,4]))),ylim=c(yaxis.min, yaxis.max),srt=30, col=colplot, font=colfont, cex=1, adj=0)
	        if(colplot == "navy") strbuf=text(xcenter[k],ypos,labels=substitute(paste(probeid),list(probeid=as.character(probeINFO[k,1]))),ylim=c(yaxis.min, yaxis.max),srt=30, col=colplot, font=colfont, cex=1, adj=0)
            if(colplot == "peru") strbuf=text(xcenter[k],ypos,labels=substitute(paste(probeid),list(probeid=as.character(probeINFO[k,1]))),ylim=c(yaxis.min, yaxis.max),srt=30, col=colplot, font=colfont, cex=1, adj=0)
  	        # plot the lines
            # 1st step
            xstart = xcbuf[k]
            ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
            if( num_anno_probe > 1 ) {
                if(adjflag[k]==1) {
                    xstart = (xcbuf[k] + xcenter[k])/2
                    segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
                }
            }
            segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
            # 2nd step
            xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*anno_dist;
            segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
        
    }
    
    # plot the threshold
    # eSMR threshold
    yebuf = -log10(as.numeric(esmr_thresh)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh)); dev_anno = yMAX/9;
    

    strbuf = paste("pESMR = ",esmr_thresh, sep="")
    segments(xmin, yebuf, xmax, yebuf, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno
    } else {
        ythrespos=yebuf+dev_anno
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=0.9,font=3);
    
    if (rmDNAm == F) {
        strbuf = paste("pMSMR = ",msmr_thresh, sep="")
        segments(xmin, ymbuf, xmax, ymbuf, col="navy", lty=3, lwd=1);
        if(yebuf<ymbuf) {
            ythrespos=ymbuf+dev_anno
        } else {
            ythrespos=ymbuf-dev_anno
        }
        text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=0.9,font=3);
    }

}


omicSMREffectPlot = function(data=SMRData, exposure_probe,cisWindow=2000, pointsize=20)
{
    if(is.na(exposure_probe))
    {
        print("ERROR: please inpute the exposure probe name.");
        quit();
    }
    # parameters for plot
    pch_top = 24; pch_cis = 21; pch_trans = 22
    col_top = "red"; col_cis = "Navy"; col_trans = "green"
    cex_coeff=3/4 * pointsize/15;
    # Extract the mprobe for plot
    snbuf = which(as.character(data$meQTL[,1])==exposure_probe)
    if(length(snbuf)==0) {
        print(paste("ERROR: no eQTL infomation found for probe",data$probeID,"!",sep=""));
        quit();
    }
    mplotData = data$meQTL[snbuf,]
    idx=which(is.na(mplotData[,5]))
    if(length(idx)>0) mplotData=mplotData[-idx,]
    if(dim(mplotData)[1]==0){
        print("ERROR: this mprobe is not significant with eprobe!");
        quit();
    }


    # Extract the probe for plot
    snbuf = which(as.character(data$eQTL[,1])==data$eprobeID)
    if(length(snbuf)==0) {
        print(paste("ERROR: no eQTL infomation found for probe",data$probeID,"!",sep=""));
        quit();
    }
    plotData = data$eQTL[snbuf,]
   
    
    # SNPs in common
    snpbuf = Reduce(intersect, list(as.character(plotData[,2]), mplotData[,2]))
    plotData = plotData[match(snpbuf, as.character(plotData[,2])),]
    mplotData = mplotData[match(snpbuf, as.character(mplotData[,2])),]
    # Effect size
    snplist = as.character(plotData[,2])
    bZX = as.numeric(as.character(mplotData[,3]));
    seZX = as.numeric(as.character(mplotData[,4]));
    snpCorr=as.numeric(as.character(mplotData[,5]));
    bZY = as.numeric(as.character(plotData[,3]));
    seZY = as.numeric(as.character(plotData[,4]));
    # Limit
    xmin =  min(bZX - seZX, na.rm=T)
    xmax =  max(bZX + seZX, na.rm=T)
    ymin =  min(bZY - seZY, na.rm=T)
    ymax =  max(bZY + seZY, na.rm=T)
    
    if(xmin>0) xmin = -xmax/2
    if(xmax<0) xmax = -xmin/2
    if(ymin>0) ymin = -ymax/2
    if(ymax<0) ymax = -ymin/2
    
    # Plots
    par(mar=c(5,6.5,5,2), xpd=FALSE)
    # Start to plot
    nsnps = dim(plotData)[1]
       # Plot the cis-eQTL
    zZX = bZX/seZX;
    zZY= bZY/seZY;
    maxid = which.max(zZX^2)
    maxsnp = snplist[maxid]
    snpCorr = snpCorr^2;
    for( k in 1 : nsnps ) {
        # effect sizes
        colbuf = rgb(0, 0, 128/255, snpCorr[k])
        colcir = col_cis;
        cex = 1
        plot(bZX[k], bZY[k], pch=pch_cis, col=colcir, bg=colbuf,
        bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
        par(new=TRUE)
        plot(bZX[k], bZY[k], pch=20, col=colcir, bg=colbuf,
        bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        cex=0.1, xlab="", ylab="", xaxt="n", yaxt="n")
        par(new=TRUE)
    }
    
    # standard error
    # cis eQTL
    for( k in 1 : nsnps ) {
        colcir = col_cis
        segments(bZX[k]-seZX[k], bZY[k], bZX[k]+seZX[k], bZY[k],
        col=colcir, lwd=0.5+snpCorr[k])
        segments(bZX[k], bZY[k]-seZY[k], bZX[k], bZY[k]+seZY[k],
        col=colcir, lwd=0.5+snpCorr[k])
    }
    
    # line
    colline = rgb(244/255,164/255,96/255,1)
    bXY = bZY[maxid]/bZX[maxid]
    abline(0, bXY, col=colline, lwd=2, lty=2)
    
    # plot effect size of the top SNP
    colbuf = "white"
    colcir = col_top
    cex=2.3
    par(new=TRUE)
    plot(bZX[maxid], bZY[maxid], pch=pch_top, col=colcir, bg=colbuf,
    bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
    cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
    colbuf = col_top
    colcir = col_top
    cex = 1
    par(new=TRUE)
    plot(bZX[maxid], bZY[maxid], pch=pch_top, col=colcir, bg=colbuf,
    bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
    cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
    
    # se of the top SNP
    colcir = rgb(1,0,0)
    segments(bZX[maxid]-seZX[maxid], bZY[maxid],
    bZX[maxid]+seZX[maxid], bZY[maxid],
    col=colcir, lwd=1.5)
    segments(bZX[maxid], bZY[maxid]-seZY[maxid],
    bZX[maxid], bZY[maxid]+seZY[maxid],
    col=colcir, lwd=1.5)
    
    
    # plot the axis
    # x axis
    devbuf = (xmax - xmin)/5
    if(xmax!=0 & xmin!=0) {
        numbuf = min(abs(xmin), abs(xmax))
        if( devbuf > numbuf ) devbuf = numbuf
    }
    numbuf = as.numeric()
    if( xmin < 0 ) numbuf = c(numbuf, -seq(0, abs(xmin), devbuf))
    if( xmax > 0 ) numbuf = c(numbuf, seq(0, xmax, devbuf))
    axis(1, at=numbuf, labels=round(numbuf,2), las=1, cex.axis=axis)
    xmid = (xmax+xmin)/2
    mtext(paste("meQTL effect sizes of ",exposure_probe,sep=""), side=1, at=xmid, line=3, cex=cex_coeff)
    # y axis
    devbuf = (ymax - ymin)/5
    if(ymax!=0 & ymin!=0) {
        numbuf = min(abs(ymin), abs(ymax))
        if( devbuf > numbuf ) devbuf = numbuf
    }
    numbuf = as.numeric()
    if( ymin < 0 ) numbuf = c(numbuf, -seq(0, abs(ymin), devbuf))
    if( ymax > 0 ) numbuf = c(numbuf, seq(0,ymax,devbuf))
    axis(2, at=numbuf, labels=round(numbuf,3), las=1, cex.axis=axis)
    ymid = (ymax + ymin)/2
    mtext(paste("eQTL effect sizes of ",data$eprobeID,sep=""), side=2, at=ymid, line=4.5, cex=cex_coeff)
    
    idx=which(data$eSMR[,1]==data$eprobeID)
    if(length(idx)==1) {
        eGeneID=as.character(data$eSMR[idx,4])
    } else {
        print("ERROR: plot file is not correct!");
        quit();
    }
    idx=which(data$mSMR[,1]==exposure_probe)
    if(length(idx)==1) {
        mGeneID=as.character(data$mSMR[idx,4])
    } else {
        print("ERROR: plot file is not correct!");
        quit();
    }
    mainstr1 = substitute(paste(eprobeid, " (", italic(egene), ")", sep=""),list(eprobeid=as.character(data$eprobeID), egene=as.character(eGeneID)))
    mainstr2 = substitute(paste(mprobeid, " (", italic(mgene), ")", sep=""), list(mprobeid=as.character(exposure_probe), mgene=as.character(mGeneID)))
    mtext(mainstr1, side=3, at=xmin, adj=0, line=2.5, cex=cex_coeff)
    mtext(mainstr2, side=3, at=xmin, adj=0, line=0.5, cex=cex_coeff)
    # Plot legend
    lstr = c("top cis-meQTL", "cis-meQTL")
    col_led = c(col_top, col_cis); pch_led = c(pch_top, pch_cis)
   
    
    if(bXY>0) {
        legend("topleft", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
    } else {
        legend("topright", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
    }
}

omicMultiTraitsSMRLocusPlot = function(data1=SMRData1,trait_name1, data2=SMRData2,trait_name2,esmr_thresh1,esmr_thresh2,msmr_thresh1,msmr_thresh2,m2esmr_thresh,esmr_heidi1=0.0499,esmr_heidi2=0.0499,msmr_heidi1=0.0499,msmr_heidi2=0.0499,m2esmr_heidi=0.05,window=500,pointsize=20,max_plot_mprobe=8,epi_plot=FALSE)
{
    ######################## epigenome parts (YW) ############################
    if(epi_plot==TRUE)
    {
        library(data.table)
        dir="/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/"
        samplefile=read.table("/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/legend_YW.txt",head=F,sep="\t")
        colfile=read.delim("/afm01/UQ/Q3003/Q0286_uqywu16/uqywu16/chromstate/annotation_25_imputed12marks.txt",head=T,sep="\t")
        collegend=c("TssA","Prom","Tx","TxWk","TxEn","EnhA","EnhW","DNase","ZNF/Rpts","Het","PromP","PromBiv","ReprPC","Quies")
        type1=c("IMR90","ESC","iPSC","ES-deriv.","Blood & T cell","HSC & B cell","Mesench.","Myosat.","Epithelial","Neurosph.","Thymus","Brain","Adipose","Muscle","Heart","Sm. Muscle","Digestive","Other","ENCODE")
        type2=c("","ESC","iPSC","ES-deriv","Blood & T-cell","HSC & B-cell","Mesenchymal","","Epithelial","","","Brain","","Muscle","Heart","","Digestive","Other","ENCODE")
        colorset=c("red","brown","slateblue","royalblue","green3","green4","lightcoral","darkorange3","darkorange","gold","gold3","goldenrod4","indianred","sienna","lightsalmon3","hotpink","lightpink","grey66","black")
        samples=samplefile[,1]
        rectplot<-function(data,colfile,ymin,ymax){
            for(i in 1:25){
                tmp=subset(data,data$V4==i)
                if(dim(tmp)[1]>0){
                    red=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][1]);green=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][2]);blue=as.numeric(strsplit(as.character(colfile[i,5]),",")[[1]][3]);
                    rect(as.numeric(tmp$V2)/1e6,ymin,as.numeric(tmp$V3)/1e6,ymax,col=rgb(red,green,blue,max=255),border=NA)
                }
            }
        }
        
    }
    ###################### ################################################
    
    if(is.na(esmr_thresh1))
    {
        print("ERROR: please inpute the p-value threshold 1 of your eSMR test.");
        quit();
    }
    if(is.na(esmr_thresh2))
    {
        print("ERROR: please inpute the p-value threshold 2 of your eSMR test.");
        quit();
    }
    if(is.na(msmr_thresh1))
    {
        print("ERROR: please inpute the p-value threshold 1 of your mSMR test.");
        quit();
    }
    if(is.na(msmr_thresh2))
    {
        print("ERROR: please inpute the p-value threshold 2 of your mSMR test.");
        quit();
    }
    if(is.na(m2esmr_thresh))
    {
        print("ERROR: please inpute the p-value threshold of your me2eSMR test.");
        quit();
    }
    if(data1$eprobeID!=data2$eprobeID)
    {
        print("ERROR: inconsistent target gene expression probe.");
        quit();
    }
    
    cex_coeff=3/4 * pointsize/15;
    if(length(which(is.na(data1$eSMR[,3])))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    if(length(which(is.na(data1$mSMR[,3])))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    idx=match(data1$eprobeID,data1$eSMR[,1]);
    #data$mprobeID from m2eSMR result
    idx2=match(data1$mprobeID,data1$mSMR[,1]);
    if(length(idx)==0 | length(idx2)==0){
        print("ERROR: Plot file is not generated correctly, can't find target probe!");
        quit();
    }
    plot_start=data1$eSMR[idx,3]-window*1000
    if(plot_start<0) plot_start=0
    plot_end=data1$eSMR[idx,3]+window*1000;
    
    idx=which(data1$eSMR[,3]>=plot_start & data1$eSMR[,3]<=plot_end)
    data1$eSMR=data1$eSMR[idx,]
    idx=which(data1$mSMR[,3]>=plot_start & data1$mSMR[,3]<=plot_end)
    data1$mSMR=data1$mSMR[idx,]
    idx=which(data1$me2e[,3]>=plot_start & data1$me2e[,3]<=plot_end)
    data1$me2e=data1$me2e[idx,]
    idx=match(data1$GWAS[,1],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$GWAS=data1$GWAS[idx,]
    idx=match(data1$eQTL[,2],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$eQTL=data1$eQTL[idx,]
    idx=match(data1$meQTL[,2],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$meQTL=data1$meQTL[idx,]
    
    data1$mSMR=data1$mSMR[order(data1$mSMR[,8]),]
    
    idx=which(data2$eSMR[,3]>=plot_start & data2$eSMR[,3]<=plot_end)
    data2$eSMR=data2$eSMR[idx,]
    idx=which(data2$mSMR[,3]>=plot_start & data2$mSMR[,3]<=plot_end)
    data2$mSMR=data2$mSMR[idx,]
    idx=which(data2$me2e[,3]>=plot_start & data2$me2e[,3]<=plot_end)
    data2$me2e=data2$me2e[idx,]
    idx=match(data2$GWAS[,1],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$GWAS=data2$GWAS[idx,]
    idx=match(data2$eQTL[,2],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$eQTL=data2$eQTL[idx,]
    idx=match(data2$meQTL[,2],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$meQTL=data2$meQTL[idx,]
    
    data2$mSMR=data2$mSMR[order(data2$mSMR[,8]),]
    
    idx=which(data1$Gene[,2]>=plot_start & data1$Gene[,3]<=plot_end )
    data1$Gene=data1$Gene[idx,]
    idx=which(data2$Gene[,2]>=plot_start & data2$Gene[,3]<=plot_end )
    data2$Gene=data2$Gene[idx,]
    
    #start to plot
    esmrindx1 = which(data1$eSMR[,8] <= esmr_thresh1)
    eheidiindx1 = which((data1$eSMR[,8] <= esmr_thresh1) & (data1$eSMR[,9] >= esmr_heidi1))
    esmrprobes1 = NA; eheidiprobes1 = NA;
    if(length(esmrindx1)>0) { esmrprobes1 =  as.character(data1$eSMR[esmrindx1,1]) }
    if(length(eheidiindx1)>0) { eheidiprobes1 = as.character(data1$eSMR[eheidiindx1,1]) }
    
    esmrindx2 = which(data2$eSMR[,8] <= esmr_thresh2)
    eheidiindx2 = which((data2$eSMR[,8] <= esmr_thresh2) & (data2$eSMR[,9] >= esmr_heidi2))
    esmrprobes2 = NA; eheidiprobes2 = NA;
    if(length(esmrindx2)>0) { esmrprobes2 =  as.character(data2$eSMR[esmrindx2,1]) }
    if(length(eheidiindx2)>0) { eheidiprobes2 = as.character(data2$eSMR[eheidiindx2,1]) }
    
    msmrindx1 = which(data1$mSMR[,8] <= msmr_thresh1)
    mheidiindx1 = which((data1$mSMR[,8] <= msmr_thresh1) & (data1$mSMR[,9] >= msmr_heidi1))
    msmrprobes1 = NA; mheidiprobes1 = NA;
    if(length(msmrindx1)>0) { msmrprobes1 =  as.character(data1$mSMR[msmrindx1,1]) }
    if(length(mheidiindx1)>0) { mheidiprobes1 = as.character(data1$mSMR[mheidiindx1,1]) }
    
    msmrindx2 = which(data2$mSMR[,8] <= msmr_thresh2)
    mheidiindx2 = which((data2$mSMR[,8] <= msmr_thresh2) & (data2$mSMR[,9] >= msmr_heidi2))
    msmrprobes2 = NA; mheidiprobes2 = NA;
    if(length(msmrindx2)>0) { msmrprobes2 =  as.character(data2$mSMR[msmrindx2,1]) }
    if(length(mheidiindx2)>0) { mheidiprobes2 = as.character(data2$mSMR[mheidiindx2,1]) }
    
    
    m2esmrindx1 = which(data1$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx1 = which((data1$me2e[,8] <= m2esmr_thresh) & (data1$me2e[,9] >= m2esmr_heidi))
    m2esmrprobes1 = NA; m2eheidiprobes1 = NA;
    if(length(m2esmrindx1)>0) { m2esmrprobes1 =  as.character(data1$me2e[m2esmrindx1,1]) }
    if(length(m2eheidiindx1)>0) { m2eheidiprobes1 = as.character(data1$me2e[m2eheidiindx1,1]) }
    
    m2esmrindx2 = which(data2$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx2 = which((data2$me2e[,8] <= m2esmr_thresh) & (data2$me2e[,9] >= m2esmr_heidi))
    m2esmrprobes2 = NA; m2eheidiprobes2 = NA;
    if(length(m2esmrindx2)>0) { m2esmrprobes2 =  as.character(data2$me2e[m2esmrindx2,1]) }
    if(length(m2eheidiindx2)>0) { m2eheidiprobes2 = as.character(data2$me2e[m2eheidiindx2,1]) }
    
    
    eprobePLOT=intersect(eheidiprobes1,eheidiprobes2)
    if(length(eprobePLOT)==0)
    {
        print("ERROR: No probe is significant in eSMR test.");
        quit();
    }
    neprobePLOT = length(eprobePLOT)
    mprobePLOT=intersect(union(mheidiprobes1,mheidiprobes2),m2eheidiprobes1) # m2eheidiprobes1 should identify with m2eheidiprobes2
    mprobehome=c()
    for(k in 1:length(mprobePLOT))
    {
        idx1=which(mheidiprobes1==mprobePLOT[k])
        idx2=which(mheidiprobes2==mprobePLOT[k])
        if(length(idx1)>0 & length(idx2)>0) mprobehome=c(mprobehome,3)
        if(length(idx1)>0 & length(idx2)==0) mprobehome=c(mprobehome,1)
        if(length(idx1)==0 & length(idx2)>0) mprobehome=c(mprobehome,2)
    }
    nmprobePLOT = length(mprobePLOT)
    num_mprb_plot=length(mprobePLOT)
    
    ################ select the common mehtylation probes ####################
    if (num_mprb_plot > max_plot_mprobe) {
        commonindx=which(mprobehome==3)
        if (length(commonindx)>=max_plot_mprobe) {selectindx=commonindx[1:max_plot_mprobe]}
        if (length(commonindx)>0 & length(commonindx)<max_plot_mprobe) {
            leftnum=max_plot_mprobe-length(commonindx)
            selectindx=c(commonindx,sample((1:num_mprb_plot)[-commonindx],leftnum))
        }
        if (length(commonindx)==0) {selectindx=sample((1:num_mprb_plot),max_plot_mprobe)}
        nmprobePLOT=num_mprb_plot=max_plot_mprobe
        mprobePLOT=mprobePLOT[selectindx]
        mprobehome=mprobehome[selectindx]
    }
    
    
    idx=which(is.na(data1$GWAS[,2]) | is.na(data1$GWAS[,3]))
    if(length(idx)>0) data1$GWAS=data1$GWAS[-idx,]
    pZY1=-log10(pchisq((data1$GWAS[,2]/data1$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY1)))!=0) {
        maxpZY1 = max(pZY1[-which(is.infinite(pZY1))])
        if (maxpZY1 > 320) {pZY1[which(is.infinite(pZY1))] = maxpZY1 + 5
            } else {pZY1[which(is.infinite(pZY1))] = 325}
    }
    
    idx=which(is.na(data2$GWAS[,2]) | is.na(data2$GWAS[,3]))
    if(length(idx)>0) data2$GWAS=data2$GWAS[-idx,]
    pZY2=-log10(pchisq((data2$GWAS[,2]/data2$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY2)))!=0) {
        maxpZY2 = max(pZY2[-which(is.infinite(pZY2))])
        if (maxpZY2 > 320) {pZY2[which(is.infinite(pZY2))] = maxpZY2 + 5
            } else {pZY2[which(is.infinite(pZY2))] = 325}
    }

    idx=match(data1$eprobeID,data1$eSMR[,1]);
    if(length(idx)>0){
        chrPLOT = data1$eSMR[idx,2]
    }else{
        print("ERROR: Plot file is not generated correctly, please report this bug!");
        quit();
    }
    
    idx=which(is.na(data1$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO1=data1$eSMR[-idx,];
    }else{
        eprobeINFO1=data1$eSMR;
    }
    idx=which(is.na(data2$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO2=data2$eSMR[-idx,];
    }else{
        eprobeINFO2=data2$eSMR;
    }
    
    idx=which(is.na(data1$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO1=data1$mSMR[-idx,];
    }else{
        mprobeINFO1=data1$mSMR;
    }
    idx=which(is.na(data2$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO2=data2$mSMR[-idx,];
    }else{
        mprobeINFO2=data2$mSMR;
    }
    
    idx=which(is.na(eprobeINFO1[,5]) | is.na(eprobeINFO1[,6]));
    idx2=which(is.na(eprobeINFO1[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO1[idx,5]=eprobeINFO1[idx,3]-7500;
    eprobeINFO1[idx,6]=eprobeINFO1[idx,3]+7500;
    eprobeINFO1[,8]=-log10(eprobeINFO1[,8]);
    eprobeINFO1[,3]=eprobeINFO1[,3]/1e6;
    eprobeINFO1[,5]=eprobeINFO1[,5]/1e6;
    eprobeINFO1[,6]=eprobeINFO1[,6]/1e6;
    
    idx=which(is.na(eprobeINFO2[,5]) | is.na(eprobeINFO2[,6]));
    idx2=which(is.na(eprobeINFO2[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO2[idx,5]=eprobeINFO2[idx,3]-7500;
    eprobeINFO2[idx,6]=eprobeINFO2[idx,3]+7500;
    eprobeINFO2[,8]=-log10(eprobeINFO2[,8]);
    eprobeINFO2[,3]=eprobeINFO2[,3]/1e6;
    eprobeINFO2[,5]=eprobeINFO2[,5]/1e6;
    eprobeINFO2[,6]=eprobeINFO2[,6]/1e6;
    
    
    idx=which(is.na(mprobeINFO1[,5]) | is.na(mprobeINFO1[,6]));
    idx2=which(is.na(mprobeINFO1[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO1[idx,5]=mprobeINFO1[idx,3]-7500;
    mprobeINFO1[idx,6]=mprobeINFO1[idx,3]+7500;
    mprobeINFO1[,8]=-log10(mprobeINFO1[,8]);
    mprobeINFO1[,3]=mprobeINFO1[,3]/1e6;
    mprobeINFO1[,5]=mprobeINFO1[,5]/1e6;
    mprobeINFO1[,6]=mprobeINFO1[,6]/1e6;
    
    idx=which(is.na(mprobeINFO2[,5]) | is.na(mprobeINFO2[,6]));
    idx2=which(is.na(mprobeINFO2[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO2[idx,5]=mprobeINFO2[idx,3]-7500;
    mprobeINFO2[idx,6]=mprobeINFO2[idx,3]+7500;
    mprobeINFO2[,8]=-log10(mprobeINFO2[,8]);
    mprobeINFO2[,3]=mprobeINFO2[,3]/1e6;
    mprobeINFO2[,5]=mprobeINFO2[,5]/1e6;
    mprobeINFO2[,6]=mprobeINFO2[,6]/1e6;
    
    epXY1=eprobeINFO1[,8];
    mpXY1=mprobeINFO1[,8];
    yMAX1 = ceiling(max(c(pZY1, epXY1,mpXY1), na.rm=T)) + 1;
    epXY2=eprobeINFO2[,8];
    mpXY2=mprobeINFO2[,8];
    yMAX2 = ceiling(max(c(pZY2, epXY2,mpXY2), na.rm=T)) + 1;
    yMAX=yMAX1
    
    #glist=cbind(eprobeINFO[,2],eprobeINFO[,5:6],as.character(eprobeINFO[,4]),eprobeINFO[,7]);#old
    glist=data1$Gene;
    glist[,2]=glist[,2]/1e6;
    glist[,3]=glist[,3]/1e6;
    colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
    idx=which(is.na(glist[,2]) | is.na(glist[,3]));
    if(length(idx>0)) glist=glist[-idx,];
    generow = GeneRowNum(glist);
    num_row = max(as.numeric(generow$ROW));
    offset_map = ceiling(yMAX1);
    
    offset_map=max(offset_map,num_row*5.5)
    
    offset_probe = yMAX1 / 2.5;
    
    offset_eqtl = ceiling(yMAX1 / 2.5) + 0.5;
    dev_axis = 0.1*yMAX1;
    if(dev_axis<1.5) dev_axis = 1.5;
    
    if(epi_plot==TRUE) {
        yaxis.min = -offset_map - offset_eqtl*neprobePLOT - offset_eqtl*num_mprb_plot- dev_axis*(neprobePLOT)- dev_axis*(num_mprb_plot)-yMAX-2*dev_axis-offset_eqtl*8-offset_eqtl/5;
    } else {
        yaxis.min = -offset_map - offset_eqtl*(neprobePLOT) - offset_eqtl*num_mprb_plot- dev_axis*(neprobePLOT)- dev_axis*(num_mprb_plot)-yMAX-2*dev_axis;
    }
    yaxis.max = yMAX + ceiling(offset_probe) + 1;
    
    # scales of x-axis
    idx=match(data1$GWAS[,1],data1$SNP[,1]);
    gwasBP = as.numeric(data1$SNP[idx,3])/1e6;
    min.pos = min(gwasBP);
    max.pos = max(gwasBP);
    start = min(as.numeric(glist[,2]));
    end = max(as.numeric(glist[,3]));
    bp = c(min.pos, max.pos, start, end);
    xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
    xmax=xmax+(xmax-xmin)*0.15 #extend
    ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
    xlab = paste("Chromosome", chrPLOT, "Mb");
    # plot GWAS p value
    par(mar=c(5,5,3,2), xpd=TRUE)
    plot(gwasBP, pZY1, yaxt="n", xaxt="n",bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    axis(1, at=round(seq(xmin,xmax,(xmax-xmin)/5),2))
    # y1 axis
    devbuf1 = yMAX/4
    axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
    mtext(ylab, side=2, line=3, at=(abs(yMAX-yMAX2)*3/4), cex=cex_coeff);
    text(xmin, max(pZY1,na.rm=T) , label=trait_name1,col="blue", cex=1, adj=0)
    
    axis.start = 0;
    axis.start = axis.start-yMAX-2*dev_axis;
    pZY2buf=pZY2/max(pZY2)*yMAX+axis.start;
    par(new=TRUE)
    idx=match(data2$GWAS[,1],data2$SNP[,1]);
    gwasBP2 = as.numeric(data2$SNP[idx,3])/1e6;
    plot(gwasBP2, pZY2buf, yaxt="n",xaxt="n",bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    devbuf1 = yMAX/4; devbuf2 = yMAX2/4
    axis(2, at=seq(axis.start,(axis.start+yMAX),devbuf1), labels=round(seq(0,yMAX2,devbuf1),0), las=1, cex.axis=axis);
    segments(xmin, axis.start+yMAX+dev_axis/2, xmax, axis.start+yMAX+dev_axis/2,
    col="dim grey", lty="24", lwd=1)
    text(xmin, axis.start+yMAX-1 , label=trait_name2,col="darkcyan", cex=1, adj=0)
    
    
    eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
    axis.down = offset_eqtl + dev_axis;
    for( k in 1 : neprobePLOT ) {
        axis.start = axis.start - axis.down
        eqtlinfobuf = data1$eQTL[which(data1$eQTL[,1]==eprobePLOT[k]),]
        if(dim(eqtlinfobuf)[1]==0) next;
        pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
        if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
        }
        if(length(which(eheidiprobes1==eprobePLOT[k]))==0) {
            col_eqtl = "navy"
        } else col_eqtl = "maroon"
        eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
        pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
        idx=match(eqtlinfobuf[,2],data1$SNP[,1]);
        eqtlbp = as.numeric(data1$SNP[idx,3])/1e6;
        probegene = unique(as.character(data1$eSMR[which(data1$eSMR[,1]==eprobePLOT[k]),4]))
        par(new=TRUE)
        pchbuf = 4;
        #if(k%%2==0) pchbuf = 20;
        plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the eQTLs
        colme2e="black"
        if(eprobePLOT[k]==data1$eprobeID) colme2e="darkorchid2"
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=eprobePLOT[k], geneid=probegene)),col=colme2e, cex=1, adj=0)
        if(k==1) labpos=axis.start+offset_eqtl-dev_axis/2;
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,eqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    ypos = labpos;
    mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff,adj=1)
    
    #plot meQTL
    meqtl.lab = expression(paste("-", log[10], "(", italic(P), " meQTL)", sep=""));
    for( k in 1 : nmprobePLOT) {
        axis.start = axis.start - axis.down
        homeid=mprobehome[k]
        if(homeid==2) meqtlinfobuf = data2$meQTL[which(data2$meQTL[,1]==mprobePLOT[k]),];
        if(homeid!=2) meqtlinfobuf = data1$meQTL[which(data1$meQTL[,1]==mprobePLOT[k]),]
        if(dim(meqtlinfobuf)[1]==0) next;
        
        pvalbuf=-log10(pchisq((meqtlinfobuf[,3]/meqtlinfobuf[,4])^2,1,lower.tail=F));
        if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
        }
        if(length(which(msmrprobes1==mprobePLOT[k]))==0) {
            col_eqtl = "navy"
        } else col_eqtl = "maroon"
        meqtl.min = 0; meqtl.max = ceiling(max(pvalbuf))
        pvalbuf = pvalbuf/meqtl.max * offset_eqtl + axis.start
        if(homeid==2) {
            idx=match(meqtlinfobuf[,2],data2$SNP[,1]);
            meqtlbp = as.numeric(data2$SNP[idx,3])/1e6;
            mprobegene = unique(as.character(data2$mSMR[which(data2$mSMR[,1]==mprobePLOT[k]),4]))
        } else {
            idx=match(meqtlinfobuf[,2],data1$SNP[,1]);
            meqtlbp = as.numeric(data1$SNP[idx,3])/1e6;
            mprobegene = unique(as.character(data1$mSMR[which(data1$mSMR[,1]==mprobePLOT[k]),4]))
        }
        
        par(new=TRUE)
        pchbuf = 20;
        #if(k%%2==0) pchbuf = 4;
        plot(meqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the meQTLs
        if(homeid==3) colme2e="red"
        if(homeid==1) colme2e="blue"
        if(homeid==2) colme2e="darkcyan"
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=mprobePLOT[k], geneid=mprobegene)),col=colme2e, cex=1, adj=0)
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = meqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,meqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    ypos = (axis.start - dev_axis-axis.down)*3/4
    mtext(meqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)
    
    ########### plot the 127 chromtin states samples #################
    if(epi_plot==TRUE)
    {
        axis.start=axis.start-dev_axis
        epi.start=color.start=axis.start
        regionchr=paste("chr",chrPLOT,sep="")
        unit=offset_eqtl*8/130
        dev=0.02*unit
        axis.down=unit+dev
        for(i in 1:length(samples)){
            filename=paste(as.character(samples[i]),"_25_imputed12marks_stateno.bed.gz",sep="")
            file=fread(paste("gunzip -c ",dir,filename,sep=""),head=F)
            file=as.data.frame(file[-dim(file)[1],])
            regionindx=which(file$V1==regionchr & as.numeric(file$V2)>=plot_start & as.numeric(file$V3)<=plot_end)
            plotfile=file[regionindx,]
            par(new=T)
            rectplot(plotfile,colfile,axis.start-unit,axis.start)
            lines(c(plot_start/1e6,plot_end/1e6),c(axis.start,axis.start),col="dim grey", lty="24", lwd=1)
            lines(c(plot_start/1e6,plot_end/1e6),c(axis.start-unit,axis.start-unit),col="dim grey", lty="24", lwd=1)
            axis.start=axis.start-axis.down
        }
        
        ########################## annotation of cell types ###################
        labwidth=(plot_end/1e6-plot_start/1e6)/60
        for(i in 1:length(unique(samplefile[,3]))){
            dup=length(which(samplefile[,3]==i))
            epi.down=epi.start-dup*(unit+dev)
            par(new=T)
            rect(plot_start/1e6-labwidth,epi.start,plot_start/1e6,epi.down,col=colorset[i],border=NA)
            #par(new=T)
            #rect(plot_end/1e6,epi.start,plot_end/1e6+labwidth,epi.down,col=colorset[i],border=NA)
            epiylab=type2[i]
            labpos=(epi.down+epi.start)/2
            mtext(epiylab, side=2, line=0, at=labpos,cex=axis*0.8,las=1);
            epi.start=epi.down-dev
        }
        
        ########################## color legend ###################
        labwidth=(plot_end/1e6-plot_start/1e6)/60
        unit=(offset_eqtl*8)/length(collegend)
        legendcol=unique(colfile[,5])
        for(i in c(1:7,14,8:13)) {
            color.down=color.start-unit
            red=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][1]);green=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][2]);blue=as.numeric(strsplit(as.character(legendcol[i]),",")[[1]][3]);
            par(new=T)
            rect(plot_end/1e6+4*labwidth,color.start,plot_end/1e6+5*labwidth,color.down,col=rgb(red,green,blue,max=255),border=NA)
            colorylab=collegend[i]
            labpos=(color.down+color.start)/2
            #mtext(collegend[i], side=4, line=-5, at=labpos,cex=0.8,las=1);
            text(plot_end/1e6+6*labwidth, labpos, labels=collegend[i], adj=0, cex=axis*0.8);
            color.start=color.down-dev
        }
        
        
    }
    
    ######################## gene label ######################
    # plot p value of bTG
    # all the probes
    
    num_gene = dim(generow)[1]
    dist = offset_map/num_row
    for( k in 1 : num_row ) {
        generowbuf = generow[which(as.numeric(generow[,5])==k),]
        xstart = as.numeric(generowbuf[,3])
        xend = as.numeric(generowbuf[,4])
        snbuf = which(xend-xstart< 1e-3)
        if(length(snbuf)>0) {
            xstart[snbuf] = xstart[snbuf] - 0.0025
            xend[snbuf] = xend[snbuf] + 0.0025
        }
        xcenter = (xstart+xend)/2
        xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
        num_genebuf = dim(generowbuf)[1]
        for( l in 1 : num_genebuf ) {
            ofs=0.3
            if(l%%2==0) ofs=-0.8
            m = num_row - k
            ypos = m*dist + yaxis.min
            code = 1
            if(generowbuf[l,2]=="+") code = 2;
            arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
            col=colors()[75], lwd=1)
            movebuf = as.numeric(generowbuf[l,6])*genemove
            text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.8)
        }
    }
    
    # plot the probes
    
    eidx1=match(eprobePLOT,eprobeINFO1[,1])
    eprobein1=eprobeINFO1[eidx1,]
    eidx2=match(eprobePLOT,eprobeINFO2[,1])
    eprobein2=eprobeINFO2[eidx2,]
    
    midx1=match(mprobePLOT,mprobeINFO1[,1])
    idxx=which(is.na(midx1))
    if(length(idxx)>0) {
        mprobein1=mprobeINFO1[midx1[-idxx],]
    } else {
        mprobein1=mprobeINFO1[midx1,]
    }
    
    midx2=match(mprobePLOT,mprobeINFO2[,1])
    idxx=which(is.na(midx2))
    if(length(idxx)>0) {
        mprobein2=mprobeINFO2[midx2[-idxx],]
    } else {
        mprobein2=mprobeINFO2[midx2,]
    }
    
    for( k in 1 : dim(mprobein1)[1])  {
        hitflag=FALSE
        if(length(which(mheidiprobes1==mprobein1[k,1]))>0 & length(which(msmrprobes1==mprobein1[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=21
        } else if(length(which(msmrprobes1==mprobein1[k,1]))>0){
            colplot = "maroon"; colfont=2; pchbuf=1
        } else {
            colplot = "navy"; colfont=1; pchbuf=1
        }
        if( as.numeric(mprobein1[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        plot_probe(mprobein1, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
        
    }
    
    for( k in 1 : dim(mprobein2)[1])  {
        hitflag=FALSE
        if(length(which(mheidiprobes2==mprobein2[k,1]))>0 & length(which(msmrprobes2==mprobein2[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=21
        } else if (length(which(msmrprobes2==mprobein2[k,1]))>0){
            colplot = "maroon"; colfont=2; pchbuf=1
        } else {
            colplot = "navy"; colfont=1; pchbuf=1
        }
        if( as.numeric(mprobein2[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        xcenter = as.numeric(mprobein2[k,3])
        pvalbuf = as.numeric(mprobein2[k,8])
        pvalbuf_=pvalbuf-yMAX2-2*dev_axis;
        strbuf = mprobein2[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }
        
    }
    
    for( k in 1 : dim(eprobein1)[1])  {
        hitflag=FALSE
        if(length(which(eheidiprobes1==eprobein1[k,1]))>0 & length(which(esmrprobes1==eprobein1[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else if (length(which(esmrprobes1==eprobein1[k,1]))>0){
            colplot = "maroon"; colfont=2; pchbuf=5
        } else {
            colplot = "navy"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobein1[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        plot_probe(eprobein1, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
    }
    
    for( k in 1 : dim(eprobein2)[1])  {
        hitflag=FALSE
        if(length(which(eheidiprobes2==eprobein2[k,1]))>0 & length(which(esmrprobes2==eprobein2[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else if(length(which(esmrprobes2==eprobein2[k,1]))>0){
            colplot = "maroon"; colfont=2; pchbuf=5
        } else {
            colplot = "navy"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobein2[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        
        
        xcenter = as.numeric(eprobein2[k,3])
        pvalbuf = as.numeric(eprobein2[k,8])
        pvalbuf_=pvalbuf-yMAX2-2*dev_axis;
        strbuf = eprobein2[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }
        
    }
    
    
    #anno probe
    probeINFO=rbind(eprobein1)
    num_anno_probe=dim(probeINFO)[1]
    xcenter = as.numeric(probeINFO[,3])
    xcbuf = xcenter
    xcenter = spread.labs(xcenter[1:num_anno_probe], mindiff=0.08, maxiter=1000, min = xmin, max = xmax-1)
    # adjust the line position
    
    adjflag = rep(0, num_anno_probe)
    if(num_anno_probe>1) {
        dbuf = c(0, xcbuf[1:(num_anno_probe-1)])
        mflag = as.numeric(abs(xcbuf[1:(num_anno_probe)] - dbuf) < 0.01)
        adjflag = as.numeric( mflag | c(mflag[2:num_anno_probe],0) )
    }
    
    for( k in 1 : num_anno_probe)  {
        
        # annotate the probes
        
        ypos = 1.02*yMAX
        strbuf =
        text(xcenter[k], ypos,
        labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),
        list(probeid=as.character(probeINFO[k,1]),
        genename=as.character(probeINFO[k,4]))),
        ylim=c(yaxis.min, yaxis.max),
        srt=30, col=colplot, font=colfont, cex=1, adj=0)
        # plot the lines
        # 1st step
        xstart = xcbuf[k]
        ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
        
        if(adjflag[k]==1) {
            xstart = (xcbuf[k] + xcenter[k])/2
            segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
        }
        
        segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
        # 2nd step
        xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*1.01;
        segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
    }
    
    # plot the threshold
    # eSMR threshold
    yebuf = -log10(as.numeric(esmr_thresh1)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh1)); dev_anno = yMAX/9;
    
    strbuf = paste("pESMR = ",esmr_thresh1, sep="")
    segments(xmin, yebuf, xmax, yebuf, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno
    } else {
        ythrespos=yebuf+dev_anno
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
    
    
    strbuf = paste("pMSMR = ",msmr_thresh1, sep="")
    segments(xmin, ymbuf, xmax, ymbuf, col="navy", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=ymbuf+dev_anno
    } else {
        ythrespos=ymbuf-dev_anno
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=axis,font=3);
    
    # plot the threshold
    # eSMR threshold
    yebuf = -log10(as.numeric(esmr_thresh2)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh2)); dev_anno = yMAX/9;
    
    strbuf = paste("pESMR = ",esmr_thresh2, sep="")
    segments(xmin, yebuf-yMAX2-2*dev_axis, xmax, yebuf-yMAX2-2*dev_axis, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno-yMAX2-2*dev_axis
    } else {
        ythrespos=yebuf+dev_anno-yMAX2-2*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
    
    
    strbuf = paste("pMSMR = ",msmr_thresh2, sep="")
    segments(xmin, ymbuf-yMAX2-2*dev_axis, xmax, ymbuf-yMAX2-2*dev_axis, col="navy", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=ymbuf+dev_anno-yMAX2-2*dev_axis
    } else {
        ythrespos=ymbuf-dev_anno-yMAX2-2*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=axis,font=3);
    
}

omicMultiTraitsSMRLocusPlot3G = function(data1=SMRData1,trait_name1,data2=SMRData2,data3=SMRData3,trait_name2,trait_name3,esmr_thresh1,esmr_thresh2,esmr_thresh3,msmr_thresh1,msmr_thresh2,msmr_thresh3,esmr_heidi1,esmr_heidi2,esmr_heidi3,msmr_heidi1,msmr_heidi2,msmr_heidi3,m2esmr_heidi,max_plot_mprobe=10,m2esmr_thresh=2.9e-7,window=500,pointsize=20)
{
    if(is.na(esmr_thresh1))
    {
        print("ERROR: please inpute the p-value threshold 1 of your eSMR test.");
        quit();
    }
    if(is.na(esmr_thresh2))
    {
        print("ERROR: please inpute the p-value threshold 2 of your eSMR test.");
        quit();
    }
    if(is.na(esmr_thresh3))
    {
        print("ERROR: please inpute the p-value threshold 2 of your eSMR test.");
        quit();
    }
    if(is.na(msmr_thresh1))
    {
        print("ERROR: please inpute the p-value threshold 1 of your mSMR test.");
        quit();
    }
    if(is.na(msmr_thresh2))
    {
        print("ERROR: please inpute the p-value threshold 2 of your mSMR test.");
        quit();
    }
    if(is.na(msmr_thresh3))
    {
        print("ERROR: please inpute the p-value threshold 2 of your mSMR test.");
        quit();
    }
    if(is.na(m2esmr_thresh))
    {
        print("ERROR: please inpute the p-value threshold of your me2eSMR test.");
        quit();
    }
    if(data1$eprobeID!=data2$eprobeID || data1$eprobeID!=data3$eprobeID)
    {
        print("ERROR: inconsistent target gene expression probe.");
        quit();
    }

    cex_coeff=3/4 * pointsize/15;
    if(length(which(is.na(data1$eSMR[,3])))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    if(length(which(is.na(data1$mSMR[,3])))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    idx=match(data1$eprobeID,data1$eSMR[,1]);
    #data$mprobeID from m2eSMR result
    idx2=match(data1$mprobeID,data1$mSMR[,1]);
    if(length(idx)==0 | length(idx2)==0){
        print("ERROR: Plot file is not generated correctly, can't find target probe!");
        quit();
    }
    plot_start=data1$eSMR[idx,3]-window*1000
    if(plot_start<0) plot_start=0
    plot_end=data1$eSMR[idx,3]+window*1000;
    
    idx=which(data1$eSMR[,3]>=plot_start & data1$eSMR[,3]<=plot_end)
    data1$eSMR=data1$eSMR[idx,]
    idx=which(data1$mSMR[,3]>=plot_start & data1$mSMR[,3]<=plot_end)
    data1$mSMR=data1$mSMR[idx,]
    idx=which(data1$me2e[,3]>=plot_start & data1$me2e[,3]<=plot_end)
    data1$me2e=data1$me2e[idx,]
    idx=match(data1$GWAS[,1],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$GWAS=data1$GWAS[idx,]
    idx=match(data1$eQTL[,2],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$eQTL=data1$eQTL[idx,]
    idx=match(data1$meQTL[,2],data1$SNP[,1])
    tmpsnpbp=data1$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data1$meQTL=data1$meQTL[idx,]
    
    data1$mSMR=data1$mSMR[order(data1$mSMR[,8]),]
    
    idx=which(data2$eSMR[,3]>=plot_start & data2$eSMR[,3]<=plot_end)
    data2$eSMR=data2$eSMR[idx,]
    idx=which(data2$mSMR[,3]>=plot_start & data2$mSMR[,3]<=plot_end)
    data2$mSMR=data2$mSMR[idx,]
    idx=which(data2$me2e[,3]>=plot_start & data2$me2e[,3]<=plot_end)
    data2$me2e=data2$me2e[idx,]
    idx=match(data2$GWAS[,1],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$GWAS=data2$GWAS[idx,]
    idx=match(data2$eQTL[,2],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$eQTL=data2$eQTL[idx,]
    idx=match(data2$meQTL[,2],data2$SNP[,1])
    tmpsnpbp=data2$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data2$meQTL=data2$meQTL[idx,]
    
    data2$mSMR=data2$mSMR[order(data2$mSMR[,8]),]

    idx=which(data3$eSMR[,3]>=plot_start & data3$eSMR[,3]<=plot_end)
    data3$eSMR=data3$eSMR[idx,]
    idx=which(data3$mSMR[,3]>=plot_start & data3$mSMR[,3]<=plot_end)
    data3$mSMR=data3$mSMR[idx,]
    idx=which(data3$me2e[,3]>=plot_start & data3$me2e[,3]<=plot_end)
    data3$me2e=data3$me2e[idx,]
    idx=match(data3$GWAS[,1],data3$SNP[,1])
    tmpsnpbp=data3$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data3$GWAS=data3$GWAS[idx,]
    idx=match(data3$eQTL[,2],data3$SNP[,1])
    tmpsnpbp=data3$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data3$eQTL=data3$eQTL[idx,]
    idx=match(data3$meQTL[,2],data3$SNP[,1])
    tmpsnpbp=data3$SNP[idx,3]
    idx=which(tmpsnpbp>=plot_start &tmpsnpbp<=plot_end)
    data3$meQTL=data3$meQTL[idx,]
    
    data3$mSMR=data3$mSMR[order(data3$mSMR[,8]),]
    
    idx=which(data1$Gene[,2]>=plot_start & data1$Gene[,3]<=plot_end )
    data1$Gene=data1$Gene[idx,]
    idx=which(data2$Gene[,2]>=plot_start & data2$Gene[,3]<=plot_end )
    data2$Gene=data2$Gene[idx,]
    idx=which(data3$Gene[,2]>=plot_start & data3$Gene[,3]<=plot_end )
    data3$Gene=data3$Gene[idx,]

    #start to plot
    esmrindx1 = which(data1$eSMR[,8] <= esmr_thresh1)
    eheidiindx1 = which((data1$eSMR[,8] <= esmr_thresh1) & (data1$eSMR[,9] >= esmr_heidi1))
    esmrprobes1 = NA; eheidiprobes1 = NA;
    if(length(esmrindx1)>0) { esmrprobes1 =  as.character(data1$eSMR[esmrindx1,1]) }
    if(length(eheidiindx1)>0) { eheidiprobes1 = as.character(data1$eSMR[eheidiindx1,1]) }
    
    esmrindx2 = which(data2$eSMR[,8] <= esmr_thresh2)
    eheidiindx2 = which((data2$eSMR[,8] <= esmr_thresh2) & (data2$eSMR[,9] >= esmr_heidi2))
    esmrprobes2 = NA; eheidiprobes2 = NA;
    if(length(esmrindx2)>0) { esmrprobes2 =  as.character(data2$eSMR[esmrindx2,1]) }
    if(length(eheidiindx2)>0) { eheidiprobes2 = as.character(data2$eSMR[eheidiindx2,1]) }

    esmrindx3 = which(data3$eSMR[,8] <= esmr_thresh3)
    eheidiindx3 = which((data3$eSMR[,8] <= esmr_thresh3) & (data3$eSMR[,9] >= esmr_heidi3))
    esmrprobes3 = NA; eheidiprobes3 = NA;
    if(length(esmrindx3)>0) { esmrprobes3 =  as.character(data3$eSMR[esmrindx3,1]) }
    if(length(eheidiindx3)>0) { eheidiprobes3 = as.character(data3$eSMR[eheidiindx3,1]) }
    
    msmrindx1 = which(data1$mSMR[,8] <= msmr_thresh1)
    mheidiindx1 = which((data1$mSMR[,8] <= msmr_thresh1) & (data1$mSMR[,9] >= msmr_heidi1))
    msmrprobes1 = NA; mheidiprobes1 = NA;
    if(length(msmrindx1)>0) { msmrprobes1 =  as.character(data1$mSMR[msmrindx1,1]) }
    if(length(mheidiindx1)>0) { mheidiprobes1 = as.character(data1$mSMR[mheidiindx1,1]) }
    
    msmrindx2 = which(data2$mSMR[,8] <= msmr_thresh2)
    mheidiindx2 = which((data2$mSMR[,8] <= msmr_thresh2) & (data2$mSMR[,9] >= msmr_heidi2))
    msmrprobes2 = NA; mheidiprobes2 = NA;
    if(length(msmrindx2)>0) { msmrprobes2 =  as.character(data2$mSMR[msmrindx2,1]) }
    if(length(mheidiindx2)>0) { mheidiprobes2 = as.character(data2$mSMR[mheidiindx2,1]) }
    
    msmrindx3 = which(data3$mSMR[,8] <= msmr_thresh3)
    mheidiindx3 = which((data3$mSMR[,8] <= msmr_thresh3) & (data3$mSMR[,9] >= msmr_heidi3))
    msmrprobes3 = NA; mheidiprobes3 = NA;
    if(length(msmrindx3)>0) { msmrprobes3 =  as.character(data3$mSMR[msmrindx3,1]) }
    if(length(mheidiindx3)>0) { mheidiprobes3 = as.character(data3$mSMR[mheidiindx3,1]) }

    m2esmrindx1 = which(data1$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx1 = which((data1$me2e[,8] <= m2esmr_thresh) & (data1$me2e[,9] >= m2esmr_heidi))
    m2esmrprobes1 = NA; m2eheidiprobes1 = NA;
    if(length(m2esmrindx1)>0) { m2esmrprobes1 =  as.character(data1$me2e[m2esmrindx1,1]) }
    if(length(m2eheidiindx1)>0) { m2eheidiprobes1 = as.character(data1$me2e[m2eheidiindx1,1]) }
    
    m2esmrindx2 = which(data2$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx2 = which((data2$me2e[,8] <= m2esmr_thresh) & (data2$me2e[,9] >= m2esmr_heidi))
    m2esmrprobes2 = NA; m2eheidiprobes2 = NA;
    if(length(m2esmrindx2)>0) { m2esmrprobes2 =  as.character(data2$me2e[m2esmrindx2,1]) }
    if(length(m2eheidiindx2)>0) { m2eheidiprobes2 = as.character(data2$me2e[m2eheidiindx2,1]) }

    m2esmrindx3 = which(data3$me2e[,8] <= m2esmr_thresh)
    m2eheidiindx3 = which((data3$me2e[,8] <= m2esmr_thresh) & (data3$me2e[,9] >= m2esmr_heidi))
    m2esmrprobes3 = NA; m2eheidiprobes3 = NA;
    if(length(m2esmrindx3)>0) { m2esmrprobes3 =  as.character(data3$me2e[m2esmrindx3,1]) }
    if(length(m2eheidiindx3)>0) { m2eheidiprobes3 = as.character(data3$me2e[m2eheidiindx3,1]) }
    
    eprobePLOT=intersect(intersect(eheidiprobes1,eheidiprobes2),eheidiprobes3)
    if(length(eprobePLOT)==0)
    {
        print("ERROR: No probe is significant in eSMR test.");
        quit();
    }
    neprobePLOT = length(eprobePLOT)
    mprobePLOT=intersect(union(union(mheidiprobes1,mheidiprobes2),mheidiprobes3),m2eheidiprobes1) # m2eheidiprobes1 should identify with m2eheidiprobes2
    num_mprb_plot=length(mprobePLOT)
    mprobehome=c()
    for(k in 1:length(mprobePLOT))
    {
        idx1=which(mheidiprobes1==mprobePLOT[k])
        idx2=which(mheidiprobes2==mprobePLOT[k])
        idx3=which(mheidiprobes3==mprobePLOT[k])
        if(length(idx1)>0 & length(idx2)>0 & length(idx3)>0 ) mprobehome=c(mprobehome,3)
        if((length(idx1)>0 & length(idx2)>0) ||(length(idx1)>0 & length(idx3)>0)||(length(idx2)>0 & length(idx3)>0)) mprobehome=c(mprobehome,2)
        #if(length(idx1)==0 & length(idx2)>0) 
        else mprobehome=c(mprobehome,1)
    }

    ########### select plot #####################
    if (num_mprb_plot > max_plot_mprobe) {        
        commonindx=which(mprobehome==3)
        secondindx=which(mprobehome==2)
        firstindx=which(mprobehome==1)
        if (length(commonindx)>=max_plot_mprobe) {selectindx=sample(commonindx,max_plot_mprobe)}
        if (length(commonindx)<max_plot_mprobe & (length(commonindx)+length(secondindx))>=max_plot_mprobe) {selectindx=c(commonindx,sample(secondindx,(max_plot_mprobe-length(commonindx))))}
        if ((length(commonindx)+length(secondindx))<max_plot_mprobe)
        {selectindx=c(commonindx,secondindx,sample(firstindx,(max_plot_mprobe-length(commonindx)-length(secondindx))))}
    }

    mprobePLOT=mprobePLOT[selectindx]
    mprobehome=mprobehome[selectindx]

    nmprobePLOT = length(mprobePLOT)
    num_mprb_plot=length(mprobePLOT)
    
    idx=which(is.na(data1$GWAS[,2]) | is.na(data1$GWAS[,3]))
    if(length(idx)>0) data1$GWAS=data1$GWAS[-idx,]
    pZY1=-log10(pchisq((data1$GWAS[,2]/data1$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY1)))!=0) {
        maxpZY1 = max(pZY1[-which(is.infinite(pZY1))])
        if (maxpZY1 > 320) {pZY1[which(is.infinite(pZY1))] = maxpZY1 + 5
            } else {pZY1[which(is.infinite(pZY1))] = 325}
    }
    
    idx=which(is.na(data2$GWAS[,2]) | is.na(data2$GWAS[,3]))
    if(length(idx)>0) data2$GWAS=data2$GWAS[-idx,]
    pZY2=-log10(pchisq((data2$GWAS[,2]/data2$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY2)))!=0) {
        maxpZY2 = max(pZY2[-which(is.infinite(pZY2))])
        if (maxpZY2 > 320) {pZY2[which(is.infinite(pZY2))] = maxpZY2 + 5
            } else {pZY2[which(is.infinite(pZY2))] = 325}
    }
    
    idx=which(is.na(data3$GWAS[,2]) | is.na(data3$GWAS[,3]))
    if(length(idx)>0) data3$GWAS=data3$GWAS[-idx,]
    pZY3=-log10(pchisq((data3$GWAS[,2]/data3$GWAS[,3])^2,1,lower.tail=F))
    if(length(which(is.infinite(pZY3)))!=0) {
        maxpZY3 = max(pZY3[-which(is.infinite(pZY3))])
        if (maxpZY3 > 320) {pZY3[which(is.infinite(pZY3))] = maxpZY3 + 5
            } else {pZY3[which(is.infinite(pZY3))] = 325}
    }

    idx=match(data1$eprobeID,data1$eSMR[,1]);
    if(length(idx)>0){
        chrPLOT = data1$eSMR[idx,2]
    }else{
        print("ERROR: Plot file is not generated correctly, please report this bug!");
        quit();
    }
    
    idx=which(is.na(data1$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO1=data1$eSMR[-idx,];
    }else{
        eprobeINFO1=data1$eSMR;
    }
    idx=which(is.na(data2$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO2=data2$eSMR[-idx,];
    }else{
        eprobeINFO2=data2$eSMR;
    }
    idx=which(is.na(data3$eSMR[,8]) )
    if(length(idx)>0) {
        eprobeINFO3=data3$eSMR[-idx,];
    }else{
        eprobeINFO3=data3$eSMR;
    }

    idx=which(is.na(data1$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO1=data1$mSMR[-idx,];
    }else{
        mprobeINFO1=data1$mSMR;
    }
    idx=which(is.na(data2$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO2=data2$mSMR[-idx,];
    }else{
        mprobeINFO2=data2$mSMR;
    }
    idx=which(is.na(data3$mSMR[,8]) )
    if(length(idx)>0) {
        mprobeINFO3=data3$mSMR[-idx,];
    }else{
        mprobeINFO3=data3$mSMR;
    }

    idx=which(is.na(eprobeINFO1[,5]) | is.na(eprobeINFO1[,6]));
    idx2=which(is.na(eprobeINFO1[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO1[idx,5]=eprobeINFO1[idx,3]-7500;
    eprobeINFO1[idx,6]=eprobeINFO1[idx,3]+7500;
    eprobeINFO1[,8]=-log10(eprobeINFO1[,8]);
    eprobeINFO1[,3]=eprobeINFO1[,3]/1e6;
    eprobeINFO1[,5]=eprobeINFO1[,5]/1e6;
    eprobeINFO1[,6]=eprobeINFO1[,6]/1e6;
    
    idx=which(is.na(eprobeINFO2[,5]) | is.na(eprobeINFO2[,6]));
    idx2=which(is.na(eprobeINFO2[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO2[idx,5]=eprobeINFO2[idx,3]-7500;
    eprobeINFO2[idx,6]=eprobeINFO2[idx,3]+7500;
    eprobeINFO2[,8]=-log10(eprobeINFO2[,8]);
    eprobeINFO2[,3]=eprobeINFO2[,3]/1e6;
    eprobeINFO2[,5]=eprobeINFO2[,5]/1e6;
    eprobeINFO2[,6]=eprobeINFO2[,6]/1e6;
    
    idx=which(is.na(eprobeINFO3[,5]) | is.na(eprobeINFO3[,6]));
    idx2=which(is.na(eprobeINFO3[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some eprobes' physical positon is missing!");
        quit();
    }
    eprobeINFO3[idx,5]=eprobeINFO3[idx,3]-7500;
    eprobeINFO3[idx,6]=eprobeINFO3[idx,3]+7500;
    eprobeINFO3[,8]=-log10(eprobeINFO3[,8]);
    eprobeINFO3[,3]=eprobeINFO3[,3]/1e6;
    eprobeINFO3[,5]=eprobeINFO3[,5]/1e6;
    eprobeINFO3[,6]=eprobeINFO3[,6]/1e6;
    

    idx=which(is.na(mprobeINFO1[,5]) | is.na(mprobeINFO1[,6]));
    idx2=which(is.na(mprobeINFO1[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO1[idx,5]=mprobeINFO1[idx,3]-7500;
    mprobeINFO1[idx,6]=mprobeINFO1[idx,3]+7500;
    mprobeINFO1[,8]=-log10(mprobeINFO1[,8]);
    mprobeINFO1[,3]=mprobeINFO1[,3]/1e6;
    mprobeINFO1[,5]=mprobeINFO1[,5]/1e6;
    mprobeINFO1[,6]=mprobeINFO1[,6]/1e6;
    
    idx=which(is.na(mprobeINFO2[,5]) | is.na(mprobeINFO2[,6]));
    idx2=which(is.na(mprobeINFO2[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO2[idx,5]=mprobeINFO2[idx,3]-7500;
    mprobeINFO2[idx,6]=mprobeINFO2[idx,3]+7500;
    mprobeINFO2[,8]=-log10(mprobeINFO2[,8]);
    mprobeINFO2[,3]=mprobeINFO2[,3]/1e6;
    mprobeINFO2[,5]=mprobeINFO2[,5]/1e6;
    mprobeINFO2[,6]=mprobeINFO2[,6]/1e6;

    idx=which(is.na(mprobeINFO3[,5]) | is.na(mprobeINFO3[,6]));
    idx2=which(is.na(mprobeINFO3[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some mprobes' physical positon is missing!");
        quit();
    }
    mprobeINFO3[idx,5]=mprobeINFO3[idx,3]-7500;
    mprobeINFO3[idx,6]=mprobeINFO3[idx,3]+7500;
    mprobeINFO3[,8]=-log10(mprobeINFO3[,8]);
    mprobeINFO3[,3]=mprobeINFO3[,3]/1e6;
    mprobeINFO3[,5]=mprobeINFO3[,5]/1e6;
    mprobeINFO3[,6]=mprobeINFO3[,6]/1e6;

    epXY1=eprobeINFO1[,8];
    mpXY1=mprobeINFO1[,8];
    yMAX1 = ceiling(max(c(pZY1, epXY1,mpXY1), na.rm=T)) + 1;
    epXY2=eprobeINFO2[,8];
    mpXY2=mprobeINFO2[,8];
    yMAX2 = ceiling(max(c(pZY2, epXY2,mpXY2), na.rm=T)) + 1;
    epXY3=eprobeINFO3[,8];
    mpXY3=mprobeINFO3[,8];
    yMAX3 = ceiling(max(c(pZY3, epXY3,mpXY3), na.rm=T)) + 1;
    yMAX=yMAX1
    
    #glist=cbind(eprobeINFO[,2],eprobeINFO[,5:6],as.character(eprobeINFO[,4]),eprobeINFO[,7]);#old
    glist=data1$Gene;
    glist[,2]=glist[,2]/1e6;
    glist[,3]=glist[,3]/1e6;
    colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
    idx=which(is.na(glist[,2]) | is.na(glist[,3]));
    if(length(idx>0)) glist=glist[-idx,];
    generow = GeneRowNum(glist);
    num_row = max(as.numeric(generow$ROW));
    offset_map = ceiling(yMAX1);
    
    offset_map=max(offset_map,num_row*5.5)
    
    offset_probe = yMAX1 / 2.5;
    
    offset_eqtl = ceiling(yMAX1 / 2.5) + 0.5;
    dev_axis = 0.1*yMAX;
    if(dev_axis<1.5) dev_axis = 1.5;
    yaxis.min = -offset_map - offset_eqtl*(neprobePLOT) - offset_eqtl*num_mprb_plot- dev_axis*(neprobePLOT-1)- dev_axis*(num_mprb_plot-1)-yMAX*3+3*dev_axis;
    
    yaxis.max = yMAX + ceiling(offset_probe) + 1;
    
    # scales of x-axis
    idx=match(data1$GWAS[,1],data1$SNP[,1]);
    gwasBP = as.numeric(data1$SNP[idx,3])/1e6;
    min.pos = min(gwasBP);
    max.pos = max(gwasBP);
    start = min(as.numeric(glist[,2]));
    end = max(as.numeric(glist[,3]));
    bp = c(min.pos, max.pos, start, end);
    xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
    #xmax=xmax+(xmax-xmin)*0.1 #extend
    ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
    xlab = paste("Chromosome", chrPLOT, "Mb");
    # plot GWAS p value
    par(mar=c(8,8,4,2), xpd=TRUE)
    plot(gwasBP, pZY1, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    # y1 axis
    devbuf1 = yMAX/4
    axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
    mtext(ylab, side=2, line=3, at=((yMAX1-yMAX2-yMAX3)*3/4), cex=cex_coeff);
    text(xmin, max(pZY1,na.rm=T) , label=trait_name1,col="blue", cex=1, adj=0)
    
    axis.start = 0;
    axis.start = axis.start-yMAX-2*dev_axis;
    pZY2buf=pZY2/max(pZY2)*yMAX+axis.start;
    par(new=TRUE)
    idx=match(data2$GWAS[,1],data2$SNP[,1]);
    gwasBP2 = as.numeric(data2$SNP[idx,3])/1e6;
    plot(gwasBP2, pZY2buf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    devbuf1 = yMAX/4 ;devbuf2 = yMAX2/4; 
    axis(2, at=seq(axis.start,(axis.start+yMAX),devbuf1), labels=round(seq(0,yMAX2,devbuf2),0), las=1, cex.axis=axis);
    segments(xmin, axis.start+yMAX+dev_axis/2, xmax, axis.start+yMAX+dev_axis/2,
    col="dim grey", lty="24", lwd=1)
    text(xmin, axis.start+yMAX-dev_axis, label=trait_name2,col="darkcyan", cex=1, adj=0)
    
    axis.start = axis.start-yMAX-1*dev_axis-0.8;
    pZY3buf=pZY3/max(pZY3)*yMAX+axis.start;
    par(new=TRUE)
    idx=match(data3$GWAS[,1],data3$SNP[,1]);
    gwasBP3 = as.numeric(data3$SNP[idx,3])/1e6;
    plot(gwasBP3, pZY3buf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    devbuf1 = yMAX/4; devbuf2 = yMAX3/4
    axis(2, at=seq(axis.start,(axis.start+yMAX),devbuf1), labels=round(seq(0,yMAX3,devbuf2),0), las=1, cex.axis=axis);
    segments(xmin, axis.start+yMAX+dev_axis/2, xmax, axis.start+yMAX+dev_axis/2,
    col="dim grey", lty="24", lwd=1)
    text(xmin, axis.start+yMAX-dev_axis, label=trait_name3,col="darkolivegreen4", cex=1, adj=0)
    
    eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
    axis.down = offset_eqtl + dev_axis;
    for( k in 1 : neprobePLOT ) {
        axis.start = axis.start - axis.down
        eqtlinfobuf = data1$eQTL[which(data1$eQTL[,1]==eprobePLOT[k]),]
        if(dim(eqtlinfobuf)[1]==0) next;
        pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
        if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
        }
        if(length(which(eheidiprobes1==eprobePLOT[k]))==0) {
            col_eqtl = "maroon"
        } else col_eqtl = "maroon"
        eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
        pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
        idx=match(eqtlinfobuf[,2],data1$SNP[,1]);
        eqtlbp = as.numeric(data1$SNP[idx,3])/1e6;
        probegene = unique(as.character(data1$eSMR[which(data1$eSMR[,1]==eprobePLOT[k]),4]))
        par(new=TRUE)
        pchbuf = 4;
        #if(k%%2==0) pchbuf = 20;
        plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the eQTLs
        colme2e="black"
        if(eprobePLOT[k]==data1$eprobeID) colme2e="darkorchid2"
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=eprobePLOT[k], geneid=probegene)),col=colme2e, cex=1, adj=0)
        if(k==1) labpos=axis.start+offset_eqtl-dev_axis/2;
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,eqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    ypos = labpos;
    mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff,adj=1)
   
    #plot meQTL
    meqtl.lab = expression(paste("-", log[10], "(", italic(P), " meQTL)", sep=""));
    for( k in 1 : nmprobePLOT) {
        axis.start = axis.start - axis.down
        homeid=mprobehome[k]
        if(homeid==2) meqtlinfobuf = data2$meQTL[which(data2$meQTL[,1]==mprobePLOT[k]),];
        if(homeid!=2) meqtlinfobuf = data1$meQTL[which(data1$meQTL[,1]==mprobePLOT[k]),]
        if(dim(meqtlinfobuf)[1]==0) next;
        pvalbuf=-log10(pchisq((meqtlinfobuf[,3]/meqtlinfobuf[,4])^2,1,lower.tail=F));
        if(length(which(is.infinite(pvalbuf)))!=0) {
            maxpvalbuf = max(pvalbuf[-which(is.infinite(pvalbuf))])
            if (maxpvalbuf > 320) {pvalbuf[which(is.infinite(pvalbuf))] = maxpvalbuf + 5
                } else {pvalbuf[which(is.infinite(pvalbuf))] = 325}
        }
        col_eqtl = "navy"
        meqtl.min = 0; meqtl.max = ceiling(max(pvalbuf))
        pvalbuf = pvalbuf/meqtl.max * offset_eqtl + axis.start
        if(homeid==2) {
            idx=match(meqtlinfobuf[,2],data2$SNP[,1]);
             meqtlbp = as.numeric(data2$SNP[idx,3])/1e6;
             mprobegene = unique(as.character(data2$mSMR[which(data2$mSMR[,1]==mprobePLOT[k]),4]))
        } else {
            idx=match(meqtlinfobuf[,2],data1$SNP[,1]);
            meqtlbp = as.numeric(data1$SNP[idx,3])/1e6;
            mprobegene = unique(as.character(data1$mSMR[which(data1$mSMR[,1]==mprobePLOT[k]),4]))
        }
        
        par(new=TRUE)
        pchbuf = 20;
        #if(k%%2==0) pchbuf = 4;
        plot(meqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the meQTLs
        if(homeid==3) colme2e="red"
        if(homeid==1) colme2e="chocolate"
        if(homeid==2) colme2e="darkorange"
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=mprobePLOT[k], geneid=mprobegene)),col=colme2e, cex=1, adj=0)
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = meqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,meqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    ypos = (axis.start - dev_axis-axis.down)*3/4
    mtext(meqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)
    
    
    # plot p value of bTG
    # all the probes
    
    num_gene = dim(generow)[1]
    dist = offset_map/num_row
    for( k in 1 : num_row ) {
        generowbuf = generow[which(as.numeric(generow[,5])==k),]
        xstart = as.numeric(generowbuf[,3])
        xend = as.numeric(generowbuf[,4])
        snbuf = which(xend-xstart< 1e-3)
        if(length(snbuf)>0) {
            xstart[snbuf] = xstart[snbuf] - 0.0025
            xend[snbuf] = xend[snbuf] + 0.0025
        }
        xcenter = (xstart+xend)/2
        xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
        num_genebuf = dim(generowbuf)[1]
        for( l in 1 : num_genebuf ) {
            ofs=0.3
            if(l%%2==0) ofs=-0.8
            m = num_row - k
            ypos = m*dist + yaxis.min
            code = 1
            if(generowbuf[l,2]=="+") code = 2;
            arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
            col=colors()[75], lwd=1)
            movebuf = as.numeric(generowbuf[l,6])*genemove
            text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.6)
        }
    }
    
    # plot the probes
    
    eidx1=match(eprobePLOT,eprobeINFO1[,1])
    eprobein1=eprobeINFO1[eidx1,]
    eidx2=match(eprobePLOT,eprobeINFO2[,1])
    eprobein2=eprobeINFO2[eidx2,]
    eidx3=match(eprobePLOT,eprobeINFO3[,1])
    eprobein3=eprobeINFO3[eidx3,]
    
    midx1=match(mprobePLOT,mprobeINFO1[,1])
    idxx=which(is.na(midx1))
    if(length(idxx)>0) {
        mprobein1=mprobeINFO1[midx1[-idxx],]
    } else {
        mprobein1=mprobeINFO1[midx1,]
    }
    
    midx2=match(mprobePLOT,mprobeINFO2[,1])
    idxx=which(is.na(midx2))
    if(length(idxx)>0) {
        mprobein2=mprobeINFO2[midx2[-idxx],]
    } else {
        mprobein2=mprobeINFO2[midx2,]
    }

    midx3=match(mprobePLOT,mprobeINFO3[,1])
    idxx=which(is.na(midx3))
    if(length(idxx)>0) {
        mprobein3=mprobeINFO3[midx3[-idxx],]
    } else {
        mprobein3=mprobeINFO3[midx3,]
    }
    
    for( k in 1 : dim(mprobein1)[1])  {
        hitflag=FALSE
        if(length(which(mheidiprobes1==mprobein1[k,1]))>0) {
            hitflag=TRUE
            colplot = "navy"; colfont=2; pchbuf=21
        } else {
            colplot = "navy"; colfont=1; pchbuf=1
        }
        if( as.numeric(mprobein1[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        plot_probe(mprobein1, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
        
    }
    
    for( k in 1 : dim(mprobein2)[1])  {
        hitflag=FALSE
        if(length(which(mheidiprobes2==mprobein2[k,1]))>0) {
            hitflag=TRUE
            colplot = "navy"; colfont=2; pchbuf=21
        } else {
            colplot = "navy"; colfont=1; pchbuf=1
        }
        if( as.numeric(mprobein2[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        xcenter = as.numeric(mprobein2[k,3])
        pvalbuf = as.numeric(mprobein2[k,8])
        pvalbuf_=pvalbuf-yMAX2-2*dev_axis;
        strbuf = mprobein2[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }

    }

    for( k in 1 : dim(mprobein3)[1])  {
        hitflag=FALSE
        if(length(which(mheidiprobes3==mprobein3[k,1]))>0) {
            hitflag=TRUE
            colplot = "navy"; colfont=2; pchbuf=21
        } else {
            colplot = "navy"; colfont=1; pchbuf=1
        }
        if( as.numeric(mprobein3[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        
        xcenter = as.numeric(mprobein3[k,3])
        pvalbuf = as.numeric(mprobein3[k,8])
        pvalbuf_=pvalbuf-yMAX3-yMAX2-2*dev_axis;
        strbuf = mprobein3[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }

    }

    for( k in 1 : dim(eprobein1)[1])  {
        hitflag=FALSE
        if(length(which(eheidiprobes1==eprobein1[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else {
            colplot = "maroon"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobein1[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        plot_probe(eprobein1, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
    }
    
    for( k in 1 : dim(eprobein2)[1])  {
        hitflag=FALSE
        if(length(which(eheidiprobes2==eprobein2[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else {
            colplot = "maroon"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobein2[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        xcenter = as.numeric(eprobein2[k,3])
        pvalbuf = as.numeric(eprobein2[k,8])
        pvalbuf_=pvalbuf-yMAX2-2*dev_axis;
        strbuf = eprobein2[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }

    }
    
    for( k in 1 : dim(eprobein3)[1])  {
        hitflag=FALSE
        if(length(which(eheidiprobes3==eprobein3[k,1]))>0) {
            hitflag=TRUE
            colplot = "maroon"; colfont=2; pchbuf=23
        } else {
            colplot = "maroon"; colfont=1; pchbuf=5
        }
        if( as.numeric(eprobein3[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        
        
        xcenter = as.numeric(eprobein3[k,3])
        pvalbuf = as.numeric(eprobein3[k,8])
        pvalbuf_=pvalbuf-yMAX3-yMAX2-2*dev_axis;
        strbuf = eprobein2[k,1]
        par(new=TRUE)
        if(hitflag==TRUE) {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        } else {
            plot(xcenter, pvalbuf_, ylim=c(yaxis.min,yaxis.max),  xlim=c(xmin,xmax),cex.axis=axis,
            xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
        }

    }
    
    #anno probe
    probeINFO=rbind(eprobein1)
    num_anno_probe=dim(probeINFO)[1]
    xcenter = as.numeric(probeINFO[,3])
    xcbuf = xcenter
    xcenter = spread.labs(xcenter[1:num_anno_probe], mindiff=0.08, maxiter=1000, min = xmin, max = xmax-1)
    # adjust the line position
    
    adjflag = rep(0, num_anno_probe)
    if(num_anno_probe>1) {
        dbuf = c(0, xcbuf[1:(num_anno_probe-1)])
        mflag = as.numeric(abs(xcbuf[1:(num_anno_probe)] - dbuf) < 0.01)
        adjflag = as.numeric( mflag | c(mflag[2:num_anno_probe],0) )
    }
    
    for( k in 1 : num_anno_probe)  {
        
        # annotate the probes
      
            ypos = 1.02*yMAX
            strbuf =
            text(xcenter[k], ypos,
            labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),
            list(probeid=as.character(probeINFO[k,1]),
            genename=as.character(probeINFO[k,4]))),
            ylim=c(yaxis.min, yaxis.max),
            srt=30, col=colplot, font=colfont, cex=1, adj=0)
            # plot the lines
            # 1st step
            xstart = xcbuf[k]
            ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
       
                if(adjflag[k]==1) {
                    xstart = (xcbuf[k] + xcenter[k])/2
                    segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
                }
            
            segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
            # 2nd step
            xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*1.01;
            segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
        }
   
    
    # plot the threshold
    # eSMR threshold
    yebuf = -log10(as.numeric(esmr_thresh1)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh1)); dev_anno = yMAX/9;
    
    strbuf = paste("pESMR = ",esmr_thresh1, sep="")
    segments(xmin, yebuf, xmax, yebuf, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno
    } else {
        ythrespos=yebuf+dev_anno
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
    
    
    strbuf = paste("pMSMR = ",msmr_thresh1, sep="")
    segments(xmin, ymbuf, xmax, ymbuf, col="navy", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=ymbuf+dev_anno
    } else {
        ythrespos=ymbuf-dev_anno
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=axis,font=3);
    
    # plot the threshold
    # eSMR threshold
    yebuf = -log10(as.numeric(esmr_thresh2)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh2)); dev_anno = yMAX/9;
    
    strbuf = paste("pESMR = ",esmr_thresh2, sep="")
    segments(xmin, yebuf-yMAX2-2*dev_axis, xmax, yebuf-yMAX2-2*dev_axis, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno-yMAX2-2*dev_axis
    } else {
        ythrespos=yebuf+dev_anno-yMAX2-2*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
    
    
    strbuf = paste("pMSMR = ",msmr_thresh2, sep="")
    segments(xmin, ymbuf-yMAX2-2*dev_axis, xmax, ymbuf-yMAX2-2*dev_axis, col="navy", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=ymbuf+dev_anno-yMAX2-2*dev_axis
    } else {
        ythrespos=ymbuf-dev_anno-yMAX2-2*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=axis,font=3);

    yebuf = -log10(as.numeric(esmr_thresh3)); dev_anno = yMAX/9;
    # mSMR threshold
    ymbuf = -log10(as.numeric(msmr_thresh3)); dev_anno = yMAX/9;
    
    strbuf = paste("pESMR = ",esmr_thresh3, sep="")
    segments(xmin, yebuf-yMAX3-yMAX2-3*dev_axis, xmax, yebuf-yMAX3-yMAX2-3*dev_axis, col="maroon", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=yebuf-dev_anno-yMAX3-yMAX2-3*dev_axis
    } else {
        ythrespos=yebuf+dev_anno-yMAX3-yMAX2-3*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
    
    
    strbuf = paste("pMSMR = ",msmr_thresh3, sep="")
    segments(xmin, ymbuf-yMAX3-yMAX2-3*dev_axis, xmax, ymbuf-yMAX3-yMAX2-3*dev_axis, col="navy", lty=2, lwd=1);
    if(yebuf<ymbuf) {
        ythrespos=ymbuf+dev_anno-yMAX3-yMAX2-3*dev_axis
    } else {
        ythrespos=ymbuf-dev_anno-yMAX3-yMAX2-3*dev_axis
    }
    text(xmax, ythrespos, labels=strbuf, adj=1, col="navy", cex=axis,font=3);
}



