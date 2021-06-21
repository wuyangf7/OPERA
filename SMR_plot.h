//
//  SMR_plot.h
//  SMR_CPP
//
//  Created by Futao Zhang on 29/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_plot__
#define __SMR_CPP__SMR_plot__

#include "SMR_data.h"
#include "SMR_data_p1.h"

namespace SMRDATA
{
    void plot_newheidi(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName, double threshpsmrest,bool new_het_mtd,double threshphet, double ld_min,bool sampleoverlap, double pmecs, int minCor);
    
    void plot_triple(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName,char* meqtlFileName, double maf,char* indilstName, char* snplstName,double p_hetero,double ld_top,int m_hetero , int opt_hetero,char* indilst2remove, char* snplst2exclde, double p_smr, char* refSNP, int cis_itvl, char* prbname, int prbWind,bool prbwindFlag, int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName, double pthres_me2esmr,double threshpsmrest,bool new_het_mtd,double threshphet, bool opt, double ld_min, bool sampleoverlap, double pmecs, int minCor,char* targetsnpproblstName);
    void count_trans(char* outFileName,char* eqtlFileName, double p_smr, int cis_itvl);
    void count_cis(char* outFileName,char* eqtlFileName, double p_smr, int cis_itvl);
}

#endif /* defined(__SMR_CPP__SMR_plot__) */
