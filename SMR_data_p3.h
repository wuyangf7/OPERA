//
//  SMR_data_p3.h
//  SMR_CPP
//
//  Created by Futao Zhang on 15/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_data_p3__
#define __SMR_CPP__SMR_data_p3__

#include "SMR_data.h"

namespace SMRDATA
{
    typedef struct{
        vector<int> pid;
        string besdpath;
    } F2Prb;
    
    void combineBesd(char* eqtlsmaslstName, char* outFileName,bool save_dense_flag, int cis_itvl, int trans_itvl, float transThres, float restThres, bool genouni,int addn);
    void make_sparse_besd(char* eqtlFileName, char* outFileName, int cis_itvl, int trans_itvl, float transThres, float restThres,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* snplstName,char* problstName, char* snplst2exclde, char* problst2exclde, bool qcflag, int qcmtd, int z_thresh,bool extract_cis_only,char* prbseqregion,double ptech, double pinsnp,double pexsnp,int addn);
    void diff(char* eqtlFileName1,char* eqtlFileName2);
}

#endif /* defined(__SMR_CPP__SMR_data_p3__) */
