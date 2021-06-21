//
//  bfile.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#include "bfile.hpp"

namespace SMRDATA
{
    int getMaxNum(bInfo* bdata,int ldWind, vector<uint64_t> &cols)
    {
        int n=0, loopj=0, preldnum=1;
        long window=ldWind*1000;
        int maxldnum=0;
        cols.resize(bdata->_include.size()+1);
        for(int i=0;i<bdata->_include.size();i++)
        {
             int chri=bdata->_chr[bdata->_include[i]];
            int bpi=bdata->_bp[bdata->_include[i]];
            int ldnum=0;
            for(int j=loopj+1;j<bdata->_include.size();j++)
            {
                int chrj=bdata->_chr[bdata->_include[j]];
                int bpj=bdata->_bp[bdata->_include[j]];
                if(chri==chrj && abs(bpj-bpi)<=window)
                {
                    ldnum++;
                    loopj=j;
                }
                else
                {
                    break;
                }
            }
            ldnum += preldnum-1;
            preldnum=ldnum;
            cols[i+1]=cols[i]+ldnum;
            if(ldnum>maxldnum) maxldnum=ldnum;
        }
        ++maxldnum;
        while((maxldnum>>=1) != 0) n++;
        maxldnum=2<<n;
        return maxldnum;
    }
    void initX(bInfo* bdata,MatrixXd &X, long snpnum)
    {
        vector<uint32_t> snpids(snpnum);
        for(int i=0;i<snpnum;i++) snpids[i]=i;
        make_XMat(bdata, snpids,X,true);
    }
    void ld_report(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr,double maf, bool ldr, bool ldr2, int ldWind)
    {
        bInfo bdata;
        MatrixXd X;
        vector<uint64_t> cols;
        vector<float> lds;
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(chr) extract_snp(&bdata, chr);
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0) filter_snp_maf(&bdata, maf);
        int m=getMaxNum(&bdata,ldWind, cols);
        if(m==1)
        {
            printf("No SNP pair included in %d Kb.\n",ldWind);
            exit(EXIT_FAILURE);
        }
        FILE* outfile=fopen(outFileName, "w");
        if (!(outfile)) {
            printf("Error: Failed to open file %s.\n",outFileName);
            exit(EXIT_FAILURE);
        }
        
        long window=ldWind*1000;
        lds.resize(cols[bdata._include.size()]);
        for(int i=0;i<bdata._include.size();i++)
        {
            if(X.size()==0) initX(&bdata, X,m);
            int start = (i & (m-1));
            VectorXd x=X.col(start);
            int chri=bdata._chr[bdata._include[i]];
            int bpi=bdata._bp[bdata._include[i]];
            string rsi=bdata._snp_name[bdata._include[i]];
        //    clock_t begin_time = clock();
            long tmpcont=0;
            for(int j=1;j<m && i+j<bdata._include.size();j++)
            {
                int chrj=bdata._chr[bdata._include[i+j]];
                int bpj=bdata._bp[bdata._include[i+j]];
                string rsj=bdata._snp_name[bdata._include[i+j]];
                int cur=((start+j) & (m-1));
                if(chri==chrj && abs(bpj-bpi)<=window)
                {
                    VectorXd y=X.col(cur);
                    double ldv=cor(x,y,true);
                    lds[tmpcont]=ldv;
                  //  string tmpstr = atos(chri)+'\t'+rsi+ '\t'+atos(bpi)+'\t'+atos(chrj)+'\t'+rsj+ '\t'+atos(bpj)+'\t'+atos(ldv)+'\n';
                 //   fputs(tmpstr.c_str(),outfile);
                    tmpcont++;
                }
            }
//            VectorXd ldv;
//            ld_calc_o2m(ldv, start, X, true);
//            for(int j=1;j<m && i+j<bdata._include.size();j++)
//            {
//                int chrj=bdata._chr[bdata._include[i+j]];
//                int bpj=bdata._bp[bdata._include[i+j]];
//                string rsj=bdata._snp_name[bdata._include[i+j]];
//                int cur=((start+j) & (m-1));
//                if(chri==chrj && abs(bpj-bpi)<=window)
//                {
//                    string tmpstr = atos(chri)+'\t'+rsi+ '\t'+atos(bpi)+'\t'+atos(chrj)+'\t'+rsj+ '\t'+atos(bpj)+'\t'+atos(ldv(cur))+'\n';
//                    fputs(tmpstr.c_str(),outfile);
//                }
//            }
         //   printf(" cost2: %f ms.\n",float( clock () - begin_time ) /  1000);
            if(i+m<bdata._include.size())
            {
                makex_eigenVector(&bdata, i+m, x, false, true);
                X.col(start)=x;
            }
        }
        fclose(outfile);

    }
}
