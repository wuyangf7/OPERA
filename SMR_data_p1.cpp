//
//  SMR_data_p1.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 10/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_data_p1.h"
namespace SMRDATA
{
    void get_top_sets(eqtlInfo* eqtlinfo, vector<string> &prbIds, vector<float> &beta, vector<float> &se, vector<string> &rs, float thres)
    {
        if(eqtlinfo->_rowid.empty())
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                double top_zsqr=0, top_beta=0, top_se=0;
                string snprs="";
                for (int j = 0; j<eqtlinfo->_esi_include.size(); j++)
                {
                    if (fabs(eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]] + 9) > 1e-6)
                    {
                        int snpchr=eqtlinfo->_esi_chr[eqtlinfo->_esi_include[j]];
                        string rstmp=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[j]];
                        if(snpchr==probechr)
                        {
                            float bxz=eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float sexz=eqtlinfo->_sexz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float zxz=bxz/sexz;
                            zxz*=zxz;
                            if(zxz-top_zsqr>1e-8) {
                                top_zsqr=zxz;
                                top_beta=bxz;
                                top_se=sexz;
                                snprs=rstmp;
                            }
                        }
                    }
                }
                double pval=pchisq(top_zsqr, 1);
                if(pval<thres)
                {
                    prbIds.push_back(probeId);
                    beta.push_back(top_beta);
                    se.push_back(top_se);
                    rs.push_back(snprs);
                }
            }
            
        }
        else
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                double top_zsqr=0, top_beta=0, top_se=0;
                string snprs="";
                
                uint64_t beta_start=eqtlinfo->_cols[eqtlinfo->_include[i]<<1];
                uint64_t se_start=eqtlinfo->_cols[1+(eqtlinfo->_include[i]<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=eqtlinfo->_rowid[beta_start+j];
                    int snpchr=eqtlinfo->_esi_chr[ge_rowid];
                    string snptmp=eqtlinfo->_esi_rs[ge_rowid];
                    if(snpchr==probechr )
                    {
                        float bxz=eqtlinfo->_val[beta_start+j];
                        float sexz=eqtlinfo->_val[se_start+j];
                        float zxz=bxz/sexz;
                        zxz*=zxz;
                        if(zxz-top_zsqr>1e-8) {
                            top_zsqr=zxz;
                            top_beta=bxz;
                            top_se=sexz;
                            snprs=snptmp;
                        }
                    }
                }
                
                double pval=pchisq(top_zsqr, 1);
                if(pval<thres)
                {
                    prbIds.push_back(probeId);
                    beta.push_back(top_beta);
                    se.push_back(top_se);
                    rs.push_back(snprs);
                }
                
            }
        }

    }
    
    void get_thres_sets(eqtlInfo* eqtlinfo, vector<string> &prbIds, vector<float> &beta, vector<float> &se, vector<string> &rs, float thres)
    {
        if(eqtlinfo->_rowid.empty())
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                for (int j = 0; j<eqtlinfo->_esi_include.size(); j++)
                {
                    if (fabs(eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]] + 9) > 1e-6)
                    {
                        int snpchr=eqtlinfo->_esi_chr[eqtlinfo->_esi_include[j]];
                        string rstmp=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[j]];
                        if(snpchr==probechr)
                        {
                            float bxz=eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float sexz=eqtlinfo->_sexz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float zxz=bxz/sexz;
                            zxz*=zxz;
                            double pval=pchisq(zxz, 1);
                            if(pval<thres) {
                                prbIds.push_back(probeId);
                                beta.push_back(bxz);
                                se.push_back(sexz);
                                rs.push_back(rstmp);
                            }
                        }
                    }
                }
            }
            
        }
        else
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                
                uint64_t beta_start=eqtlinfo->_cols[eqtlinfo->_include[i]<<1];
                uint64_t se_start=eqtlinfo->_cols[1+(eqtlinfo->_include[i]<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=eqtlinfo->_rowid[beta_start+j];
                    int snpchr=eqtlinfo->_esi_chr[ge_rowid];
                    string snptmp=eqtlinfo->_esi_rs[ge_rowid];
                    if(snpchr==probechr )
                    {
                        float bxz=eqtlinfo->_val[beta_start+j];
                        float sexz=eqtlinfo->_val[se_start+j];
                        float zxz=bxz/sexz;
                        zxz*=zxz;
                        double pval=pchisq(zxz, 1);
                        if(pval<thres) {
                            prbIds.push_back(probeId);
                            beta.push_back(bxz);
                            se.push_back(sexz);
                            rs.push_back(snptmp);
                        }
                    }
                }
            }
        }
        
    }

    
    long est_n(vector<float> &beta,vector<float> &se)
    {
        long n;
        VectorXd Y(se.size());
        MatrixXd X(beta.size(),2);
        for(int i=0;i<se.size();i++)
        {
            Y(i)=se[i]*se[i];
            X(i,0)=1.0;
            X(i,1)=beta[i]*beta[i]-Y(i);
        }
        
        MatrixXd XtX_i=(X.transpose()*X).inverse();
        VectorXd w_hat=XtX_i*X.transpose()*Y;
        n=-1/w_hat[1];
        return n;
    }
    void est_effect_splsize(char* eqtlsmaslstName, char* eqtlFileName, char* snplstName,char* problstName,char* snplst2exclde, char* problst2exclde,float thres)
    {
        vector<string> prbIds;
        vector<float> beta;
        vector<float> se;
        vector<string> rs;
        
        vector<string> smasNames;
        if(eqtlsmaslstName!=NULL)
        {
            read_smaslist(smasNames, string(eqtlsmaslstName));
            if(smasNames.size()==0) throw("No eqtl summary file list in [ "+ string(eqtlsmaslstName)  +" ]");
            for(int ii=0;ii<smasNames.size();ii++)
            {
                eqtlInfo eqtlinfo;
                read_esifile(&eqtlinfo, smasNames[ii]+".esi");
                if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
                if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
                read_epifile(&eqtlinfo, smasNames[ii]+".epi");
                if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
                if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
                read_besdfile(&eqtlinfo, smasNames[ii]+".besd");
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    printf("No data included from %s in the analysis.\n",smasNames[ii].c_str());
                    exit(EXIT_FAILURE);
                }
                //get_top_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
                get_thres_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
            
            }

        } else if(eqtlFileName!=NULL)
        {
            eqtlInfo eqtlinfo;
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s in the analysis.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            //get_top_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
            get_thres_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
           
        }else {
            throw("Error: please input eQTL summary data files list or eQTL summary data by the flag --beqtl-summaries or --beqtl-summary.");
        }
        
        cout<<prbIds.size()<<" eQTLs are inluded to estimate the effective population size."<<endl;
        ///////
        
        string tstfile = "tst.top.bse.txt";
        ofstream tst(tstfile.c_str());
        if (!tst) throw ("Error: can not open the fam file " + tstfile + " to save!");
        
        tst << "Probe" <<'\t'<<"SNP"<<'\t'<< "b" <<'\t' << "se"<<'\n';
        
        for (int i = 0;i <beta.size(); i++) {
            tst<<prbIds[i]<<'\t'<<rs[i]<<'\t'<<beta[i]<<'\t'<<se[i]<< '\n';
        }
        
        tst.close();
        
        //////
        
        long n=est_n(beta,se);
        cout<<"The estimated effective population size is: "<<n<<endl;
    }
    
    void read_esmr(eSMRrlt* erlt, string smrFileName)
    {
        cout << "Reading SMR result information from [" + smrFileName + "]." << endl;
       
        ifstream flptr;
        flptr.open(smrFileName.c_str());
        if (!flptr) throw ("Error: can not open the file [" + smrFileName + "] to read.");
        
       
        char tbuf[MAX_LINE_SIZE];
        int lineNum(0);
       
        flptr.getline(tbuf,MAX_LINE_SIZE);// the header
        while(!flptr.eof())
        {
            string tmpStr;
            flptr.getline(tbuf,MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                istringstream iss(tbuf);
                iss>>tmpStr;
                erlt->_Expo_id.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Expo_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Expo_gene.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Expo_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Outco_id.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Outco_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Outco_gene.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Outco_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_snp_rs.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_snp_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_snp_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                to_upper(tmpStr);
                erlt->_snp_a1.push_back(tmpStr.c_str());
                iss>>tmpStr;
                to_upper(tmpStr);
                erlt->_snp_a2.push_back(tmpStr.c_str());
                iss>>tmpStr;
                erlt->_snp_frq.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr; // b_outcome
                iss>>tmpStr; // se_outcome
                iss>>tmpStr; // p_outcome
                iss>>tmpStr; // b_exposure
                iss>>tmpStr; // se_exposure
                iss>>tmpStr; // p_exposure
                iss>>tmpStr;
                erlt->_b.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_se.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_p_smr.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                if(!tmpStr.compare("NA")) erlt->_p_heidi.push_back(-9);
                else  erlt->_p_heidi.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                if(!tmpStr.compare("NA")) erlt->_nsnp.push_back(-9);
                else  erlt->_nsnp.push_back(atoi(tmpStr.c_str()));
                erlt->_include.push_back(lineNum);
                lineNum++;
            }
        }
        erlt->lineNum=lineNum;
        flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + smrFileName + "]." << endl;
        
    }
    
    
    void read_probevarfile(eqtlInfo* eqtlinfo, string vpFileName)
    {
        cout << "Reading variance information from [" + vpFileName + "]." << endl;
        
        ifstream flptr;
        flptr.open(vpFileName.c_str());
        if (!flptr) throw ("Error: can not open the file [" + vpFileName + "] to read.");
        
        
        char tbuf[MAX_LINE_SIZE];
        int lineNum(0);
        vector<string> tmp_prbid;
        vector<double> tmp_var;
        
        while(!flptr.eof())
        {
            string tmpStr;
            flptr.getline(tbuf,MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                istringstream iss(tbuf);
                iss>>tmpStr;
                tmp_prbid.push_back(tmpStr.c_str());
                iss>>tmpStr;
                tmp_var.push_back(atof(tmpStr.c_str()));
                lineNum++;
            }
        }
        flptr.close();
        
        
        eqtlinfo->_epi_var.resize(eqtlinfo->_epi_prbID.size());
        vector<int> idx;
        match_only(eqtlinfo->_epi_prbID, tmp_prbid, idx);
        if(idx.size()!=eqtlinfo->_epi_prbID.size())
        {
            cout<<"Some Probes in summary data are not in variance data!"<<endl;
            exit(1);
        }
        for(int i=0;i<idx.size();i++)
        {
            eqtlinfo->_epi_var[i]=tmp_var[idx[i]];
        }
        
        cout << idx.size() << " probes variance info to be included from [" + vpFileName + "]." << endl;

    }
    

    void iternal_test(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,int cis_itvl,char* smrFileName)
    {
        setNbThreads(thread_num);
        
        eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        double threshold= chi_val(1,p_hetero);
        cis_itvl=cis_itvl*1000;
       
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe here presumes the probe list should contain both probes of exposure dataset and probes of outcome dataset.\n If you want to only extract probes from one dataset please include these probes in the file and all the probes of the other dataset as well.\n"<<endl;
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&etrait, problstName);
        if(problst2exclde != NULL) exclude_prob(&etrait, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&etrait, oproblstName);
        if(oproblst2exclde != NULL) exclude_prob(&etrait, oproblst2exclde);
        
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(&bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
            allele_check(&bdata, &etrait, &esdata);
            // if no snp left after check
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, &etrait, &esdata);
            }
            

        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName2);
            exit(EXIT_FAILURE);
        }
        
        
        
        
        for(int i=0;i<etrait._include.size();i++)
            etrait._probe_name_map.insert(pair<string, int>(etrait._epi_prbID[etrait._include[i]],etrait._include[i]));
        
        for(int i=0;i<esdata._include.size();i++)
            esdata._probe_name_map.insert(pair<string, int>(esdata._epi_prbID[esdata._include[i]],esdata._include[i]));
        
        cout<<endl<<"Performing interanl test..... "<<endl;
        vector<long> out_probid;
        vector<string> bxy;
        vector<string> sexy;
        vector<string> pxy;
        vector<string> bgwas;
        vector<string> segwas;
        vector<string> beqtl;
        vector<string> seeqtl;
        vector<string> pgwas;
        vector<string> peqtl;
        vector<string> rsid;
        vector<string> rschr;
        vector<string> rsbp;
        vector<string> rsfreq;
        
        vector<string> rsa1;
        vector<string> rsa2;
        vector<string> prb1;
        vector<string> nsnp_test1;
        
        vector<string> etrait_id;
        vector<string> etrait_bp;
        vector<string> etrait_gene;
        vector<string> etrait_chr;
        
        eSMRrlt erlt;
        read_esmr(&erlt,smrFileName);
        
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        for( int lid=0;lid<erlt.lineNum;lid++)
        {
            progr1=1.0*lid/erlt.lineNum;
            if(progr1-progr0-0.05>1e-6 || lid+1==erlt.lineNum)
            {
                if(lid+1==erlt.lineNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }

            
            string traitname=erlt._Outco_id[lid];
            
            map<string, int>::iterator iter;
            iter =  etrait._probe_name_map.find(traitname);
            if (iter == etrait._probe_name_map.end()) throw("Could not find probe "+traitname+" in "+eqtlFileName+" dataset, please check!");
            uint32_t ii = iter->second;
            
            gwasData gdata;
            gdata.allele_1.resize(etrait._esi_include.size());
            gdata.allele_2.resize(etrait._esi_include.size());
            gdata.byz.resize(etrait._esi_include.size());
            gdata.seyz.resize(etrait._esi_include.size());
            gdata.freq.resize(etrait._esi_include.size());
            gdata.pvalue.resize(etrait._esi_include.size());
            gdata.splSize.resize(etrait._esi_include.size());
            
            
          
           // cout<<"\nPerforming analysis of eTrait [ "+traitname+" ]..."<<endl;
            gdata._include.clear();
            gdata.snpName.clear();
            int count=0;
            if(etrait._rowid.empty())
            {
                for (int j = 0; j<etrait._esi_include.size(); j++)
                {
                    if (fabs(etrait._bxz[ii][etrait._esi_include[j]] + 9) > 1e-6)
                    {
                        gdata.byz[count]=etrait._bxz[ii][etrait._esi_include[j]];
                        gdata.seyz[count]=etrait._sexz[ii][etrait._esi_include[j]];
                        gdata.snpName.push_back(etrait._esi_rs[etrait._esi_include[j]]);
                        gdata.allele_1[count]=etrait._esi_allele1[etrait._esi_include[j]];
                        gdata.allele_2[count]=etrait._esi_allele2[etrait._esi_include[j]];
                        gdata._include.push_back(etrait._esi_include[j]); // row id selected
                        count++;
                    }
                }
            }
            else
            {
                uint64_t beta_start=etrait._cols[ii<<1];
                uint64_t se_start=etrait._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=etrait._rowid[beta_start+j];
                    if(binary_search(etrait._esi_include.begin(), etrait._esi_include.end(), ge_rowid))
                    {
                        gdata.byz[count]=etrait._val[beta_start+j];
                        gdata.seyz[count]=etrait._val[se_start+j];
                        gdata.snpName.push_back(etrait._esi_rs[ge_rowid]);
                        gdata.allele_1[count]=etrait._esi_allele1[ge_rowid];
                        gdata.allele_2[count]=etrait._esi_allele2[ge_rowid];
                        gdata._include.push_back(ge_rowid);
                        count++;
                    }
                    
                }
            }
            if(gdata.snpName.size()< m_hetero)
            {
                cout<<gdata.snpNum<<" common SNPs (less than parameter m_hetero: "+atos(m_hetero)+" ) are included from eTrait [ "+traitname+" ] summary."<<endl;
                continue;
            }
            gdata.snpNum=gdata.snpName.size();
            //cout<<gdata.snpNum<<" common SNPs are included from eTrait [ "+traitname+" ] summary."<<endl;
            
            vector<double> bxz;
            vector<double> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            vector<int> snpchrom;
            
            vector<string> allele1;
            vector<string> allele2;
            vector<uint32_t> bpsnp;
            vector<double> freq;
            
            vector<double> byz;
            vector<double> seyz;
            VectorXd zsxz;
            
            vector<int> sn_ids;
            
            VectorXd _byz;
            VectorXd _seyz;
            VectorXd _bxz;
            VectorXd _sexz;
            VectorXd _zsxz;
            MatrixXd _X;
            
            VectorXd ld_v;
            MatrixXd _LD_heidi;
            MatrixXd _X_heidi;
            
           
            string expo_name=erlt._Expo_id[lid];
            string refSNP=erlt._snp_rs[lid];
            map<string, int>::iterator iter1;
            iter1 =  esdata._probe_name_map.find(expo_name);
            if (iter1 == esdata._probe_name_map.end())
            {
                cout<<"Warning: Could not find probe "+expo_name+" in "+eqtlFileName2+" dataset, please check!"<<endl;
                continue;
            }
            int i = iter1->second;           
                
                //extract info from eqtl summary and gwas summary
                bxz.clear();
                sexz.clear();
                curId.clear(); // is the idxes of bfile._include not the values of
                eName.clear();
                snpchrom.clear();
                byz.clear();
                seyz.clear();
                allele1.clear();
                allele2.clear();
                bpsnp.clear();
                freq.clear();
                long maxid =-9;
                int probebp=esdata._epi_bp[i];
                int probechr=esdata._epi_chr[i];
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == etrait._esi_include.size()
                    {
                        if (fabs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            int snpbp=esdata._esi_bp[j];
                            int snpchr=esdata._esi_chr[j];
                            
                            int etrait_rid=etrait._esi_include[j];
                            long pos=find(gdata._include.begin(), gdata._include.end(), etrait_rid)-gdata._include.begin();
                            
                            if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl && pos!=gdata._include.size())
                            {
                                bxz.push_back(esdata._bxz[i][j]);
                                sexz.push_back(esdata._sexz[i][j]);
                                byz.push_back(gdata.byz[pos]);
                                seyz.push_back(gdata.seyz[pos]);
                                curId.push_back(j);
                                eName.push_back(esdata._esi_rs[j]);
                                snpchrom.push_back(esdata._esi_chr[j]);
                                if(esdata._esi_rs[j]==refSNP) maxid=(eName.size()-1);
                                freq.push_back(bdata._mu[bdata._include[j]]/2);
                                allele1.push_back(esdata._esi_allele1[j]);
                                allele2.push_back(esdata._esi_allele2[j]);
                                bpsnp.push_back(esdata._esi_bp[j]);
                            }
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata._rowid[beta_start+j];
                        int snpbp=esdata._esi_bp[ge_rowid];
                        int snpchr=esdata._esi_chr[ge_rowid];
                        
                        int etrait_rid=etrait._esi_include[ge_rowid];
                        long pos=find(gdata._include.begin(), gdata._include.end(), etrait_rid)-gdata._include.begin();
                        
                        if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl && pos!=gdata._include.size())
                        {
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[pos]);
                            seyz.push_back(gdata.seyz[pos]);
                            curId.push_back(ge_rowid);
                            eName.push_back(esdata._esi_rs[ge_rowid]);
                            snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                            if(esdata._esi_rs[ge_rowid]==refSNP) maxid=(eName.size()-1);
                            allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            freq.push_back(bdata._mu[bdata._include[ge_rowid]]/2);
                        }
                    }
                }
                if(maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (bxz.size() == 0) continue;
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                double chisqyz = byz[maxid] / seyz[maxid];
                double pyz_val = pchisq(chisqyz*chisqyz, 1);
                
                out_probid.push_back(i);
                bxy.push_back(dtosf(bxy_val));
                sexy.push_back(dtosf(sexy_val));
                pxy.push_back(dtos(pxy_val));
                bgwas.push_back(dtosf(byz[maxid]));
                segwas.push_back(dtosf(seyz[maxid]));
                beqtl.push_back(dtosf(bxz[maxid]));
                seeqtl.push_back(dtosf(sexz[maxid]));
                pgwas.push_back(dtos(pyz_val));
                peqtl.push_back(dtos(pxz_val));
                rsid.push_back(eName[maxid]);
                rschr.push_back(atos(snpchrom[maxid]));
                rsbp.push_back(itos(bpsnp[maxid]));
                rsfreq.push_back(dtosf(freq[maxid]));
                
                rsa1.push_back(allele1[maxid]);
                rsa2.push_back(allele2[maxid]);
            
            etrait_id.push_back(etrait._epi_prbID[ii]);
            etrait_chr.push_back(atos(etrait._epi_chr[ii]));
            etrait_gene.push_back(etrait._epi_gene[ii]);
            etrait_bp.push_back(atos(etrait._epi_bp[ii]));
            
                    //extract info from reference
                    make_XMat(&bdata,curId, _X); //_X: one row one individual, one column one SNP
                    ld_calc_o2m(ld_v,maxid,_X);
                    
            
                    
                    sn_ids.clear(); //increase order
                    if(fabs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                    else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                    
                    if(sn_ids.size() < m_hetero)
                    {
                        prb1.push_back(string("NA"));
                        nsnp_test1.push_back( string("NA"));
                        continue;
                    }
            
                    _byz.resize(sn_ids.size());
                    _seyz.resize(sn_ids.size());
                    _bxz.resize(sn_ids.size());
                    _sexz.resize(sn_ids.size());
                    _zsxz.resize(sn_ids.size());
                    _X_heidi.resize(_X.rows(), sn_ids.size());
                    
                    #pragma omp parallel for
                    for(int j=0;j<sn_ids.size();j++)
                    {
                        _byz[j]=byz[sn_ids[j]];
                        _seyz[j]=seyz[sn_ids[j]];
                        _bxz[j]=bxz[sn_ids[j]];
                        _sexz[j]=sexz[sn_ids[j]];
                        _zsxz[j]=zsxz[sn_ids[j]];
                        _X_heidi.col(j)=_X.col(sn_ids[j]);
                    }
                    _X.resize(0,0);
                    cor_calc(_LD_heidi, _X_heidi);
                    
            
                    _X_heidi.resize(0,0);
                    
                    long nsnp = sn_ids.size();
                    double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
                    
                    prb1.push_back(dtos(pdev));
                    if(nsnp>0) nsnp_test1.push_back(itos(nsnp));
                    else nsnp_test1.push_back(string("NA"));
            
            free_gwas_data( &gdata);
            
        }
        
        if(bxy.size()>0)
        {
            string smrfile = string(outFileName)+".smr";
            ofstream smr(smrfile.c_str());
            if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
            
            smr << "Expo_ID" <<'\t'<< "Expo_Chr" <<'\t' << "Expo_Gene"  << '\t' << "Expo_bp" << '\t'<< "Outco_ID" <<'\t'<< "Outco_Chr" <<'\t' << "Outco_Gene"  << '\t' << "Outco_bp" << '\t'<< "SNP"<< '\t' << "SNP_Chr"<< '\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"Freq"<<'\t'<<"b_Outco"<<'\t'<<"se_Outco"<<'\t'<< "p_Outco" << '\t'<<"b_Expo"<<'\t'<<"se_Expo"<<'\t'<< "p_Expo" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HEIDI"<< "\t" << "nsnp" << '\n';
            
            for (int i = 0;i <bxy.size(); i++) {
                smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<etrait_id[i]<<'\t'<<etrait_chr[i]<<'\t'<<etrait_gene[i]<<'\t'<<etrait_bp[i]<<'\t'<<rsid[i]<<'\t'<<rschr[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<rsfreq[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
            }
            cout<<"Internal Test finished.\nInternal Test results of "<<bxy.size()<<" probes have been saved in the file [" + smrfile + "]."<<endl;
            smr.close();
            
        }else cout<<"Internal Test finished.\nInternal Test results of "<<bxy.size()<<" probes have been saved."<<endl;
        


    }
    
    void make_cojo(char* outFileName,char* eqtlFileName, char* snplstName,char* snplst2exclde, char* problstName, char* problst2exclde, char* genelistName, bool bFlag)
    {
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
            if(genelistName != NULL) extract_prob_by_gene(&eqtlinfo, genelistName);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s in the analysis.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            
            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.");
        
        vector<int> out_esi_id;
        vector<int> out_epi_id;
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;
        if(eqtlinfo._valNum==0)
        {
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                for(uint32_t j=0;j<eqtlinfo._snpNum;j++)
                {
                    double beta=eqtlinfo._bxz[i][j];
                    double se=eqtlinfo._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                        out_esi_id.push_back(j);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                   
                }
            }
        }
        else
        {
            if(eqtlinfo._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.");
            }
            
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                uint64_t proid=eqtlinfo._include[i];
                uint64_t pos=eqtlinfo._cols[proid<<1];
                uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=eqtlinfo._val[pos+j];
                    double se=eqtlinfo._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    
                        out_esi_id.push_back(eqtlinfo._rowid[pos+j]);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                   
                }
            }
        }
        
        string smrfile = string(outFileName)+".raw";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        
        smr << "ProbeID"<<'\t' << "SNP" <<'\t'<< "Chr" <<'\t' << "A1" << '\t'<< "A2"<< '\t' << "b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
        
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" items have been saved in the plaint text file [" + smrfile + "]."<<endl;
        
    }
    /*
    void standardization(char* outFileName, char* eqtlFileName,bool bFlag,char* freqName, char* vpFileName)
    {
        eqtlInfo esdata;
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(freqName==NULL && vpFileName==NULL) throw("Error: please input feq data or variance data for standardisation by the flag --freq or --probe-var.");
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        read_epifile(&esdata, string(eqtlFileName)+".epi");
         read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        if(vpFileName!=NULL)
        {
            
            read_probevarfile(&esdata, string(vpFileName));
            for(int i=0;i<esdata._probNum;i++)
            {
                double prbvar_sqrt=sqrt(esdata._epi_var[i]);
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)
                    {
                        if (fabs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            
                            float bxz=esdata._bxz[i][j];
                            float sexz=esdata._sexz[i][j];
                            esdata._bxz[i][j]=bxz/prbvar_sqrt;
                            esdata._sexz[i][j]=sexz/prbvar_sqrt;
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(uint64_t j=0;j<numsnps;j++)
                    {
                        float bxz=esdata._val[beta_start+j];
                        float sexz=esdata._val[se_start+j];
                        esdata._val[beta_start+j]=bxz/prbvar_sqrt;
                        esdata._val[se_start+j]=sexz/prbvar_sqrt;
                    }
                }
                
            }

        }
        else if(freqName!=NULL)
        {
            
            update_freq(&esdata, string(freqName));
            for(int i=0;i<esdata._probNum;i++)
            {
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)
                    {
                        if (fabs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            
                            float bxz=esdata._bxz[i][j];
                            float sexz=esdata._sexz[i][j];
                            float p=esdata._esi_freq[j];
                            float z=bxz/sexz;
                            float b=z/sqrt(2*p*(1-p)*(n+z*z));
                            float se=1/sqrt(2*p*(1-p)*(n+z*z));
                            esdata._bxz[i][j]=b;
                            esdata._sexz[i][j]=se;
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(uint64_t j=0;j<numsnps;j++)
                    {
                        uint64_t ge_rowid=esdata._rowid[beta_start+j];
                        float bxz=esdata._val[beta_start+j];
                        float sexz=esdata._val[se_start+j];
                        float p=esdata._esi_freq[ge_rowid];
                        float z=bxz/sexz;
                        float b=z/sqrt(2*p*(1-p)*(n+z*z));
                        float se=1/sqrt(2*p*(1-p)*(n+z*z));
                        esdata._val[beta_start+j]=b;
                        esdata._val[se_start+j]=se;
                    }
                }
                
            }

        }
       write_besd(outFileName, &esdata);
    }
   */

    void lookup(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName,char* genelistName, double plookup,bool bFlag, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl)
    {
        string logstr;
        int flag4chr=0;
        if(chr!=0) flag4chr++;
        if(prbchr!=0 || snpchr!=0) flag4chr++;
        if(flag4chr==2)
        {
            logstr="WARNING: --chr is not surpposed to use together with --probe-chr or --snp-chr. --chr will be disabled.\n";
            chr=0;
            fputs(logstr.c_str(), stdout);
        }       
        
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
           
            
           read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s in the analysis.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }

            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
        
        vector<int> out_esi_id;
        vector<int> out_epi_id;
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;
        if(eqtlinfo._valNum==0)
        {
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                for(uint32_t j=0;j<eqtlinfo._snpNum;j++)
                {
                    double beta=eqtlinfo._bxz[i][j];
                    double se=eqtlinfo._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=plookup)
                    {
                        out_esi_id.push_back(j);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                    }
                }
            }
        }
        else
        {
            if(eqtlinfo._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.\n");
            }
            
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                uint64_t proid=eqtlinfo._include[i];
                uint64_t pos=eqtlinfo._cols[proid<<1];
                uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=eqtlinfo._val[pos+j];
                    double se=eqtlinfo._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=plookup)
                    {
                        out_esi_id.push_back(eqtlinfo._rowid[pos+j]);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                    }
                }
            }
        }
        
        string smrfile = string(outFileName)+".txt";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
        
        smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' <<"Freq"<<'\t'<< "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"Orientation"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<((eqtlinfo._esi_freq[out_esi_id[i]]+9>1e-6)?atos(eqtlinfo._esi_freq[out_esi_id[i]]):"NA")<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_orien[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
        
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" SNPs have been saved in the file [" + smrfile + "]."<<endl;
        
    }
    
    void rm_unmatched_snp(gwasData* gdata, eqtlInfo* esdata)
    {
        cout<<"Checking the consistency of SNP alleles between GWAS summary data and eQTL summary data..."<<endl;
        // get the common SNPs
        vector<string> slctSNPs;
        vector<int> id_unmatched_esd;
        vector<int> id_unmatched_gwas;
        vector<int> gdId;
        vector<int> edId;
        edId.clear();
        long pre_size=esdata->_esi_include.size();
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(essnp, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }
        
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        id_unmatched_esd.clear();
        id_unmatched_gwas.clear();
        for (int i = 0; i<edId.size(); i++)
        {
            string ga1, ga2, ea1, ea2, rs1, rs2;
            rs1 = gdata->snpName[gdId[i]];
            rs2 = esdata->_esi_rs[edId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ga1 && ea2 == ga2)
            {
                
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                
            }else {
                id_unmatched_esd.push_back(edId[i]);
                 id_unmatched_gwas.push_back(gdId[i]);
            }
        }
        if(id_unmatched_esd.size()>0)
        {
            vector<int> tmpvec;
            esdata->_esi_include.swap(tmpvec);
            set_complement(id_unmatched_esd,tmpvec, esdata->_esi_include);
        }
        if(id_unmatched_gwas.size()>0)
        {
            vector<int> tmpvec;
            gdata->_include.swap(tmpvec);
            set_complement(id_unmatched_gwas,tmpvec, gdata->_include);
        }
        cout<<id_unmatched_esd.size()<<" SNPs failed in allele check and have been excluded. Total "<<gdata->_include.size()<<" SNPs left in GWAS summary dataset and "<<esdata->_esi_include.size()<<" SNPs left in eQTL summary dtaset."<<endl;
       
    }
    
    void rm_unmatched_snp(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        vector<int> id_unmatched_esd;
        vector<int> id_unmatched_etrait;
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP allelesbetween 2 eQTL summary datasets. ";
        cout<<logstr<<endl;
        vector<string> etsnp;
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum || etrait->_esi_include.size()<etrait->_snpNum)
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            for(int i=0;i<etrait->_esi_include.size();i++) etsnp.push_back(etrait->_esi_rs[etrait->_esi_include[i]]);
            match_only(essnp, etsnp, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etsnp[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, etrait->_esi_rs, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etrait->_esi_rs[gdId[i]]);
        }
        
        
        //alleles check
        match(slctSNPs, esdata->_esi_rs, edId);
        gdId.clear();
        match(slctSNPs, etrait->_esi_rs, gdId);
        cmmnSNPs.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string ta1, ta2, ea1, ea2;
            
            ta1 = etrait->_esi_allele1[gdId[i]];
            ta2 = etrait->_esi_allele2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ta1 && ea2 == ta2)
            {
            }
            else if(ea1 == ta2 && ea2 == ta1)
            {
                if(etrait->_val.size()>0)
                {
                    int count=0;
                    for(int j=0;j<etrait->_rowid.size();j++)
                    {
                        if(etrait->_rowid[j]==gdId[i])
                        {
                            count++;
                            if(count & 1)
                                etrait->_val[j]=-etrait->_val[j];
                        }
                    }
                }
                else
                {
                    for(int j=0;j<etrait->_include.size();j++) if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 ) etrait->_bxz[j][gdId[i]]=-etrait->_bxz[j][gdId[i]];
                }
                
            }
            else
            {
                id_unmatched_esd.push_back(edId[i]);
                id_unmatched_etrait.push_back(gdId[i]);
            }
        }
       
        if(id_unmatched_esd.size()>0)
        {
            vector<int> tmpvec;
            esdata->_esi_include.swap(tmpvec);
            set_complement(id_unmatched_esd,tmpvec, esdata->_esi_include);
        }
        if(id_unmatched_etrait.size()>0)
        {
            vector<int> tmpvec;
            etrait->_esi_include.swap(tmpvec);
            set_complement(id_unmatched_etrait,tmpvec, etrait->_esi_include);
        }

        cout<<id_unmatched_esd.size()<<" SNPs failed in allele check and have been excluded. Total "<<etrait->_esi_include.size()<<" SNPs left in one eQTL summary dataset and "<<esdata->_esi_include.size()<<" SNPs left in another eQTL summary dtaset."<<endl;
        
    }

    void read_gene_anno(char* geneAnnoName,vector<int> &chr, vector<string> &genename,vector<int> &start,vector<int> &end)
    {
       
        ifstream flptr;
        flptr.open(geneAnnoName);
        if (!flptr) throw ("Error: can not open the file [" + string(geneAnnoName) + "] to read.");
      
        cout << "Reading gene annotation information from [" + string(geneAnnoName) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        
        while(!flptr.eof() )
        {
            string tmpStr;
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                istringstream iss(buf);
                iss>>tmpStr; //chr
                if(tmpStr=="X" || tmpStr=="x") chr.push_back(23);
                else if(tmpStr=="Y" || tmpStr=="y") chr.push_back(24);
                else chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // start
                start.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // end
                end.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // gene
                genename.push_back(tmpStr.c_str());
                lineNum++;
                
            }
        }
        cout << lineNum << " gene annotation infomation to be included from [" + string(geneAnnoName) + "]." << endl;
        flptr.close();

    }
    void read_gene_anno_strand(char* geneAnnoName,vector<int> &chr, vector<string> &genename,vector<int> &start,vector<int> &end, vector<string> &strand)
    {
        
        ifstream flptr;
        flptr.open(geneAnnoName);
        if (!flptr) throw ("Error: can not open the file [" + string(geneAnnoName) + "] to read.");
        
        cout << "Reading gene annotation information from [" + string(geneAnnoName) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        vector<string> vs_buf;
        int col_num = -9;
     
        while(!flptr.eof() )
        {
            string tmpStr;
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int tmpnum = split_string(buf, vs_buf, ", \t\n");
                if(tmpnum!=4 && tmpnum!=5) {
                    printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(col_num==-9)
                {
                    col_num=tmpnum;
                } else {
                    if(col_num!=tmpnum)
                    {
                        printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                start.push_back(atoi(vs_buf[1].c_str()));
                end.push_back(atoi(vs_buf[2].c_str()));
                genename.push_back(vs_buf[3].c_str());
                if(col_num==5) strand.push_back(vs_buf[4].c_str());
                lineNum++;
                
            }
        }
        cout << lineNum << " gene annotation infomation to be included from [" + string(geneAnnoName) + "]." << endl;
        flptr.close();
        
    }


    void plot(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName)
    {
        
        setNbThreads(thread_num);
        if(!prbwindFlag) throw("Error: please input probe window by the flag --probe-wind.");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data  by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data  by the flag --eqtl-summary.");
        if(geneAnnoName==NULL) throw("Error: please input gene annotation file by the flag --gene-list.");
        
        vector<int> gene_anno_chr;
        vector<string> gene_anno_genename;
        vector<int> gene_anno_start;
        vector<int> gene_anno_end;
        read_gene_anno(geneAnnoName,gene_anno_chr, gene_anno_genename,gene_anno_start,gene_anno_end);
        map<string, int> gene_anno_map;
        map<string, int>::iterator iter;
        for(int i=0;i<gene_anno_genename.size();i++) gene_anno_map.insert(pair<string,int>(gene_anno_genename[i], i));
        
        eqtlInfo prbhead;
        read_epifile(&prbhead, string(eqtlFileName)+".epi");
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_prob(&prbhead, problstName);
            if(genelistName != NULL) extract_prob_by_gene(&prbhead, genelistName);
        }
        else if(prbname!=NULL)
        {
            extract_eqtl_single_probe(&prbhead, prbname);
        }
  
        cout<<"\nThere would create "<<prbhead._include.size()<<" plot files..."<<endl;
        if(prbhead._include.size()>1024) cout<<"WARNING: Too many files!!! We strongly recommend using --probe or --extract-probe."<<endl;
        
        for(int plotid=0;plotid<prbhead._include.size();plotid++)
        {
            string plotprbname=prbhead._epi_prbID[prbhead._include[plotid]];
          
            gwasData gdata_;
            eqtlInfo esdata_;
            bInfo bdata;
            gwasData gdata;
            eqtlInfo esdata;
            double threshold= chi_val(1,p_hetero);
            bool heidiFlag=false;
            if(refSNP!=NULL) heidiFlag=true;
            
            
            
            cout<<"\nRetrieving SMR results for plot..."<<endl;
            read_gwas_data( &gdata, gwasFileName);
            read_esifile(&esdata, string(eqtlFileName)+".esi");
            esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, plotprbname.c_str());
            if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
            if(!heidioffFlag)
            {
                read_famfile(&bdata, string(bFileName)+".fam");
                if(indilstName != NULL) keep_indi(&bdata,indilstName);
                if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
                read_bimfile(&bdata, string(bFileName)+".bim");
                if(snplstName != NULL) extract_snp(&bdata, snplstName);
                if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
                allele_check(&bdata, &gdata, &esdata);
                read_bedfile(&bdata, string(bFileName)+".bed");
                if (bdata._mu.empty()) calcu_mu(&bdata);
                if (maf > 0)
                {
                    filter_snp_maf(&bdata, maf);
                    update_geIndx(&bdata, &gdata, &esdata);
                }
                
            }else
            {
                allele_check(&gdata, &esdata);
            }
            
            update_gwas(&gdata);
            cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&esdata, string(eqtlFileName)+".epi");
            extract_prob(&esdata, plotprbname, prbWind);
           
            read_besdfile(&esdata, string(eqtlFileName)+".besd");
            if(esdata._rowid.empty() && esdata._bxz.empty())
            {
                printf("No data included from %s in the analysis.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            
            
            unsigned int probNum = esdata._probNum;
            
            vector<double> bxz;
            vector<double> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            vector<int> snpchrom;
            
            vector<string> allele1;
            vector<string> allele2;
            vector<uint32_t> bpsnp;
            vector<double> freq;
            
            vector<double> byz;
            vector<double> seyz;
            vector<double> pyz;
            VectorXd zsxz;
            
            vector<int> sn_ids;
            
            VectorXd _byz;
            VectorXd _seyz;
            VectorXd _bxz;
            VectorXd _sexz;
            VectorXd _zsxz;
            MatrixXd _X;
            MatrixXd _LD;
            VectorXd ld_v;
            MatrixXd _LD_heidi;
            MatrixXd _X_heidi;
            int cis_itvl_bp=cis_itvl*1000;
            
            //for plot
            string plotdir="";
            string plotnm=outFileName;
            for(long j=strlen(outFileName)-1;j>=0;j--)
                if(outFileName[j]=='/')
                {
                    plotdir=string(outFileName).substr(0,j+1);
                    plotnm=string(outFileName).substr(j+1,strlen(outFileName));
                    break;
                }
            if(plotdir=="") plotdir="./";
            
            plotdir=string(plotdir)+"plot";
            struct stat st = {0};
            if (stat(plotdir.c_str(), &st) == -1) {
#if defined _WIN64 || defined _WIN32
                _mkdir(plotdir.c_str());
#else
                mkdir(plotdir.c_str(), 0755);
#endif
            }
            string plot_path= string(plotdir)+"/"+plotnm+"."+plotprbname+".txt";
            
            //probe info
            vector<string> epi_out_id;
            vector<int> epi_out_chr;
            vector<int> epi_out_bp;
            vector<string> epi_out_gene;
            vector<char> epi_out_orien;
            vector<double> epi_out_pheidi;
            vector<double> epi_out_psmr;
            vector<int> epi_out_start;
            vector<int> epi_out_end;
            
            //eqtl ld info
            vector<int> ldprbid;
            vector<string> ldprb;
            vector<int> ldnperprb;
            vector<string> ldrs;
            vector<double> outld;
            ldnperprb.push_back(0);
            
            
            long idx=find(esdata._epi_prbID.begin(), esdata._epi_prbID.end(), plotprbname)-esdata._epi_prbID.begin();
            if(idx==esdata._epi_prbID.size())
            {
                string logstr="ERROR: Can't find probe "+string(plotprbname)+".\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            int minBP=esdata._epi_bp[idx]-prbWind*1000>0?(esdata._epi_bp[idx]-prbWind*1000):0;
            int maxBP=esdata._epi_bp[idx]+prbWind*1000;
            int refprbchr=esdata._epi_chr[idx];
            
            for(int i=0;i<probNum;i++)
            {
                
                //extract info from eqtl summary and gwas summary
                bxz.clear();
                sexz.clear();
                curId.clear(); // is the idxes of bfile._include not the values of
                eName.clear();
                snpchrom.clear();
                byz.clear();
                seyz.clear();
                pyz.clear();
                allele1.clear();
                allele2.clear();
                bpsnp.clear();
                freq.clear();
                long maxid =-9;
                int probebp=esdata._epi_bp[i];
                int probechr=esdata._epi_chr[i];
                string probenm=esdata._epi_prbID[i];
                int tmpminBP=probebp;
                int tmpmaxBP=probebp;
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<bdata._include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                    {
                        if (fabs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            int snpbp=esdata._esi_bp[j];
                            int snpchr=esdata._esi_chr[j];
                            if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl_bp)
                            {
                                tmpminBP=snpbp<tmpminBP?snpbp:tmpminBP;
                                tmpmaxBP=snpbp>tmpmaxBP?snpbp:tmpmaxBP;
                                bxz.push_back(esdata._bxz[i][j]);
                                sexz.push_back(esdata._sexz[i][j]);
                                byz.push_back(gdata.byz[j]);
                                seyz.push_back(gdata.seyz[j]);
                                pyz.push_back(gdata.pvalue[j]);
                                curId.push_back(j);
                                eName.push_back(esdata._esi_rs[j]);
                                snpchrom.push_back(esdata._esi_chr[j]);
                                if(heidiFlag && esdata._esi_rs[j]==string(refSNP)) maxid=(eName.size()-1);
                                if(!heidioffFlag) //if heidi off , bfile is not necessary to read.
                                {
                                    freq.push_back(bdata._mu[bdata._include[j]]/2);
                                }
                                allele1.push_back(esdata._esi_allele1[j]);
                                allele2.push_back(esdata._esi_allele2[j]);
                                bpsnp.push_back(esdata._esi_bp[j]);
                            }
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata._rowid[beta_start+j];
                        int snpbp=esdata._esi_bp[ge_rowid];
                        int snpchr=esdata._esi_chr[ge_rowid];
                        
                        if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl_bp)
                        {
                            tmpminBP=snpbp<tmpminBP?snpbp:tmpminBP;
                            tmpmaxBP=snpbp>tmpmaxBP?snpbp:tmpmaxBP;
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[ge_rowid]);
                            seyz.push_back(gdata.seyz[ge_rowid]);
                            pyz.push_back(gdata.pvalue[ge_rowid]);
                            curId.push_back(ge_rowid);
                            eName.push_back(esdata._esi_rs[ge_rowid]);
                            snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                            if(heidiFlag && esdata._esi_rs[ge_rowid]==string(refSNP)) maxid=(eName.size()-1);
                            allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            if(!heidioffFlag){
                                freq.push_back(bdata._mu[bdata._include[ge_rowid]]/2);
                            }
                        }
                    }
                }
                if(heidiFlag && maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (bxz.size() == 0) continue;
                
                epi_out_id.push_back(esdata._epi_prbID[i]);
                epi_out_chr.push_back(probechr);
                epi_out_bp.push_back(probebp);
                epi_out_gene.push_back(esdata._epi_gene[i]);
                epi_out_orien.push_back(esdata._epi_orien[i]);
                epi_out_pheidi.push_back(-9);
                epi_out_psmr.push_back(-9);
                iter = gene_anno_map.find((esdata._epi_gene[i]));
                if (iter != gene_anno_map.end())
                                          {
                                              epi_out_start.push_back(gene_anno_start[iter->second]);
                                              epi_out_end.push_back(gene_anno_end[iter->second]);
                                          } else {
                                              epi_out_start.push_back(-9);
                                              epi_out_end.push_back(-9);
                                          }
               
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(!heidiFlag) maxid=max_abs_id(zsxz);
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                if(!heidiFlag && pxz_val>p_smr) continue;
                
                
                if(tmpminBP<minBP) minBP=tmpminBP;
                if(tmpmaxBP>maxBP) maxBP=tmpmaxBP;
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                epi_out_psmr[i]=pxy_val;
                
                make_XMat(&bdata,curId, _X);
                ld_calc_o2m(ld_v,maxid,_X);
                
                // out ld
                ldprbid.push_back(i);
                ldprb.push_back(probenm);
                ldnperprb.push_back((int)curId.size()+ldnperprb[ldnperprb.size()-1]);
                for(int jj=0;jj<curId.size();jj++){
                    ldrs.push_back(esdata._esi_rs[curId[jj]]);
                    outld.push_back(ld_v(jj));
                }
                
                
                sn_ids.clear(); //increase order
                if(fabs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                
                if(sn_ids.size() < m_hetero)   continue;
                
                
                _byz.resize(sn_ids.size());
                _seyz.resize(sn_ids.size());
                _bxz.resize(sn_ids.size());
                _sexz.resize(sn_ids.size());
                _zsxz.resize(sn_ids.size());
                _X_heidi.resize(_X.rows(), sn_ids.size());
                
#pragma omp parallel for
                for(int j=0;j<sn_ids.size();j++)
                {
                    _byz[j]=byz[sn_ids[j]];
                    _seyz[j]=seyz[sn_ids[j]];
                    _bxz[j]=bxz[sn_ids[j]];
                    _sexz[j]=sexz[sn_ids[j]];
                    _zsxz[j]=zsxz[sn_ids[j]];
                    _X_heidi.col(j)=_X.col(sn_ids[j]);
                }
                _X.resize(0,0);
                cor_calc(_LD_heidi, _X_heidi);
                
                
                
                _X_heidi.resize(0,0);
                
                long nsnp = sn_ids.size();
                double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
                epi_out_pheidi[i]=pdev;
            }
            
            if(ldprbid.size()==0)
            {
                string logstr="No SMR results fetched.\n";
                fputs(logstr.c_str(),stdout);
                continue;
            }
            
            cout<<"\nRetrieving GWAS summary information and eQTL summary information for plot..."<<endl;
            read_gwas_data( &gdata_, gwasFileName);
            
            read_epifile(&esdata_, string(eqtlFileName)+".epi");
            extract_prob(&esdata_, plotprbname, prbWind);
            
            
            //get eQTL info in the region [minBP,maxBP]
            read_esifile(&esdata_, string(eqtlFileName)+".esi");
            vector<int> newIcld;
            for(int i=0;i<esdata_._esi_include.size();i++)
            {
                int tmpint=esdata_._esi_include[i];
                if( esdata_._esi_chr[tmpint]==refprbchr &&esdata_._esi_bp[tmpint]>=minBP && esdata_._esi_bp[tmpint]<=maxBP) newIcld.push_back(tmpint);
            }
            esdata_._esi_include.clear();
            esdata_._esi_include=newIcld;
            cout << esdata_._esi_include.size() << " SNPs are extracted from SNP BP: " +atos(minBP)+" bp to SNP BP: " + atos(maxBP) + " bp of eQTL summary dataset." << endl;
            //get RSs in the region [minBP,maxBP] to select gwas summary
            read_bimfile(&bdata, string(bFileName)+".bim");
            vector<string> bsnprs;
            vector<int> bsnpbp;
            for(int i=0;i<bdata._snp_num;i++)
                if(bdata._chr[i]==refprbchr && bdata._bp[i]<=maxBP && bdata._bp[i]>=minBP)
                {
                    bsnprs.push_back(bdata._snp_name[i]);
                    bsnpbp.push_back(bdata._bp[i]);
                }
            
            vector<int> idx1;
            gdata_.snpBp.clear();
            gdata_.snpBp.resize(gdata_.snpName.size());
            gdata_._include.clear();
            match(gdata_.snpName,bsnprs,idx1); //get gwas info
            for(int i=0;i<gdata_.snpName.size();i++)
            {
                if(idx1[i]!=-9)
                {
                    gdata_._include.push_back(i);
                    gdata_.snpBp[i]=bsnpbp[idx1[i]];
                }
                
            }
            cout << gdata_._include.size() << " SNPs are extracted from SNP BP: " +atos(minBP)+" bp to SNP BP: " + atos(maxBP) + " bp of GWAS summary dataset." << endl;
            
            update_gwas(&gdata_);
            rm_unmatched_snp(&gdata_, &esdata_); //allele check at the same time
            update_gwas(&gdata_);
            
            read_besdfile(&esdata_, string(eqtlFileName)+".besd");
            if(esdata_._rowid.empty() && esdata_._bxz.empty())
            {
                printf("No data included from %s in the analysis.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            //SNP info, the union of gwas and eqtl
            vector<string> out_rs;
            vector<int> out_chr;
            vector<int> out_bp;
            vector<string> out_a1;
            vector<string> out_a2;
            map<string,int> snp_name_map;
            long mapsize=0;
            for(int i=0;i<esdata_._esi_include.size();i++)
            {
                snp_name_map.insert(pair<string, int>(esdata_._esi_rs[esdata_._esi_include[i]], mapsize));
                if (mapsize < snp_name_map.size()) {
                    out_rs.push_back(esdata_._esi_rs[esdata_._esi_include[i]]);
                    out_chr.push_back(esdata_._esi_chr[esdata_._esi_include[i]]);
                    out_bp.push_back(esdata_._esi_bp[esdata_._esi_include[i]]);
                    out_a1.push_back(esdata_._esi_allele1[esdata_._esi_include[i]]);
                    out_a2.push_back(esdata_._esi_allele2[esdata_._esi_include[i]]);
                    mapsize=snp_name_map.size();
                }
                
            }
            for(int i=0;i<gdata_._include.size();i++)
            {
                snp_name_map.insert(pair<string, int>(gdata_.snpName[gdata_._include[i]], mapsize));
                if (mapsize < snp_name_map.size()) {
                    out_rs.push_back(gdata_.snpName[gdata_._include[i]]);
                    out_chr.push_back(refprbchr);
                    out_bp.push_back(gdata_.snpBp[gdata_._include[i]]);
                    out_a1.push_back(gdata_.allele_1[gdata_._include[i]]);
                    out_a2.push_back(gdata_.allele_2[gdata_._include[i]]);
                    mapsize=snp_name_map.size();
                }
            }
            vector<int> bprank;
            getRank_norep(out_bp, bprank);
            vector<string> tmptmpstr;
            tmptmpstr.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpstr[bprank[i]]=out_rs[i];
            out_rs.swap(tmptmpstr);
            vector<int> tmptmpint;
            tmptmpint.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpint[bprank[i]]=out_bp[i];
            out_bp.swap(tmptmpint);
            vector<string> tmptmpchar;
            tmptmpchar.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpchar[bprank[i]]=out_a1[i];
            out_a1.swap(tmptmpchar);
            tmptmpchar.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpchar[bprank[i]]=out_a2[i];
            out_a2.swap(tmptmpchar);
            
            //gwas info
            vector<string> gwas_rs;
            vector<float> gwas_be;
            vector<float> gwas_se;
            for(int i=0;i<gdata_.snpNum;i++)
            {
                gwas_rs.push_back(gdata_.snpName[i]);
                gwas_be.push_back(gdata_.byz[i]);
                gwas_se.push_back(gdata_.seyz[i]);
            }
            //eqtl info
            vector<int> out_esi_id;
            vector<int> out_epi_id;
            vector<string> out_epi_name;
            vector<float> out_beta;
            vector<float> out_se;
            vector<double> out_pval;
            if(esdata_._valNum==0)
            {
                for(uint32_t ii=0;ii<ldprbid.size();ii++)
                {
                    int i=ldprbid[ii];
                    for(uint32_t j=0;j<esdata_._snpNum;j++)
                    {
                        double beta=esdata_._bxz[i][j];
                        double se=esdata_._sexz[i][j];
                        if(ABS(se+9)<1e-6) continue;
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=1)
                        {
                            out_esi_id.push_back(j);
                            out_epi_id.push_back(i);
                            out_epi_name.push_back(esdata_._epi_prbID[i]);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }
            else
            {
                if(esdata_._val.size()==0)
                {
                    throw ("Error: No data extracted from the input, please check.\n");
                }
                
                for(uint32_t ii=0;ii<ldprbid.size();ii++)
                {
                    int i=ldprbid[ii];
                    uint64_t proid=esdata_._include[i];
                    uint64_t pos=esdata_._cols[proid<<1];
                    uint64_t pos1=esdata_._cols[(proid<<1)+1];
                    uint64_t num=pos1-pos;
                    for(int j=0;j<num;j++)
                    {
                        double beta=esdata_._val[pos+j];
                        double se=esdata_._val[pos+j+num];
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=1)
                        {
                            out_esi_id.push_back(esdata_._rowid[pos+j]);
                            out_epi_id.push_back(i);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }
            
            vector<double> out_esi_ld;
            vector<string> out_esi_rs;
            vector<int> stend;
            stend.push_back(0);
            int curprid=out_epi_id[0];
            for(int i=0;i<out_esi_id.size();i++)
            {
                if(out_epi_id[i]!=curprid)
                {
                    stend.push_back(i);
                    curprid=out_epi_id[i];
                }
                out_esi_rs.push_back(esdata_._esi_rs[out_esi_id[i]]);
                out_esi_ld.push_back(-9);
            }
            stend.push_back((int)out_esi_id.size());
            
            
            for(uint32_t ii=0;ii<ldprbid.size();ii++)
            {
                vector<string> outrs;
                for(int j=stend[ii];j<stend[ii+1];j++) outrs.push_back(out_esi_rs[j]);
                vector<string> ldsnp;
                for(int j=ldnperprb[ii];j<ldnperprb[ii+1];j++) ldsnp.push_back(ldrs[j]);
                vector<int> idx;
                match(ldsnp,outrs,idx);
                for(int j=0;j<idx.size();j++) out_esi_ld[stend[ii]+idx[j]]=outld[ldnperprb[ii]+j];
            }
            cout<<"Total "<<out_esi_id.size()<<" eQTLs for "<<ldprbid.size()<<" probes are extracted."<<endl;
            free_gwas_data( &gdata_);
            
            
            
            
            FILE* plotfile=NULL;
            plotfile = fopen(plot_path.c_str(), "w");
            if (!(plotfile)) {
                printf("Open error %s\n", plot_path.c_str());
                exit(1);
            }
            
            string outstr="$probe "+atos(epi_out_id.size())+" "+plotprbname+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<epi_out_id.size();i++)
            {
                outstr=epi_out_id[i]+' '+atos(epi_out_chr[i])+' '+atos(epi_out_bp[i])+' '+epi_out_gene[i]+' '+(epi_out_start[i]==-9?"NA":(epi_out_start[i]==23?"X":(epi_out_start[i]==24?"Y":(atos(epi_out_start[i])))))+' '+(epi_out_end[i]==-9?"NA":(epi_out_end[i]==23?"X":(epi_out_end[i]==24?"Y":(atos(epi_out_end[i])))))+' '+epi_out_orien[i]+' '+(epi_out_psmr[i]+9<1e-6?"NA":dtos(epi_out_psmr[i]))+' '+(epi_out_pheidi[i]+9<1e-6?"NA":dtos(epi_out_pheidi[i]))+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$SNP "+atos(out_rs.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<out_rs.size();i++)
            {
                outstr=out_rs[i]+' '+atos(out_chr[i])+' '+atos(out_bp[i])+' '+out_a1[i]+' '+out_a2[i]+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$GWAS "+atos(gwas_rs.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<gwas_rs.size();i++)
            {
                outstr=gwas_rs[i]+' '+atos(gwas_be[i])+' '+atos(gwas_se[i])+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$eQTL "+atos(ldprbid.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<stend.size()-1;i++)
            {
                outstr=esdata_._epi_prbID[out_epi_id[stend[i]]] +" "+atos(stend[i+1]-stend[i])+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int j=stend[i];j<stend[i+1];j++)
                {
                    outstr=out_esi_rs[j]+' '+atos(out_beta[j])+' '+atos(out_se[j])+' '+(out_esi_ld[j]+9<1e-6?"NA":atos(out_esi_ld[j]))+'\n';
                    if(fputs_checked(outstr.c_str(),plotfile))
                    {
                        printf("ERROR: in writing file %s .\n", plot_path.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
            }
            fclose(plotfile);
            free_gwas_data( &gdata);
            cout<<"Information for plot has been saved in "<<plot_path<<"."<<endl;
            
        }
    }
    
    
    void read_geneAnno(string gAnno_file, vector<string> &gene_name, vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2) {
        ifstream in_gAnno(gAnno_file.c_str());
        if (!in_gAnno) throw ("Error: can not open the file [" + gAnno_file + "] to read.");
        cout << "Reading physical positions of the genes from [" + gAnno_file + "]." << endl;
        string str_buf;
        vector<string> vs_buf;
        while (getline(in_gAnno, str_buf)) {
            if (StrFunc::split_string(str_buf, vs_buf) != 4) throw ("Error: in line \"" + str_buf + "\".");
            gene_chr.push_back(atoi(vs_buf[0].c_str()));
            gene_bp1.push_back(atoi(vs_buf[1].c_str()));
            gene_bp2.push_back(atoi(vs_buf[2].c_str()));
            gene_name.push_back(vs_buf[3]);
        }
        in_gAnno.close();
        cout << "Physical positions of " << gene_name.size() << " genes have been include." << endl;
    }
    
    void sbat_read_snpset(bInfo* bdata, char* snpset_file, vector<string> &set_name,  vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2, vector< vector<string> > &snpset)
    {
        ifstream in_snpset(snpset_file);
        if (!in_snpset) throw ("Error: can not open the file [" + string(snpset_file) + "] to read.");
        cout << "\nReading SNP sets from [" + string(snpset_file) + "]." << endl;
        string str_buf;
        vector<string> vs_buf, snpset_buf, snp_name;
        map<string, int>::iterator iter;
        int i = 0;
        map<int, int> chr_map;
        bool warning=false;
        while (in_snpset>>str_buf) {
            if(str_buf!="END" && str_buf!="end") vs_buf.push_back(str_buf);
            else{
                if(vs_buf.empty()) continue;
                int bp1=INT_MAX;
                int bp2=0;
                int chr=0;
                chr_map.clear();
                for(i = 1; i < vs_buf.size(); i++){
                    iter = bdata->_snp_name_map.find(vs_buf[i]);
                    if (iter != bdata->_snp_name_map.end()){
                        int snpbp=bdata->_bp[iter->second];
                        bp1=bp1<snpbp?bp1:snpbp;
                        bp2=bp2>snpbp?bp2:snpbp;
                        chr=snpbp=bdata->_chr[iter->second];
                        chr_map.insert(pair<int,int>(chr,i));
                        if(chr_map.size()>1)
                        {
                            printf("ERROR: SNPs from different chromosomes found in SNP set: %s.\n", vs_buf[0].c_str());
                            exit(EXIT_FAILURE);
                        }
                        snpset_buf.push_back(vs_buf[i]);
                        snp_name.push_back(vs_buf[i]);
                    }
                }
                
                if(snpset_buf.size()>0)
                {
                    if(bp2-bp1>1000000 && !warning)
                    {
                        printf("WARNING: The size of at least one SNP set is over 1Mb.\n");
                        warning=true;
                    }
                    
                    set_name.push_back(vs_buf[0]);
                    gene_chr.push_back(chr);
                    gene_bp1.push_back(bp1);
                    gene_bp2.push_back(bp2);
                    snpset.push_back(snpset_buf);
                }
                vs_buf.clear();
                snpset_buf.clear();
            }
        }
        in_snpset.close();
        snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());        
        cout << snp_name.size() << " SNPs in " << snpset.size() << " sets have been matched and included." << endl;
    }

   
    int maxabsid(vector<double> &zsxz, vector<int> &ids)
    {
        
        if(ids.size()==0)
        {
            printf("ERROR: no id values found!\n");
            exit(EXIT_FAILURE);
        }
        int id=ids[0];
        double tmpVal, cmpVal=fabs(zsxz[ids[0]]);
        for( int i=1;i<ids.size();i++)
        {
            tmpVal=fabs(zsxz[ids[i]]);
            if( cmpVal-tmpVal < 1e-6)
            {
                cmpVal=tmpVal;
                id=ids[i];
            }
        }
        return(id);
    }
    

    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1,vector<double> &zxz4smr) {
        
        vector<int> remain,select;
        for(int i=0;i<zxz4smr.size();i++) remain.push_back(i);
        while(!remain.empty())
        {
            
            int maxid=maxabsid(zxz4smr,remain);
            select.push_back(maxid);
            vector<int> new_remain;
            for(int i=0;i<remain.size();i++)
            {
                if(remain[i]!= maxid )
                {
                    if( fabs(R(maxid,remain[i])) > R_cutoff){
                        rm_ID1.push_back(remain[i]);
                    }else {
                        new_remain.push_back(remain[i]);
                    }
                }
            }
            remain.swap(new_remain);
        }
         stable_sort(rm_ID1.begin(), rm_ID1.end());
    }
    
    void sbat_calcu_lambda(MatrixXd &X, VectorXd &eigenval,VectorXd &eigenvalxy, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx,vector<double> &zxz4smr, vector<double> &zyz4smr)
    {
        int m = snp_count;
       
        vector<int> rm_ID1;
        double R_cutoff = sbat_ld_cutoff;
        int qi = 0; //alternate index
        
        MatrixXd C;
        cor_calc(C, X);
        if (sbat_ld_cutoff < 1) rm_cor_sbat(C, R_cutoff, m, rm_ID1,zxz4smr);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = (int)sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            MatrixXd D(sub_indx.size(),sub_indx.size());
            for (int i = 0 ; i < sub_indx.size() ; i++) {
                for (int j = 0 ; j < sub_indx.size() ; j++) {
                    D(i,j) = C(sub_indx[i],sub_indx[j]);
                }
            }
            C = D;
        }
        SelfAdjointEigenSolver<MatrixXd> saes(C);
        eigenval = saes.eigenvalues().cast<double>();
         /* save out ld*/
        /*
        cout<<C.rows()<<"*"<<C.cols()<<endl;
        FILE* bsefile=NULL;
        string ldfilename=atos(xh)+".ld.txt";
        bsefile = fopen(ldfilename.c_str(), "w");
        if (!(bsefile)) {
            printf("Open error %s\n", ldfilename.c_str());
            exit(1);
        }
        string outstr="";
        for(int j=0;j<C.rows();j++)
        {
            for(int k=0;k<C.cols()-1;k++) outstr+=atos(C(j,k))+'\t';
            outstr+=atos(C(j,C.rows()-1))+'\n';
        }
        
        if(fputs_checked(outstr.c_str(),bsefile))
        {
            printf("ERROR: in writing file %s .\n", ldfilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(bsefile);
       
        cout<<eigenval<<endl;
        */
         /* end of saving*/
        
        VectorXd zs_xz(sub_indx.size()),zs_yz(sub_indx.size());
        for (int i = 0 ; i < sub_indx.size() ; i++){
            zs_xz(i)=zxz4smr[sub_indx[i]];
            zs_yz(i)=zyz4smr[sub_indx[i]];
        }
        MatrixXd xtx=zs_xz*zs_xz.transpose();
        MatrixXd yty=zs_yz*zs_yz.transpose();
        MatrixXd zz=yty.array()/xtx.array();
        MatrixXd numerator=C.array()*(xtx+yty).array()-zz.array();
        VectorXd zsq=zs_xz.array()*zs_xz.array()+zs_yz.array()*zs_yz.array();
        MatrixXd denominator=(zsq*zsq.transpose()).array().sqrt();
        C=numerator.array()/denominator.array();
        SelfAdjointEigenSolver<MatrixXd> saesxy(C);
        eigenvalxy = saesxy.eigenvalues().cast<double>();
    }
    
    double heidi_test(bInfo* bdata,vector<double> &slct_zsxz,vector<uint32_t> &slctId, long slct_maxid,double ld_top, double threshold, int m_hetero, vector<double> &slct_byz, vector<double> &slct_seyz,vector<double> &slct_bxz,vector<double> &slct_sexz, long &nsnp )
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        VectorXd tmp_zsxz(slct_zsxz.size());
        for(int j=0;j<slct_zsxz.size();j++) tmp_zsxz(j)=slct_zsxz[j];
        
        make_XMat(bdata,slctId, _X);
        ld_calc_o2m(ld_v,slct_maxid,_X);
        if(fabs(ld_top-1)<1e-6) get_square_idxes(sn_ids,tmp_zsxz,threshold);
        else get_square_ldpruning_idxes(sn_ids,tmp_zsxz,threshold,ld_v, slct_maxid,ld_top);
        if(sn_ids.size() < m_hetero) return -9;
        
        VectorXd _byz,_seyz, _bxz,_sexz,_zsxz;
        MatrixXd _X_heidi, _LD_heidi;
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
        _X_heidi.resize(_X.rows(), sn_ids.size());
        
#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            _byz[j]=slct_byz[sn_ids[j]];
            _seyz[j]=slct_seyz[sn_ids[j]];
            _bxz[j]=slct_bxz[sn_ids[j]];
            _sexz[j]=slct_sexz[sn_ids[j]];
            _zsxz[j]=slct_zsxz[sn_ids[j]];
            _X_heidi.col(j)=_X.col(sn_ids[j]);
        }
        _X.resize(0,0);
        cor_calc(_LD_heidi, _X_heidi);
        
        _X_heidi.resize(0,0);
        
        nsnp = sn_ids.size();
        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
        
        return pdev;
    }

    int smr_setbased_test(bInfo* bdata, vector<uint32_t> &slctId, vector<double> &slct_bxz,vector<double> &slct_sexz,vector<double> &slct_byz,vector<double> &slct_seyz, double p_smr, double ld_top,double &set_pval_smr, double &set_pval_gwas,double &set_pval_eqtl, vector<string> &snp4msmr)
    {
        /* step4: Filter out the SNPs with p-smr threshold and ld-pruning */
        printf("Conducting multi-SNP SMR test...\n");
        vector<uint32_t> Id4smr;
        vector<double> bxz4smr;
        vector<double> sexz4smr;
        vector<double> byz4smr;
        vector<double> seyz4smr;
        vector<double> zxz4smr;
        vector<double> zyz4smr;
        double z_smr=sqrt(qchisq(p_smr,1));
        for(int j=0;j<slctId.size();j++)
        {
            double ztmp=fabs(slct_bxz[j]/slct_sexz[j]);
            if(ztmp>=z_smr)
            {
                Id4smr.push_back(slctId[j]);
                bxz4smr.push_back(slct_bxz[j]);
                sexz4smr.push_back(slct_sexz[j]);
                zxz4smr.push_back(slct_bxz[j]/slct_sexz[j]);
                byz4smr.push_back(slct_byz[j]);
                seyz4smr.push_back(slct_seyz[j]);
                zyz4smr.push_back(slct_byz[j]/slct_seyz[j]);
            }
        }
        int snp_count=(int)Id4smr.size();
        printf("%ld SNPs passed the p-value threshold %6.2e and %ld SNPs are excluded.\n",Id4smr.size(), p_smr,slctId.size()-Id4smr.size());
        if(snp_count==0) return -9;
        /* step5: multiple-SNP SMR test */
        
        VectorXd eigenval;
        VectorXd eigenvalxy;
        vector<int> sub_indx;
        MatrixXd _X;
        make_XMat(bdata,Id4smr, _X); //_X: one row one individual, one column one SNP
        double sbat_ld_cutoff=sqrt(ld_top);
        sbat_calcu_lambda(_X, eigenval, eigenvalxy, snp_count,  sbat_ld_cutoff, sub_indx, zxz4smr, zyz4smr); //the index of slectId, snp_count can chage here
        printf("%ld SNPs passed LD-square threshold %6.2f and %ld SNPs are excluded.\n",sub_indx.size(), ld_top,Id4smr.size()-sub_indx.size());
        vector<double> zsxysq_slct(sub_indx.size());
        double chisq_zy=0;
        double chisq_zx=0;
        for(int j=0;j<sub_indx.size();j++)
        {
            double z1=byz4smr[sub_indx[j]]/seyz4smr[sub_indx[j]];
            double z2=bxz4smr[sub_indx[j]]/sexz4smr[sub_indx[j]];
            chisq_zy += z1*z1;
            chisq_zx += z2*z2;
            zsxysq_slct[j]=z1*z1*z2*z2/(z1*z1+z2*z2);
        }
        double chisq_o = 0;
        for (int j = 0; j < sub_indx.size(); j++)  chisq_o += zsxysq_slct[j];
        for (int j = 0; j < sub_indx.size(); j++) snp4msmr.push_back(bdata->_snp_name[bdata->_include[Id4smr[sub_indx[j]]]]);
        
        /* save out byz,seyz,bxz,sexz*/
        /*
        FILE* bsefile=NULL;
        string bsefilename=atos(xh)+".bse.txt";
        bsefile = fopen(bsefilename.c_str(), "w");
        if (!(bsefile)) {
            printf("Open error %s\n", bsefilename.c_str());
            exit(1);
        }
         string outstr="";
        for(int j=0;j<sub_indx.size();j++)
        {
             outstr+=atos(byz4smr[sub_indx[j]])+'\t'+atos(seyz4smr[sub_indx[j]])+'\t'+atos(bxz4smr[sub_indx[j]])+'\t'+atos(sexz4smr[sub_indx[j]])+'\n';
        }
        
        if(fputs_checked(outstr.c_str(),bsefile))
        {
            printf("ERROR: in writing file %s .\n", bsefilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(bsefile);
        xh++;
         */
        /* end of saving*/
        printf("%ld SNPs are included in the multi-SNP SMR test.\n",sub_indx.size());
        if(sub_indx.size() == 1)
        {
            set_pval_smr = pchisq(chisq_o, 1.0);
            set_pval_gwas= pchisq(chisq_zy, 1.0);
            set_pval_eqtl= pchisq(chisq_zx, 1.0);
            
        }
        else {
            set_pval_smr = chisq_o<1e-8?1:pchisqsum(chisq_o, eigenval);
            set_pval_gwas= chisq_zy<1e-8?1:pchisqsum(chisq_zy, eigenval);
            set_pval_eqtl= chisq_zx<1e-8?1:pchisqsum(chisq_zx, eigenval);
        }
        return snp_count;
    }
    void ssmr_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, double ld_min,int opt_hetero, int expanWind,bool sampleoverlap, double pmecs, int minCor, double ld_top_multi)
    {
        
        vector<string> set_name;
        vector< vector<string> > snpset;
        vector<int> gene_chr,gene_bp1,gene_bp2;
        double thresh_heidi= chi_val(1,p_hetero),theta=0;
        
        cout<<endl<<"Performing multi-SNP based SMR analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        cis_itvl=cis_itvl*1000;
        if(expanWind!=-9) {
            expanWind=expanWind*1000;
            if(expanWind>cis_itvl) {
                cis_itvl=expanWind;
                printf("Cis-region window size is extended to %dKb.\n",cis_itvl);
            }
        }
        
        FILE *smr=NULL, *glst=NULL, *setlst=NULL;
        long write_count=0;
        string smrfile="", setlstfile="", genelstfile="", outstr="";
        if(outFileName!=NULL)
        {
            setlstfile = string(outFileName)+".snps4msmr.list";
            setlst = fopen(setlstfile.c_str(), "w");
            if (!(setlst)) {
                printf("Open error %s\n", setlstfile.c_str());
                exit(1);
            }
            
            genelstfile = string(outFileName)+".prbregion4msmr.list";
            glst = fopen(genelstfile.c_str(), "w");
            if (!(glst)) {
                printf("Open error %s\n", genelstfile.c_str());
                exit(1);
            }

            smrfile = string(outFileName)+".msmr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("Open error %s\n", smrfile.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_SMR_multi\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
     
        map<string, int>::iterator iter;
        SMRWK smrwk;
            for(int ii=0;ii<esdata->_include.size();ii++)
            {
                
                progr1=1.0*ii/esdata->_include.size();
                if(progr1-progr0-0.05>1e-6 || ii+1==esdata->_include.size())
                {
                    if(ii+1==esdata->_include.size()) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
                int i=esdata->_include[ii];
                int probebp=esdata->_epi_bp[i];
                int probechr=esdata->_epi_chr[i];
                string probename=esdata->_epi_prbID[i];
                string probegene=esdata->_epi_gene[i];
                char probeorien=esdata->_epi_orien[i];
                init_smr_wk(&smrwk);
                smrwk.cur_prbidx=i;
                /* step1: get cis-eQTLs */
                printf("\nInitiating the workspace of probe %s for multi-SNP SMR analysis....\n",probename.c_str());
                long maxid =fill_smr_wk(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, heidioffFlag);
                if(refSNP!=NULL && maxid==-9) {
                    printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                    continue;
                } //ref heidi SNP is not in selected SNPs
                if (smrwk.bxz.size() == 0) {
                    printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                    continue;
                }
                printf("%ld SNPs are included from the cis-region of the probe %s.\n",smrwk.bxz.size(),probename.c_str());
                //now if you sepcify reference SNP, maxid point to this SNP, otherwise maxid is -9
                /* step2: get top-SNP */
                printf("Checking the top-SNP in the region....\n");
                Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
                Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
                if(sampleoverlap)
                {
                    printf("Estimating the correlation ...\n");
                    double z2mecs=qchisq(pmecs,1);
                    double zmecs=sqrt(z2mecs);
                    vector<double> zxz, zyz;
                    for(int k=0;k<smrwk.bxz.size();k++)
                    {
                        double z1=smrwk.bxz[k]/smrwk.sexz[k];
                        double z2=smrwk.byz[k]/smrwk.seyz[k];
                        if(fabs(z1)<zmecs && fabs(z2)<zmecs )
                        {
                            zxz.push_back(z1);
                            zyz.push_back(z2);
                        }
                        
                    }
                    if(zxz.size()< minCor){
                        printf("WARNING: less than %d common SNPs obtained from the cis-region of probe %s at a p-value threshold %5.2e.\n ", minCor,probename.c_str(), pmecs);
                        printf("probe %s is skipped.\n ", probename.c_str());
                        continue;
                    }
                    else
                    {
                        theta=cor(zxz,zyz);
                        printf("The estimated correlation is %f.\n",theta);
                    }
                }
                VectorXd zsxz;
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(refSNP==NULL) maxid=max_abs_id(zsxz); // now maxid point to the sig eQTL SNP or ref SNP in the new datastruct(not the raw).
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                //double computing, consistency should be checked
                for(int tid=0;tid<zsxz.size();tid++) {
                    if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                    {
                        printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                        exit(EXIT_FAILURE);
                    }
                }
                string topsnpname=smrwk.rs[maxid];
                printf("The top SNP of probe %s is %s with p-value %e.\n", probename.c_str(), topsnpname.c_str(),pxz_val);
                if(refSNP==NULL && pxz_val>p_smr){
                    printf("WARNING: no SNP passed the p-value threshold %e for Multiple-SNP SMR analysis for probe %s.\n", p_smr, probename.c_str());
                    continue;
                } else {
                    printf("Conducting multi-SNP SMR and HEIDI test for probe %s...\n", probename.c_str());
                }
                //cout<<maxid<<":"<<topsnpname<<":"<<esdata._esi_rs[smrwk.curId[maxid]]<<":"<<bdata._snp_name[bdata._include[smrwk.curId[maxid]]]<<":"<<gdata.snpName[smrwk.curId[maxid]]<<endl;
                /* step3: extract SNPs around the --set-wind around sig (or ref) SNP */
                if(expanWind!=-9) printf("Extracting SNPs in a specified window around top-SNP/ref-SNP....\n");
                else printf("Extracting SNPs in the cis-region....\n");
                vector<uint32_t> slctId;
                vector<int> slct_bpsnp,slct_snpchr;
                vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_zsxz,slct_zxz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
                vector<string> slct_snpName,slct_a1, slct_a2;
                long slct_maxid=-9;
                if(expanWind!=-9){
                    for(int j=0;j<zsxz.size();j++)
                    {
                        int maxid_bp=smrwk.bpsnp[maxid];
                        int tmplower=maxid_bp-expanWind>0?maxid_bp-expanWind:0;
                        int tmpupper=maxid_bp+expanWind;
                        if(smrwk.bpsnp[j]>=tmplower && smrwk.bpsnp[j]<=tmpupper)
                        {
                            slctId.push_back(smrwk.curId[j]); // for get X
                            slct_bxz.push_back(smrwk.bxz[j]);
                            slct_sexz.push_back(smrwk.sexz[j]);
                            slct_byz.push_back(smrwk.byz[j]);
                            slct_seyz.push_back(smrwk.seyz[j]);
                            slct_zsxz.push_back(zsxz(j));
                            slct_snpName.push_back(smrwk.rs[j]);
                            if(j==maxid) slct_maxid=slctId.size()-1;
                            slct_zxz.push_back(smrwk.zxz[j]);
                            slct_pyz.push_back(smrwk.pyz[j]);
                            slct_bpsnp.push_back(smrwk.bpsnp[j]);
                            slct_snpchr.push_back(smrwk.snpchrom[j]);
                            slct_a1.push_back(smrwk.allele1[j]);
                            slct_a2.push_back(smrwk.allele2[j]);
                            slct_freq.push_back(smrwk.freq[j]);
                        }
                    }
                }else {
                    slctId.swap(smrwk.curId);
                    slct_bxz.swap(smrwk.bxz);
                    slct_sexz.swap(smrwk.sexz);
                    slct_byz.swap(smrwk.byz);
                    slct_seyz.swap(smrwk.seyz);
                    slct_snpName.swap(smrwk.rs);
                    slct_maxid=maxid;
                    for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
                    slct_zxz.swap(smrwk.zxz);
                    slct_pyz.swap(smrwk.pyz);
                    slct_bpsnp.swap(smrwk.bpsnp);
                    slct_snpchr.swap(smrwk.snpchrom);
                    slct_a1.swap(smrwk.allele1);
                    slct_a2.swap(smrwk.allele2);
                    slct_freq.swap(smrwk.freq);
                }
                
                int out_raw_id = slctId[slct_maxid];
                double bxz_max = slct_bxz[slct_maxid];
                double sexz_max = slct_sexz[slct_maxid];
                double byz_max = slct_byz[slct_maxid];
                double seyz_max = slct_seyz[slct_maxid];
                double bxy_max = byz_max / bxz_max;
                double sexy_max = sqrt((seyz_max * seyz_max * bxz_max * bxz_max + sexz_max * sexz_max * byz_max * byz_max) / (bxz_max * bxz_max * bxz_max * bxz_max));
                double chisqxy = bxy_max*bxy_max / (sexy_max*sexy_max);
                double zsxz_max = bxz_max / sexz_max;
                double zsyz_max = byz_max / seyz_max;
                double pxz_max = pchisq(zsxz_max * zsxz_max, 1);
                double pyz_max = pchisq(zsyz_max * zsyz_max, 1);
                double pxy_max = pchisq(chisqxy, 1);
                // cout<<slct_maxid<<":"<<slct_snpName[slct_maxid]<<endl;
                printf("%ld SNPs are included into SMR and HEIDI test.\n",slctId.size());
                
                double set_pval_smr=-9;
                double set_pval_gwas=-9;
                double set_pval_eqtl=-9;
                vector<string> snp4msmr;
                int snp_count=smr_setbased_test(bdata, slctId, slct_bxz,slct_sexz,slct_byz,slct_seyz, p_smr, ld_top_multi,set_pval_smr, set_pval_gwas,set_pval_eqtl,snp4msmr);
                if(snp_count==-9) continue;
                
                if(outFileName!=NULL) {
                    /* output snp set list*/
                    string setstr=probename+'\n';
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    for(int j=0;j<snp4msmr.size();j++)
                    {
                        setstr=snp4msmr[j]+'\n';
                        if(fputs_checked(setstr.c_str(),setlst))
                        {
                            printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                    setstr="end\n";
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    /* end of output */
                    
                    /* output gene list */
                    int Lbp= *min_element(slct_bpsnp.begin(),slct_bpsnp.end());
                    int Rbp= *max_element(slct_bpsnp.begin(),slct_bpsnp.end());
                    string gstr=atos(probechr) + "\t" + atos(Lbp) + "\t" + atos(Rbp)+"\t"+ probename + "\n";
                    if(fputs_checked(gstr.c_str(),glst))
                    {
                        printf("ERROR: in writing file %s .\n", genelstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    /* end of output */
                    

                }
                /* step6: HEIDI test */
                long nsnp=-9;
                double pdev=-9;
                smrwk.curId.swap(slctId);
                smrwk.bxz.swap(slct_bxz);
                smrwk.sexz.swap(slct_sexz);
                smrwk.byz.swap(slct_byz);
                smrwk.seyz.swap(slct_seyz);
                smrwk.rs.swap(slct_snpName);
                smrwk.zxz.swap(slct_zsxz);
                smrwk.allele1.swap(slct_a1);
                smrwk.allele2.swap(slct_a2);
                smrwk.pyz.swap(slct_pyz);
                smrwk.bpsnp.swap(slct_bpsnp);
                smrwk.snpchrom.swap(slct_snpchr);
                smrwk.freq.swap(slct_freq);
                
                if(!heidioffFlag) {
                    printf("Conducting HEIDI test...\n");
                    pdev= heidi_test_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,ld_min,opt_hetero,sampleoverlap, theta);
                    printf("HEIDI test complete.\n");
                }
                if(outFileName!=NULL)
                {
                    outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + atos(esdata->_esi_chr[out_raw_id]) + '\t' + atos(esdata->_esi_bp[out_raw_id]) + '\t' + esdata->_esi_allele1[out_raw_id] + '\t' + esdata->_esi_allele2[out_raw_id] + '\t' + atos(bdata->_mu[bdata->_include[out_raw_id]] / 2) + '\t';
                    outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t';
                    outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t';
                    outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    SMRRLT currlt;
                    currlt.ProbeID=probename;
                    currlt.ProbeChr=probechr;
                    currlt.Gene=probegene;
                    currlt.Probe_bp=probebp;
                    currlt.Orien=probeorien;
                    currlt.SNP=topsnpname;
                    currlt.SNP_Chr=esdata->_esi_chr[out_raw_id];
                    currlt.SNP_bp=esdata->_esi_bp[out_raw_id];
                    currlt.A1=esdata->_esi_allele1[out_raw_id];
                    currlt.A2=esdata->_esi_allele2[out_raw_id];
                    currlt.Freq=bdata->_mu[bdata->_include[out_raw_id]] / 2;
                    currlt.b_GWAS=byz_max;
                    currlt.se_GWAS=seyz_max;
                    currlt.p_GWAS=pyz_max;
                    currlt.b_eQTL=bxz_max;
                    currlt.se_eQTL=sexz_max;
                    currlt.p_eQTL=pxz_val;
                    currlt.b_SMR=bxy_max;
                    currlt.se_SMR=sexy_max;
                    currlt.p_SMR=pxy_max;
                    currlt.p_SSMR=set_pval_smr;
                    currlt.p_HET=pdev;
                    currlt.nsnp=(int)nsnp;
                    smrrlts.push_back(currlt);
                }
            }
        if(outFileName!=NULL)
        {
            cout<<"\nMultiple-SNP SMR and HEIDI analyses completed.\nSMR and heterogeneity analysis results of "<<write_count<<" sets have been saved in the file [" + smrfile + "]."<<endl;
            cout<<"SNP sets included in multi-SNP SMR have been saved in the file [" + setlstfile + "]."<<endl;
            fclose(smr);
            fclose(setlst);
            fclose(glst);
        }
        free_gwas_data(gdata);
        
    }
    void smr_multipleSNP(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind, double ld_min,double threshpsmrest, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi,double afthresh,double percenthresh)
    {
        double theta=0;
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        double threshold= chi_val(1,p_hetero);
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &gdata, &esdata);
        
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &gdata, &esdata);
        }
        if(forcefrqck)
        {
            double prop= freq_check(&bdata, &gdata, &esdata,afthresh,percenthresh);
            if(prop > percenthresh)
            {
                printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                exit(EXIT_FAILURE);
            }
        }
        update_gwas(&gdata);
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        /**/
        //vector<SMRRLT> smrrlts;
        //ssmr_heidi_func(smrrlts,  outFileName, &bdata,&gdata,&esdata,  cis_itvl,  heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest, ld_min,opt_hetero,expanWind,sampleoverlap,pmecs, minCor,ld_top_multi);
        /**/
        
        vector<string> set_name;
        vector< vector<string> > snpset;
        vector<int> gene_chr,gene_bp1,gene_bp2;
        if(setlstName!=NULL) sbat_read_snpset(&bdata,setlstName,set_name,gene_chr, gene_bp1,gene_bp2, snpset );
        else if(geneAnnoFileName!=NULL) read_geneAnno(geneAnnoFileName, set_name, gene_chr, gene_bp1, gene_bp2);
        
      
        unsigned int probNum = esdata._probNum;
        
        
        cout<<endl<<"Performing multi-SNP based SMR analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);

        cis_itvl=cis_itvl*1000;
        if(expanWind!=-9) {
            expanWind=expanWind*1000;
            if(expanWind>cis_itvl) {
                cis_itvl=expanWind;
                printf("Cis-region window size is extended to %dKb.\n",cis_itvl);
            }
        }
        
        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string genelstfile = string(outFileName)+".prbregion4msmr.list";
        FILE* glst=NULL;
        glst = fopen(genelstfile.c_str(), "w");
        if (!(glst)) {
            printf("Open error %s\n", genelstfile.c_str());
            exit(1);
        }

        string smrfile = string(outFileName)+".msmr";
        FILE* smr=NULL;
        smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("Open error %s\n", smrfile.c_str());
            exit(1);
        }
        
        //string outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tp_GWAS_multi\tb_eQTL\tse_eQTL\tp_eQTL\tp_eQTL_multi\tb_SMR\tse_SMR\tp_SMR\tnsnp_msmr\tp_SMR_multi\tp_HET\tnsnp_heidi\n";
        
         string outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_SMR_multi\tp_HEIDI\tnsnp_HEIDI\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }
        long write_count=0;
        map<string, int>::iterator iter;
        SMRWK smrwk;
        if(set_name.size()>0)
        {
            //gene list based or set list based
            for(int ii=0;ii<set_name.size();ii++)
            {
                progr1=1.0*ii/set_name.size();
                if(progr1-progr0-0.05>1e-6 || ii+1==set_name.size())
                {
                    if(ii+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
                
                int probebp=-9;
                int probechr=gene_chr[ii];
                string probename=set_name[ii];
                string probegene="";
                
                printf("\nInitiating the workspace of probe %s for multi-SNP SMR analysis....\n",set_name[ii].c_str());
                init_smr_wk(&smrwk);
                iter=esdata._probe_name_map.find(set_name[ii]);
                if(iter==esdata._probe_name_map.end())
                {
                    printf("%s is not found in BESD file.\n",set_name[ii].c_str());
                    continue;
                } else {
                    smrwk.cur_prbidx=iter->second;
                    if(gene_chr[ii] != esdata._epi_chr[iter->second])
                    {
                        printf("Error: Inconsistency of probe chromosome of probe %s in the BESD file and the set file.\n",set_name[ii].c_str());
                        printf("Currently this software only supports multi-SNP SMR analysis for the cis-region.\n");
                        exit(EXIT_FAILURE);
                    } else {
                        smrwk.cur_chr=gene_chr[ii];
                        probebp = esdata._epi_bp[iter->second];
                        probegene = esdata._epi_gene[iter->second];
                    }
                    
                }
                
                
                int lowerbp=gene_bp1[ii];
                int upperbp=gene_bp2[ii];
                long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refSNP,lowerbp, upperbp, heidioffFlag);
                if(refSNP!=NULL && maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (smrwk.bxz.size() == 0) continue;
                
                vector<uint32_t> slctId;
                vector<int> slct_bpsnp,slct_snpchr;
                vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_zsxz,slct_zxz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
                vector<string> slct_snpName,slct_a1, slct_a2;
                long slct_maxid;
                if(snpset.size()==0)
                {
                    // gene list
                    slctId.swap(smrwk.curId);
                    slct_bxz.swap(smrwk.bxz);
                    slct_sexz.swap(smrwk.sexz);
                    slct_byz.swap(smrwk.byz);
                    slct_seyz.swap(smrwk.seyz);
                    slct_snpName.swap(smrwk.rs);
                    slct_bpsnp.swap(smrwk.bpsnp);
                    slct_snpchr.swap(smrwk.snpchrom);
                    slct_a1.swap(smrwk.allele1);
                    slct_a2.swap(smrwk.allele2);
                    slct_freq.swap(smrwk.freq);
                    slct_zxz.swap(smrwk.zxz);
                    slct_pyz.swap(smrwk.pyz);
                    
                }else {
                    // set list
                    vector<int> matchidx;
                    match_only(snpset[ii], smrwk.rs, matchidx);
                    if(refSNP!=NULL && find(matchidx.begin(),matchidx.end(),maxid)==matchidx.end()) continue;
                    
                    for(int j=0;j<matchidx.size();j++)
                    {
                        slctId.push_back(smrwk.curId[matchidx[j]]);
                        slct_bxz.push_back(smrwk.bxz[matchidx[j]]);
                        slct_sexz.push_back(smrwk.sexz[matchidx[j]]);
                        slct_byz.push_back(smrwk.byz[matchidx[j]]);
                        slct_seyz.push_back(smrwk.seyz[matchidx[j]]);
                        slct_snpName.push_back(smrwk.rs[matchidx[j]]);
                        slct_pyz.push_back(smrwk.pyz[matchidx[j]]);
                        slct_zxz.push_back(smrwk.zxz[matchidx[j]]);
                        slct_bpsnp.push_back(smrwk.bpsnp[matchidx[j]]);
                        slct_snpchr.push_back(smrwk.snpchrom[matchidx[j]]);
                        slct_a1.push_back(smrwk.allele1[matchidx[j]]);
                        slct_a2.push_back(smrwk.allele2[matchidx[j]]);
                        slct_freq.push_back(smrwk.freq[matchidx[j]]);
                    }
                }
                
                Map<VectorXd> ei_bxz(&slct_bxz[0],slct_bxz.size());
                Map<VectorXd> ei_sexz(&slct_sexz[0],slct_sexz.size());
                VectorXd zsxz;
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(refSNP==NULL) maxid=max_abs_id(zsxz);
                string topsnpname=slct_snpName[maxid];
                slct_maxid=maxid;
                for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
                
                int out_raw_id = slctId[slct_maxid];
                double bxz_max = slct_bxz[slct_maxid];
                double sexz_max = slct_sexz[slct_maxid];
                double byz_max = slct_byz[slct_maxid];
                double seyz_max = slct_seyz[slct_maxid];
                double bxy_max = byz_max / bxz_max;
                double sexy_max = sqrt((seyz_max * seyz_max * bxz_max * bxz_max + sexz_max * sexz_max * byz_max * byz_max) / (bxz_max * bxz_max * bxz_max * bxz_max));
                double chisqxy = bxy_max*bxy_max / (sexy_max*sexy_max);
                double zsxz_max = bxz_max / sexz_max;
                double zsyz_max = byz_max / seyz_max;
                double pxz_max = pchisq(zsxz_max * zsxz_max, 1);
                double pyz_max = pchisq(zsyz_max * zsyz_max, 1);
                double pxy_max = pchisq(chisqxy, 1);
                
                double set_pval_smr=-9;
                double set_pval_gwas=-9;
                double set_pval_eqtl=-9;
                vector<string> snp4msmr;
                int snp_count=smr_setbased_test(&bdata, slctId, slct_bxz,slct_sexz,slct_byz,slct_seyz, p_smr, ld_top_multi,set_pval_smr, set_pval_gwas,set_pval_eqtl, snp4msmr);
                if(snp_count==-9) continue;
                
                // output snp set list of MSMR test
                string setstr=probename+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int j=0;j<snp4msmr.size();j++)
                {
                    setstr=snp4msmr[j]+'\n';
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                setstr="end\n";
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                // end of output
                
                // step6: HEIDI test
                
                long nsnp=-9;
                double pdev=-9;
                //if(!heidioffFlag)  pdev= heidi_test(&bdata,slct_zsxz,slctId, slct_maxid, ld_top,  threshold,  m_hetero, slct_byz, slct_seyz,slct_bxz,slct_sexz, nsnp );
                smrwk.curId.swap(slctId);
                smrwk.bxz.swap(slct_bxz);
                smrwk.sexz.swap(slct_sexz);
                smrwk.byz.swap(slct_byz);
                smrwk.seyz.swap(slct_seyz);
                smrwk.rs.swap(slct_snpName);
                smrwk.allele1.swap(slct_a1);
                smrwk.allele2.swap(slct_a2);
                smrwk.pyz.swap(slct_pyz);
                smrwk.zxz.swap(slct_zxz);
                smrwk.bpsnp.swap(slct_bpsnp);
                smrwk.snpchrom.swap(slct_snpchr);
                smrwk.freq.swap(slct_freq);
                
                if(!heidioffFlag) {
                    printf("Conducting HEIDI test...\n");
                    pdev= heidi_test_new(&bdata,&smrwk, ld_top,  threshold,  m_hetero, nsnp,ld_min,opt_hetero,sampleoverlap, theta);
                }
                
                
                outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(bdata._mu[bdata._include[out_raw_id]] / 2) + '\t';
                //outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t' + dtos(set_pval_gwas) + '\t';
                //outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t' + dtos(set_pval_eqtl) + '\t';
                //outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + atos(snp_count) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t';
                outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t';
                outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr))
                {
                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                    exit(EXIT_FAILURE);
                }
                write_count++;
                
            }
        } else {
            // top-SNP (/ ref-SNP) based
            for(int i=0;i<probNum;i++)
            {
                
                progr1=1.0*i/probNum;
                if(progr1-progr0-0.05>1e-6 || i+1==probNum)
                {
                    if(i+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
               }
                
                int probebp=esdata._epi_bp[i];
                int probechr=esdata._epi_chr[i];
                string probename=esdata._epi_prbID[i];
                string probegene=esdata._epi_gene[i];
                init_smr_wk(&smrwk);
                smrwk.cur_prbidx=i;
                // step1: get cis-eQTLs
                printf("\nInitiating the workspace of probe %s for multi-SNP SMR analysis....\n",probename.c_str());
                long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refSNP, cis_itvl, heidioffFlag);
                if(refSNP!=NULL && maxid==-9) {
                    printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                    continue;
                } //ref heidi SNP is not in selected SNPs
                if (smrwk.bxz.size() == 0) {
                    printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                    continue;
                }
                printf("%ld SNPs are included from the cis-region of the probe %s.\n",smrwk.bxz.size(),probename.c_str());
                //now if you sepcify reference SNP, maxid point to this SNP, otherwise maxid is -9
                // step2: get top-SNP
                printf("Checking the top-SNP in the region....\n");
                Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
                Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
                VectorXd zsxz;
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(refSNP==NULL) maxid=max_abs_id(zsxz); // now maxid point to the sig eQTL SNP or ref SNP in the new datastruct(not the raw).
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                //double computing, consistency should be checked
                for(int tid=0;tid<zsxz.size();tid++) {
                    if(fabs(zsxz(tid)-smrwk.zxz[tid])>1e-3)
                    {
                        printf("ERROR: zxz not consistent %f vs %f. please report this.\n",zsxz(tid),smrwk.zxz[tid]);
                        exit(EXIT_FAILURE);
                    }
                }
                string topsnpname=smrwk.rs[maxid];
                printf("The top SNP of probe %s is %s with p-value %e.\n", probename.c_str(), topsnpname.c_str(),pxz_val);
                if(refSNP==NULL && pxz_val>p_smr){
                    printf("WARNING: no SNP passed the p-value threshold %e for Multiple-SNP SMR analysis for probe %s.\n", p_smr, probename.c_str());
                    continue;
                } else {
                    printf("Conducting multi-SNP SMR and HEIDI test for probe %s...\n", probename.c_str());
                }
                //cout<<maxid<<":"<<topsnpname<<":"<<esdata._esi_rs[smrwk.curId[maxid]]<<":"<<bdata._snp_name[bdata._include[smrwk.curId[maxid]]]<<":"<<gdata.snpName[smrwk.curId[maxid]]<<endl;
                // step3: extract SNPs around the --set-wind around sig (or ref) SNP
                if(expanWind!=-9) printf("Extracting SNPs in a specified window around top-SNP/ref-SNP....\n");
                else printf("Extracting SNPs in the cis-region....\n");
                vector<uint32_t> slctId;
                vector<int> slct_bpsnp,slct_snpchr;
                vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_zsxz,slct_zxz, slct_pyz,slct_freq; //slct_zsxz,slct_zxz would be removed one of them
                vector<string> slct_snpName,slct_a1, slct_a2;
                long slct_maxid=-9;
                if(expanWind!=-9){
                    for(int j=0;j<zsxz.size();j++)
                    {
                        int maxid_bp=smrwk.bpsnp[maxid];
                        int tmplower=maxid_bp-expanWind>0?maxid_bp-expanWind:0;
                        int tmpupper=maxid_bp+expanWind;
                        if(smrwk.bpsnp[j]>=tmplower && smrwk.bpsnp[j]<=tmpupper)
                        {
                            slctId.push_back(smrwk.curId[j]); // for get X
                            slct_bxz.push_back(smrwk.bxz[j]);
                            slct_sexz.push_back(smrwk.sexz[j]);
                            slct_byz.push_back(smrwk.byz[j]);
                            slct_seyz.push_back(smrwk.seyz[j]);
                            slct_zsxz.push_back(zsxz(j));
                            slct_snpName.push_back(smrwk.rs[j]);
                            if(j==maxid) slct_maxid=slctId.size()-1;
                            slct_zxz.push_back(smrwk.zxz[j]);
                            slct_pyz.push_back(smrwk.pyz[j]);
                            slct_bpsnp.push_back(smrwk.bpsnp[j]);
                            slct_snpchr.push_back(smrwk.snpchrom[j]);
                            slct_a1.push_back(smrwk.allele1[j]);
                            slct_a2.push_back(smrwk.allele2[j]);
                            slct_freq.push_back(smrwk.freq[j]);
                        }
                    }
                }else {
                    slctId.swap(smrwk.curId);
                    slct_bxz.swap(smrwk.bxz);
                    slct_sexz.swap(smrwk.sexz);
                    slct_byz.swap(smrwk.byz);
                    slct_seyz.swap(smrwk.seyz);
                    slct_snpName.swap(smrwk.rs);
                    slct_maxid=maxid;
                    for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
                    slct_zxz.swap(smrwk.zxz);
                    slct_pyz.swap(smrwk.pyz);
                    slct_bpsnp.swap(smrwk.bpsnp);
                    slct_snpchr.swap(smrwk.snpchrom);
                    slct_a1.swap(smrwk.allele1);
                    slct_a2.swap(smrwk.allele2);
                    slct_freq.swap(smrwk.freq);
                }
                
                int out_raw_id = slctId[slct_maxid];
                double bxz_max = slct_bxz[slct_maxid];
                double sexz_max = slct_sexz[slct_maxid];
                double byz_max = slct_byz[slct_maxid];
                double seyz_max = slct_seyz[slct_maxid];
                double bxy_max = byz_max / bxz_max;
                double sexy_max = sqrt((seyz_max * seyz_max * bxz_max * bxz_max + sexz_max * sexz_max * byz_max * byz_max) / (bxz_max * bxz_max * bxz_max * bxz_max));
                double chisqxy = bxy_max*bxy_max / (sexy_max*sexy_max);
                double zsxz_max = bxz_max / sexz_max;
                double zsyz_max = byz_max / seyz_max;
                double pxz_max = pchisq(zsxz_max * zsxz_max, 1);
                double pyz_max = pchisq(zsyz_max * zsyz_max, 1);
                double pxy_max = pchisq(chisqxy, 1);
                // cout<<slct_maxid<<":"<<slct_snpName[slct_maxid]<<endl;
                printf("%ld SNPs are included into SMR and HEIDI test.\n",slctId.size());
                
                double set_pval_smr=-9;
                double set_pval_gwas=-9;
                double set_pval_eqtl=-9;
                vector<string> snp4msmr;
                int snp_count=smr_setbased_test(&bdata, slctId, slct_bxz,slct_sexz,slct_byz,slct_seyz, p_smr, ld_top,set_pval_smr, set_pval_gwas,set_pval_eqtl,snp4msmr);
                if(snp_count==-9) continue;
                // output snp set list
                string setstr=probename+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int j=0;j<snp4msmr.size();j++)
                {
                    setstr=snp4msmr[j]+'\n';
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                setstr="end\n";
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                // end of output
                
                // output gene list
                int Lbp= *min_element(slct_bpsnp.begin(),slct_bpsnp.end());
                int Rbp= *max_element(slct_bpsnp.begin(),slct_bpsnp.end());
                string gstr=atos(probechr) + "\t" + atos(Lbp) + "\t" + atos(Rbp)+"\t"+ probename + "\n";
                if(fputs_checked(gstr.c_str(),glst))
                {
                    printf("ERROR: in writing file %s .\n", genelstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                // end of output
                
                // step6: HEIDI test
                long nsnp=-9;
                double pdev=-9;
                smrwk.curId.swap(slctId);
                smrwk.bxz.swap(slct_bxz);
                smrwk.sexz.swap(slct_sexz);
                smrwk.byz.swap(slct_byz);
                smrwk.seyz.swap(slct_seyz);
                smrwk.rs.swap(slct_snpName);
                smrwk.zxz.swap(slct_zsxz);
                smrwk.allele1.swap(slct_a1);
                smrwk.allele2.swap(slct_a2);
                smrwk.pyz.swap(slct_pyz);
                smrwk.bpsnp.swap(slct_bpsnp);
                smrwk.snpchrom.swap(slct_snpchr);
                smrwk.freq.swap(slct_freq);
                
                if(!heidioffFlag) {
                    printf("Conducting HEIDI test...\n");
                    pdev= heidi_test_new(&bdata,&smrwk, ld_top,  threshold,  m_hetero, nsnp,ld_min,opt_hetero,sampleoverlap, theta);
                    printf("HEIDI test complete.\n");
                }
                
                outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(bdata._mu[bdata._include[out_raw_id]] / 2) + '\t';
                //outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t' + dtos(set_pval_gwas) + '\t';
                //outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t' + dtos(set_pval_eqtl) + '\t';
                //outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + atos(snp_count) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t';
                outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t';
                outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr))
                {
                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                    exit(EXIT_FAILURE);
                }
                write_count++;
                
            }



        }

        cout<<"\nMultiple-SNP SMR and HEIDI analyses completed.\nSMR and heterogeneity analysis results of "<<write_count<<" sets have been saved in the file [" + smrfile + "]."<<endl;
        cout<<"SNP sets included in multi-SNP SMR have been saved in the file [" + setlstfile + "]."<<endl;
        fclose(smr);
        fclose(setlst);
        fclose(glst);
        free_gwas_data( &gdata);
        
    }
    void getMetaEsi(vector<eqtlInfo> &eqtls, vector<string> &besds)
    {
        vector<string> alfairs;
        map<string, int>::iterator iter;
        long f2r=besds.size();
        map<string, int> allel_map;
        for (int i = 0; i < f2r; i++)
        {
            eqtlInfo etmp;
            string eifile = besds[i]+".esi";
            read_esifile(&etmp, eifile);
            eifile = besds[i]+".epi";
            read_epifile(&etmp, eifile);
            eqtls.push_back(etmp);
        }
        printf("\nPerforming Allele checking. This step could be a little long....\n");
        printf("We set the first BESD as the baseline to conduct allele check.\n");
        vector<string> commrs;
        vector<string> commprb;
        vector< vector<int> > epi_include;
        vector< vector<int> > esi_include;
        vector< vector<bool> > reverse;
        for(int i=0; i<eqtls[0]._snpNum;i++)
        {
            string rs=eqtls[0]._esi_rs[i];
            allel_map.clear();
            allel_map.insert(pair<string,int>(eqtls[0]._esi_allele1[i],0));
            allel_map.insert(pair<string,int>(eqtls[0]._esi_allele2[i],1));
            bool hitall=false;
            vector<int> icldtmp;
            vector<bool> revtmp;
            icldtmp.push_back(i);
            revtmp.push_back(false);
            for(int j=1;j<eqtls.size();j++)
            {
                iter = eqtls[j]._snp_name_map.find(rs);
                if (iter != eqtls[j]._snp_name_map.end()) {
                    int id=iter->second;
                    allel_map.insert(pair<string,int>(eqtls[j]._esi_allele1[id],allel_map.size()));
                    allel_map.insert(pair<string,int>(eqtls[j]._esi_allele2[id],allel_map.size()));
                    if(allel_map.size()>2) {
                        alfairs.push_back(rs);
                        hitall=false;
                        break;
                    } else{
                        hitall=true;
                        icldtmp.push_back(id);
                        if(eqtls[j]._esi_allele1[id]==eqtls[0]._esi_allele1[i] && eqtls[j]._esi_allele2[id]==eqtls[0]._esi_allele2[i])  revtmp.push_back(false);
                        else if(eqtls[j]._esi_allele1[id]==eqtls[0]._esi_allele2[i] && eqtls[j]._esi_allele2[id]==eqtls[0]._esi_allele1[i])  revtmp.push_back(true);
                        else printf("Can't happen!\n");
                    }
                } else {
                    hitall=false;
                    break;
                }
            }
            if(hitall) {
                commrs.push_back(rs);
                esi_include.push_back(icldtmp);
                reverse.push_back(revtmp);
            }
        }
        printf("%ld common SNPs are included from %ld datasets.\n", commrs.size(), eqtls.size());
        
        printf("\nPerforming common probe selection....\n");
        for(int i=0; i<eqtls[0]._probNum;i++)
        {
            string prb=eqtls[0]._epi_prbID[i];
            bool hitall=false;
            vector<int> icldtmp;
            icldtmp.push_back(i);
            for(int j=1;j<eqtls.size();j++)
            {
                iter = eqtls[j]._probe_name_map.find(prb);
                if (iter != eqtls[j]._probe_name_map.end()) {
                    int id=iter->second;
                    hitall=true;
                    icldtmp.push_back(id);
                } else {
                    hitall=false;
                    break;
                }
            }
            if(hitall) {
                commprb.push_back(prb);
                epi_include.push_back(icldtmp);
            }
        }
        
        printf("%ld common probes are included from %ld datasets.\n", commprb.size(), eqtls.size());
        for(int i=0;i<eqtls.size();i++){
            eqtls[i]._esi_include.clear();
            eqtls[i]._snp_name_map.clear();
            for(int j=0;j<esi_include.size();j++)
            {
                int idx=esi_include[j][i];
                eqtls[i]._esi_include.push_back(idx);
                eqtls[i]._snp_name_map.insert(pair<string,int>(eqtls[i]._esi_rs[idx],idx));
                if(reverse[j][i]) {
                    eqtls[i]._esi_gd[idx]=1;
                    //cout<<eqtls[i]._esi_rs[idx]<<":"<<commrs[j]<<endl;
                }
            }
        }
        
        for(int i=0;i<eqtls.size();i++){
            eqtls[i]._include.clear();
            eqtls[i]._probe_name_map.clear();
            for(int j=0;j<epi_include.size();j++)
            {
                int idx=epi_include[j][i];
                eqtls[i]._include.push_back(idx);
                eqtls[i]._probe_name_map.insert(pair<string,int>(eqtls[i]._epi_prbID[idx],idx));
            }
        }

        
    }
    void extract_one_probe(eqtlInfo* esdata,int prbidx, vector<double> &bxz, vector<double> &sexz, vector<string> &rs,vector<int> &curId)
    {
        //here prbidx should be the index of esdata->_epi_probeID not the index of esdata->_include
        int i=prbidx;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size(); j++) // esdata->_esi_include.size() should equal esdata->_snp_num
            {
                if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    bxz.push_back(esdata->_bxz[i][j]);
                    sexz.push_back(esdata->_sexz[i][j]);
                    curId.push_back(j);
                    rs.push_back(esdata->_esi_rs[j]);
                }
            }
            
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                    bxz.push_back(esdata->_val[beta_start+j]);
                    sexz.push_back(esdata->_val[se_start+j]);
                    curId.push_back(ge_rowid); //save snp id of the raw
                    rs.push_back(esdata->_esi_rs[ge_rowid]);
            }
        }
    }
    void getCommonPerProbe(vector<eqtlInfo> &eqtls, int i,vector< vector<double> > &slct_beta,vector< vector<double> > &slct_se,vector<int> &slct_idx)
    {
        vector< vector<double> > beta, se;
        vector< vector<int> > ids;
        vector< map<int,int> > idmaps;
        map<int, int>::iterator iter;
        for(int j=0;j<eqtls.size();j++)
        {
            vector<double> bxz,sexz;
            vector<string> rs;
            vector<int> curId;
            extract_one_probe(&eqtls[j],i,bxz,sexz,rs,curId);
            string tmstr="th";
            if(!j) tmstr="st";
            else if(j==1) tmstr="nd";
            else if (j==2) tmstr="rd";
            else tmstr="th";
            printf("%ld SNPs are included from the %d%s BESD dataset.\n", rs.size(), j+1,tmstr.c_str());
            for(int k=0;k<curId.size();k++)
            {
                if(eqtls[j]._esi_gd[curId[k]]) {
                    bxz[k]*=-1;
                   // printf("Beta of SNP %s has been changed to %f because of allele check. \n",rs[k].c_str(),bxz[k]);
                }
            }
            map<int,int> tmpmap;
            for(int k=0;k<rs.size();k++) tmpmap.insert(pair<int,int>(curId[k],k));
            idmaps.push_back(tmpmap);
            beta.push_back(bxz);
            se.push_back(sexz);
            ids.push_back(curId);
        }
        slct_beta.resize(eqtls.size());
        slct_se.resize(eqtls.size());
        for(int j=0;j<ids[0].size();j++)
        {
            int snpid=ids[0][j];
            bool hitall=false;
            vector<double> tbeta;
            vector<double> tse;
            tbeta.push_back(beta[0][j]);
            tse.push_back(se[0][j]);
            for( int k=1;k<idmaps.size();k++)
            {
                iter = idmaps[k].find(snpid);
                if (iter != idmaps[k].end()) {
                    hitall=true;
                    tbeta.push_back(beta[k][iter->second]);
                    tse.push_back(se[k][iter->second]);
                } else {
                    hitall= false;
                    break;
                }
            }
            if(hitall) {
                for(int k=0;k<eqtls.size();k++)
                {
                    slct_beta[k].push_back(tbeta[k]);
                    slct_se[k].push_back(tse[k]);
                }
                slct_idx.push_back(snpid);
            }
        }
    }
    void getV(MatrixXd &V, vector< vector<double> > &slct_beta,vector< vector<double> > &slct_se)
    {
        vector<double> SD,MEAN;
        vector<double> bi_bj;
        for(int i=0;i<slct_se.size();i++) {
            bi_bj.clear();
            for(int j=0;j<slct_se[i].size();j++) bi_bj.push_back(slct_se[i][j]*slct_se[i][j]);
            MEAN.push_back(mean(bi_bj));
            SD.push_back(sqrt(var(slct_se[i])));
        }
        for(int i=0;i<slct_beta.size();i++)
        {
            for(int j=i;j<slct_beta.size();j++) {
                bi_bj.clear();
                for(int k=0;k<slct_beta[i].size();k++)
                    bi_bj.push_back(slct_beta[i][k]-slct_beta[j][k]);
                V(i,j)=V(j,i)=(MEAN[i]+MEAN[j]+cov(bi_bj,slct_beta[j])-cov(bi_bj,slct_beta[i]))/(2*sqrt(MEAN[i]*MEAN[j]));
            }
        }
    }
    
    void inverse_V(MatrixXd &Vi, bool &determinant_zero)
    {
        SelfAdjointEigenSolver<MatrixXd> eigensolver(Vi);
        VectorXd eval = eigensolver.eigenvalues();
        for(int i=0;i<eval.size();i++)
        {
            if(fabs(eval(i))<1e-6) {
                determinant_zero=true;
                eval(i)=0;
            } else {
                eval(i) = 1.0 / eval(i);
            }
        }
        Vi = eigensolver.eigenvectors() * DiagonalMatrix<double, Dynamic, Dynamic>(eval) * eigensolver.eigenvectors().transpose();
    }
    void write_epi(string outFileName, eqtlInfo* esdata)
    {
        printf("\nGenerating the .epi file...\n");
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the EPI file " + epifile + " to save!");
        for (int j = 0;j <esdata->_include.size(); j++) {
            epi<<((esdata->_epi_chr[esdata->_include[j]]==-9)?"NA":atos(esdata->_epi_chr[esdata->_include[j]]))<<'\t'<<esdata->_epi_prbID[esdata->_include[j]]<<'\t'<<esdata->_epi_gd[esdata->_include[j]]<<'\t'<<((esdata->_epi_bp[esdata->_include[j]]==-9)?"NA":atos(esdata->_epi_bp[esdata->_include[j]]))<<'\t'<<esdata->_epi_gene[esdata->_include[j]]<<'\t'<<((esdata->_epi_orien[esdata->_include[j]]=='N')?"NA":atos(esdata->_epi_orien[esdata->_include[j]]))<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",esdata->_include.size(),epifile.c_str());
        
    }
    void write_esi(string outFileName, eqtlInfo* esdata)
    {
        printf("\nGenerating the .esi file...\n");
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the ESI file to save!");
        for (int j = 0;j <esdata->_esi_include.size(); j++) {
            esi<<((esdata->_esi_chr[esdata->_esi_include[j]]==-9)?"NA":atos(esdata->_esi_chr[esdata->_esi_include[j]])) <<'\t'<<esdata->_esi_rs[esdata->_esi_include[j]]<<'\t'<<esdata->_esi_gd[esdata->_esi_include[j]]<<'\t'<<((esdata->_esi_bp[esdata->_esi_include[j]]==-9)?"NA":atos(esdata->_esi_bp[esdata->_esi_include[j]]))<<'\t'<<esdata->_esi_allele1[esdata->_esi_include[j]]<<'\t'<<esdata->_esi_allele2[esdata->_esi_include[j]]<<'\t'<<(esdata->_esi_freq[esdata->_esi_include[j]]+9>1e-6?atos(esdata->_esi_freq[esdata->_esi_include[j]]):"NA")<<'\n';
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",esdata->_esi_include.size(),esifile.c_str());
    }
    void write_sbesd3(char* outFileName,vector<uint64_t> &cols, vector<uint32_t> &rowids, vector<float> &val)
    {
        printf("\nGenerating the .besd file...\n");
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        uint32_t filetype=SPARSE_FILE_TYPE_3F;
        fwrite (&filetype,sizeof(uint32_t), 1, smr1);
        
        uint64_t valNum=val.size();
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        printf("eQTL summary statistics have been saved in binary file %s.\n", outFileName);

    }
     void meta_nooverlap_func(vector< vector<double> > &slct_beta,vector< vector<double> > &slct_se,vector<int> &slct_idx)
    {
        for(int j=0;j<slct_idx.size();j++)
        {
            double numerator=0.0;
            double deno=0.0;
            for(int k=0;k<slct_beta.size();k++)
            {
                double tmp2=slct_se[k][j]*slct_se[k][j];
                deno+=1/tmp2;
                numerator+=slct_beta[k][j]/tmp2;
            }
            slct_beta[0][j]=numerator/deno;
            slct_se[0][j]=1/sqrt(deno);
        }

    }
    void meta_overlap_func(vector< vector<double> > &slct_beta,vector< vector<double> > &slct_se,vector<int> &slct_idx,bool detailout, vector<int> &noninvertible, vector<int> &negativedeno)
    {
        
         MatrixXd V(slct_beta.size(),slct_beta.size());
         getV(V, slct_beta,slct_se);
        if(detailout){
            FILE* tmpfile=fopen("est_cor.txt","w");
            if(!tmpfile)
            {
                printf("error open file.\n");
                exit(EXIT_FAILURE);
            }
            for(int t=0;t<V.cols();t++)
            {
                string str="";
                for(int tt=0;tt<V.rows();tt++)
                {
                    str+=atos(V(tt,t))+'\t';
                }
                str+='\n';
                fputs(str.c_str(),tmpfile);
            }
            
            fclose(tmpfile);

        }
         
         for(int j=0;j<slct_idx.size();j++)
         {
             VectorXd sev(slct_se.size());
             for(int k=0;k<slct_se.size();k++) sev(k)=slct_se[k][j];
             MatrixXd W=sev*sev.transpose();
             W=W.array()*V.array();
             bool determinant_zero=false;
             inverse_V(W,determinant_zero);
             if(determinant_zero) noninvertible.push_back(slct_idx[j]);
             double deno=W.sum();
             if(deno<=0) {
                 negativedeno.push_back(slct_idx[j]);
                 slct_beta[0][j]=0;
                 slct_se[0][j]=-9;
             } else {
                 VectorXd colsum=W.colwise().sum();
                 double numerator=0.0;
                 for(int k=0;k<slct_beta.size();k++) numerator+=colsum(k)*slct_beta[k][j];
                 slct_beta[0][j]=numerator/deno;
                 slct_se[0][j]=1/sqrt(deno);
             }
             
         }
        
    }
    void combine_epi(vector<smr_probeinfo> &probeinfo, vector<string> &besds, vector<uint64_t> &nprb)
    {
        long counter = 0;
        map<string, int> prb_map;
        map<string, int> prbbp_map;
        map<string, int>::iterator iter;
        nprb.clear();
        char inputname[FNAMESIZE];
        for (int i = 0; i < besds.size(); i++)
        {
            eqtlInfo etmp;
            memcpy(inputname,besds[i].c_str(),besds[i].length()+1);
            char* suffix=inputname+besds[i].length();
            memcpy(suffix,".epi",5);
            read_epifile(&etmp, inputname);
            nprb.push_back(etmp._probNum);
            for (int j = 0; j<etmp._probNum; j++)
            {
                string crsbpstr=etmp._epi_prbID[j]+":"+atos(etmp._epi_bp[j]);
                prb_map.insert(pair<string, int>(etmp._epi_prbID[j].c_str(), counter));
                prbbp_map.insert(pair<string, int>(crsbpstr.c_str(), counter));
                if(prb_map.size() != prbbp_map.size())
                {
                    printf("ERROR: inconsistent position for the probe %s  in different .epi files. Please check.\n", etmp._epi_prbID[j].c_str()) ;
                    exit(EXIT_FAILURE);
                }
                
                if (counter < prb_map.size())
                {
                    smr_probeinfo probinfotmp;
                    counter=prb_map.size();
                    probinfotmp.probechr=etmp._epi_chr[j];
                    strcpy2(&probinfotmp.probeId, etmp._epi_prbID[j]);
                    probinfotmp.bp=etmp._epi_bp[j];
                    probinfotmp.gd=etmp._epi_gd[j];
                    strcpy2(&probinfotmp.genename, etmp._epi_gene[j]);
                    probinfotmp.orien=etmp._epi_orien[j];
                    probinfotmp.bfilepath=NULL;
                    probinfotmp.esdpath=NULL;
                    probinfotmp.ptr=new int[besds.size()];
                    for(int k=0;k<besds.size();k++){
                        if(i==k){
                            probinfotmp.ptr[k]=j;
                        } else {
                            probinfotmp.ptr[k]=-9;
                        }
                    }
                    probeinfo.push_back(probinfotmp);
                    
                } else {
                    iter=prb_map.find(etmp._epi_prbID[j]);
                    if(iter!=prb_map.end())
                    {
                        probeinfo[iter->second].ptr[i]=j; //probeinfo with prb_map
                    }
                    else
                    {
                        printf("ERROR: This would never happen. please help to report this bug.\n") ;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        printf("Total %ld probes to be included from %ld epi files.\n",probeinfo.size(),besds.size());
    }
    void combine_esi(vector<smr_snpinfo> &snpinfo, vector<string> &besds, vector<uint64_t> &nsnp)
    {
        long ex_counter = 0;
        snpinfo.clear();
        vector<string> ex_snp;
        map<string, int> in_map;
        map<string, int> ex_map;
        map<string, int>::iterator iter;
        nsnp.clear();
        char inputname[FNAMESIZE];
        printf("\nPerforming Allele checking. This step could be a little long....\n");
        for (int i = 0; i < besds.size(); i++)
        {
            eqtlInfo etmp;
            memcpy(inputname,besds[i].c_str(),besds[i].length()+1);
            char* suffix=inputname+besds[i].length();
            memcpy(suffix,".esi",5);
            read_esifile(&etmp, inputname);
            nsnp.push_back(etmp._snpNum);
            for (int j = 0; j<etmp._snpNum; j++)
            {
                if(ex_map.size()>0) {
                    iter=ex_map.find(etmp._esi_rs[j]);
                    if(iter!=ex_map.end()) continue;
                }
                
                iter=in_map.find(etmp._esi_rs[j]);
                if(iter==in_map.end())
                {
                    in_map.insert(pair<string, int>(etmp._esi_rs[j].c_str(), snpinfo.size()));//the second is snpinfo id
                    
                    smr_snpinfo snpinfotmp;
                    snpinfotmp.snpchr=etmp._esi_chr[j];
                    strcpy2(&snpinfotmp.snprs, etmp._esi_rs[j]);
                    snpinfotmp.bp=etmp._esi_bp[j];
                    snpinfotmp.gd=etmp._esi_gd[j];
                    strcpy2(&snpinfotmp.a1, etmp._esi_allele1[j]);
                    strcpy2(&snpinfotmp.a2, etmp._esi_allele2[j]);
                    snpinfotmp.freq=etmp._esi_freq[j];
                    
                    snpinfotmp.rstr=new int[besds.size()];
                    snpinfotmp.revs=new bool[besds.size()];
                    for(int k=0;k<besds.size();k++){
                        if(i==k){
                            snpinfotmp.rstr[k]=j;
                            snpinfotmp.revs[k]=false;
                        } else {
                            snpinfotmp.rstr[k]=-9;
                            snpinfotmp.revs[k]=false;
                        }
                    }
                    snpinfo.push_back(snpinfotmp);
                } else {
                    
                    if((snpinfo[iter->second].snpchr != etmp._esi_chr[j]) ||(snpinfo[iter->second].bp != etmp._esi_bp[j]))
                    {
                        printf("ERROR: inconsistent chromosome or position for the SNP %s in different .epi files. Please check.\n", etmp._esi_rs[j].c_str()) ;
                        exit(EXIT_FAILURE);
                    }
                    string a1=etmp._esi_allele1[j];
                    string a2=etmp._esi_allele2[j];
                    string a3=snpinfo[iter->second].a1;
                    string a4=snpinfo[iter->second].a2;
                    if(a1==a3 && a2==a4) {
                        snpinfo[iter->second].rstr[i]=j;
                    }
                    else if(a1==a4 && a2==a3 ){
                        snpinfo[iter->second].rstr[i]=j;
                        snpinfo[iter->second].revs[i]=true;
                    }
                    else {
                        ex_snp.push_back(etmp._esi_rs[j]);
                        ex_map.insert(pair<string,int>(etmp._esi_rs[j],ex_counter));
                        ex_counter=ex_map.size();
                        in_map.erase(iter->first);
                    }
                }
            }
        }
        if(in_map.size() + ex_map.size() != snpinfo.size()){
            printf("ERROR: bugs found in allele check. please report.\n") ;
            exit(EXIT_FAILURE);
        }
        long ttl_snp_common=snpinfo.size();
        string failName="failed.snp.list";
        FILE* failfptr=fopen(failName.c_str(),"w");
        if(failfptr==NULL)
        {
            printf("ERROR: failed in open file %s.\n",failName.c_str()) ;
            exit(EXIT_FAILURE);
        }
        vector<smr_snpinfo> snpinfo_adj;
        snpinfo_adj.resize(in_map.size());
        int ids=0;
        for(int i=0; i<snpinfo.size();i++)
        {
            iter=in_map.find(snpinfo[i].snprs);
            if(iter!=in_map.end()) {
                snpinfo_adj[ids++] = snpinfo[i];
            } else {
                string snpstr=string(snpinfo[i].snprs) + '\n';
                if(fputs_checked(snpstr.c_str(),failfptr))
                {
                    printf("ERROR: in writing file %s .\n", failName.c_str());
                    exit(EXIT_FAILURE);
                }
                if(snpinfo[i].a1) free2(&snpinfo[i].a1);
                if(snpinfo[i].a2) free2(&snpinfo[i].a2);
                if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
                if(snpinfo[i].rstr) free2(&snpinfo[i].rstr);
                if(snpinfo[i].revs) free2(&snpinfo[i].revs);
            }
        }
        snpinfo.swap(snpinfo_adj);
        snpinfo_adj.clear();
        fclose(failfptr);
        
        printf("Total %ld SNPs to be included from %ld esi files. %ld SNPs failed in allele check. %ld SNPs included in analysis.\n",ttl_snp_common,besds.size(),ex_snp.size(),in_map.size());
        printf("%ld SNPs that failed in allele check were saved in file %s.\n",ex_snp.size(),failName.c_str());
    }
    int comp_epi(const void *a,const void *b){ return (((*(smr_probeinfo *)a).probechr>(*(smr_probeinfo *)b).probechr) || ( ((*(smr_probeinfo *)a).probechr ==(*(smr_probeinfo *)b).probechr) && ((*(smr_probeinfo *)a).bp > (*(smr_probeinfo *)b).bp) ))?1:-1; }
    void get_BesdHeaders(char* besdFileName, vector<int> &headers)//
    {
        headers.resize(RESERVEDUNITS);
        FILE* besd=fopen(besdFileName,"rb");
        if(besd==NULL) {
            exit(EXIT_FAILURE);
        }
        printf("Reading eQTL summary data from %s. \n",besdFileName);
        if(fread(&headers[0], sizeof(int),RESERVEDUNITS, besd)<1)
        {
            printf("ERROR: File %s read failed!\n", besdFileName);
            exit(EXIT_FAILURE);
        }
        fclose(besd);
    }
    void check_besds_format( vector<string> &besds, vector<int> &format, vector<int> &smpsize) {
        
        char inputname[FNAMESIZE];
        format.clear();
        smpsize.clear();
        vector<int> headers;
        for(int i=0;i<besds.size();i++)
        {
            string tmpstr=besds[i]+".besd";
            memcpy(inputname,tmpstr.c_str(),tmpstr.length()+1);
            get_BesdHeaders(inputname, headers);
            format.push_back(headers[0]);
            if(headers[0]==DENSE_FILE_TYPE_1 || headers[0]==SPARSE_FILE_TYPE_3F) smpsize.push_back(-9);
            else smpsize.push_back(headers[1]);
        }
    }
    void extract_prb_sparse(FILE* fptr, uint64_t pid, uint64_t probnum,vector<uint32_t> &row_ids, vector<float> &betases)
    {
        if(pid>probnum) {
            printf("ERROR: probe index %llu is larger than the totoal probe number %llu.\n", pid, probnum);
            exit(EXIT_FAILURE);
        }
        row_ids.clear();
        betases.clear();
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        if(indicator==SPARSE_FILE_TYPE_3F || indicator==SPARSE_FILE_TYPE_3)
        {
            int infoLen=sizeof(uint32_t);
            if(indicator==SPARSE_FILE_TYPE_3)
            {
                infoLen=RESERVEDUNITS*sizeof(int);
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    printf("The sample size is %d.\n",ss);
                }
                delete[] indicators;
            }
            
            uint64_t colNum=(probnum<<1)+1;
            uint64_t valnum=readuint64(fptr);
            fseek(fptr,(pid<<1)*sizeof(uint64_t),SEEK_CUR);
            uint64_t betaStart=readuint64(fptr);
            uint64_t seStart=readuint64(fptr);
            long num=seStart-betaStart;
            if(num>0)
            {
                row_ids.resize(num);
                betases.resize(2*num);
                uint64_t rowSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=infoLen + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valnum*sizeof(uint32_t);
                fseek(fptr, rowSTART+betaStart*sizeof(uint32_t), SEEK_SET);
                fread(&row_ids[0], sizeof(uint32_t),num,fptr);
                fseek(fptr,valSTART+betaStart*sizeof(float),SEEK_SET);
                fread(&betases[0],sizeof(float), 2*num,fptr);
                
            }
        }
        else
        {
            printf("ERROR: OSCA format. please use OSCA to do this analysis.\n");
            exit(EXIT_FAILURE);
        }
    }
    void extract_prb_dense(FILE* fptr,  uint64_t pid, uint64_t epinum,uint64_t esinum, vector<float> &betases)
    {
        if(epinum==0 || esinum==0) {
            printf("ERROR: .epi file or .esi file is empty. please check.\n");
            exit(EXIT_FAILURE);
        }
        if(pid>epinum) {
            printf("ERROR: probe index %llu is larger than the totoal probe number %llu.\n", pid, epinum);
            exit(EXIT_FAILURE);
        }
        fseek(fptr,0L,SEEK_SET);
        uint32_t indicator=readuint32(fptr);
        if(indicator==DENSE_FILE_TYPE_1 || indicator==DENSE_FILE_TYPE_3) {
            if(indicator==DENSE_FILE_TYPE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                fread(indicators, sizeof(int),(RESERVEDUNITS-1), fptr);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    printf("The sample size is %d.\n",ss);
                }
                delete[] indicators;
            }
            int infoLen=sizeof(uint32_t);
            if(indicator==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);
            
            betases.resize(2*esinum);
            fseek(fptr,((pid*esinum)<<3)+infoLen, SEEK_SET);
            fread(&betases[0], sizeof(float),2*esinum,fptr);
        }
        else
        {
            printf("ERROR: OSCA format. please use OSCA to do this analysis.\n");
            exit(EXIT_FAILURE);
        }
    }

    void pcc(MatrixXd &PCC, float* buffer_beta,float* buffer_se,long snpnum, long cohortnum, double pmecs)
    {
        //pearson correlation with pairwise.complete.obs
        double zmecs=qchisq(pmecs,1);
        vector<double> beta1,beta2;
        vector<int> pairsNoCor1,pairsNoCor2;
        double sumcor=0.0;
        int pairHasCorNUm=0;
        for( int i=0;i<cohortnum;i++)
        for(int j=i+1;j<cohortnum;j++)
        {
            
            beta1.clear();
            beta2.clear();
            for(int k=0;k<snpnum;k++)
            {
                double sei=buffer_se[i*snpnum+k];
                double betai=buffer_beta[i*snpnum+k];
                double sej=buffer_se[j*snpnum+k];
                double betaj=buffer_beta[j*snpnum+k];
                if(fabs(sei+9)>1e-6 && fabs(sej+9)>1e-6) {
                    double zi=betai/sei;
                    double zj=betaj/sej;
                    zi*=zi;
                    zj*=zj;
                    if(zi < zmecs && zj < zmecs)
                    {
                        
                        beta1.push_back(betai);
                        beta2.push_back(betaj);
                    }
                }
            }
            if(beta1.size()<1) {
                printf("WARNING: %ld SNP in common between cohort %i and cohort %d (cohort number stats form 0).\n",beta1.size(),i,j);
                printf("The correlation value of cohort %d and cohort %d would be imputed with the mean of all the correlation values excluding the diagnoal.\n",i,j);
                pairsNoCor1.push_back(i);
                pairsNoCor2.push_back(j);
                PCC(i,j)=PCC(j,i)=0;
            } else {
                double corrtmp=cor(beta1,beta2);
                sumcor +=corrtmp;
                PCC(i,j)=PCC(j,i)=corrtmp;
                pairHasCorNUm++;
            }
        }
        if(pairsNoCor1.size()>0) {
            printf("WARNING: %ld cohort pairs didn't get enough common SNPs to calcualte the correlation.\n",pairsNoCor1.size());
            if(pairHasCorNUm==0) {
                printf("ERROR: Every pair of cohort has not enough common SNPs to calcualte the correlation.\n");
                exit(EXIT_FAILURE);
            }
            double corMean=sumcor/pairHasCorNUm;
            printf("WARNING: These missing correlation values are imputed with the mean %f.\n",corMean);
            for(int i=0;i<pairsNoCor1.size();i++)
            {
                int p1=pairsNoCor1[i];
                int p2=pairsNoCor2[i];
                PCC(p1,p2)=PCC(p2,p1)=corMean;
            }
        }
        for( int i=0;i<cohortnum;i++) PCC(i,i)=1;
    }
    void subMatrix_symm(MatrixXd &to,MatrixXd &from,vector<int> &idx)
    {
        long num = idx.size();
        if(num>from.cols())
        {
            printf("ERROR: to extract sub-matrix. the size of index %ld is larger than the matrix dimension (%ld,%ld).\n",num, from.rows(),from.cols());
            exit(EXIT_FAILURE);
        } else if(num==0) {
            printf("ERROR: the size of index is 0.\n");
            exit(EXIT_FAILURE);
        } else {
            to.resize(num,num);
            for(int i=0;i<num;i++)
            for(int j=i;j<num;j++)
            to(i,j)=to(j,i)=from(idx[i],idx[j]);
        }
    }
    void mecs_per_prob(float* buffer_beta,float* buffer_se, long snpnum, long cohortnum,double pmecs,vector<int> &noninvertible, vector<int> &negativedeno)
    {
        
        MatrixXd Corr(cohortnum,cohortnum);
        printf("Estimate the cohort correlation using the beta values of pair-wised common SNPs.\n");
        printf("We exclude the significant common SNPs with a p-value threshold %e to estimate the correlation matrix.\n",pmecs);
        pcc(Corr,buffer_beta,buffer_se,snpnum,cohortnum,pmecs);
        //cout<<Corr<<endl;
#pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            vector<double> ses, betas;
            vector<int> keep;
            //MatrixXd Corr_work = Corr;
            MatrixXd Corr_work;
            int nmiss=0, miss=0;
            for(int k=0;k<cohortnum;k++)
            {
                double se=buffer_se[k*snpnum+j];
                double beta=buffer_beta[k*snpnum+j];
                buffer_se[k*snpnum+j]=-9;
                if(fabs(se+9)>1e-6){
                    ses.push_back(se);
                    betas.push_back(beta);
                    keep.push_back(k);
                    nmiss++;
                } else {
                    //removeRow(Corr_work, k-miss);
                    //removeColumn(Corr_work, k-miss);
                    miss++;
                }
            }
            
            if(nmiss==1)
            {
                buffer_beta[j]=betas[0];
                buffer_se[j]=ses[0];
            }
            else if(nmiss>1)
            {
                if(nmiss<cohortnum) subMatrix_symm(Corr_work, Corr, keep);
                else if(nmiss==cohortnum) Corr_work = Corr;
                else {
                    printf("Can't happen. I can guarantee!\n");
                }
                VectorXd sev(ses.size());
                for(int k=0;k<ses.size();k++) sev(k)=ses[k];
                MatrixXd W=sev*sev.transpose();
                W=W.array()*Corr_work.array();
                bool determinant_zero=false;
                inverse_V(W,determinant_zero);
                if(determinant_zero) noninvertible.push_back(j);
                double deno=W.sum();
                if(deno<=0) {
                    negativedeno.push_back(j);
                } else {
                    VectorXd colsum=W.colwise().sum();
                    double numerator=0.0;
                    for(int k=0;k<betas.size();k++) numerator+=colsum(k)*betas[k];
                    buffer_beta[j]=numerator/deno;
                    buffer_se[j]=1/sqrt(deno);
                }
                
            }
        }
    }
    void meta_per_prob(float* buffer_beta,float* buffer_se, long snpnum, long cohortnum)
    {
        #pragma omp parallel for
        for(int j=0;j<snpnum;j++)
        {
            double numerator=0.0;
            double deno=0.0;
            int nmiss=0;
            for(int k=0;k<cohortnum;k++)
            {
                double se=buffer_se[k*snpnum+j];
                double beta=buffer_beta[k*snpnum+j];
                buffer_se[k*snpnum+j]=-9;
                if(fabs(se+9)>1e-6){
                    double tmp2=se*se;
                    deno+=1/tmp2;
                    numerator+=beta/tmp2;
                    nmiss++;
                }
            }
            if(nmiss>0)
            {
                buffer_beta[j]=numerator/deno;
                buffer_se[j]=1/sqrt(deno);
            }
        }
    }
    void meta(char* besdlistFileName, char* outFileName, int meta_mth, double pthresh, bool cis_flag, int cis_itvl)
    {
        cis_flag=true; // for later update. !!!
        printf("NOTE: Only the information in the cis-region would be used.\n");
        
        string analysisType="";
        if(meta_mth) analysisType="MeCS";
        else analysisType="Meta";
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(besdlistFileName, besds,"eQTL summary file names");
        if(besds.size()<=1) {
            printf("Less than 2 BESD files list in %s.\n",besdlistFileName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
        printf("Checking the BESD format...\n");
        vector<int> format, smpsize;
        check_besds_format(besds, format, smpsize);
        int label=-1;
        for(int i=0;i<format.size();i++)
        {
            if(format[i]==DENSE_FILE_TYPE_1 || format[i]==DENSE_FILE_TYPE_3 ) {
                if(label==-1) {
                    label=0;
                } else if(label==1) {
                    label=2;
                    break;
                }
                
            } else if (format[i]==SPARSE_FILE_TYPE_3F || format[i]==SPARSE_FILE_TYPE_3 ) {
                if(label==-1) {
                    label=1;
                } else if(label==0) {
                    label=2;
                    break;
                }
            } else {
                printf("Some BESDs are from old sparse format. please use SMR to re-make it.\n");
                exit(EXIT_FAILURE);
            }
        }
        
        if(meta_mth) label=1; // for later update. !!!!
        
        combine_epi(probeinfo, besds,nprb);
        combine_esi(snpinfo, besds,nsnp);
        if(probeinfo.size()==0)
        {
            printf("ERROR: No probe to be included!\n");
            exit(EXIT_FAILURE);
        }
        smr_probeinfo* epiptr=&probeinfo[0];
        qsort(epiptr,probeinfo.size(),sizeof(smr_probeinfo),comp_epi);
        if(snpinfo.size()==0)
        {
            printf("ERROR: No SNP to be included!\n");
            exit(EXIT_FAILURE);
        }
        smr_snpinfo* esiptr=&snpinfo[0];
        qsort(esiptr,snpinfo.size(),sizeof(smr_snpinfo),comp_esi);
        
        long besdNum=besds.size();
        long metaPrbNum=probeinfo.size();
        long metaSNPnum=snpinfo.size();
        
        printf("\nGenerating epi file...\n");
       
        string epiName=string(outFileName)+".epi";
         FILE* efile=fopen( epiName.c_str(),"w");
        if(efile==NULL) exit(EXIT_FAILURE);
        for(int i=0;i<probeinfo.size();i++)
        {
            string chrstr;
            if(probeinfo[i].probechr==23) chrstr="X";
            else if(probeinfo[i].probechr==24) chrstr="Y";
            else chrstr=atosm(probeinfo[i].probechr);
            
            string str=chrstr+'\t'+probeinfo[i].probeId+'\t'+atos(0)+'\t'+atosm(probeinfo[i].bp)+'\t'+probeinfo[i].genename+'\t'+(probeinfo[i].orien=='*'?"NA":atos(probeinfo[i].orien))+'\n';
            if(fputs_checked(str.c_str(),efile)) 
            {
                printf("ERROR: in writing file %s .\n", epiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        printf("%ld probes have been saved in the file %s .\n", probeinfo.size(), epiName.c_str());
        
        
        printf("\nGenerating esi file...\n");
        string esiName=string(outFileName)+".esi";
        efile=fopen(esiName.c_str(),"w");
        if(efile==NULL) exit(EXIT_FAILURE);
        for(int i=0;i<snpinfo.size();i++)
        {
            string chrstr;
            if(snpinfo[i].snpchr==23) chrstr="X";
            else if(snpinfo[i].snpchr==24) chrstr="Y";
            else chrstr=atosm(snpinfo[i].snpchr);
            string str=chrstr+'\t'+snpinfo[i].snprs+'\t'+atos(0)+'\t'+atosm(snpinfo[i].bp)+'\t'+snpinfo[i].a1+'\t'+snpinfo[i].a2+'\t'+(fabs(snpinfo[i].freq+9)>1e-6?atos(snpinfo[i].freq):"NA")+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                printf("ERROR: in writing file %s .\n", esiName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        fclose(efile);
        printf("%ld SNPs have been saved in the file %s .\n", snpinfo.size(), esiName.c_str());
        
        lookup.resize(besdNum);
        for(int i=0;i<besdNum;i++) lookup[i].resize(nsnp[i]);
        for(int i=0;i<besdNum;i++)
        for(int j=0;j<lookup[i].size();j++)
        lookup[i][j]=-9;
        for(int i=0;i<metaSNPnum;i++)
        {
            for(int j=0;j<besdNum;j++)
            {
                int tmpval=snpinfo[i].rstr[j];
                if(tmpval>=0) {
                    if(tmpval >=nsnp[j])
                    {
                        printf("ERROR: bug found in snpinfo. Please report.\n");
                        exit(EXIT_FAILURE);
                    }
                    lookup[j][tmpval]=i;
                }
            }
        }
        
        printf("\nPerforming %s analysis and save the result in BESD file....\n",analysisType.c_str());
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        vector<uint64_t> cols;
        vector<uint32_t> rowids;
        vector<float> val;
        string besdName=string(outFileName)+".besd";
        efile=fopen( besdName.c_str(),"wb");
        if(efile==NULL) exit(EXIT_FAILURE);
        
        uint32_t filetype=SPARSE_FILE_TYPE_3;
        cols.resize((metaPrbNum<<1)+1);
        cols[0]=0;
        
        vector<int> ten_ints(RESERVEDUNITS);
        ten_ints[0]=filetype;
        ten_ints[1]=-9;
        ten_ints[2]=(int)snpinfo.size();
        ten_ints[3]=(int)probeinfo.size();
        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
        fwrite(&ten_ints[0],sizeof(int), RESERVEDUNITS,efile);
        
        
        float* buffer_beta=(float *)malloc(sizeof(float)*besdNum*metaSNPnum);
        if (buffer_beta == NULL) {
            printf("Memory buffer for beta values error.\n");
            exit(EXIT_FAILURE);
        } //probe major
        
        float* buffer_se=(float *)malloc(sizeof(float)*besdNum*metaSNPnum);
        if (buffer_se == NULL) {
            printf("Memory buffer for SEs error.\n");
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<besdNum*metaSNPnum;i++) buffer_se[i]=-9;
        vector<string> noninvtb_prbs;
        vector<string> nega_prbs;
        for(int i=0;i<metaPrbNum;i++)
        {
            printf("%3.0f%%\r", 100.0*i/(metaPrbNum));
            fflush(stdout);
            printf("processing with probe %s...\n",probeinfo[i].probeId);
            vector<float> betases;
            vector<uint32_t> row_ids;
            int probebp=probeinfo[i].bp;
            int probechr=probeinfo[i].probechr;
            long cohortnum=0;
            for(int j=0;j<besds.size();j++)
            {
                int pid=probeinfo[i].ptr[j];
                betases.clear();
                row_ids.clear();
                if(pid>=0)
                {
                    if(format[j]==SPARSE_FILE_TYPE_3 || format[j]==SPARSE_FILE_TYPE_3F )
                    {
                        extract_prb_sparse(fptrs[j], (uint64_t)pid, nprb[j],row_ids, betases);
                        long num=row_ids.size();
                        if(num==0) {
                            printf("WARNING: empty probe %s found in the BESD file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                            continue;
                        }
                        long alignnum=0;
                        if( cis_flag) {
                            printf("Extract the cis-region of probe %s in the BESD file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        } else {
                            printf("Extract the eQTLs of probe %s in the BESD file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        }
                        for(int k=0;k<num;k++)
                        {
                            //align each cohort to the buffer. we don't have -9 of se from sparse.
                            int idx=lookup[j][row_ids[k]];
                            if(idx>=0)
                            {
                                int snpchr=snpinfo[idx].snpchr;
                                int snpbp=snpinfo[idx].bp;
                                if(cis_flag) {
                                    //extract cis if mecs
                                    int startpos=(probebp-cis_itvl*1000)>0?(probebp-cis_itvl*1000):0;
                                    int endpos=probebp+cis_itvl*1000;
                                    if(snpchr==probechr && snpbp >= startpos && snpbp <= endpos) {
                                        if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                        else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                        buffer_se[cohortnum*metaSNPnum+idx]=betases[k+num];
                                        alignnum++;
                                    }
                                } else {
                                    if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                    else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                    buffer_se[cohortnum*metaSNPnum+idx]=betases[k+num];
                                    alignnum++;
                                }
                            }
                        }
                        if(alignnum==0) {
                            printf("No eQTLs extracted from probe %s in the BESD file %s for %s analysis.\n",probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                            continue;
                        } else {
                            printf("%ld eQTLs extracted from probe %s in the BESD file %s for %s analysis.\n",alignnum,probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                        }
                    }
                    else if(format[j]==DENSE_FILE_TYPE_1 || format[j]==DENSE_FILE_TYPE_3)
                    {
                        extract_prb_dense(fptrs[j], (uint64_t)pid, nprb[j], nsnp[j], betases);
                        printf("%llu eQTLs extracted from probe %s in the BESD file %s.\n",nsnp[j],probeinfo[i].probeId,besds[j].c_str());
                        long alignnum=0;
                        if( cis_flag) {
                            printf("Extract the cis-region of probe %s in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        } else {
                            printf("Extract the eQTLs of probe %s in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                        }
                        for(int k=0;k<nsnp[j];k++)
                        {
                            int idx=lookup[j][k];
                            if(idx>=0)
                            {
                                int snpchr=snpinfo[idx].snpchr;
                                int snpbp=snpinfo[idx].bp;
                                if(cis_flag)
                                {
                                    //extract cis if mecs
                                    int startpos=(probebp-cis_itvl*1000)>0?(probebp-cis_itvl*1000):0;
                                    int endpos=probebp+cis_itvl*1000;
                                    if(snpchr==probechr && snpbp >= startpos && snpbp <= endpos) {
                                        if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                        else buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                        buffer_se[cohortnum*metaSNPnum+idx]=betases[k+nsnp[j]];
                                        alignnum++;
                                    }
                                } else
                                {
                                    if(snpinfo[idx].revs[j]) buffer_beta[cohortnum*metaSNPnum+idx]=-1.0*betases[k];
                                    else  buffer_beta[cohortnum*metaSNPnum+idx]=betases[k];
                                    buffer_se[cohortnum*metaSNPnum+idx]=betases[k+nsnp[j]];
                                }
                            }
                        }
                        if(alignnum==0) {
                            printf("no eQTLs extracted from probe %s in the BESD file %s for %s analysis.\n",probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                            continue;
                        } else {
                            printf("%ld eQTLs extracted from probe %s in the BESD file %s for %s analysis.\n",alignnum,probeinfo[i].probeId,besds[j].c_str(),analysisType.c_str());
                        }
                    }
                    cohortnum++;
                }
                else {
                    printf("probe %s is not in the file %s.\n",probeinfo[i].probeId,besds[j].c_str());
                }
            }
            if(cohortnum==0) {
                printf("No information of probe %s is included from any cohort for %s analysis.\n\n",probeinfo[i].probeId,analysisType.c_str());
            } else if(cohortnum==1) {
                printf("The information of probe %s is included from %ld / %ld cohorts for %s analysis.\n",probeinfo[i].probeId,cohortnum,besds.size(),analysisType.c_str());
                printf("The information of the probe %s would be saved in the result.\n\n",probeinfo[i].probeId);
                
            } else if(cohortnum>1) {
                printf("The information of probe %s is included from %ld / %ld cohorts for %s analysis.\n",probeinfo[i].probeId,cohortnum,besds.size(),analysisType.c_str());
                if(meta_mth){
                    printf("Performing %s analysis of probe %s...\n",analysisType.c_str(), probeinfo[i].probeId);
                    vector<int> noninvertible, negativedeno;
                    mecs_per_prob( buffer_beta, buffer_se, metaSNPnum, cohortnum, pthresh,noninvertible,negativedeno);
                    if(noninvertible.size()>0) {
                        printf("%ld SNPs of probe %s have non-invertible S matrix.\n",noninvertible.size(), probeinfo[i].probeId);
                        noninvtb_prbs.push_back(probeinfo[i].probeId);
                    }
                    if(negativedeno.size()>0) {
                        printf("%ld SNPs of probe %s have negative 1'inv(S)1 .\n",negativedeno.size(), probeinfo[i].probeId);
                        nega_prbs.push_back(probeinfo[i].probeId);
                    }
                    printf("end of %s analysis of probe %s.\n\n",analysisType.c_str(), probeinfo[i].probeId);
                }
                else {
                    printf("Performing %s analysis of probe %s...\n",analysisType.c_str(), probeinfo[i].probeId);
                    meta_per_prob(buffer_beta,buffer_se, metaSNPnum, cohortnum);
                    printf("End of %s analysis of probe %s.\n\n",analysisType.c_str(), probeinfo[i].probeId);
                }
            }
            
            
            vector<float> tmpse;
            vector<uint32_t> tmprid;
            for(int k=0;k<metaSNPnum;k++)
            {
                if(fabs(buffer_se[k]+9)>1e-6)
                {
                    val.push_back(buffer_beta[k]);
                    rowids.push_back(k);
                    tmpse.push_back(buffer_se[k]);
                    tmprid.push_back(k);
                    buffer_se[k]=-9;
                }
            }
            for(int k=0;k<tmpse.size();k++)
            {
                val.push_back(tmpse[k]);
                rowids.push_back(tmprid[k]);
            }
            uint64_t real_num=tmpse.size();
            cols[(i<<1)+1]=real_num+cols[i<<1];
            cols[i+1<<1]=(real_num<<1)+cols[i<<1];
        }
        
        
        
        uint64_t valNum=val.size();
        fwrite(&valNum,sizeof(uint64_t),1, efile) ;
        fwrite(&cols[0],sizeof(uint64_t),cols.size(), efile);
        fwrite(&rowids[0],sizeof(uint32_t),rowids.size(), efile);
        fwrite(&val[0],sizeof(float),val.size(), efile);
        
        if(noninvtb_prbs.size()>0)
        {
            printf("\nWARNING: %ld probes have at least one eQTL whose S is non-invertible.\n",noninvtb_prbs.size());
            string filename=string(outFileName)+".non-invertible.probe.list";
            FILE* tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                printf("error open file.\n");
                exit(EXIT_FAILURE);
            }
            for(int t=0;t<noninvtb_prbs.size();t++)
            {
                string str=noninvtb_prbs[t]+'\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
            printf("These probes are saved in file %s.\n",filename.c_str());
        }
        
        if(nega_prbs.size()>0)
        {
            printf("\nWARNING: %ld probes have at least one eQTL whose 1'inv(S)1 is negative .\n",nega_prbs.size());
            printf("WARNING: That means we can't get SE by doing square-root.\n");
            printf("WARNING: In term of such case we set effect size as 0 and SE as missing (-9).\n");
            string filename=string(outFileName)+".negative.probe.list";
            FILE* tmpfile=fopen(filename.c_str(),"w");
            if(!tmpfile)
            {
                printf("error open file.\n");
                exit(EXIT_FAILURE);
            }
            for(int t=0;t<nega_prbs.size();t++)
            {
                string str=nega_prbs[t]+'\n';
                fputs(str.c_str(),tmpfile);
            }
            fclose(tmpfile);
            printf("These probes are saved in file %s.\n",filename.c_str());
        }
        
        
        for(int i=0;i<probeinfo.size();i++)
        {
            if(probeinfo[i].genename) free2(&probeinfo[i].genename);
            if(probeinfo[i].probeId) free2(&probeinfo[i].probeId);
            if(probeinfo[i].ptr) free2(&probeinfo[i].ptr);
            if(probeinfo[i].bfilepath) free2(&probeinfo[i].bfilepath);
            if(probeinfo[i].esdpath) free2(&probeinfo[i].esdpath);
        }
        for(int i=0;i<snpinfo.size();i++)
        {
            if(snpinfo[i].a1) free2(&snpinfo[i].a1);
            if(snpinfo[i].a2) free2(&snpinfo[i].a2);
            if(snpinfo[i].snprs) free2(&snpinfo[i].snprs);
            if(snpinfo[i].rstr) free2(&snpinfo[i].rstr);
            if(snpinfo[i].revs) free2(&snpinfo[i].revs);
        }
        for(int i=0;i<besds.size();i++)
        {
            fclose(fptrs[i]);
        }
        fclose(efile);
        free(fptrs);
        free(buffer_se);
        free(buffer_beta);
        
        printf("\nThe eQTL infomation of %ld probes and %ld SNPs has been in binary file %s.\n",metaPrbNum,metaSNPnum,besdName.c_str());
    }
    
    void read_epi4u(eqtlInfo* eqtlinfo, string epifile)
    {
        ifstream epi(epifile.c_str());
        if (!epi) throw ("ERROR: can not open the file [" + epifile + "] to read.");
        cout << "Reading eQTL probe information from [" + epifile + "]." << endl;
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        while(!epi.eof())
        {
            epi.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        eqtlinfo->_probNum=lineNum;
        cout << eqtlinfo->_probNum << " Probes to be included from [" + epifile + "]." << endl;
        
        eqtlinfo->_epi_chr.resize(lineNum);
        eqtlinfo->_epi_prbID.resize(lineNum);
        eqtlinfo->_epi_gd.resize(lineNum);
        eqtlinfo->_epi_bp.resize(lineNum);
        eqtlinfo->_epi_gene.resize(lineNum);
        eqtlinfo->_epi_orien.resize(lineNum);
        eqtlinfo->_include.resize(lineNum);
        epi.clear(ios::goodbit);
        epi.seekg (0, ios::beg);
        for(int i=0;i<lineNum;i++)
        {
            string tmpStr;
            epi.getline(buf,MAX_LINE_SIZE);
            istringstream iss(buf);
            iss>>tmpStr;
            if(tmpStr=="NA" || tmpStr=="na" || tmpStr=="-9") {
                printf("chromosome is \"NA\" in row %d.\n", i+1);
                eqtlinfo->_epi_chr[i]=-9;
            } else {
                int tmpchr;
                if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                else tmpchr=atoi(tmpStr.c_str());
                eqtlinfo->_epi_chr[i]=tmpchr;
            }
            iss>>tmpStr;
            eqtlinfo->_include[i]=i;
            if(eqtlinfo->_probe_name_map.find(tmpStr) != eqtlinfo->_probe_name_map.end()){
                cout << "Warning: Duplicated probe ID \"" + tmpStr + "\" ";
                stringstream ss;
                ss << tmpStr << "_" << i + 1;
                tmpStr = ss.str();
                cout<<"has been changed to \"" + tmpStr + "\".\n";
            }
            eqtlinfo->_probe_name_map.insert(pair<string, int>(tmpStr, i));
            
            
            eqtlinfo->_epi_prbID[i]=tmpStr;
            iss>>tmpStr;
            eqtlinfo->_epi_gd[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            if(tmpStr=="NA" || tmpStr=="na" || tmpStr=="0") {
                printf("probe BP is \"NA\" in row %d.\n", i+1);
                eqtlinfo->_epi_bp[i]=-9;
            }
            else {
                eqtlinfo->_epi_bp[i]=atoi(tmpStr.c_str());
            }
            iss>>tmpStr;
            eqtlinfo->_epi_gene[i]=tmpStr.c_str();
            iss>>tmpStr;
            eqtlinfo->_epi_orien[i]=tmpStr.c_str()[0];
        }
        
        epi.close();
    }    
    void read_esi4u(eqtlInfo* eqtlinfo, string esifile)
    {
        ifstream esi(esifile.c_str());
        if (!esi) throw ("ERROR: can not open the file [" + esifile + "] to read.");
        cout << "Reading eQTL SNP information from [" + esifile + "]." << endl;
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        eqtlinfo->_esi_freq.clear();
        
        char buf[MAX_LINE_SIZE];
        vector<string> vs_buf;
        int lineNum(0);
        bool ptrnullfrq=false;
        while(!esi.eof())
        {
            esi.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=6 && col_num!=7) {
                    printf("ERROR: the number of columns is incorrect in row %d!\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na" || vs_buf[0]=="0"){
                    printf(" chromosome is \"NA\" in row %d.\n",lineNum+1);
                    eqtlinfo->_esi_chr.push_back(-9);
                } else {
                    int tmpchr;
                    if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                    else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                    else tmpchr=atoi(vs_buf[0].c_str());
                    eqtlinfo->_esi_chr.push_back(tmpchr);
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(eqtlinfo->_snp_name_map.find(vs_buf[1]) != eqtlinfo->_snp_name_map.end())
                {
                    cout << "WARNING: Duplicated SNP ID \"" + vs_buf[1] + "\" ";
                    stringstream ss;
                    ss << vs_buf[1] << "_" << lineNum + 1;
                    vs_buf[1] = ss.str();
                    cout<<"has been changed to \"" + vs_buf[1] + "\".\n";
                }
                eqtlinfo->_snp_name_map.insert(pair<string, int>(vs_buf[1], lineNum));
                eqtlinfo->_esi_rs.push_back(vs_buf[1]);
                eqtlinfo->_esi_gd.push_back(atoi(vs_buf[2].c_str()));
                if(vs_buf[3]=="NA" || vs_buf[3]=="na"){
                    printf("SNP BP is \'NA\' in row %d.\n", lineNum+1);
                    eqtlinfo->_esi_bp.push_back(-9);
                } else {
                    eqtlinfo->_esi_bp.push_back(atoi(vs_buf[3].c_str()));
                }
                
                if(vs_buf[4]=="NA" || vs_buf[4]=="na") printf("WARNING: allele1 is \"NA\" in row %d.\n", lineNum+1);
                to_upper(vs_buf[4]);
                eqtlinfo->_esi_allele1.push_back(vs_buf[4]);
                if(vs_buf[5]=="NA" || vs_buf[5]=="na") printf("WARNING: allele2 is \"NA\" in row %d.\n", lineNum+1);
                to_upper(vs_buf[5]);
                eqtlinfo->_esi_allele2.push_back(vs_buf[5]);
                if(col_num==7)
                {
                    if(vs_buf[6]=="NA" || vs_buf[6]=="na"){
                        if(!ptrnullfrq){
                            printf("WARNING: frequency is \"NA\" in one or more rows.\n");
                            ptrnullfrq=true;
                        }
                        eqtlinfo->_esi_freq.push_back(-9);
                    } else {
                        eqtlinfo->_esi_freq.push_back(atof(vs_buf[6].c_str()));
                    }
                } else {
                    eqtlinfo->_esi_freq.push_back(-9);
                }
                eqtlinfo->_esi_include.push_back(lineNum);
                lineNum++;
            }
        }
        eqtlinfo->_snpNum=lineNum;
        cout << eqtlinfo->_snpNum << " SNPs to be included from [" + esifile + "]." << endl;
        esi.close();
    }
    void update_epifile(char* eqtlFileName,char* s_epiFileName)
    {
        eqtlInfo eqtlinfo;
        if(eqtlFileName==NULL)
        {
            printf("Error: please input the BESD file by the option --beqtl-summary.\n");
            exit(EXIT_FAILURE);
        }
        if(s_epiFileName==NULL)
        {
            printf("Error: please input the probe information file by the option --update-epi.\n");
            exit(EXIT_FAILURE);
        }
        string fname=string(eqtlFileName)+".epi";
        read_epi4u(&eqtlinfo,fname);
        fname=string(eqtlFileName)+".bak";
        write_epi(fname, &eqtlinfo);
        
        FILE* epifile=fopen(s_epiFileName,"r");
        if (!(epifile)) {
            printf("Open error %s\n", s_epiFileName);
            exit(EXIT_FAILURE);
        }
        char Tbuf[MAX_LINE_SIZE];
        map<string, int>::iterator iter;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int hit =0;
        while(fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(strlist.size()>6)
            {
                printf("WARNING: line %u has more than 6 items.\n", line_idx);
            }
            if(strlist.size()<6)
            {
                printf("ERROR: line %u has less than 6 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            iter=eqtlinfo._probe_name_map.find(strlist[1]);
            if(iter!=eqtlinfo._probe_name_map.end())
            {
                int probeidx=iter->second;
                if(strlist[0]=="X" || strlist[0]=="x") eqtlinfo._epi_chr[probeidx]=23;
                else if(strlist[0]=="Y" || strlist[0]=="y") eqtlinfo._epi_chr[probeidx]=24;
                else if (atoi(strlist[0].c_str())==0 ) {
                    printf("ERROR: unrecongized chromomose found:\n");
                    printf("%s\n",Tbuf);
                    printf("ERROR: please use a different number to keep this probe:\n");
                    exit(EXIT_FAILURE);
                } else if (atoi(strlist[0].c_str())>24) {
                    printf("WARNING: abnormal chromomose found:\n");
                    printf("%s\n",Tbuf);
                    eqtlinfo._epi_chr[probeidx]=atoi(strlist[0].c_str());
                } else eqtlinfo._epi_chr[probeidx]=atoi(strlist[0].c_str());
                if(strlist[3]=="NA" | strlist[3]=="na") {
                    printf("ERROR: probe BP is missing:\n");
                    printf("%s\n",Tbuf);
                    exit(EXIT_FAILURE);
                }
                eqtlinfo._epi_bp[probeidx]=atoi(strlist[3].c_str());
                if(strlist[4]=="NA" | strlist[4]=="na") {
                    printf("WARNING: Gene id is missing:\n");
                    printf("%s\n",Tbuf);
                }
                eqtlinfo._epi_gene[probeidx]=strlist[4].c_str();
                if(strlist[5]=="NA" | strlist[5]=="na") {
                    printf("WARNING: Gene strand is missing:\n");
                    printf("%s\n",Tbuf);
                    eqtlinfo._epi_orien[probeidx]='*';
                } else eqtlinfo._epi_orien[probeidx]=strlist[5][0];
                hit++;
            }
            
            line_idx++;
        }
        fclose(epifile);
        write_epi(string(eqtlFileName),&eqtlinfo);
        printf("%d of %llu probes are updated.\n",hit, eqtlinfo._probNum);
    }
    void update_esifile(char* eqtlFileName,char* s_esiFileName)
    {
        eqtlInfo eqtlinfo;
        if(eqtlFileName==NULL)
        {
            printf("Error: please input the BESD file by the option --beqtl-summary.\n");
            exit(EXIT_FAILURE);
        }
        if(s_esiFileName==NULL)
        {
            printf("Error: please input the probe information file by the option --update-epi.\n");
            exit(EXIT_FAILURE);
        }
        string fname=string(eqtlFileName)+".esi";
        read_esi4u(&eqtlinfo,fname);
        fname=string(eqtlFileName)+".bak";
        write_esi(fname, &eqtlinfo);
        
        FILE* esifile=fopen(s_esiFileName,"r");
        if (!(esifile)) {
            printf("Open error %s\n", s_esiFileName);
            exit(EXIT_FAILURE);
        }
        char Tbuf[MAX_LINE_SIZE];
        map<string, int>::iterator iter;
        vector<string> strlist;
        uint32_t line_idx = 0;
        int hit =0;
        while(fgets(Tbuf, MAX_LINE_SIZE, esifile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(strlist.size()>7)
            {
                printf("WARNING: line %u has more than 7 items.\n", line_idx);
            }
            if(strlist.size()<6)
            {
                printf("ERROR: line %u has less than 6 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            iter=eqtlinfo._snp_name_map.find(strlist[1]);
            if(iter!=eqtlinfo._snp_name_map.end())
            {
                int snpidx=iter->second;
                if(strlist[0]=="X" || strlist[0]=="x") eqtlinfo._esi_chr[snpidx]=23;
                else if(strlist[0]=="Y" || strlist[0]=="y") eqtlinfo._esi_chr[snpidx]=24;
                else if (atoi(strlist[0].c_str())==0 ) {
                    printf("ERROR: unrecongized chromomose found:\n");
                    printf("%s\n",Tbuf);
                    printf("ERROR: please use a different number to keep this SNP:\n");
                    exit(EXIT_FAILURE);
                } else if (atoi(strlist[0].c_str())>24) {
                    printf("WARNING: abnormal chromomose found:\n");
                    printf("%s\n",Tbuf);
                    eqtlinfo._esi_chr[snpidx]=atoi(strlist[0].c_str());
                } else eqtlinfo._esi_chr[snpidx]=atoi(strlist[0].c_str());
                if(strlist[3]=="NA" | strlist[3]=="na") {
                    printf("ERROR: SNP BP is missing:\n");
                    printf("%s\n",Tbuf);
                    exit(EXIT_FAILURE);
                }
                eqtlinfo._esi_bp[snpidx]=atoi(strlist[3].c_str());
                if(strlist[4]=="NA" | strlist[4]=="na") {
                    printf("WARNING: the effect allele is missing:\n");
                    printf("%s\n",Tbuf);
                }
                eqtlinfo._esi_allele1[snpidx]=strlist[4].c_str();
                if(strlist[5]=="NA" | strlist[5]=="na") {
                    printf("WARNING: the other allele is missing:\n");
                    printf("%s\n",Tbuf);
                }
                eqtlinfo._esi_allele2[snpidx]=strlist[5].c_str();
                if(strlist.size()==7)
                {
                    if(strlist[6]=="NA" || strlist[6]=="na"){
                        printf("WARNING: frequency is \"NA\" in one or more rows.\n");
                    } else {
                        eqtlinfo._esi_freq[snpidx]=atof(strlist[6].c_str());
                    }
                }
                hit++;
            }
            line_idx++;
        }
        fclose(esifile);
        write_esi(string(eqtlFileName),&eqtlinfo);
        printf("%d of %llu SNPs are updated.\n",hit, eqtlinfo._snpNum);
    }

    
}
