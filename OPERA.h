//
//  OPERA.h
//
//  Created by Yang Wu on 06/08/22.
//  Copyright (c) 2022 Yang Wu. All rights reserved.
//

#ifndef __OPERA_CPP__OPERA__
#define __OPERA_CPP__OPERA__

#include "SMR_data.h"

namespace SMRDATA
{
    typedef struct {
        vector<int> _chr;
        vector<string> _snp_name;
        vector<int> _bp;
        vector<int> _include;
    } lociData;
    typedef struct{
        int cur_chr;
        int cur_prbidx;
        vector<double> bxz, sexz,freq,zxz;
        vector<vector<double>> byz,seyz,pyz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
    } MTSMRWK;
    typedef struct{
        int cur_chr;
        int cur_prbidx;
        vector<uint32_t> splSize;
        vector<vector<double>> bxz,sexz,freq,zxz;
        vector<double> byz,seyz,pyz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
    } MTSMRWKEXP;

    // OPERA Yang wu 06/08/2022 //
    void ppafile_summary(char* ppaFilename,char* GWAScojosnplstName,char* numFilename,char* piFilename, string priorstr, char* outFileName, double thresh_PP, double thresh_heidi); // take the mean of combinations when marginal ppa > 0.9;
    void ppafile_summary_old1(char* ppaFilename, char* numFilename,char* piFilename,string priorstr,char* outFileName, double thresh_PP, double thresh_heidi); // take the mean of combinations for all marginal ppa;
    void multiexposurepi_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* eqtlsmaslstName, char* gwasFileName, double maf, string sigmastr, double sigma_def, double alpha_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,bool jointsmrflag,int cis_itvl,int piWind,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor,char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh);
    void multiexposure_jointsmr_old(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* rhoFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr,double sigma_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP,double thresh_PP_out,double thresh_smr,double thresh_gwas,double thresh_heidi,char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, bool printcombppaflag, bool printsmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh);
    void multiexposure_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* rhoFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr,double sigma_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP,double thresh_PP_out,double thresh_smr,double thresh_gwas,double thresh_heidi,char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, bool printcombppaflag, bool printsmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh);
    void multioutcomesmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* eqtlFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor,char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh);
    void read_GWAS_cojo_snplist(lociData* ldata, char* lociFileName);
    void read_pifile(string piFilename, vector<string> &priorsplit);
    void read_varfile(string varFilename, vector<string> &sigmasplit);
    void read_rhofile(string rhoFilename, MatrixXd &rho);
    void read_numfile(string numFilename, vector<long> &numsum);
    void cis_xqtl_probe_include_only(eqtlInfo* eqtlinfo, double p_smr, int cis_itvl, string eqtlFileName);
    void read_prb_cojo_snplist(string snpprblistfile, vector<string> &prblist,  map<string, vector<string>> &prb_snp);
    void extract_ldata_by_chr(lociData* lociData, int prbchr);
    void extract_prob_opera(eqtlInfo* eqtlinfo,string problstName);
    void update_ldata(lociData* ldata);
    void e2gconvert(eqtlInfo* etrait, gwasData* gdata, int &ii);
    void e2econvert(eqtlInfo* etrait, eqtlInfo* esdata);
    void e2econvert_old(eqtlInfo* etrait, eqtlInfo* esdata); // Not used; failed to improve
    void e2econvert_old2(eqtlInfo* etrait, eqtlInfo* esdata); // Not used; failed to improve
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata);
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, eqtlInfo* esdata);
    void allele_check_multi(vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata);
    void allele_check_multi(vector<eqtlInfo> &etrait, eqtlInfo* esdata);
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);
    void allele_check_multi_opt_old(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);
    void allele_check_multi_opt_old2(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);
    void allele_check_multi_opt_old3(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);//relative new version
    void allele_check_multi_opt(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);// parallel jobs,least switch
    void allele_check_multi_opt_test(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);
    void allele_check_multi(vector<eqtlInfo> &etrait, gwasData* gdata);
    void allele_compare(string a1, string a2, string s1,string s2, int &Id, int &flip);    
    void smr_heidi_func_para(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_het_mtd, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg);
    void multi_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, vector<gwasData> &gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_heidi_func_so(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_joint_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, vector<vector<string>> &prbcojolist, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void run_joint_effect_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy);
    void run_joint_effect_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp);
    void run_joint_effect_eigenVar_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp);
    void run_joint_effect_shrink_NoLDprune_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp);
    void multi_joint_smr_func_old(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_joint_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_joint_smr_func_v2(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, eqtlInfo &esdata, double ngwas, vector<string> &cojolist, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void multi_cond_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor);
    void init_smr_wk_mlt(MTSMRWK* smrwk);
    void init_smr_wk_mlt(MTSMRWKEXP* smrwk);
    void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
    long fill_smr_wk_new(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_include(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt(bInfo* bdata,vector<gwasData> &gdata,eqtlInfo* esdata,MTSMRWK* smrwk, const char* refSNP, int lowerbp,int upperbp,bool heidioffFlag);
    long fill_smr_wk_mlt(bInfo* bdata,vector<gwasData> &gdata,eqtlInfo* esdata,MTSMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt_old(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt_old2(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt_old3(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag);
    long fill_smr_wk_mlt_include(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag);
    double multi_heidi_test_new(bInfo* bdata,MTSMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, vector<double> theta);
    double multi_heidi_test_new(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, vector<double> theta);
    double multi_heidi_test_new_v2(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, vector<double> theta);
    double multi_heidi_test_new_v2_so(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, double ld_min, int opt_hetero, bool sampleoverlap, MatrixXd theta);
    void multi_heidi_pruning(bInfo* bdata,MTSMRWKEXP* smrwk,int t,double ldr2_top, double threshold, int m_hetero ,double ld_min,int opt_hetero,vector<string> &selectSNPs);
    void multi_heidi_pruning_v2(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero ,double ld_min,int opt_hetero,vector<string> &slctsnps);
    float bxy_mltheter(vector<VectorXd> &_byz, VectorXd &_bxz, vector<VectorXd> &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num, vector<double> theta);
    float bxy_mltheter(VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz, vector<VectorXd> &_sexz, vector<VectorXd> &_zsxz, MatrixXd &_LD_heidi, long* snp_num, vector<double> theta);
    float bxy_mltheter_v2_so(VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz, vector<VectorXd> &_sexz, vector<VectorXd> &_zsxz, MatrixXd &_LD_heidi, long* snp_num, MatrixXd theta);
    void est_cov_dxy(MatrixXd &covbxy,VectorXd &vardev,int maxid,vector<VectorXd> &_byz, VectorXd &_bxz, vector<VectorXd> &_seyz,VectorXd &_sexz, MatrixXd &_LD_heidi, vector<double> theta);
    void est_cov_dxy(MatrixXd &covdxy,VectorXd &vardev,vector<int> maxid, VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz,vector<VectorXd> &_sexz, MatrixXd &_LD_heidi, vector<double> theta);
    void est_cov_dxy_so(MatrixXd &covdxy,VectorXd &vardev,vector<int> maxid, VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz,vector<VectorXd> &_sexz, MatrixXd &_LD_heidi, MatrixXd theta);
    void combn_marg_pi1(long expoNum, vector<vector<int>> &combins, vector<vector<int>> &combmarg, vector<vector<int>> &idxmarg, vector<double> &pi1, vector<double> prior);
    void stage2_output_file_format(string &smrfile, FILE* &smr, string &resfile, FILE* &res, string &propfile, FILE* &prop, vector<string> &ppafile, vector<FILE*> &ppa, char* gwasFileName, char* GWAScojosnplstName, bool printcombppaflag, bool printsmrflag, vector<vector<int>> &combmarg);

    // rename function to aviod conflict with SMR_data.cpp
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata);
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, eqtlInfo* esdata);
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata);
    void update_snidx_opera(MTSMRWK* smrwk, vector<int> &sn_ids,int max_snp_slct, string forwhat);
    void update_snidx_opera(MTSMRWKEXP* smrwk, vector<int> &sn_ids,int max_snp_slct, string forwhat);
    void update_snidx_opera(MTSMRWKEXP* smrwk, int t, vector<int> &sn_ids,int max_snp_slct, string forwhat);
    void extract_smrwk_opera_new(SMRWK* smrwk,vector<int> &sn_ids,SMRWK* smrwk2);
    void extract_smrwk_opera(MTSMRWK* smrwk,vector<int> &sn_ids,MTSMRWK* smrwk2);
    void extract_smrwk_opera(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MTSMRWKEXP* smrwk2);
    void update_smrwk_x_opera(MTSMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X);
    void update_smrwk_x_opera(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MatrixXd &X);    
    int max_zsmr_id_opera(MTSMRWKEXP *smrwk , double p_smr);
    void extract_snp_by_chr(bInfo* bdata, int chr);
    void exclude_eqtl_snp_opera(eqtlInfo* eqtlinfo, string snplstName);
    void smr_heidi_func_opera(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_het_mtd, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg);
    
    // read mbfile //
    void read_multi_bimfiles(bInfo* bdata, vector<string> multi_bfiles, map<string, string> &snp_name_per_chr);
    void read_multi_famfiles(bInfo* bdata, vector<string> multi_bfiles);
    void read_multi_bedfiles(bInfo* bdata, vector<string> multi_bfiles, map<string, string> &snp_name_per_chr);
    void read_single_bimfile(bInfo* bdata, string bimfile, bool msg_flag);
    void read_single_famfile(bInfo* bdata, string famfile, bool msg_flag);
    void read_single_bedfile(bInfo* bdata, string bedfile, vector<pair<int,int>> rsnp, vector<int> rindi, bool msg_flag);
    void update_include(bInfo* bdata, bInfo* bdatatmp, int file_indx, map<string, string> &snp_name_per_chr);
    void update_keep(bInfo* bdata, bInfo* bdatatmp, string famfile);
    void update_fam(bInfo* bdata,vector<int> &rindi);
    void update_bim(bInfo* bdata,vector<int> &rsnp);
    void update_id_chr_map(map<string, string> &chr_map, map<string, int> id_map);
    void retrieve_snp(map<string,string> snp_chr_map, map<string,int> snp_id_map, vector<vector<pair<int,int>>> &rsnp, int nbfiles);    
    // end //
}
#endif /* defined(__SRM_CPP__SMR_data__) */