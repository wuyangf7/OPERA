//  SMR_data.cpp
//  SRM_CPP
//
//  Created by Yang Wu on 04/06/18.
//  Copyright (c) 2021 Yang Wu. All rights reserved.
//

#include "SMR_data.h"
#include "stat.hpp"
#include "SMR_data_p1.h"

namespace SMRDATA
{
   
    int file_read_check(ifstream* in_file, const char* filename)
    {
        in_file->open(filename);
        if(!*in_file)  return 0;
        return 1;
    }
    void progress_print(float progress)
    {
        int barWidth = 70;
        cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();    
    }

    void read_famfile(bInfo* bdata, string famfile) {
		bdata->_autosome_num = 22;
        ifstream Fam(famfile.c_str());
        if (!Fam) throw ("Error: can not open the file [" + famfile + "] to read.");
        cout << "Reading PLINK FAM file from [" + famfile + "]." << endl;
        
        int i = 0;
        string str_buf;
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        while (Fam) {
            Fam >> str_buf;
            if (Fam.eof()) break;
            bdata->_fid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_pid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_fa_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_mo_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_sex.push_back(atoi(str_buf.c_str()));
            Fam >> str_buf;
            bdata->_pheno.push_back(atoi(str_buf.c_str()));
        }
        Fam.clear();
        Fam.close();
        bdata->_indi_num = bdata->_fid.size();
        cout << bdata->_indi_num << " individuals to be included from [" + famfile + "]." << endl;
        
        // Initialize _keep
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        int size = 0;
        for (int i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
            if (size == bdata->_id_map.size()) throw ("Error: Duplicate individual ID found: \"" + bdata->_fid[i] + "\t" + bdata->_pid[i] + "\".");
            size = bdata->_id_map.size();
        }
    }
    
    void update_bim(bInfo* bdata,vector<int> &rsnp) {
        int i = 0;
        
        //update bim information
        vector<int> chr_buf, bp_buf;
        vector<string> a1_buf, a2_buf, ref_A_buf, other_A_buf;
        vector<string> snp_name_buf;
        vector<double> genet_dst_buf, impRsq_buf;
        for (i = 0; i < bdata->_snp_num; i++) {
            if (!rsnp[i]) continue;
            chr_buf.push_back(bdata->_chr[i]);
            snp_name_buf.push_back(bdata->_snp_name[i]);
            genet_dst_buf.push_back(bdata->_genet_dst[i]);
            bp_buf.push_back(bdata->_bp[i]);
            a1_buf.push_back(bdata->_allele1[i]);
            a2_buf.push_back(bdata->_allele2[i]);
            ref_A_buf.push_back(bdata->_ref_A[i]);
            other_A_buf.push_back(bdata->_other_A[i]);
            if(bdata->_impRsq.size()>0) impRsq_buf.push_back(bdata->_impRsq[i]);
        }
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        bdata->_ref_A.clear();
        bdata->_other_A.clear();
        bdata->_impRsq.clear();
        bdata->_chr = chr_buf;
        bdata->_snp_name = snp_name_buf;
        bdata->_genet_dst = genet_dst_buf;
        bdata->_bp = bp_buf;
        bdata->_allele1 = a1_buf;
        bdata->_allele2 = a2_buf;
        bdata->_ref_A = ref_A_buf;
        bdata->_other_A = other_A_buf;
        bdata->_impRsq=impRsq_buf;
        bdata->_snp_num = bdata->_chr.size();
        bdata->_include.clear();
        bdata-> _include.resize(bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (i = 0; i < bdata->_snp_num; i++) {
            bdata->_include[i] = i;
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }
    
    void update_fam(bInfo* bdata,vector<int> &rindi) {
        //update fam information
        int i = 0;
        vector<string> fid_buf, pid_buf, fa_id_buf, mo_id_buf;
        vector<int> sex_buf;
        vector<double> pheno_buf;
        for (i = 0; i < bdata->_indi_num; i++) {
            if (!rindi[i]) continue;
            fid_buf.push_back(bdata->_fid[i]);
            pid_buf.push_back(bdata->_pid[i]);
            fa_id_buf.push_back(bdata->_fa_id[i]);
            mo_id_buf.push_back(bdata->_mo_id[i]);
            sex_buf.push_back(bdata->_sex[i]);
            pheno_buf.push_back(bdata->_pheno[i]);
        }
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        bdata->_fid = fid_buf;
        bdata->_pid = pid_buf;
        bdata->_fa_id = fa_id_buf;
        bdata->_mo_id = mo_id_buf;
        bdata->_sex = sex_buf;
        bdata->_pheno = pheno_buf;
        
        bdata->_indi_num = bdata->_fid.size();
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        for (i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
        }
    }
    
	void update_epi(eqtlInfo* eqtlinfo)
	{
		eqtlinfo->_probNum = eqtlinfo->_include.size();

        vector<int> chr_buf, gd_buf,bp_buf, start, end;
		vector<string> prbID_buf, gene_buf;
		vector<char> orien_buf;
		for (int i = 0; i < eqtlinfo->_probNum; i++)
		{
			chr_buf.push_back(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]);
			gd_buf.push_back(eqtlinfo->_epi_gd[eqtlinfo->_include[i]]);
			bp_buf.push_back(eqtlinfo->_epi_bp[eqtlinfo->_include[i]]);
			prbID_buf.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
			gene_buf.push_back(eqtlinfo->_epi_gene[eqtlinfo->_include[i]]);
			orien_buf.push_back(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]);
            if(eqtlinfo->_epi_start.size()>0) start.push_back(eqtlinfo->_epi_start[eqtlinfo->_include[i]]);
            if(eqtlinfo->_epi_end.size()>0) start.push_back(eqtlinfo->_epi_end[eqtlinfo->_include[i]]);
		}
		eqtlinfo->_epi_chr.clear();
		eqtlinfo->_epi_gd.clear();
		eqtlinfo->_epi_bp.clear();
		eqtlinfo->_epi_prbID.clear();
		eqtlinfo->_epi_gene.clear();
		eqtlinfo->_epi_orien.clear();
        eqtlinfo->_epi_start.clear();
        eqtlinfo->_epi_end.clear();
		eqtlinfo->_epi_chr.swap(chr_buf);
		eqtlinfo->_epi_gd.swap(gd_buf);
		eqtlinfo->_epi_bp.swap(bp_buf);
		eqtlinfo->_epi_prbID.swap(prbID_buf);
		eqtlinfo->_epi_gene.swap(gene_buf);
		eqtlinfo->_epi_orien.swap(orien_buf);
        eqtlinfo->_epi_start.swap(start);
        eqtlinfo->_epi_end.swap(end);
        
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        for (int i = 0; i < eqtlinfo->_probNum; i++)
        {
            eqtlinfo->_include.push_back(i);
            eqtlinfo->_probe_name_map.insert(pair<string,int>(eqtlinfo->_epi_prbID[i],i));
        }
	}

	void update_esi(eqtlInfo* eqtlinfo)
	{
		eqtlinfo->_snpNum = eqtlinfo->_esi_include.size();	

		vector<int> chr_buf, gd_buf, bp_buf;
		vector<string> rs_buf;
		vector<string> allele1_buf, allele2_buf;
        vector<float> freq_buf;
		for (int i = 0; i < eqtlinfo->_snpNum; i++)
		{
			chr_buf.push_back(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
			gd_buf.push_back(eqtlinfo->_esi_gd[eqtlinfo->_esi_include[i]]);
			bp_buf.push_back(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]]);
			rs_buf.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
			allele1_buf.push_back(eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]);
			allele2_buf.push_back(eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]);
            freq_buf.push_back(eqtlinfo->_esi_freq[eqtlinfo->_esi_include[i]]);
		}
		eqtlinfo->_esi_chr.clear();
		eqtlinfo->_esi_gd.clear();
		eqtlinfo->_esi_bp.clear();
		eqtlinfo->_esi_rs.clear();
		eqtlinfo->_esi_allele1.clear();
		eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_freq.clear();
		eqtlinfo->_esi_chr.swap(chr_buf);
		eqtlinfo->_esi_gd.swap(gd_buf);
		eqtlinfo->_esi_bp.swap(bp_buf);
		eqtlinfo->_esi_rs.swap(rs_buf);
		eqtlinfo->_esi_allele1.swap(allele1_buf);
		eqtlinfo->_esi_allele2.swap(allele2_buf);
        eqtlinfo->_esi_freq.swap(freq_buf);
        
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        for (int i = 0; i < eqtlinfo->_snpNum; i++)
        {
            eqtlinfo->_esi_include.push_back(i);
            eqtlinfo->_snp_name_map.insert(pair<string, int>(eqtlinfo->_esi_rs[i], i));
        }
	}

    void read_bimfile(bInfo* bdata,string bimfile) {
        // Read bim file: recombination rate is defined between SNP i and SNP i-1
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        ifstream Bim(bimfile.c_str());
        if (!Bim) throw ("Error: can not open the file [" + bimfile + "] to read.");
        cout << "Reading PLINK BIM file from [" + bimfile + "]." << endl;
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        while (Bim) {
            Bim >> ibuf;
            if (Bim.eof()) break;
            bdata->_chr.push_back(ibuf);
            Bim >> str_buf;
            bdata->_snp_name.push_back(str_buf);
            Bim >> dbuf;
            bdata->_genet_dst.push_back(dbuf);
            Bim >> ibuf;
            bdata->_bp.push_back(ibuf);
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele1.push_back(cbuf.c_str());
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str());
        }
        Bim.close();
        bdata->_snp_num = bdata->_chr.size();
        bdata->_ref_A = bdata->_allele1;
        bdata->_other_A = bdata->_allele2;
        cout << bdata->_snp_num << " SNPs to be included from [" + bimfile + "]." << endl;
        
        // Initialize _include
         bdata->_include.clear();
         bdata->_include.resize( bdata->_snp_num);
         bdata->_snp_name_map.clear();
        
        for (int i = 0; i <  bdata->_snp_num; i++) {
            
            bdata->_include[i] = i;
            if( bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()){
                cout << "Warning: Duplicated SNP ID \"" + bdata->_snp_name[i] + "\" ";
                stringstream ss;
                ss << bdata->_snp_name[i] << "_" << i + 1;
                bdata->_snp_name[i] = ss.str();
                cout<<"has been changed to \"" + bdata->_snp_name[i] + "\".\n";
            }
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }
    // some code are adopted from PLINK with modificationsm
    void read_bedfile(bInfo* bdata, string bedfile)
    {
        int i = 0, j = 0, k = 0;
        
        // Flag for reading individuals and SNPs
        vector<int> rindi, rsnp;
        //get_rindi
        rindi.clear();
        rindi.resize(bdata->_indi_num);
        for (int i = 0; i < bdata->_indi_num; i++) {
            if (bdata->_id_map.find(bdata->_fid[i] + ":" + bdata->_pid[i]) != bdata->_id_map.end()) rindi[i] = 1;
            else rindi[i] = 0;
        }
        //get_rsnp
        rsnp.clear();
        rsnp.resize(bdata->_snp_num);
        for (int i = 0; i < bdata->_snp_num; i++) {
            if (bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()) rsnp[i] = 1;
            else rsnp[i] = 0;
        }
        
        if (bdata->_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        if (bdata->_keep.size() == 0) throw ("Error: No individual is retained for analysis.");
        
        // Read bed file
        char ch[1];
        bitset<8> b;
        bdata->_snp_1.resize(bdata->_include.size());
        bdata->_snp_2.resize(bdata->_include.size());
        for (i = 0; i < bdata->_include.size(); i++) {
			bdata->_snp_1[i].reserve(bdata->_keep.size());
			bdata->_snp_2[i].reserve(bdata->_keep.size());
        }
        fstream BIT(bedfile.c_str(), ios::in | ios::binary);
        if (!BIT) throw ("Error: can not open the file [" + bedfile + "] to read.");
        cout << "Reading PLINK BED file from [" + bedfile + "] in SNP-major format ..." << endl;
        for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
        int snp_indx = 0, indi_indx = 0;
        for (j = 0, snp_indx = 0; j < bdata->_snp_num; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
            if (!rsnp[j]) {
                for (i = 0; i < bdata->_indi_num; i += 4) BIT.read(ch, 1);
                continue;
            }
            for (i = 0, indi_indx = 0; i < bdata->_indi_num;) {
                BIT.read(ch, 1);
                if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
                b = ch[0];                
                k = 0;
                while (k < 7 && i < bdata->_indi_num) { // change code: 11 for AA; 00 for BB;
                    if (!rindi[i]) k += 2;
                    else {						
                        bdata->_snp_2[snp_indx][indi_indx] = (!b[k++]);
                        bdata->_snp_1[snp_indx][indi_indx] = (!b[k++]);
                        indi_indx++;
                    }
                    i++;
                }
            }
            if (snp_indx == bdata->_include.size()) break;
            snp_indx++;
        }
        BIT.clear();
        BIT.close();
        cout << "Genotype data for " << bdata->_keep.size() << " individuals and " << bdata->_include.size() << " SNPs to be included from [" + bedfile + "]." << endl;
        
        update_fam(bdata, rindi);
        update_bim(bdata, rsnp);
    }

    // read mbfile, modified from GCTA, Yang WU 29/11/2021
     void read_multi_bimfiles(bInfo* bdata, vector<string> multi_bfiles, map<string, string> &snp_name_per_chr)
     {
         cout << "Reading the PLINK BIM files ..." << endl;
         int i=0, nbfiles = multi_bfiles.size();
         string bimfile = "";

         bdata->_chr.clear();
         bdata->_snp_name.clear();
         bdata->_genet_dst.clear();
         bdata->_bp.clear();
         bdata->_allele1.clear();
         bdata->_allele2.clear();

         snp_name_per_chr.clear();
         for( i=0; i<nbfiles; i++ ) {
             bInfo bdatatmp;
             bimfile = multi_bfiles[i]+".bim";
             read_single_bimfile(&bdatatmp, bimfile, false);
             update_include(bdata, &bdatatmp, i, snp_name_per_chr);
         }

         // Initialize _include
         bdata->_snp_num = bdata->_snp_name.size();
         bdata->_include.clear(); bdata->_include.resize(bdata->_snp_num);
         bdata->_snp_name_map.clear();
         for (int i = 0; i <  bdata->_snp_num; i++) {
             bdata->_include[i] = i;
             bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
         }
         bdata->_ref_A = bdata->_allele1; bdata->_other_A = bdata->_allele2;
         printf("%ld SNPs to be included from PLINK BIM files.\n", bdata->_snp_num);
     }

     void update_include(bInfo* bdata, bInfo* bdatatmp, int file_indx, map<string, string> &snp_name_per_chr)
     {
         for (int i = 0; i < bdatatmp->_snp_num; i++) {
             if( snp_name_per_chr.find(bdatatmp->_snp_name[i]) != snp_name_per_chr.end()) {
                 cout << "Warning: Duplicated SNP ID \"" + bdatatmp->_snp_name[i] + "\" ";
                 stringstream ss;
                 ss << bdatatmp->_snp_name[i] << "_" << i + 1;
                 bdatatmp->_snp_name[i] = ss.str();
                 cout<<"has been changed to \"" + bdatatmp->_snp_name[i] + "\".\n";
             }
             snp_name_per_chr.insert(pair<string, string>(bdatatmp->_snp_name[i], to_string(file_indx)+":"+to_string(i)));
         }
         // Add new SNPs
         bdata->_chr.insert(bdata->_chr.end(), bdatatmp->_chr.begin(), bdatatmp->_chr.end());
         bdata->_snp_name.insert(bdata->_snp_name.end(), bdatatmp->_snp_name.begin(), bdatatmp->_snp_name.end());
         bdata->_genet_dst.insert(bdata->_genet_dst.end(), bdatatmp->_genet_dst.begin(), bdatatmp->_genet_dst.end());
         bdata->_bp.insert(bdata->_bp.end(), bdatatmp->_bp.begin(), bdatatmp->_bp.end());
         bdata->_allele1.insert(bdata->_allele1.end(), bdatatmp->_allele1.begin(), bdatatmp->_allele1.end());
         bdata->_allele2.insert(bdata->_allele2.end(), bdatatmp->_allele2.begin(), bdatatmp->_allele2.end());
     }

     void read_single_bimfile(bInfo* bdata,string bimfile, bool msg_flag)
     {
         // Read bim file: recombination rate is defined between SNP i and SNP i-1
         int ibuf = 0;
         string cbuf = "0";
         double dbuf = 0.0;
         string str_buf;
         ifstream Bim(bimfile.c_str());
         if (!Bim) throw ("Error: can not open the file [" + bimfile + "] to read.");
         if(msg_flag) cout << "Reading PLINK BIM file from [" + bimfile + "]." << endl;
         bdata->_chr.clear();
         bdata->_snp_name.clear();
         bdata->_genet_dst.clear();
         bdata->_bp.clear();
         bdata->_allele1.clear();
         bdata->_allele2.clear();
         while (Bim) {
             Bim >> ibuf;
             if (Bim.eof()) break;
             bdata->_chr.push_back(ibuf);
             Bim >> str_buf;
             bdata->_snp_name.push_back(str_buf);
             Bim >> dbuf;
             bdata->_genet_dst.push_back(dbuf);
             Bim >> ibuf;
             bdata->_bp.push_back(ibuf);
             Bim >> cbuf;
             StrFunc::to_upper(cbuf);
             bdata->_allele1.push_back(cbuf.c_str());
             Bim >> cbuf;
             StrFunc::to_upper(cbuf);
             bdata->_allele2.push_back(cbuf.c_str());
         }
         Bim.close();
         bdata->_snp_num = bdata->_chr.size();
         if(msg_flag) cout << bdata->_snp_num << " SNPs to be included from [" + bimfile + "]." << endl;
     }

     void read_multi_famfiles(bInfo* bdata, vector<string> multi_bfiles)
     {
         cout << "Reading the PLINK FAM files ..." << endl;
         int i=0, nbfiles = multi_bfiles.size();
         string famfile = "";

         bdata->_fid.clear();
         bdata->_pid.clear();
         bdata->_fa_id.clear();
         bdata->_mo_id.clear();
         bdata->_sex.clear();
         bdata->_pheno.clear();
         for( i=0; i<nbfiles; i++ ) {
             bInfo bdatatmp;
             famfile = multi_bfiles[i]+".fam";
             read_single_famfile(&bdatatmp, famfile, false);
             update_keep(bdata, &bdatatmp, famfile);
         }

         // Sample size
         bdata->_indi_num = bdata->_fid.size();
         printf("%ld individuals have been included from the PLINK FAM files.\n", bdata->_indi_num);
     }

     void update_keep(bInfo* bdata, bInfo* bdatatmp, string famfile)
     {
         int indx = 0; string indi_str = "";
         bdatatmp->_indi_num = bdatatmp->_fid.size();
         // Initial sample size
         if(bdata->_fid.size() == 0) {
             bdata->_fid = bdatatmp->_fid;
             bdata->_pid = bdatatmp->_pid;
             bdata->_fa_id = bdatatmp->_fa_id;
             bdata->_mo_id = bdatatmp->_mo_id;
             bdata->_sex = bdatatmp->_sex;
             bdata->_pheno = bdatatmp->_pheno;
             // Initialize the id map per chr
             int size=0;
             bdata->_keep.clear(); bdata->_id_map.clear();
             for(int i = 0; i < bdatatmp->_indi_num; i++) {
                 bdata->_keep.push_back(i);
                 bdata->_id_map.insert(pair<string,int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
                 if (size == bdata->_id_map.size()) throw ("Error: Duplicate individual ID found: \"" + bdata->_fid[i] + "\t" + bdata->_pid[i] + "\".");
                 size = bdata->_id_map.size();
             }
         } else {
             // Add new individuals
             for(int i = 0; i < bdatatmp->_indi_num; i++) {
                 // Search conflicted information of individuals
                 indi_str = bdatatmp->_fid[i] + ":" + bdatatmp->_pid[i];
                 map<string,int>::iterator iter = bdata->_id_map.find(indi_str);

                 if(iter!=bdata->_id_map.end()) {
                     // already existed
                     indx = iter->second;
                     if(bdatatmp->_fa_id[i] != bdata->_fa_id[indx]) throw("Inconsistent paternal ID found, " + bdatatmp->_fid[i] + " " + bdatatmp->_pid[i] + ", from [" + famfile + "].");
                     if(bdatatmp->_mo_id[i] != bdata->_mo_id[indx]) throw("Inconsistent maternal ID found, " + bdatatmp->_fid[i] + " " + bdatatmp->_pid[i] + ", from [" + famfile + "].");
                     if(bdatatmp->_sex[i] != bdata->_sex[indx]) throw("Inconsistent gender found, " + bdatatmp->_fid[i] + " " + bdatatmp->_pid[i] + ", from [" + famfile + "].");
                     if(bdatatmp->_pheno[i] != bdata->_pheno[indx]) throw("Inconsistent phenotype found, " + bdatatmp->_fid[i] + " " + bdatatmp->_pid[i] + ", from [" + famfile + "].");
                     if(i!=indx) throw("Inconsistent orders of individuals found from [" + famfile + "]. Please make sure that the orders of individuals are the same across the fam files.");
                 } else {
                     // not existed
                     throw("Unexpected individual ID found, " + bdatatmp->_fid[i] + " " + bdatatmp->_pid[i] + ", from [" + famfile + "].");
                 }
             }
         }
     }

     void read_single_famfile(bInfo* bdata, string famfile, bool msg_flag)
     {
         bdata->_autosome_num = 22;
         ifstream Fam(famfile.c_str());
         if (!Fam) throw ("Error: can not open the file [" + famfile + "] to read.");
         if(msg_flag) cout << "Reading PLINK FAM file from [" + famfile + "]." << endl;
        
         int i = 0;
         string str_buf;
         bdata->_fid.clear();
         bdata->_pid.clear();
         bdata->_fa_id.clear();
         bdata->_mo_id.clear();
         bdata->_sex.clear();
         bdata->_pheno.clear();
         while (Fam) {
             Fam >> str_buf;
             if (Fam.eof()) break;
             bdata->_fid.push_back(str_buf);
             Fam >> str_buf;
             bdata->_pid.push_back(str_buf);
             Fam >> str_buf;
             bdata->_fa_id.push_back(str_buf);
             Fam >> str_buf;
             bdata->_mo_id.push_back(str_buf);
             Fam >> str_buf;
             bdata->_sex.push_back(atoi(str_buf.c_str()));
             Fam >> str_buf;
             bdata->_pheno.push_back(atoi(str_buf.c_str()));
         }
         Fam.clear();
         Fam.close();
         if(msg_flag) {
             bdata->_indi_num = bdata->_fid.size();
             cout << bdata->_indi_num << " individuals to be included from [" + famfile + "]." << endl;
         }
     }

     void read_multi_bedfiles(bInfo* bdata, vector<string> multi_bfiles, map<string, string> &snp_name_per_chr)
     {
         int i=0, nbfiles = multi_bfiles.size();
         string bedfile = "";
         vector<vector<pair<int, int>>> rsnp;
         vector<int> rindi_flag, rsnp_flag;

         if (bdata->_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
         if (bdata->_keep.size() == 0) throw ("Error: No individual is retained for analysis.");

         cout << "Reading the PLINK BIM files ..." << endl;
         // Initialize the matrix
         bdata->_snp_1.resize(bdata->_include.size());
         bdata->_snp_2.resize(bdata->_include.size());
         for (i = 0; i < bdata->_include.size(); i++) {
             bdata->_snp_1[i].reserve(bdata->_keep.size());
             bdata->_snp_2[i].reserve(bdata->_keep.size());
         }

         // Update the map to retrieve individuals and SNPs
         update_id_chr_map(snp_name_per_chr, bdata->_snp_name_map);

         // Flag for reading individuals and SNPs
         //get_rindi_flag
         rindi_flag.clear();
         rindi_flag.resize(bdata->_indi_num);
         for (int i = 0; i < bdata->_indi_num; i++) {
             if (bdata->_id_map.find(bdata->_fid[i] + ":" + bdata->_pid[i]) != bdata->_id_map.end()) rindi_flag[i] = 1;
             else rindi_flag[i] = 0;
         }
         //get_rsnp_flag
         rsnp_flag.clear();
         rsnp_flag.resize(bdata->_snp_num);
         for (int i = 0; i < bdata->_snp_num; i++) {
             if (bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()) rsnp_flag[i] = 1;
             else rsnp_flag[i] = 0;
         }
         update_fam(bdata, rindi_flag);
         update_bim(bdata, rsnp_flag);

         retrieve_snp(snp_name_per_chr, bdata->_snp_name_map, rsnp, nbfiles);

         // Read the coded genotypes
         for(i=0; i<nbfiles; i++) {
             if(rsnp[i].size()==0) {
                 cout<<"Skip reading " + multi_bfiles[i] + ".bed, no SNPs retained on this chromosome."<<endl;
                 continue;
             } else {
                 bedfile = multi_bfiles[i] + ".bed";
                 read_single_bedfile(bdata, bedfile, rsnp[i], rindi_flag, false);
             }
         }
         cout << "Genotype data for " + to_string(bdata->_keep.size()) + " individuals and " + to_string(bdata->_include.size()) + " SNPs have been included." <<endl;
         
     }

     void update_id_chr_map(map<string, string> &chr_map, map<string, int> id_map)
     {
         int i = 0;
         map<string, string> chr_map_buf(chr_map);
         map<string, int>::iterator iter1;
         map<string, string>::iterator iter2;

         for(iter1=id_map.begin(); iter1!=id_map.end(); iter1++) chr_map_buf.erase(iter1->first);
         for(iter2=chr_map_buf.begin(); iter2!=chr_map_buf.end(); iter2++) chr_map.erase(iter2->first);
     }

     void retrieve_snp(map<string,string> snp_chr_map, map<string,int> snp_id_map, vector<vector<pair<int,int>>> &rsnp, int nbfiles)
     {
         int i = 0, j=0, snp_indx = 0;
         string snp_indx_str = "";
         vector<string> vs_buf;
         map<string,string>::iterator iter1;
         map<string,int>::iterator iter2;

         rsnp.clear(); rsnp.resize(nbfiles);

         for(iter1=snp_chr_map.begin(), iter2=snp_id_map.begin(); iter1 != snp_chr_map.end(); iter1++, iter2++) {
             vs_buf.clear();
             snp_indx_str = iter1->second;
             StrFunc::split_string(snp_indx_str, vs_buf, ":");
             snp_indx = iter2->second;
             rsnp[atoi(vs_buf[0].c_str())].push_back(make_pair(atoi(vs_buf[1].c_str()), snp_indx));
         }

         for(i=0; i<nbfiles; i++) stable_sort(rsnp[i].begin(), rsnp[i].end());
     }

     void read_single_bedfile(bInfo* bdata, string bedfile, vector<pair<int,int>> rsnp, vector<int> rindi, bool msg_flag)
     {
         int i = 0, j = 0, k = 0, t1=0, nsnp_chr = rsnp.size(), nindi_chr = rindi.size();
        
         // Read bed file
         char ch[1];
         bitset<8> b;
         fstream BIT(bedfile.c_str(), ios::in | ios::binary);
         if (!BIT) throw ("Error: can not open the file [" + bedfile + "] to read.");
         if(msg_flag) cout << "Reading PLINK BED file from [" + bedfile + "] in SNP-major format ..." << endl;

         for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
         int snp_indx = 0, indi_indx = 0;
         for (j = 0, t1 = 0, snp_indx = 0; t1 < nsnp_chr; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
             if (j!=rsnp[t1].first) {
                 for (i = 0; i < nindi_chr; i += 4) BIT.read(ch, 1);
                 continue;
             }
             snp_indx = rsnp[t1].second;
             for (i = 0, indi_indx = 0; i < nindi_chr;) {
                 BIT.read(ch, 1);
                 if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
                 b = ch[0];
                 k = 0;
                 while (k < 7 && i < nindi_chr) { // change code: 11 for AA; 00 for BB;
                     if (!rindi[i]) k += 2;
                     else {
                         bdata->_snp_2[snp_indx][indi_indx] = (!b[k++]);
                         bdata->_snp_1[snp_indx][indi_indx] = (!b[k++]);
                         indi_indx++;
                     }
                     i++;
                 }
             }
             t1++;
         }
         BIT.clear();
         BIT.close();
         if(msg_flag) cout << "Genotype data for " << bdata->_keep.size() << " individuals and " << bdata->_include.size() << " SNPs to be included from [" + bedfile + "]." << endl;
     }
    // end
    
    void read_gwas_data(gwasData* gdata, char* gwasFileName)
    {
        bool warnnullfreq=false;
        ifstream gwasFile;
        if(!file_read_check(&gwasFile, gwasFileName))
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     gwasFileName, strerror (errno));
            exit (EXIT_FAILURE);
        }
        cout << "Reading GWAS summary data from [" + string(gwasFileName) + "]." << endl;
        gdata->_include.clear();
        gdata->_snp_name_map.clear();
        gdata->snpName.clear();
        gdata->snpBp.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        gwasFile.getline(buf,MAX_LINE_SIZE);// the header
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",gwasFileName);
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: %s should have headers that start with \"snp\".\n", gwasFileName);
            exit(EXIT_FAILURE);
        }
        while(!gwasFile.eof())
        {
            gwasFile.getline(buf,MAX_LINE_SIZE);
            
            if(buf[0]!='\0')
            {
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=8) {
                    printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(gdata->_snp_name_map.find(vs_buf[0]) != gdata->_snp_name_map.end()){
                    cout << "WARNING: Duplicated SNP ID \"" + vs_buf[0] + "\" ";
                    stringstream ss;
                    ss << vs_buf[0] << "_" << lineNum + 1;
                    vs_buf[0] = ss.str();
                    cout<<"has been changed to \"" + vs_buf[0] + "\".\n";
                }
                gdata->_snp_name_map.insert(pair<string, int>(vs_buf[0], lineNum));
                gdata->snpName.push_back(vs_buf[0]);
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: allele1 is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[1]);
                gdata->allele_1.push_back(vs_buf[1]);
                
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: allele2 is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[2]);
                gdata->allele_2.push_back(vs_buf[2]);
                
                if(vs_buf[3]=="NA" || vs_buf[3]=="na")
                {
                   
                        if(!warnnullfreq){
                            warnnullfreq=true;
                            printf("WARNING: frequency is \'NA\' in one or more rows.\n");
                        }
                        gdata->freq.push_back(-9);
                    
                }
                else {
                    gdata->freq.push_back(atof(vs_buf[3].c_str()));
                }
               
                
                if(vs_buf[4]=="NA" || vs_buf[4]=="na"){
                    printf("WARNING: effect size is \'NA\' in row %d.\n", lineNum+2);
                    gdata->byz.push_back(0);
                } else {
                    gdata->byz.push_back(atof(vs_buf[4].c_str()));
                }
                if(vs_buf[5]=="NA" || vs_buf[5]=="na"){
                    printf("WARNING: standard error is \'NA\' in row %d.\n", lineNum+2);
                    gdata->seyz.push_back(-9);
                } else {
                    gdata->seyz.push_back(atof(vs_buf[5].c_str()));
                }

                gdata->pvalue.push_back(atof(vs_buf[6].c_str()));
                gdata->splSize.push_back(atof(vs_buf[7].c_str()));
                gdata->_include.push_back(lineNum);
                lineNum++;
            }
        }
        gdata->snpNum=gdata->_include.size();
       cout <<"GWAS summary data of "<<gdata->snpNum << " SNPs to be included from [" + string(gwasFileName) + "]." << endl;
        gwasFile.close();
    }
	/*
    void update_include_map(eqtlInfo* eqtlinfo)
    {
        eqtlinfo->_incld_id_map.clear();
        long size=0;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            eqtlinfo->_incld_id_map.insert(pair<int,int>(eqtlinfo->_esi_include[i],i));
            if (size == eqtlinfo->_incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
            size = eqtlinfo->_incld_id_map.size();
        }
        
    }
	*/
    void read_esifile(eqtlInfo* eqtlinfo, string esifile, bool prtscr)
    {
        ifstream esi(esifile.c_str());
        if (!esi) throw ("ERROR: can not open the file [" + esifile + "] to read.");
        if(prtscr) cout << "Reading xQTL SNP information from [" + esifile + "]." << endl;
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
                    printf("ERROR: chromosome is \"NA\" in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                eqtlinfo->_esi_chr.push_back(tmpchr);
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(eqtlinfo->_snp_name_map.find(vs_buf[1]) != eqtlinfo->_snp_name_map.end()){
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
                    printf("ERROR: SNP BP is \'NA\' in row %d.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                eqtlinfo->_esi_bp.push_back(atoi(vs_buf[3].c_str()));
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
                        //printf("WARNING: frequency column missing.\n");
                        eqtlinfo->_esi_freq.push_back(-9);
                    }
                
                
                eqtlinfo->_esi_include.push_back(lineNum);
                lineNum++;
            }
        }
        eqtlinfo->_snpNum=lineNum;
        if(prtscr) cout << eqtlinfo->_snpNum << " SNPs to be included from [" + esifile + "]." << endl;
        esi.close();
    }
    void read_epifile(eqtlInfo* eqtlinfo, string epifile, bool prtscr)
    {
        ifstream epi(epifile.c_str());
        if (!epi) throw ("ERROR: can not open the file [" + epifile + "] to read.");
        if(prtscr) cout << "Reading xQTL probe information from [" + epifile + "]." << endl;
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
        if(prtscr) cout << eqtlinfo->_probNum << " Probes to be included from [" + epifile + "]." << endl;
        
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
                printf("ERROR: chromosome is \"NA\" in row %d.\n", i+1);
                exit(EXIT_FAILURE);
            }
            int tmpchr;
            if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
            else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
            else tmpchr=atoi(tmpStr.c_str());
            eqtlinfo->_epi_chr[i]=tmpchr;
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
                printf("ERROR: probe BP is \"NA\" in row %d.\n", i+1);
                exit(EXIT_FAILURE);
            }
            eqtlinfo->_epi_bp[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_epi_gene[i]=tmpStr.c_str();
            iss>>tmpStr;
            eqtlinfo->_epi_orien[i]=tmpStr.c_str()[0];
            
        }
        
        epi.close();
    }    
    int shown(string besdfile)
    {
        string fname=besdfile+".besd";
        FILE* besd=fopen(fname.c_str(), "rb");
        if(!besd)
        {
            printf ( "ERROR: Couldn't open file %s\n", fname.c_str());
            exit (EXIT_FAILURE);
        }

        printf("Reading sample size from %s \n",fname.c_str());
        uint32_t indicator;
        int ss=-9;
        if(fread(&indicator, sizeof(uint32_t),1, besd)!=1)
        {
            printf("ERROR: File %s read failed!\n", fname.c_str());
            exit (EXIT_FAILURE);
        }
        if(indicator==SPARSE_FILE_TYPE_3 || indicator==DENSE_FILE_TYPE_3)
        {
            if(fread(&ss, sizeof(int),1, besd)!=1)
            {
                printf("ERROR: File %s read failed!\n", fname.c_str());
                exit (EXIT_FAILURE);
            }
            if(ss==-9)
            {
                printf("The sample size is missing. You may use --add-n to add it to the BESD file.\n");
            }
            else {
                printf("The sample size is %d\n",ss);
            }
        } else {
            printf("This file is in the old BESD format which doesn't contain the information of sample size.\n");
        }
        fclose(besd);
        return ss;
    }
    void read_besdfile(eqtlInfo* eqtlinfo, string besdfile, bool prtscr)
    {
        if (eqtlinfo->_include.size() == 0) throw ("Error: No probe is retained for analysis.");
        if (eqtlinfo->_esi_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        
        eqtlinfo->_cols.clear();
        eqtlinfo->_rowid.clear();
        eqtlinfo->_val.clear();
        eqtlinfo->_valNum = 0;
        eqtlinfo->_bxz.clear();
        eqtlinfo->_sexz.clear();
        
        // the fastest way is using malloc and memcpy
        char SIGN[sizeof(uint64_t)+8];
        ifstream besd(besdfile.c_str(), ios::in|ios::binary);
        if(!besd)
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     besdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        if(prtscr)  cout << "Reading xQTL summary data from [" + besdfile + "]." << endl;
        
        besd.read(SIGN, 4);
        uint32_t gflag = *(uint32_t *)SIGN;
        /*
        if(gflag==4){
            // clear datastruct for sparse befor read dense
            eqtlinfo->_cols.clear();
            eqtlinfo->_rowid.clear();
            eqtlinfo->_val.clear();
            eqtlinfo->_valNum = 0;
            uint64_t memsize2use=eqtlinfo->_include.size()*eqtlinfo->_esi_include.size()*2*sizeof(float);
            if(memsize2use>0x200000000) printf("WARNING: %llu GB should be allocated for your besd file.\n",memsize2use>>30);
            eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
            eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
            for(unsigned int i=0;i<eqtlinfo->_include.size();i++)
            {
                eqtlinfo->_bxz[i].resize(eqtlinfo->_esi_include.size());
                eqtlinfo->_sexz[i].resize(eqtlinfo->_esi_include.size());
            }
            char* buffer;
            buffer = (char*) malloc (sizeof(char)*eqtlinfo->_snpNum<<3);
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            float* ft;
            float* se_ptr;
           
                for(int j=0;j<eqtlinfo->_esi_include.size();j++)
                {
                    unsigned long sid=eqtlinfo->_esi_include[j];
                    besd.seekg(((sid<<1)*eqtlinfo->_probNum+1)<<2);
                    memset(buffer,0,sizeof(char)*eqtlinfo->_probNum<<3);
                    besd.read(buffer,eqtlinfo->_probNum<<3);
                    ft=(float *)buffer;
                    for (int i = 0; i<eqtlinfo->_include.size(); i++) eqtlinfo->_bxz[i][j] = *(ft + eqtlinfo->_include[i]);
                    se_ptr = ft + eqtlinfo->_probNum;
                    for (int i = 0; i<eqtlinfo->_include.size(); i++) eqtlinfo->_sexz[i][j] = *(se_ptr + eqtlinfo->_include[i]);
                }
                if(prtscr)  std::cout << "eQTL summary-level statistics of " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_esi_include.size() << " SNPs to be included from [" + besdfile + "]." << endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            
                free(buffer);

        }
         */
        
        if(gflag == 0x40000000){
            // clear datastruct for dense befor read sparse
            cout<<"This is an old file format. Please use --make-besd to update the file format."<<endl;
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            
            uint64_t colNum=(eqtlinfo->_probNum<<1)+1;
            uint64_t valNum;
            uint64_t lSize;
            char* buffer;
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            besd.read(SIGN, sizeof(uint64_t));
            valNum=*(uint64_t *)SIGN;
            if( lSize - (sizeof(float) + sizeof(uint64_t) + (colNum+valNum)*sizeof(uint32_t) + valNum*sizeof(float)) != 0)
            {
                printf("The file size is %llu",lSize);
                printf(" %zu + %zu + %lld + %lld + %lld \n",sizeof(float),sizeof(uint64_t), colNum*sizeof(uint32_t),valNum*sizeof(uint32_t), valNum*sizeof(float));
                printf("ERROR: failed in binary file check.\n");
                exit(EXIT_FAILURE);
            }
            
            
            buffer = (char*) malloc (sizeof(char)*(lSize));
            if (buffer == NULL) {fputs ("Memory error.\n",stderr); exit (1);}
            besd.read(buffer,lSize);
            if (besd.gcount()+sizeof(float) + sizeof(uint64_t) != lSize) {fputs ("Reading error",stderr); exit (2);}
            
            
            uint32_t* ptr;
            ptr=(uint32_t *)buffer;
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=*ptr;
                uint32_t* row_ptr;
                row_ptr=ptr+colNum;
                float* val_ptr;
                val_ptr=(float*)(row_ptr+valNum);
                
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
                    size = _incld_id_map.size();
                }

                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
                    uint32_t pos=*(ptr+(pid<<1));
                    uint32_t pos1=*(ptr+(pid<<1)+1);
                    uint32_t num=pos1-pos;
                    uint32_t real_num=0;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=*(row_ptr+pos+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                        
                       // long sid=find(eqtlinfo->_esi_include.begin(),eqtlinfo->_esi_include.end(),rid)-eqtlinfo->_esi_include.begin(); //slow
                      //  if(sid<eqtlinfo->_esi_include.size())
                      //  {
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                            real_num++;
                        }
                       
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                if(prtscr)  cout<<"xQTL summary data of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                eqtlinfo->_cols.resize(colNum);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) eqtlinfo->_cols[i]=*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=*ptr++;
                float* val_ptr=(float*)ptr;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
               if(prtscr)  cout<<"xQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        else if(gflag == DENSE_FILE_TYPE_1 || gflag == DENSE_FILE_TYPE_3)
        {
            // clear datastruct for sparse befor read dense
            if(gflag==DENSE_FILE_TYPE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                besd.read(indicators,length);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    printf("The sample size is %d.\n",ss);
                }
                if(*tmp++!=eqtlinfo->_snpNum)
                {
                    printf("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                if(*tmp++!=eqtlinfo->_probNum)
                {
                    printf("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                delete[] indicators;
            }
            int infoLen=sizeof(uint32_t);
            if(gflag==DENSE_FILE_TYPE_3) infoLen=RESERVEDUNITS*sizeof(int);
            eqtlinfo->_cols.clear();
            eqtlinfo->_rowid.clear();
            eqtlinfo->_val.clear();
            eqtlinfo->_valNum = 0;
            uint64_t memsize2use=eqtlinfo->_include.size()*eqtlinfo->_esi_include.size()*2*sizeof(float);
            if(memsize2use>0x200000000) printf("WARNING: %llu GB should be allocated for your besd file.\n",memsize2use>>30);
            eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
            eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
            for(unsigned int i=0;i<eqtlinfo->_include.size();i++)
            {
                eqtlinfo->_bxz[i].resize(eqtlinfo->_esi_include.size());
                eqtlinfo->_sexz[i].resize(eqtlinfo->_esi_include.size());
            }
            char* buffer;
            buffer = (char*) malloc (sizeof(char)*eqtlinfo->_snpNum<<3);
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            float* ft;
            float* se_ptr;
            if (eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)  //means with the parameter --extract-probe. This also can read all the probes, but currently I don't think it is good for too many I/Os.
            {
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    unsigned long pid=eqtlinfo->_include[i];
                    besd.seekg(((pid*eqtlinfo->_snpNum)<<3)+infoLen);
                    memset(buffer,0,sizeof(char)*eqtlinfo->_snpNum<<3);
                    besd.read(buffer,eqtlinfo->_snpNum<<3);
                    ft=(float *)buffer;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_bxz[i][j] = *(ft + eqtlinfo->_esi_include[j]);
                    se_ptr = ft + eqtlinfo->_snpNum;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_sexz[i][j] = *(se_ptr + eqtlinfo->_esi_include[j]);
                }
                if(prtscr)  std::cout << "xQTL summary-level statistics of " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_esi_include.size() << " SNPs to be included from [" + besdfile + "]." << endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                //without --extract-probe, read with less I/O. and need not to update epi.
                //read with static buffer. If dynamic buffer, 2GB per I/O can be more efficient.
                /*
                 unsigned long long count=0;
                 while(!besd.eof())
                 {
                 besd.read(buf,MAX_LINE_NUM);
                 unsigned long Bread=besd.gcount();
                 buf[Bread]='\0';
                 char* ptr=buf;
                 //while(*ptr != '\0') //can not use this, too many 0x00 in buf
                 while(Bread)
                 {
                 unsigned long pid=count/eqtlinfo->_snpNum;
                 unsigned long sid=count++%eqtlinfo->_snpNum;
                 ft=(float *)ptr;
                 if(pid&1) eqtlinfo->_sexz[pid>>1][sid]=*ft;
                 else eqtlinfo->_bxz[pid>>1][sid]=*ft;
                 ptr+=4;
                 Bread-=4;
                 }
                 }
                 cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                 */
                
                /*
                 
                 // besd.seekg(0,besd.end);
                 // uint64_t lSize = besd.tellg();
                 //  lSize-=4;
                 //  besd.seekg(4); // same as besd.seekg(4, besd.beg);
                 uint64_t alread=0;
                 
                 char* buff;
                 uint64_t buffszie=0x40000000;
                 buff = (char*) malloc (sizeof(char)*buffszie);
                 if (buff == NULL) {fputs ("Memory error",stderr); exit (1);}
                 memset(buff,0,sizeof(char)*buffszie);
                 uint64_t count=0;
                 while(!besd.eof())
                 {
                 besd.read(buff,buffszie);
                 unsigned long Bread=besd.gcount();
                 alread+=Bread;
                 char* ptr=buff;
                 while(Bread)
                 {
                 unsigned long pid=count/eqtlinfo->_snpNum;
                 unsigned long sid=count++%eqtlinfo->_snpNum;
                 ft=(float *)ptr;
                 if(pid&1) eqtlinfo->_sexz[pid>>1][sid]=*ft;
                 else eqtlinfo->_bxz[pid>>1][sid]=*ft;
                 ptr+=4;
                 Bread-=4;
                 }
                 cout<<alread<<":"<<(alread>>30)<<"GB "<<endl;
                 }
                 cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besd
                 file + "]." <<endl;
                 free(buff);
                 */
                
                /*
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {fputs ("Memory error",stderr); exit (1);}
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;  //should be even number
                probonce>>=1;
                probonce<<=1;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!besd.eof())
                {
                    besd.read(buff,readsize);
                    unsigned long Bread=besd.gcount();
                    float* rptr=(float *)buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],rptr,perbeta);
                        rptr+=eqtlinfo->_snpNum;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],rptr,perbeta);
                        rptr+=eqtlinfo->_snpNum;
                        Bread-=(perbeta<<1);
                    }
                     cout<<probcount<<" done! "<<endl;
                }
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                free(buff);
                 */
                
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {fputs ("Memory error when reading dense BESD file.",stderr); exit (1);}
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;  //should be even number
                probonce>>=1;
                probonce<<=1;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!besd.eof())
                {
                    besd.read(buff,readsize);
                    uint64_t Bread=besd.gcount();
                    char* rptr=buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],rptr,perbeta);
                        rptr+=perbeta;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],rptr,perbeta);
                        rptr+=perbeta;
                        Bread-=(perbeta<<1);
                    }
                    printf("Reading... %3.0f%%\r", 100.0*probcount/eqtlinfo->_probNum);
                    fflush(stdout);
                }
                if(prtscr) cout<<"\neQTL summary data of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                free(buff);

            }
            free(buffer);
        }
        
        else if (gflag == 0x3f800000 )
        {
            cout<<"This is an old file format. Please use --make-besd to update the file format."<<endl;
            // clear datastruct for dense befor read sparse
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            
            uint64_t colNum=eqtlinfo->_probNum<<1;
            uint64_t valNum;
            uint64_t lSize;
            char* buffer;
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            besd.read(SIGN, 4);
            valNum=(uint32_t)*(float *)SIGN; // int to float then float to int back can lose pricision. hence this clause and bellow are unbelievable
             if(lSize-((3+colNum+(valNum<<1))<<2) != 0)
             {
                 printf("The file size is %llu",lSize);
                 printf(" %zu + %zu + %lld + %lld + %lld \n",sizeof(float),sizeof(float), (1+colNum)*sizeof(float),valNum*sizeof(float), valNum*sizeof(float));
                 printf("ERROR: failed in binary file check.\n");
                 exit(EXIT_FAILURE);
             }
            
            valNum=((lSize>>2)-3-colNum)>>1;
            
            buffer = (char*) malloc (sizeof(char)*(lSize-8));
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            besd.read(buffer,(lSize-8));
            if (besd.gcount()+8 != lSize) {fputs ("Reading error",stderr); exit (2);}
            float* ptr;
            ptr=(float *)buffer;
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {                
                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=(uint32_t)*ptr;
                float* row_ptr;
                row_ptr=ptr+colNum+1;
                float* val_ptr;
                val_ptr=row_ptr+valNum;

				map<int, int > _incld_id_map;
				long size = 0;
				for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
				{
					_incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
					if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
					size = _incld_id_map.size();
				}

                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    unsigned long pid=eqtlinfo->_include[i];
                    uint32_t pos=(uint32_t)*(ptr+(pid<<1));
                    uint32_t pos1=(uint32_t)*(ptr+(pid<<1)+1);
                    uint32_t num=pos1-pos;
                    uint32_t real_num=0;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=(uint32_t)*(row_ptr+pos+j);
						 
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;							
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                            real_num++;
                        }
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                if(prtscr)  cout<<"xQTL summary data of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                eqtlinfo->_cols.resize(colNum+1);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                for(int i=0;i<=colNum;i++) eqtlinfo->_cols[i]=(uint32_t)*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=(uint32_t)*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*ptr++;
                eqtlinfo->_valNum = valNum;
                if(prtscr)  cout<<"xQTL summary data of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        
        else if (gflag == SPARSE_FILE_TYPE_3F || gflag == SPARSE_FILE_TYPE_3)
        {
            // clear datastruct for dense befor read sparse
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            char* buffer;
            uint64_t colNum=(eqtlinfo->_probNum<<1)+1;
            uint64_t valNum;
            uint64_t lSize;
           
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            if(gflag==SPARSE_FILE_TYPE_3)
            {
                int length=(RESERVEDUNITS-1)*sizeof(int);
                char* indicators=new char[length];
                besd.read(indicators,length);
                int* tmp=(int *)indicators;
                int ss=*tmp++;
                if(ss!=-9)
                {
                    printf("The sample size is %d.\n",ss);
                }
                if(*tmp++!=eqtlinfo->_snpNum)
                {
                    printf("ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file %s.\n", besdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                if(*tmp++!=eqtlinfo->_probNum)
                {
                    printf("ERROR: The probes in your .epi file are not in consistency with the one in .besd file %s.\n", besdfile.c_str());
                    exit(EXIT_FAILURE);
                }
                delete[] indicators;
            }

            besd.read(SIGN, sizeof(uint64_t));
            valNum=*(uint64_t *)SIGN;
            if(gflag==SPARSE_FILE_TYPE_3F) {
                if( lSize - (sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0)
                {
                    printf("The file size is %llu",lSize);
                    printf(" %zu + %zu + %lld + %lld + %lld \n",sizeof(uint32_t),sizeof(uint64_t), colNum*sizeof(uint64_t),valNum*sizeof(uint32_t), valNum*sizeof(float));
                    printf("ERROR: wrong value number. File %s is ruined.\n", besdfile.c_str());
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                if( lSize - (RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0)
                {
                    printf("The file size is %llu",lSize);
                    printf(" %zu + %zu + %lld + %lld + %lld \n",RESERVEDUNITS*sizeof(int),sizeof(uint64_t), colNum*sizeof(uint64_t),valNum*sizeof(uint32_t), valNum*sizeof(float));
                    printf("ERROR: wrong value number. File %s is ruined.\n", besdfile.c_str());
                    exit (EXIT_FAILURE);
                }
                
            }
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                
                uint64_t colsize=colNum*sizeof(uint64_t);
                buffer = (char*) malloc (sizeof(char)*(colsize));
                if (buffer == NULL) {fputs ("Memory error when reading sparse BESD file.",stderr); exit (1);}
                besd.read(buffer,colsize);
                
                uint64_t* ptr;
                ptr=(uint64_t *)buffer;

                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=*ptr;
             
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
                    size = _incld_id_map.size();
                }

                uint64_t rowSTART=0;
                uint64_t valSTART=0;
                if(gflag==SPARSE_FILE_TYPE_3F)
                {
                    rowSTART=sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                    valSTART=sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                } else {
                    rowSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                    valSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                }
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
                    uint64_t pos=*(ptr+(pid<<1)); //BETA START
                    uint64_t pos1=*(ptr+(pid<<1)+1); //SE START
                    uint64_t num=pos1-pos;
                    uint64_t real_num=0;
                    if(num==0) {
                        eqtlinfo->_cols[(i<<1)+1]=eqtlinfo->_cols[i<<1];
                        eqtlinfo->_cols[i+1<<1]=eqtlinfo->_cols[i<<1];
                        //printf("WARNING: Probe %s with no eQTL found.\n",eqtlinfo->_epi_prbID[pid].c_str());
                        continue;
                        
                    }
                    char* row_char_ptr;
                    row_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(uint32_t));
                    if (row_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
                    char* val_char_ptr;
                    val_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(float));
                    if (val_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
                    memset(row_char_ptr,0,sizeof(char)*2*num*sizeof(uint32_t));
                    memset(val_char_ptr,0,sizeof(char)*2*num*sizeof(float));
                    besd.seekg(rowSTART+pos*sizeof(uint32_t));
                    besd.read(row_char_ptr, 2*num*sizeof(uint32_t));
                    uint32_t* row_ptr=(uint32_t *)row_char_ptr;                    
                    besd.seekg(valSTART+pos*sizeof(float));
                    besd.read(val_char_ptr, 2*num*sizeof(float));
                    float* val_ptr=(float*)val_char_ptr;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=*(row_ptr+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                            
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+j));
                            real_num++;
                        }
                        
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                    free(row_char_ptr);
                    free(val_char_ptr);
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
               
                if(prtscr)  cout<<"xQTL summary data of "<<eqtlinfo->_include.size()<<" Probes to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                buffer = (char*) malloc (sizeof(char)*(lSize));
                if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
                besd.read(buffer,lSize);
                if(gflag==SPARSE_FILE_TYPE_3F)
                {
                    if (besd.gcount()+sizeof(uint32_t) + sizeof(uint64_t) != lSize)
                    {
                        printf("ERROR: reading file %s error.\n", besdfile.c_str());
                        exit (EXIT_FAILURE);
                    }
                }else {
                    if (besd.gcount()+RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) != lSize)
                    {
                        printf("ERROR: reading file %s error.\n", besdfile.c_str());
                        exit (EXIT_FAILURE);
                    }
                }
                
                uint64_t* ptr;
                ptr=(uint64_t *)buffer;

                eqtlinfo->_cols.resize(colNum);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) eqtlinfo->_cols[i]=*ptr++;
                uint32_t* ptr4B=(uint32_t *)ptr;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=*ptr4B++;
                float* val_ptr=(float*)ptr4B;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
               if(prtscr)  cout<<"xQTL summary data of "<<eqtlinfo->_probNum<<" Probes to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        else {
            cout<<"SMR doesn't support this format. Please use OSCA (http://cnsgenomics.com/software/osca) to transform it to SMR format."<<endl;
        }
        besd.close();
        /*
        if(eqtlinfo->_rowid.empty() && eqtlinfo->_bxz.empty())
        {
            printf("NO data included from eQTL summary data %s, please check.\n",besdfile.c_str()); exit (EXIT_FAILURE);
        }
         */
    }
    
    void filter_probe_null(eqtlInfo* eqtlinfo)
    {
        vector<string> nullprobes;
        cout<<"\nfiltering out the probes with no value..."<<endl;
        if(eqtlinfo->_valNum==0)
        {
            eqtlinfo->_include.clear();
            for (int i = 0; i < eqtlinfo->_probNum; i++)
            {
                bool NA_flag = true;
                for (int j = 0; j<eqtlinfo->_snpNum; j++)
                {
                    if (fabs(eqtlinfo->_sexz[i][j] + 9) > 1e-6)
                    {
                        NA_flag = false;
                        break;
                    }
                }
                if (!NA_flag) eqtlinfo->_include.push_back(i);
                else nullprobes.push_back(eqtlinfo->_epi_prbID[i]);
            }
        }
        else{
            eqtlinfo->_include.clear();
            for (int i = 0; i < eqtlinfo->_probNum; i++)
            {
                if (eqtlinfo->_cols[(i<<1)+1] > eqtlinfo->_cols[i<<1]) eqtlinfo->_include.push_back(i);
                else nullprobes.push_back(eqtlinfo->_epi_prbID[i]);
            }
        }
        eqtlinfo->_probe_name_map.clear();
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            eqtlinfo->_probe_name_map.insert(pair<string,int>(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]],eqtlinfo->_include[i]));
        }
        if(nullprobes.size()>0)
        {
            string fname="chr"+atos(eqtlinfo->_esi_chr[0])+".nullprobes.log";
            FILE* nullprobefile=fopen(fname.c_str(), "w");
            if (!(nullprobefile)) {
                printf("Error: Failed to open null probe log file.\n");
            }
            for(int i=0;i<nullprobes.size();i++)
            {
                string tmpstr=nullprobes[i]+'\n';
                fputs(tmpstr.c_str(),nullprobefile);
            }
            fclose(nullprobefile);
        }
        
        cout<<eqtlinfo->_include.size()<<" probes to be included."<<endl;
    }
    void filter_snp_null(eqtlInfo* eqtlinfo)
    {
        vector<string> nullSNPs;
        cout<<"\nfiltering out the SNPs with no value..."<<endl;
        vector<int> esi_include;
        if(eqtlinfo->_valNum==0)
        {
            for (int i = 0; i < eqtlinfo->_esi_include.size(); i++)
            {
                bool NA_flag = true;
                for (int j = 0; j<eqtlinfo->_include.size(); j++)
                {
                    if (fabs(eqtlinfo->_sexz[eqtlinfo->_include[j]][eqtlinfo->_esi_include[i]] + 9) > 1e-6)
                    {
                        NA_flag = false;
                        break;
                    }
                }
                if (!NA_flag)esi_include.push_back(i);
                else nullSNPs.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
            }
        }
        else{
            vector<uint32_t> rowid;
            rowid=eqtlinfo->_rowid;
            getUnique(rowid);
            vector<int> idx;
           
            map<uint32_t, int> id_map;
            map<uint32_t, int>::iterator iter;
            for(int i=0; i<rowid.size(); i++) id_map.insert(pair<uint32_t,int>(rowid[i], i));
            for(int i=0; i< eqtlinfo->_esi_include.size(); i++){
                iter=id_map.find(eqtlinfo->_esi_include[i]);
                if(iter==id_map.end()) idx.push_back(-9);
                else idx.push_back(iter->second);
            }
        
            for (int i = 0; i < idx.size(); i++)
            {
                if (idx[i]!=-9) esi_include.push_back(eqtlinfo->_esi_include[i]);
                else nullSNPs.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
            }
        }
        eqtlinfo->_esi_include.swap(esi_include);
        eqtlinfo->_snp_name_map.clear();
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            eqtlinfo->_snp_name_map.insert(pair<string,int>(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]],eqtlinfo->_esi_include[i]));
        }
        if(nullSNPs.size()>0)
        {
            string fname="chr"+atos(eqtlinfo->_esi_chr[0])+".nullSNPs.log";
            FILE* nullsnpfile=fopen(fname.c_str(), "w");
            if (!(nullsnpfile)) {
                printf("Error: Failed to open null probe log file.\n");
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<nullSNPs.size();i++)
            {
                string tmpstr=nullSNPs[i]+'\n';
                fputs(tmpstr.c_str(),nullsnpfile);
            }
            fclose(nullsnpfile);
        }
        
        cout<<eqtlinfo->_esi_include.size()<<" SNPs to be included."<<endl;
    }
    


    bool has_suffix(const std::string &str, const std::string &suffix)
    {
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    bool has_prefix(const std::string &str, const std::string &prefix)
    {
        return str.size() >= prefix.size() &&
        str.compare(0, prefix.size(), prefix) == 0;
    }

    void get_square_idxes(vector<int> &sn_ids,VectorXd &zsxz,double threshold)
    {
        
        for(int i=0;i<zsxz.size();i++)
        {
           if(zsxz[i]*zsxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
    }
    void get_square_idxes_sort(vector<int> &sn_ids,VectorXd &zsxz,double threshold)
    {
        vector<double> z;
        vector<int> rk;
        for(int i=0;i<zsxz.size();i++)
        {
            if(zsxz[i]*zsxz[i]-threshold>1e-6)
            {
                sn_ids.push_back(i);
                z.push_back(zsxz(i));
            }
        }
        getRank(z,rk);
        
    }
	void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, MatrixXd &LD, long maxid, double ld_top)
    {
        for(int i=0;i<zsxz.size();i++)
        {
            if(i!= maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (LD(maxid,i)*LD(maxid,i)-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                 if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
            
        }
        
    }
    void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, VectorXd &ld_v, long maxid, double ld_top)
    {
        for(int i=0;i<zsxz.size();i++)
        {
            if(i!= maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (ld_v(i)*ld_v(i)-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
        }
    }
   void est_cov_bxy(MatrixXd &covbxy, VectorXd &_zsxz,VectorXd &_bxy,VectorXd &_seyz,VectorXd &_bxz, MatrixXd &_LD_heidi)
    {
        long nsnp =_zsxz.size();
        if(nsnp>1)
        {          
           MatrixXd bxytbxy= _bxy*_bxy.transpose();
           MatrixXd zsxztzsxz= _zsxz*_zsxz.transpose();
           covbxy=_LD_heidi.array()*((_seyz*_seyz.transpose()).array()/(_bxz*_bxz.transpose()).array() + bxytbxy.array()/zsxztzsxz.array()) - bxytbxy.array()/(zsxztzsxz.array()*zsxztzsxz.array());
             
        }
    }
    void est_cov_bxy_so(MatrixXd &covbxy,VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz,VectorXd &_sexz, MatrixXd &_LD_heidi, double theta)
    {
        long nsnp =_byz.size();
        if(nsnp>1)
        {
            VectorXd ivect = VectorXd::Ones(nsnp);
            MatrixXd bexpo=_bxz*_bxz.transpose();
            MatrixXd seexpo=_sexz*_sexz.transpose();
            MatrixXd boutco=_byz*_byz.transpose();
            MatrixXd seoutco=_seyz*_seyz.transpose();
            VectorXd _bxz2=_bxz.array()*_bxz.array();
            MatrixXd cov1 = _LD_heidi.array()*(seoutco.array()/bexpo.array()) + _LD_heidi.array()*seexpo.array()*boutco.array()/(bexpo.array()*bexpo.array());
            MatrixXd cov2 = -2 * theta*_LD_heidi.array()*(_seyz*_sexz.transpose()).array()*(_byz*ivect.transpose()).array()/(_bxz2*_bxz.transpose()).array();
        
            covbxy= cov1 + cov2;
            
        }
    }
  
   float bxy_hetero3(VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num)
    {
        VectorXd _bxy;
        VectorXd _sexy;
        VectorXd dev;
        VectorXd tmp3;
        long nsnp=*snp_num;
        int maxid;        
        float pdev=-1.0;
        MatrixXd covbxy(nsnp,nsnp);
        MatrixXd vdev(nsnp-1,nsnp-1);
       
        VectorXd bxz2;
        dev.resize(nsnp-1);
        if(nsnp>1)
        {
            _bxy=_byz.array()/_bxz.array();
            bxz2=_bxz.array()*_bxz.array();
            _sexy=(_seyz.array()*_seyz.array()*bxz2.array()+_sexz.array()*_sexz.array()*_byz.array()*_byz.array())/(bxz2.array()*bxz2.array()).sqrt();
            
            maxid=max_abs_id(_zsxz);
            
            for(int j=0;j<maxid;j++) dev[j]=_bxy[maxid]-_bxy[j];
            for(int j=maxid+1;j<nsnp;j++) dev[j-1]=_bxy[maxid]-_bxy[j];			

            est_cov_bxy(covbxy, _zsxz,_bxy,_seyz,_bxz,_LD_heidi);            
            
            
            double tmp1=covbxy(maxid,maxid);           
            tmp3.resize(nsnp-1);
            for(int i=0; i<maxid; i++) tmp3[i]=covbxy(maxid,i);
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=covbxy(maxid,i);
            // vdev as tmp2
			vdev.block(0, 0, maxid, maxid) = covbxy.block(0, 0, maxid, maxid);
			vdev.block(0, maxid, maxid, nsnp - maxid - 1) = covbxy.block(0, maxid + 1, maxid, nsnp - maxid - 1);
			vdev.block(maxid, 0, nsnp - maxid - 1, maxid) = covbxy.block(maxid + 1, 0, nsnp - maxid - 1, maxid);
			vdev.block(maxid, maxid, nsnp - maxid - 1, nsnp - maxid - 1) = covbxy.block(maxid + 1, maxid + 1, nsnp - maxid - 1, nsnp - maxid - 1);            
			
			// get vdev
			VectorXd v1 = VectorXd::Zero(nsnp - 1);
			v1 = v1.array() + 1.0;
		
			vdev = tmp1 + vdev.array() - (v1*tmp3.transpose()).array() - (tmp3*v1.transpose()).array();
			for (int i = 0; i<nsnp - 1; i++)  vdev(i,i) += 1e-8; // in R code  		

            //tmp3 as vardev            
            for(int i=0; i<maxid; i++) tmp3[i]=tmp1+covbxy(i,i)-2*tmp3[i]+ 1e-8;
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy(i,i)-2*tmp3[i-1]+ 1e-8;
			
            //dev as chisq_dev
            
            for(int i=0;i<nsnp-1;i++) dev[i]=dev[i]*dev[i]/tmp3[i];
            
            double sumChisq_dev=0.0;
            for(int i=0;i<nsnp-1;i++)sumChisq_dev+=dev[i];
           
            //using covbxy to store corr_dev
			covbxy.resize(nsnp - 1, nsnp - 1);
			covbxy = vdev.array() / sqrt((vdev.diagonal()*vdev.diagonal().transpose()).array());
           
            // using Eigen Library
            
			SelfAdjointEigenSolver<MatrixXd> es(covbxy);			
            VectorXd lambda;
            lambda=es.eigenvalues();
            
            /*
             EigenSolver<MatrixXd> es(A);
             cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
             cout<< es.eigenvalues().transpose()<<endl;
             */
            
            /*
             MatrixXd D = es.pseudoEigenvalueMatrix();
             cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
            */
            
            pdev= pchisqsum(sumChisq_dev,lambda);
            *snp_num=lambda.size();
         
        }else *snp_num=-9;      
        
        return(pdev);
    }
    
    float bxy_mltheter_so(VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num, double theta)
    {
        VectorXd _bxy;
        VectorXd _sexy;
        VectorXd dev;
        VectorXd tmp3;
        long nsnp=*snp_num;
        int maxid;
        float pdev=-1.0;
        MatrixXd covbxy(nsnp,nsnp);
        MatrixXd vdev(nsnp-1,nsnp-1);
        
        VectorXd bxz2;
        dev.resize(nsnp-1);
        if(nsnp>1)
        {
            _bxy=_byz.array()/_bxz.array();
            bxz2=_bxz.array()*_bxz.array();
            _sexy=(_seyz.array()*_seyz.array()*bxz2.array()+_sexz.array()*_sexz.array()*_byz.array()*_byz.array())/(bxz2.array()*bxz2.array()).sqrt();
            
            maxid=max_abs_id(_zsxz);
            
            for(int j=0;j<maxid;j++) dev[j]=_bxy[maxid]-_bxy[j];
            for(int j=maxid+1;j<nsnp;j++) dev[j-1]=_bxy[maxid]-_bxy[j];
            
            est_cov_bxy_so(covbxy, _byz,_bxz,_seyz,_sexz,_LD_heidi,theta);
            
            double tmp1=covbxy(maxid,maxid);
            tmp3.resize(nsnp-1);
            for(int i=0; i<maxid; i++) tmp3[i]=covbxy(maxid,i);
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=covbxy(maxid,i);
            // vdev as tmp2
            vdev.block(0, 0, maxid, maxid) = covbxy.block(0, 0, maxid, maxid);
            vdev.block(0, maxid, maxid, nsnp - maxid - 1) = covbxy.block(0, maxid + 1, maxid, nsnp - maxid - 1);
            vdev.block(maxid, 0, nsnp - maxid - 1, maxid) = covbxy.block(maxid + 1, 0, nsnp - maxid - 1, maxid);
            vdev.block(maxid, maxid, nsnp - maxid - 1, nsnp - maxid - 1) = covbxy.block(maxid + 1, maxid + 1, nsnp - maxid - 1, nsnp - maxid - 1);
            
            // get vdev
            VectorXd v1 = VectorXd::Zero(nsnp - 1);
            v1 = v1.array() + 1.0;
            
            vdev = tmp1 + vdev.array() - (v1*tmp3.transpose()).array() - (tmp3*v1.transpose()).array();
            for (int i = 0; i<nsnp - 1; i++)  vdev(i,i) += 1e-8; // in R code
            
            //tmp3 as vardev
            for(int i=0; i<maxid; i++) tmp3[i]=tmp1+covbxy(i,i)-2*tmp3[i]+ 1e-8;
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy(i,i)-2*tmp3[i-1]+ 1e-8;
           
            //dev as chisq_dev
            
            for(int i=0;i<nsnp-1;i++) dev[i]=dev[i]*dev[i]/tmp3[i];
           
            double sumChisq_dev=0.0;
            for(int i=0;i<nsnp-1;i++)sumChisq_dev+=dev[i];
            
            //using covbxy to store corr_dev
            covbxy.resize(nsnp - 1, nsnp - 1);
            covbxy = vdev.array() / sqrt((vdev.diagonal()*vdev.diagonal().transpose()).array());
            covbxy.triangularView<Lower>() = covbxy.transpose(); // they make the matrix symmetrical in the Rscript
           
            // using Eigen Library
            SelfAdjointEigenSolver<MatrixXd> es(covbxy);
            VectorXd lambda;
            lambda=es.eigenvalues();
            
            /*
             EigenSolver<MatrixXd> es(A);
             cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
             cout<< es.eigenvalues().transpose()<<endl;
             */
            
            /*
             MatrixXd D = es.pseudoEigenvalueMatrix();
             cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
             */
            
            pdev= pchisqsum(sumChisq_dev,lambda);
            *snp_num=lambda.size();
            
        }else *snp_num=-9;      
        
        return(pdev);
    }

    string dtos(double value)
    {
        stringstream ss;
        ss << scientific<< value;
        // ss << fixed << setprecision(400) << __value;
        return(ss.str());
    }
    string dtosf(double value)
    {
        stringstream ss;        
        ss << fixed << value;
        return(ss.str());
    }
    string itos(int value)
    {
        stringstream ss;
        ss << value;      
        return(ss.str());
    }
    string ltos(long value)
    {
        stringstream ss;
        ss << value;
        return(ss.str());
    }
    
    void free_gwas_data(gwasData* gdata)
    {
        gdata->snpName.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
    }
            
    void mu_func(bInfo* bdata, int j, vector<double> &fac) {
        int i = 0;
		bdata->_dosage_flag = 0;
        double fcount = 0.0, f_buf = 0.0;
        if (bdata->_dosage_flag) {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                    bdata->_mu[bdata->_include[j]] += fac[i] * bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    fcount += fac[i];
                }
            }
        } else {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                    f_buf = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    if (bdata->_allele2[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) f_buf = 2.0 - f_buf;
                    bdata->_mu[bdata->_include[j]] += fac[i] * f_buf;
                    fcount += fac[i];
                }
            }
        }
        
        if (fcount > 0.0)bdata->_mu[bdata->_include[j]] /= fcount;
    }
    

    void calcu_mu(bInfo* bdata, bool ssq_flag)
    {
        int i = 0, j = 0;
        
        vector<double> auto_fac(bdata->_keep.size()), xfac(bdata->_keep.size()), fac(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++)
        {
            auto_fac[i] = 1.0;
            if (bdata->_sex[bdata->_keep[i]] == 1) xfac[i] = 0.5;
            else if (bdata->_sex[bdata->_keep[i]] == 2) xfac[i] = 1.0;
            fac[i] = 0.5;
        }
        
        cout << "Calculating allele frequencies ..." << endl;
        bdata->_mu.clear();
        bdata->_mu.resize(bdata->_snp_num);
        
        #pragma omp parallel for
        for (j = 0; j < bdata->_include.size(); j++)
        {
            if (bdata->_chr[bdata->_include[j]]<(bdata->_autosome_num + 1)) mu_func(bdata, j, auto_fac);
            else if (bdata->_chr[bdata->_include[j]] == (bdata->_autosome_num + 1)) mu_func(bdata,j, xfac);
            else mu_func(bdata, j, fac);
        }
    }

    eigenMatrix reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table=false)
    {
        int N = x.size();
        if (N != y.size() || N < 1) throw ("Error: The lengths of x and y do not match.");
        
        int i = 0;
        double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
        for (i = 0; i < N; i++) {
            x_mu += x[i];
            y_mu += y[i];
        }
        x_mu /= (double) N;
        y_mu /= (double) N;
        for (i = 0; i < N; i++) {
            d_buf = (x[i] - x_mu);
            x_var += d_buf*d_buf;
            d_buf = (y[i] - y_mu);
            y_var += d_buf*d_buf;
        }
        x_var /= (double) (N - 1.0);
        y_var /= (double) (N - 1.0);
        for (i = 0; i < N; i++) cov += (x[i] - x_mu)*(y[i] - y_mu);
        cov /= (double) (N - 1);
        double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, p = 0.0, rsq = 0.0, r = 0.0;
        if (x_var > 0.0) b = cov / x_var;
        a = y_mu - b*x_mu;
        for (i = 0; i < N; i++) {
            d_buf = y[i] - a - b * x[i];
            sse += d_buf*d_buf;
        }
        if (x_var > 0.0) {
            a_se = sqrt((sse / (N - 2.0))*(1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
            b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
        }
        if (x_var > 0.0 && y_var > 0.0) {
            r = cov / sqrt(y_var * x_var);
            rsq = r*r;
        }
        double t = 0.0;
        if (b_se > 0.0) t = fabs(b / b_se);
        p = StatFunc::t_prob(N - 2.0, t, true);
        rst.clear();
        rst.push_back(b);
        rst.push_back(b_se);
        rst.push_back(p);
        rst.push_back(rsq);
        rst.push_back(r);
        
        eigenMatrix reg_sum(3, 3);
        if (table) {
            reg_sum(2, 0) = rsq;
            reg_sum(1, 0) = b;
            reg_sum(1, 1) = b_se;
            reg_sum(1, 2) = p;
            if (a_se > 0.0) t = fabs(a / a_se);
            p = StatFunc::t_prob(N - 2.0, t, true);
            reg_sum(0, 0) = a;
            reg_sum(0, 1) = a_se;
            reg_sum(0, 2) = p;
            return (reg_sum);
        }
        return (reg_sum);
    }    

  
    bool make_XMat(bInfo* bdata, MatrixXf &X)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        
        cout << "Recoding genotypes (individual major mode) ..." << endl;
        bool have_mis = false;
        unsigned long i = 0, j = 0, n = bdata->_keep.size(), m = bdata->_include.size();
        
        X.resize(0,0);
        X.resize(n, m);
        //#pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            if (bdata->_dosage_flag) {
                for (j = 0; j < m; j++) {
                    if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                        if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) X(i,j) = bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                        else X(i,j) = 2.0 - bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    }
                    else {
                        X(i,j) = 1e6;
                        have_mis = true;
                    }
                }
                bdata->_geno_dose[i].clear();
            }
            else {
                for (j = 0; j < bdata->_include.size(); j++) {
                    if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                        if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) X(i,j) = bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]];
                        else X(i,j) = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    }
                    else {
                        X(i,j) = 1e6;
                        have_mis = true;
                    }
                }
            }
        }
        return have_mis;
    }
   
    void makex_eigenVector(bInfo* bdata,int j, VectorXd &x, bool resize, bool minus_2p)
    {
        int i = 0;
        if (resize) x.resize(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) x[i] = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                else x[i] = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
            } else x[i] = bdata->_mu[bdata->_include[j]];
            if (minus_2p) x[i] -= bdata->_mu[bdata->_include[j]];
        }
    }
    
       
    void makeptrx(bInfo* bdata,int bsnpid,int cursnpid, float* X, bool minus_2p)
    {
        int i = 0;
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[bsnpid]] == bdata->_ref_A[bdata->_include[bsnpid]]) X[cursnpid*bdata->_indi_num+i] = (bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]);
                else X[cursnpid*bdata->_indi_num+i] = 2.0 - (bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]);
            } else X[cursnpid*bdata->_indi_num+i] = bdata->_mu[bdata->_include[bsnpid]];
            if (minus_2p) X[cursnpid*bdata->_indi_num+i] -= bdata->_mu[bdata->_include[bsnpid]];
        }
    }
    
    void read_indi_list(string indi_list_file, vector<string> &indi_list)
    {
        ifstream i_indi_list(indi_list_file.c_str());
        if(!i_indi_list) throw("Error: can not open the file ["+indi_list_file+"] to read.");
        string str_buf, id_buf;
        indi_list.clear();
        while(i_indi_list){
            i_indi_list>>str_buf;
            if(i_indi_list.eof()) break;
            id_buf=str_buf+":";
            i_indi_list>>str_buf;
            id_buf+=str_buf;
            indi_list.push_back(id_buf);
            getline(i_indi_list, str_buf);
        }
        i_indi_list.close();
    }
    
    void read_msglist(string msglistfile, vector<string> &msglist, string msg)
    {
        // Read msglist file
        msglist.clear();
        string StrBuf;
        ifstream i_msglist(msglistfile.c_str());
        if(!i_msglist) throw("Error: can not open the file ["+msglistfile+"] to read.");
        cout<<"Reading a list of "<<msg<<" from ["+msglistfile+"]."<<endl;
        while(i_msglist>>StrBuf){
            msglist.push_back(StrBuf);
            getline(i_msglist, StrBuf);
        }
        i_msglist.close();
    }
    void read_snpprblist(string snpprblistfile, vector<string> &prblist,  map<string, string> &prb_snp)
    {
        // Read probe file
        prblist.clear();
        prb_snp.clear();
        FILE* rfile=fopen(snpprblistfile.c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",snpprblistfile.c_str());
            exit(EXIT_FAILURE);
        }
        printf("Reading the SNP - probe list from %s ...\n", snpprblistfile.c_str());
        char Tbuf[MAX_LINE_SIZE];
        int line_idx=0;
        vector<string> strlist;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(strlist.size()!=2)
            {
                printf("ERROR: line %d doesn't have 2 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            string tarsnp=strlist[0];
            string tarprb=strlist[1];
            prb_snp.insert(pair<string, string>(tarprb,tarsnp));
            if(prb_snp.size()==line_idx)
            {
                printf("ERROR: Duplicate probe found : %s.\n", tarprb.c_str());
                exit(EXIT_FAILURE);
            }
            prblist.push_back(tarprb);
            line_idx++;
        }
        fclose(rfile);
        printf("%d SNP - probe pairs have been read from %s.\n",line_idx,snpprblistfile.c_str());        
        
    }

    void read_prb_cojo_snplist(string snpprblistfile, vector<string> &prblist,  map<string, vector<string>> &prb_snp)
    {
        // Read probe file
        prblist.clear();
        prb_snp.clear();
        FILE* rfile=fopen(snpprblistfile.c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",snpprblistfile.c_str());
            exit(EXIT_FAILURE);
        }
        printf("Reading the independent COJO signal list for each probe from %s ...\n", snpprblistfile.c_str());
        char Tbuf[MAX_LINE_SIZE];
        int line_idx=0;
        vector<string> strlist;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            split_string(Tbuf, strlist, " \t\n");
            if(strlist.size()!=2)
            {
                printf("ERROR: line %d doesn't have 2 items.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            string tarprb=strlist[0];
            vector<string> tarsnp;
            split_string(strlist[1],tarsnp);
            prb_snp.insert(pair<string, vector<string>>(tarprb,tarsnp));
            if(prb_snp.size()==line_idx)
            {
                printf("ERROR: Duplicate probe found : %s.\n", tarprb.c_str());
                exit(EXIT_FAILURE);
            }
            prblist.push_back(tarprb);
            line_idx++;
        }
        fclose(rfile);
        printf("The independent COJO signals for %d probes have been read from %s.\n",line_idx,snpprblistfile.c_str());            
    }

    void read_GWAS_cojo_snplist(lociData* ldata, char* lociFileName)
    {
        ifstream lociFile;
        if(!file_read_check(&lociFile, lociFileName))
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     lociFileName, strerror (errno));
            exit (EXIT_FAILURE);
        }
        cout << "Reading the independent GWAS COJO signals list from [" + string(lociFileName) + "]." << endl;        
        ldata->_chr.clear();
        ldata->_snp_name.clear();
        ldata->_bp.clear();
        ldata->_include.clear();
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        lociFile.getline(buf,MAX_LINE_SIZE);// the header
        if(buf[0]=='\0')
        {
            printf("ERROR: the first row of the file %s is empty.\n",lociFileName);
            exit(EXIT_FAILURE);
        }
        vector<string> vs_buf;
        split_string(buf, vs_buf, ", \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[1]!="SNP") {
            printf("ERROR: %s should have headers that start with \"snp\".\n", lociFileName);
            exit(EXIT_FAILURE);
        }
        while(!lociFile.eof())
        {
            lociFile.getline(buf,MAX_LINE_SIZE);
            
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=3) {
                    printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("WARNING: Chromosome is \'NA\' in row %d.\n", lineNum+2);
                    ldata->_chr.push_back(0);
                } else {
                    ldata->_chr.push_back(atof(vs_buf[0].c_str()));
                }
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                ldata->_snp_name.push_back(vs_buf[1]);                                               
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("WARNING: physical position is \'NA\' in row %d.\n", lineNum+2);
                    ldata->_bp.push_back(0);
                } else {
                    ldata->_bp.push_back(atof(vs_buf[2].c_str()));
                }
                ldata->_include.push_back(lineNum);
                lineNum++;
            }
        }
        cout <<"There are "<<ldata->_include.size() << " independent GWAS COJO loci information to be included from [" + string(lociFileName) + "]." << endl;
        lociFile.close();
    }

    void read_pifile(string piFilename, vector<string> &priorsplit)
    {
        // Read prior probabilities file
        priorsplit.clear();
        FILE* rfile=fopen(piFilename.c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",piFilename.c_str());
            exit(EXIT_FAILURE);
        }
        printf("\nReading the estimated prior probabilities from %s ...\n", piFilename.c_str());
        char Tbuf[MAX_LINE_SIZE];
        int line_idx=0;
        vector<string> strlist;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(line_idx == 0 && strlist[0]!="Posteriors") 
            {
                printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check.\n", piFilename.c_str());
                exit(EXIT_FAILURE);
            }
            line_idx++;
            if(line_idx == 2) {
                for(int i=1;i<strlist.size();i++) {
                    priorsplit.push_back(strlist[i]);
                }                
            }
        }
        if(line_idx != 3)
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check.\n", piFilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(rfile);
        printf("There are %d prior probabilities included in %s.\n",priorsplit.size(),piFilename.c_str());
    }
 
    void extract_targets(eqtlInfo* eqtlinfo, string snpprblistfile, map<string, string> &prb_snp)
    {
        vector<string> prblist;
        read_snpprblist( snpprblistfile, prblist,  prb_snp);
        vector<string> raw_problist;
        for(int i=0;i<eqtlinfo->_include.size();i++) raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
        vector<string> common_probes;
        set_intersect(prblist, raw_problist, common_probes);
        eqtlinfo->_include.clear();
        StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+snpprblistfile+"]."<<endl;
        
    }
    
    void update_id_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
    {
        int i=0;
        map<string, int> id_map_buf(id_map);
        for(i=0; i<id_list.size(); i++) id_map_buf.erase(id_list[i]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) id_map.erase(iter->first);
        
        keep.clear();
        for(iter=id_map.begin(); iter!=id_map.end(); iter++) keep.push_back(iter->second);
        stable_sort(keep.begin(), keep.end());
    }
    void update_id_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
    {
        int i = 0;
        for (i = 0; i < id_list.size(); i++) id_map.erase(id_list[i]);
        
        keep.clear();
        map<string, int>::iterator iter;
        for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
        stable_sort(keep.begin(), keep.end());
    }
    

    void keep_indi(bInfo* bdata,string indi_list_file)
    {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        update_id_map_kp(indi_list, bdata->_id_map, bdata->_keep);
        cout<<bdata->_keep.size()<<" individuals are kept from ["+indi_list_file+"]."<<endl;
    }
    
    void remove_indi(bInfo* bdata, string indi_list_file) {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        int prev_size = bdata->_keep.size();
        update_id_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        cout << prev_size - bdata->_keep.size() << " individuals are removed from [" + indi_list_file + "] and there are " << bdata->_keep.size() << " individuals remaining." << endl;
    }
    void extract_region_bp(bInfo* bdata, int chr, int fromkb, int tokb)
    {
        int frombp=fromkb*1000;
        int tobp=tokb*1000;
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr && bdata->_bp[j]<=tobp && bdata->_bp[j]>=frombp) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from SNP BP: "<<fromkb<<"Kb"<<"to SNP BP: "<<tokb<<"Kb."<< endl;
    }
    void extract_snp(bInfo* bdata, int chr)
    {
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from chromosome "<<chr<< "."<<endl;
    }
    void extract_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout<<bdata->_include.size()<<" SNPs are extracted from ["+snplistfile+"]."<<endl;
    }
    void exclude_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
         string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        int prev_size = bdata->_include.size();
        update_id_map_rm(snplist, bdata->_snp_name_map, bdata->_include);
        cout << prev_size - bdata->_include.size() << " SNPs are excluded from [" + snplistfile + "] and there are " << bdata->_include.size() << " SNPs remaining." << endl;
    }
   
    void extract_eqtl_by_chr(eqtlInfo* eqtlinfo, int snpchr)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr ) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from chromosome [" + atos(snpchr) + "]." << endl;
    }
    void extract_epi_by_chr(eqtlInfo* eqtlinfo, int prbchr)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr ) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from chromosome [" + atos(prbchr) + "]." << endl;
    }
    void extract_ldata_by_chr(lociData* lociData, int prbchr)
    {
        vector<int> newIcld;
        for(int i=0;i<lociData->_include.size();i++)
        {
            int tmpint=lociData->_include[i];
            if(lociData->_chr[tmpint]==prbchr) newIcld.push_back(tmpint);
        }
        lociData->_include.clear();
        lociData->_include=newIcld;
        cout << lociData->_include.size() << " GWAS loci are extracted from chromosome [" + atos(prbchr) + "]." << endl;
    }
	void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
	{
		vector<string> snplist;
		string msg = "SNPs";
		read_msglist(snplstName, snplist, msg);
        if(eqtlinfo->_esi_include.size()==eqtlinfo->_snpNum)
        {
            eqtlinfo->_esi_include.clear();
            StrFunc::match_only(snplist, eqtlinfo->_esi_rs, eqtlinfo->_esi_include);
            stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());
        }
        else
        {
            vector<int> icld_step1;
            StrFunc::match_only(snplist, eqtlinfo->_esi_rs, icld_step1);
            
            vector<int> common_probes;
            set_intersect(icld_step1, eqtlinfo->_esi_include, common_probes);
            eqtlinfo->_esi_include.clear();
            eqtlinfo->_esi_include=common_probes;
            stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());
        }
        
		cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from [" + snplstName + "]." << endl;
	}
   
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snporprb, int Wind, string msg)
    {
        string logstr;
        
        int bp=-9;
        int chr=-9;
        if(msg=="SNP")
        {
            long idx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), snporprb)-eqtlinfo->_esi_rs.begin();
            if(idx==eqtlinfo->_esi_rs.size())
            {
                logstr="ERROR: Can't find the SNP "+snporprb+". Please check.\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            bp=eqtlinfo->_esi_bp[idx];
            chr=eqtlinfo->_esi_chr[idx];
        }
        else if(msg=="probe")
        {
            long idx=-9;
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                if(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]==snporprb) idx=eqtlinfo->_include[i];
            }
            if(idx==-9)
            {
                logstr="ERROR: Can't find the probe "+snporprb+". Please check.\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            bp=eqtlinfo->_epi_bp[idx];
            chr=eqtlinfo->_epi_chr[idx];
        }
        
        int upbound=bp+Wind*1000;
        int tmpint=bp-Wind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==chr && eqtlinfo->_esi_bp[tmpint]>=lowbound && eqtlinfo->_esi_bp[tmpint]<=upbound) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from the region: " +atos(Wind)+" Kb around [" + snporprb + "]." << endl;
    }
    void extract_eqtl_single_snp(eqtlInfo* eqtlinfo, string snprs)
    {
        string logstr;
        long idx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), snprs)-eqtlinfo->_esi_rs.begin();
        if(idx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find the SNP "+snprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include.push_back((int)idx);
        cout << snprs << " is extracted. " << endl;
    }
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string fromsnprs, string tosnprs)
    {
        string logstr;
        long fromidx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), fromsnprs)-eqtlinfo->_esi_rs.begin();
        if(fromidx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find the SNP "+fromsnprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int fromsnpbp=eqtlinfo->_esi_bp[fromidx];
        int snpchr=eqtlinfo->_esi_chr[fromidx];
        
        long toidx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), tosnprs)-eqtlinfo->_esi_rs.begin();
        if(toidx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find the SNP "+tosnprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int tosnpbp=eqtlinfo->_esi_bp[toidx];
        int tosnpchr=eqtlinfo->_esi_chr[toidx];
        if(tosnpchr != snpchr)
        {
            logstr="ERROR: SNP "+fromsnprs+" and SNP "+tosnprs +" are not one the same chromosome. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
     
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr && eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from SNP " +fromsnprs+" to SNP " + tosnprs + "." << endl;
    }
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, int chr, int fromsnpkb, int tosnpkb)
    {
        int fromsnpbp=fromsnpkb*1000;
        int tosnpbp=tosnpkb*1000;
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if( eqtlinfo->_esi_chr[tmpint]==chr &&eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from SNP BP: " +atos(fromsnpkb)+" Kb to SNP BP: " + atos(tosnpkb) + " Kb." << endl;
    }



    void exclude_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist;
        vector<string> mapstr;
         vector<int> tmp;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        int pre_num=eqtlinfo->_esi_include.size();
        mapstr.resize(pre_num);
        tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mapstr[i]=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]];
            tmp[i]=eqtlinfo->_esi_include[i];
        }
        
        eqtlinfo->_esi_include.clear();
        StrFunc::set_complement(snplist, mapstr, tmp, eqtlinfo->_esi_include); //sorted
        stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());       
        cout << pre_num-eqtlinfo->_esi_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << eqtlinfo->_esi_include.size() << " SNPs remaining." << endl;
    }
    
    void extract_gwas_snp(gwasData* gdata, string snplstName)
    {
        vector<string> snplist;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        gdata->_include.clear();
        StrFunc::match_only(snplist, gdata->snpName, gdata->_include);
        stable_sort(gdata->_include.begin(), gdata->_include.end());
        cout << gdata->_include.size() << " SNPs are extracted from [" + snplstName + "]." << endl;
    }
    
    void exclude_gwas_snp(gwasData* gdata, string snplstName)
    {
        vector<string> snplist;
        vector<string> mapstr;
         vector<int> tmp;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        int pre_num=gdata->_include.size();
        mapstr.resize(pre_num);
         tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mapstr[i]=gdata->snpName[gdata->_include[i]];
            tmp[i]=gdata->_include[i];
        }
        
        gdata->_include.clear();
        StrFunc::set_complement(snplist, mapstr, tmp,gdata->_include); //sorted
        stable_sort(gdata->_include.begin(), gdata->_include.end());
        cout << pre_num-gdata->_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << gdata->_include.size() << " SNPs remaining." << endl;
    }

    void read_epistartend(eqtlInfo* eqtlinfo,char* prbseqregion)
    {
        vector<string> problist;
        vector<int> start;
        vector<int> end;
        eqtlinfo->_epi_start.resize(eqtlinfo->_probNum);
        eqtlinfo->_epi_end.resize(eqtlinfo->_probNum);
        for(int i=0;i<eqtlinfo->_probNum;i++)
        {
            eqtlinfo->_epi_start[i]=-9;
            eqtlinfo->_epi_end[i]=-9;
        }
        char Tbuf[MAX_LINE_SIZE];
        FILE* epifile=fopen(prbseqregion,"r");
        if(!epifile) {
            printf("File %s open failed.\n",prbseqregion);
            exit(EXIT_FAILURE);
        }
        printf("Reading the information of the probe hybridization region from %s ...\n", prbseqregion);
        vector<string> strlist;
        int  line_idx = 0, hit = 0;
        map<string, int>::iterator iter;
        while(fgets(Tbuf, MAX_LINE_SIZE, epifile))
        {
            split_string(Tbuf, strlist, ", \t\n");
             if(strlist.size()<4)
             {
                 printf("ERROR: line %u has more than 5 items.\n", line_idx);
                 exit(EXIT_FAILURE);
             }
            string prbname=strlist[1];
            iter=eqtlinfo->_probe_name_map.find(prbname);
            if(iter!=eqtlinfo->_probe_name_map.end())
            {
                int sta=atoi(strlist[2].c_str());
                int en=atoi(strlist[3].c_str());
                if(sta<0)
                {
                    printf("ERROR: the start position %d should not be negative.\n", sta);
                    exit(EXIT_FAILURE);
                }
                if(en<0)
                {
                    printf("ERROR: the end position %d should not be negative.\n", en);
                    exit(EXIT_FAILURE);
                }
                if(sta>en) {
                    printf("ERROR: the start position %d should not be larger than the end position %d.\n", sta,en);
                    exit(EXIT_FAILURE);
                }
                eqtlinfo->_epi_start[iter->second]=sta;
                eqtlinfo->_epi_end[iter->second]=en;
                hit++;
            }
            line_idx++;
        }
        printf("%d rows are included from %s and %d of them are matched with %ld probes in the .epi file.\n", line_idx, prbseqregion, hit, eqtlinfo->_include.size());
        fclose(epifile);
    }

    
    void extract_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        vector<string> raw_problist;
        for(int i=0;i<eqtlinfo->_include.size();i++) raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
        vector<string> common_probes;
        set_intersect(problist, raw_problist, common_probes);
        eqtlinfo->_include.clear();
		StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+problstName+"]."<<endl;
    }
    
    void extract_prob_by_gene(eqtlInfo* eqtlinfo, string genelistName)
    {
        vector<string> genelist;
        string msg="genes";
        read_msglist(genelistName, genelist,msg);
        vector<string> raw_problist;
        vector<string> raw_genelist;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
            raw_genelist.push_back(eqtlinfo->_epi_gene[eqtlinfo->_include[i]]);
        }
        vector<int> idx;
        for(int i=0;i<genelist.size();i++)
        {
            string tmpname1=genelist[i];
            for(int j=0;j<raw_genelist.size();j++)
            {
                string tmpname2=raw_genelist[j];
                if(tmpname1==tmpname2)  idx.push_back(j);
                else
                {
                    vector<string> substrs;
                    int tmpnum=split_string(tmpname2,substrs);
                    if(tmpnum>1)
                        for(int k=0;k<tmpnum;k++)
                            if(tmpname1==substrs[k])
                            {
                                idx.push_back(j);
                                break;
                            }
                }
            }
            out2:;
        }
        vector<string> common_probes;
        for(int i=0;i<idx.size();i++) common_probes.push_back(raw_problist[idx[i]]);
        eqtlinfo->_include.clear();
        StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+genelistName+"]."<<endl;
    }
    
    void extract_prob(eqtlInfo* eqtlinfo, string prbname, int prbWind)
    {
        string logstr;
        long idx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), prbname)-eqtlinfo->_epi_prbID.begin();
        if(idx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find probe "+prbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int prbbp=eqtlinfo->_epi_bp[idx];
        int prbchr=eqtlinfo->_epi_chr[idx];
        int upbound=prbbp+prbWind*1000;
        int tmpint=prbbp-prbWind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=lowbound && eqtlinfo->_epi_bp[tmpint]<=upbound) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from the region: " +atos(prbWind)+" Kb around [" + prbname + "]." << endl;
    }
    void extract_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname, bool prtscr)
    {
        string logstr;
        int idx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), prbname)-eqtlinfo->_epi_prbID.begin();
        if(idx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find probe "+prbname+". Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }

        eqtlinfo->_include.clear();
        eqtlinfo->_include.push_back(idx);
        if(prtscr) cout << prbname << " is extracted. " << endl;
    }
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, string fromprbname, string toprbname)
    {
        string logstr;
        long fromidx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), fromprbname)-eqtlinfo->_epi_prbID.begin();
        if(fromidx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find SNP "+fromprbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int fromprbbp=eqtlinfo->_epi_bp[fromidx];
        int prbchr=eqtlinfo->_epi_chr[fromidx];
        
        long toidx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), toprbname)-eqtlinfo->_epi_prbID.begin();
        if(toidx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find SNP "+toprbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int toprbbp=eqtlinfo->_epi_bp[toidx];
        int toprbchr=eqtlinfo->_epi_chr[toidx];
        if(toprbchr != prbchr)
        {
            logstr="ERROR: probe "+fromprbname+" and probe "+toprbname +" are not on the same chromosome. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from probe : " +fromprbname+" to probe " + toprbname + "." << endl;
    }
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, int chr, int fromprbkb, int toprbkb)
    {
        int fromprbbp=fromprbkb*1000;
        int toprbbp=toprbkb*1000;
        
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if( eqtlinfo->_epi_chr[tmpint]==chr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from probe BP: " +atos(fromprbkb)+"Kb to probe BP: " + atos(toprbkb) + "Kb." << endl;
    }
    void extract_prob_by_single_gene(eqtlInfo* eqtlinfo, string genename)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if( eqtlinfo->_epi_gene[tmpint]==genename) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted that are mapping to gene: " +genename+ "." << endl;
    }

    void exclude_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        vector<string> mappro;
        vector<int> tmp;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        int pre_num=eqtlinfo->_include.size();
        mappro.resize(pre_num);
        tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mappro[i]=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
            tmp[i]=eqtlinfo->_include[i];
        }
       
        eqtlinfo->_include.clear();
        StrFunc::set_complement(problist, mappro, tmp, eqtlinfo->_include);
        
        cout<<pre_num-eqtlinfo->_include.size()<<" probes are excluded from ["+problstName+"]and there are "<<eqtlinfo->_include.size()<<" probes remaining."<<endl;
    }
    
    void exclude_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname)
    {
        string logstr;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int idx=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_prbID[idx] != prbname) newIcld.push_back(idx);
        }
       
        if(newIcld.size()==eqtlinfo->_include.size())
        {
            logstr="WARNING: Can't find probe "+prbname+" in the dataset. Nothing to exclude.\n";
            fputs(logstr.c_str(),stdout);
            
        }else
        {
            eqtlinfo->_include.clear();
            eqtlinfo->_include=newIcld;
            cout << prbname << " is excluded from the dataset. " << endl;
        }
        
    }

    void filter_snp_maf(bInfo* bdata,double maf)
    {
        if(bdata->_mu.empty()) calcu_mu(bdata);
        
        cout<<"Pruning SNPs with MAF > "<<maf<<" ..."<<endl;
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        int prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mu[iter->second]*0.5;
            if(fbuf<=maf || (1.0-fbuf)<=maf) continue;
            bdata->_snp_name_map.insert(*iter);
            bdata->_include.push_back(iter->second);
        }
        if(bdata->_include.size()==0) throw("Error: No SNP is retained for analysis.");
        else{
            stable_sort(bdata->_include.begin(), bdata->_include.end());
            cout<<"After pruning SNPs with MAF > "<<maf<<", there are "<<bdata->_include.size()<<" SNPs ("<<prev_size-bdata->_include.size()<<" SNPs with MAF < "<<maf<<")."<<endl;
        }
        
    }
    
    void allele_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).";
        cout<<logstr<<endl;
        
        vector<string> bsnp;
        vector<string> essnp;
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<bdata->_include.size();i++) bsnp.push_back(bdata->_snp_name[bdata->_include[i]]);
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(bsnp, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common between reference data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            StrFunc::match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<edId.size();i++) edId[i]=esdata->_esi_include[edId[i]];
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no SNPs in common  between reference data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        StrFunc::match(slctSNPs, bdata->_snp_name, bdId);
        StrFunc::match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                   
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=gdata->allele_1[gdId[i]];
                    gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                    gdata->allele_2[gdId[i]] = tmpch;
                    if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                    
                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }

            
        }
        
        logstr=itos(bdata->_include.size())+" SNPs are included after allele checking. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    void allele_check_opt(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        map<string, int> allel_map;
        map<string, int>::iterator iter, iter1, iter2;
        vector<int> bin,ein,gin;
        for(int i=0;i<gdata->_include.size();i++) //no map in gwas data
        {
            printf("%3.0f%%\r", 100.0*i/gdata->_include.size());
            fflush(stdout);
            int rid=gdata->_include[i];
            string grs=gdata->snpName[rid];
            allel_map.clear();
            string a1, a2, ga1, ga2, ea1, ea2;
            ga1 = gdata->allele_1[rid];
            ga2 = gdata->allele_2[rid];
            allel_map.insert(pair<string,int>(ga1,0));
            allel_map.insert(pair<string,int>(ga2,1));
            bool hitall=false;
            iter1 = esdata-> _snp_name_map.find(grs);
            if (iter1 != esdata->_snp_name_map.end()) {
                iter2 = bdata-> _snp_name_map.find(grs);
                if (iter2 != bdata->_snp_name_map.end()) {
                    ea1=esdata->_esi_allele1[iter1->second];
                    ea2=esdata->_esi_allele2[iter1->second];
                    a1=bdata->_allele1[iter2->second];
                    a2=bdata->_allele2[iter2->second];
                    allel_map.insert(pair<string,int>(ea1,allel_map.size()));
                    allel_map.insert(pair<string,int>(ea2,allel_map.size()));
                    allel_map.insert(pair<string,int>(a1,allel_map.size()));
                    allel_map.insert(pair<string,int>(a2,allel_map.size()));
                    if(allel_map.size()>2) {
                        hitall=false;
                    } else {
                        hitall=true;
                        // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
                        if(ea1 == a1 &&  ea2 == a2)
                        {
                            if( ea1 == ga1 && ea2 == ga2)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                            }
                            else if(ea1 == ga2 && ea2 == ga1)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                gdata->byz[rid]=-1.0*gdata->byz[rid];
                                
                            }
                        }
                        else if(ea1 == a2 &&  ea2 == a1)
                        {
                            
                            if( ea1 == ga1 && ea2 == ga2)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                string tmpch=bdata->_ref_A[iter2->second];
                                bdata->_ref_A[iter2->second]=bdata->_other_A[iter2->second];
                                bdata->_other_A[iter2->second]=tmpch;
                                
                            }
                            else if(ea1 == ga2 && ea2 == ga1)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                gdata->byz[rid]=-1.0*gdata->byz[rid];
                                string tmpch=bdata->_ref_A[iter2->second];
                                bdata->_ref_A[iter2->second]=bdata->_other_A[iter2->second];
                                bdata->_other_A[iter2->second]=tmpch;
                                
                            }
                        }
                        
                    }
                }
            }
        }

        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        esdata->_esi_include.swap(ein);
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    bool allele_check_(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        
        vector<string> bsnp;
        vector<string> essnp;
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<bdata->_include.size();i++) bsnp.push_back(bdata->_snp_name[bdata->_include[i]]);
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(bsnp, gdata->snpName, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            StrFunc::match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) edId[i]=esdata->_esi_include[edId[i]];
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        StrFunc::match(slctSNPs, bdata->_snp_name, bdId);
        StrFunc::match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }
            
            
        }
        
        logstr=itos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
        
        return true;
        
    }
    
    void allele_check(gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data and the eQTL summary data).\n ";
        cout<<logstr<<endl;
       
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(essnp, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }
        
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        cmmnSNPs.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string ga1, ga2, ea1, ea2;
           
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ga1 && ea2 == ga2)
            {
               
                gdata->_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata->_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                
                gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                string tmpch=gdata->allele_1[gdId[i]];
                gdata->allele_1[gdId[i]] = gdata->allele_2[gdId[i]];
                gdata->allele_2[gdId[i]] = tmpch;
                if(gdata->freq[gdId[i]]>0) gdata->freq[gdId[i]] = 1 - gdata->freq[gdId[i]];
                
            }
        }
        printf("%ld SNPs are included after allele check.\n",gdata->_include.size());
    
        
    }
    
    void allele_check(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the outcome eQTL summary data and the exposure eQTL summary data ).\n ";
        cout<<logstr<<endl;
        vector<string> etsnp;
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum || etrait->_esi_include.size()<etrait->_snpNum)
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            for(int i=0;i<etrait->_esi_include.size();i++) etsnp.push_back(etrait->_esi_rs[etrait->_esi_include[i]]);
            StrFunc::match_only(essnp, etsnp, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etsnp[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, etrait->_esi_rs, gdId);
            if(gdId.empty()) throw("Error: no SNPs in common.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etrait->_esi_rs[gdId[i]]);
        }        
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        gdId.clear();
        StrFunc::match(slctSNPs, etrait->_esi_rs, gdId);
        cmmnSNPs.clear();
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        
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
                
                etrait->_esi_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
            }
            else if(ea1 == ta2 && ea2 == ta1)
            {
                etrait->_esi_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                string tmpch=etrait->_esi_allele1[gdId[i]];
                etrait->_esi_allele1[gdId[i]] = etrait->_esi_allele2[gdId[i]];
                etrait->_esi_allele2[gdId[i]] = tmpch;
                if(etrait->_esi_freq[gdId[i]] >0) etrait->_esi_freq[gdId[i]] = 1 - etrait->_esi_freq[gdId[i]];
                
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
        }
        printf("%ld SNPs are included after allele check. \n",etrait->_esi_include.size());
    }
    
    void allele_check(gwasData* gdata1, gwasData* gdata2)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles between GWAS summary dataset and eQTL summary dataset. ";
        cout<<logstr<<endl;
        
        StrFunc::match_only(gdata2->snpName, gdata1->snpName, gdId);
        if(gdId.empty()) throw("Error: no SNPs in common.");
        
        for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata1->snpName[gdId[i]]);
        
        //alleles check
        StrFunc::match(slctSNPs, gdata2->snpName, edId);
        cmmnSNPs.clear();
        gdata1->_include.clear();
        gdata2->_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string ga1, ga2, ea1, ea2;
            
            ga1 = gdata1->allele_1[gdId[i]];
            ga2 = gdata1->allele_2[gdId[i]];
            ea1 = gdata2->allele_1[edId[i]];
            ea2 = gdata2->allele_2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ga1 && ea2 == ga2)
            {
                
                gdata1->_include.push_back(gdId[i]);
                gdata2->_include.push_back(edId[i]);
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata1->_include.push_back(gdId[i]);
                gdata2->_include.push_back(edId[i]);
                
                gdata1->byz[gdId[i]]=-gdata1->byz[gdId[i]];
            }
        }
        logstr=atos(gdata1->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
    }
    void allele_check(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).\n";
        cout<<logstr<<endl;
        
        vector<string> bsnp(bdata->_include.size());
        vector<string> etsnp(etrait->_esi_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum || etrait->_esi_include.size()<etrait->_snpNum )
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<etrait->_esi_include.size();i++)
                etsnp[i]=etrait->_esi_rs[etrait->_esi_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(bsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }else
        {
            match_only(bdata->_snp_name, etrait->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etrait->_esi_rs[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        match(slctSNPs, bdata->_snp_name, bdId);
        match(slctSNPs, etrait->_esi_rs, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            fflush(stdout);
            string a1, a2, ta1, ta2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ta1 = etrait->_esi_allele1[gdId[i]];
            ta2 = etrait->_esi_allele2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
         
            // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ta1 && ea2 == ta2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ta2 && ea2 == ta1)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=etrait->_esi_allele1[gdId[i]];
                    etrait->_esi_allele1[gdId[i]] = etrait->_esi_allele2[gdId[i]];
                    etrait->_esi_allele2[gdId[i]] = tmpch;
                    if(etrait->_esi_freq[gdId[i]] >0) etrait->_esi_freq[gdId[i]] = 1 - etrait->_esi_freq[gdId[i]];
                    
                    if(etrait->_val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<etrait->_rowid.size();j++)
                        {
                            if(etrait->_rowid[j]==gdId[i])
                            {
                                count++;
                                if(count & 1)
                                    etrait->_val[j]=-1*etrait->_val[j];
                            }
                        }
                        
                        /*
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_probNum;j++)
                        {
                            uint64_t beta_start=etrait->_cols[j<<1];
                            uint64_t se_start=etrait->_cols[1+(j<<1)];
                            uint64_t numsnps=se_start-beta_start;
                            for(int k=0;k<numsnps;k++)
                            {
                                uint32_t ge_rowid=etrait->_rowid[beta_start+k];
                                if(ge_rowid==gdId[i])
                                {
                                    etrait->_val[beta_start+k]=-1*etrait->_val[beta_start+k];
                                    break; //the rowid in a probe is unique.
                                }
                            }
                        }
                         */
                    }
                    else
                    {
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_include.size();j++)
                            if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 )
                                etrait->_bxz[j][gdId[i]]=-1*etrait->_bxz[j][gdId[i]];
                    }
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ta1 && ea2 == ta2)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ta2 && ea2 == ta1)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=etrait->_esi_allele1[gdId[i]];
                    etrait->_esi_allele1[gdId[i]] = etrait->_esi_allele2[gdId[i]];
                    etrait->_esi_allele2[gdId[i]] = tmpch;
                    if(etrait->_esi_freq[gdId[i]] >0) etrait->_esi_freq[gdId[i]] = 1 - etrait->_esi_freq[gdId[i]];
                    if(etrait->_val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<etrait->_rowid.size();j++)
                        {
                            if(etrait->_rowid[j]==gdId[i])
                            {
                                count++;
                                if(count & 1) etrait->_val[j]=-etrait->_val[j];
                            }
                        }
                        
                         
                        /*
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_probNum;j++)
                        {
                            uint64_t beta_start=etrait->_cols[j<<1];
                            uint64_t se_start=etrait->_cols[1+(j<<1)];
                            uint64_t numsnps=se_start-beta_start;
                            for(int k=0;k<numsnps;k++)
                            {
                                uint32_t ge_rowid=etrait->_rowid[beta_start+k];
                                if(ge_rowid==gdId[i])
                                {
                                    etrait->_val[beta_start+k]=-1*etrait->_val[beta_start+k];
                                    break; //the rowid in a probe is unique.
                                }
                            }
                        }
                        */

                    }
                    else
                    {
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_include.size();j++)
                            if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 )
                                etrait->_bxz[j][gdId[i]]=-1*etrait->_bxz[j][gdId[i]];
                    }
                    tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }
            
            
        }
        
        logstr=atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

     void allele_check_opt(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata)
    {/*
        string logstr="Checking the consistency of SNP alleles among GWAS summary data, eQTL summary data and reference data.\n ";
        cout<<logstr<<endl;
        map<string, int> allel_map;
        map<string, int>::iterator iter, iter1, iter2;
        vector<int> bin,ein,gin;
        for(int i=0;i<gdata->_include.size();i++) //no map in gwas data
        {
            printf("%3.0f%%\r", 100.0*i/gdata->_include.size());
            fflush(stdout);
            int rid=gdata->_include[i];
            string grs=gdata->snpName[rid];
            allel_map.clear();
            string a1, a2, ga1, ga2, ea1, ea2;
            ga1 = gdata->allele_1[rid];
            ga2 = gdata->allele_2[rid];
            allel_map.insert(pair<string,int>(ga1,0));
            allel_map.insert(pair<string,int>(ga2,1));
            bool hitall=false;
            iter1 = esdata-> _snp_name_map.find(grs);
            if (iter1 != esdata->_snp_name_map.end()) {
                iter2 = bdata-> _snp_name_map.find(grs);
                if (iter2 != bdata->_snp_name_map.end()) {
                    ea1=esdata->_esi_allele1[iter1->second];
                    ea2=esdata->_esi_allele2[iter1->second];
                    a1=bdata->_allele1[iter2->second];
                    a2=bdata->_allele2[iter2->second];
                    allel_map.insert(pair<string,int>(ea1,allel_map.size()));
                    allel_map.insert(pair<string,int>(ea2,allel_map.size()));
                    allel_map.insert(pair<string,int>(a1,allel_map.size()));
                    allel_map.insert(pair<string,int>(a2,allel_map.size()));
                    if(allel_map.size()>2) {
                        hitall=false;
                    } else {
                        hitall=true;
                        // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
                        if(ea1 == a1 &&  ea2 == a2)
                        {
                            if( ea1 == ga1 && ea2 == ga2)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                            }
                            else if(ea1 == ga2 && ea2 == ga1)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                gdata->byz[rid]=-1.0*gdata->byz[rid];
                                
                            }
                        }
                        else if(ea1 == a2 &&  ea2 == a1)
                        {
                            
                            if( ea1 == ga1 && ea2 == ga2)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                string tmpch=bdata->_ref_A[iter2->second];
                                bdata->_ref_A[iter2->second]=bdata->_other_A[iter2->second];
                                bdata->_other_A[iter2->second]=tmpch;
                                
                            }
                            else if(ea1 == ga2 && ea2 == ga1)
                            {
                                bin.push_back(iter2->second);
                                gin.push_back(rid);
                                ein.push_back(iter1->second);
                                
                                gdata->byz[rid]=-1.0*gdata->byz[rid];
                                string tmpch=bdata->_ref_A[iter2->second];
                                bdata->_ref_A[iter2->second]=bdata->_other_A[iter2->second];
                                bdata->_other_A[iter2->second]=tmpch;
                                
                            }
                        }
                        
                    }
                }
            }
        }
        
        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        esdata->_esi_include.swap(ein);
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
      */
    }
    void update_gwas(gwasData* gdata) {
        
        bool hasBP=false;
        if(gdata->snpBp.size()>0) hasBP=true;
        vector<int> snpBp;
        if(hasBP) snpBp.resize(gdata->_include.size());
        vector<string> snpName(gdata->_include.size());
        vector<string> allele_1(gdata->_include.size());
        vector<string> allele_2(gdata->_include.size());
        vector<double> freq(gdata->_include.size());
        vector<double> byz(gdata->_include.size());
        vector<double> seyz(gdata->_include.size());
        vector<double> pvalue(gdata->_include.size());
        vector<uint32_t> splSize(gdata->_include.size());
        
        gdata->snpNum=gdata->_include.size();
        for(int i=0;i<gdata->_include.size();i++ )
        {
            snpName[i]=gdata->snpName[gdata->_include[i]];
            allele_1[i]=gdata->allele_1[gdata->_include[i]];
            allele_2[i]=gdata->allele_2[gdata->_include[i]];
            freq[i]=gdata->freq[gdata->_include[i]];
            byz[i]=gdata->byz[gdata->_include[i]];
            seyz[i]=gdata->seyz[gdata->_include[i]];
            pvalue[i]=gdata->pvalue[gdata->_include[i]];
            splSize[i]=gdata->splSize[gdata->_include[i]];
            if(hasBP) snpBp[i]=gdata->snpBp[gdata->_include[i]];
        }
        
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        
        gdata->snpName.swap(snpName);
        gdata->allele_1.swap(allele_1);
        gdata->allele_2.swap(allele_2);
        gdata->freq.swap(freq);
        gdata->byz.swap(byz);
        gdata->seyz.swap(seyz);
        gdata->pvalue.swap(pvalue);
        gdata->splSize.swap(splSize);
        
        if(hasBP) gdata->snpBp=snpBp;
        gdata->_snp_name_map.clear();
        for(int i=0;i<gdata->snpNum;i++){
            gdata->_include[i]=i;
            gdata->_snp_name_map.insert(pair<string,int>(gdata->snpName[i],i));
        }
        
    }
    
    void make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool minus_2p) {
        // Eigen is column-major by default. here row of X is individual, column of X is SNP.
        uint64_t snpNum=snpids.size();
        X.resize(bdata->_keep.size(),snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            uint32_t snpid=snpids[i];
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[snpid]] == bdata->_ref_A[bdata->_include[snpid]]) X(j,i)= bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]]);
                } else X(j,i) = bdata->_mu[bdata->_include[snpid]];
                if (minus_2p) X(j,i) -= bdata->_mu[bdata->_include[snpid]];
                
            }
        }
    }
    void cor_calc(MatrixXd &LD, MatrixXd &X)
    {
        long size=X.cols();
        long n=X.rows();
        VectorXd tmpX(size);
        VectorXd tmpX2(size);
        #pragma omp parallel for
        for(int i=0;i<size;i++){
            tmpX[i]=X.col(i).sum();
            tmpX2[i]=X.col(i).dot(X.col(i));
        }
        LD.noalias()=X.transpose()*X;
        VectorXd tmpXX=(sqrt(tmpX2.array()*n-tmpX.array()*tmpX.array())).matrix();
        LD = (LD*n-tmpX*tmpX.transpose()).array()/ (tmpXX*tmpXX.transpose()).array();
        
    }
    void update_geIndx(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(gdata->_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
        }
        gdata->_include.clear();
        esdata->_esi_include.clear();
        gdata->_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
    }
    double freq_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata, double &freqthresh, double &percenthresh)
    {
        printf("Checking the consistency of allele frequency of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).\n");
        bool success=true;
        long failcount=0, snpnum=bdata->_include.size();
        vector<int> pasbid, pasgid, paseid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tReferencePanel\tGWAS\teQTL\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<bdata->_include.size();i++)
        {
            double bfreq=bdata->_mu[bdata->_include[i]]/2;
            double gfreq=gdata->freq[gdata->_include[i]];
            double efreq=esdata->_esi_freq[esdata->_esi_include[i]];
            if(gfreq<0 && efreq>0)
            {
                if(fabs(bfreq-efreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else if(gfreq>0 && efreq<0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else if (gfreq>0 && efreq>0)
            {
                if( fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-efreq)>freqthresh || fabs(gfreq-efreq)>freqthresh )
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else
            {
                pasbid.push_back(bdata->_include[i]);
                pasgid.push_back(gdata->_include[i]);
                paseid.push_back(esdata->_esi_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your GWAS summmay data and your .esi file.\n");
            exit(EXIT_FAILURE);
        }
        bdata->_include = pasbid;
        gdata->_include = pasgid;
        esdata->_esi_include = paseid;

        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
       
        fclose(tmpfile);
        return prop;
    }
    double freq_check( gwasData* gdata, eqtlInfo* esdata, double &freqthresh, double &percenthresh)
    {
        printf("Checking the consistency of allele frequency of each SNP between pairwise data sets (including the GWAS summary data and the eQTL summary data).\n");
        bool success=true;
        long failcount=0, snpnum=esdata->_esi_include.size();
        vector<int>  pasgid, paseid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tGWAS\teQTL\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<esdata->_esi_include.size();i++)
        {
            double gfreq=gdata->freq[gdata->_include[i]];
            double efreq=esdata->_esi_freq[esdata->_esi_include[i]];
            if (gfreq>0 && efreq>0)
            {
                if( fabs(gfreq-efreq)>freqthresh )
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasgid.push_back(gdata->_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else
            {
                pasgid.push_back(gdata->_include[i]);
                paseid.push_back(esdata->_esi_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your GWAS summmay data and your .esi file.\n");
            exit(EXIT_FAILURE);
        }
        gdata->_include = pasgid;
        esdata->_esi_include = paseid;
        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        
        fclose(tmpfile);
        return prop;
    }
    double freq_check(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata, double &freqthresh, double &percenthresh)
    {
        printf("Checking the consistency of allele frequency of each SNP between pairwise data sets (including the outcome eQTL summary data, the exposure eQTL summary data and the LD reference data).\n");
        bool success=true;
        long failcount=0, snpnum=bdata->_include.size();
        vector<int> pasbid, pasgid, paseid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tReferencePanel\tGWAS\teQTL\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<bdata->_include.size();i++)
        {
            double bfreq=bdata->_mu[bdata->_include[i]]/2;
            double gfreq=etrait->_esi_freq[etrait->_esi_include[i]];
            double efreq=esdata->_esi_freq[esdata->_esi_include[i]];
            if(gfreq<0 && efreq>0)
            {
                if(fabs(bfreq-efreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(etrait->_esi_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else if(gfreq>0 && efreq<0)
            {
                
                if(fabs(bfreq-gfreq)>freqthresh)
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(etrait->_esi_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else if (gfreq>0 && efreq>0)
            {
                if( fabs(bfreq-gfreq)>freqthresh || fabs(bfreq-efreq)>freqthresh || fabs(gfreq-efreq)>freqthresh )
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]]+ '\t' + atos(bfreq) + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasbid.push_back(bdata->_include[i]);
                    pasgid.push_back(etrait->_esi_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else
            {
                pasbid.push_back(bdata->_include[i]);
                pasgid.push_back(etrait->_esi_include[i]);
                paseid.push_back(esdata->_esi_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your GWAS summmay data and your .esi file.\n");
            exit(EXIT_FAILURE);
        }
        bdata->_include = pasbid;
        etrait->_esi_include = pasgid;
        esdata->_esi_include = paseid;
        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        
        fclose(tmpfile);
        return prop;
    }
    double freq_check( eqtlInfo* etrait, eqtlInfo* esdata, double &freqthresh, double &percenthresh)
    {
        printf("Checking the consistency of allele frequency of each SNP between pairwise data sets (including the outcome eQTL summary data and the exposure eQTL summary data).\n");
        bool success=true;
        long failcount=0, snpnum=esdata->_esi_include.size();
        vector<int> pasbid, pasgid, paseid;
        string fname=string(outFileName) + ".snp_failed_freq_ck.list";
        FILE* tmpfile=fopen(fname.c_str(), "w");
        if (!(tmpfile)) {
            printf("Error: Failed to open file %s.\n",fname.c_str());
            exit(EXIT_FAILURE);
        }
        string tmpstr = "SNP\tAllele\tOutcome\tExposure\n";
        fputs(tmpstr.c_str(),tmpfile);
        for(int i=0;i<esdata->_esi_include.size();i++)
        {
            double gfreq=etrait->_esi_freq[etrait->_esi_include[i]];
            double efreq=esdata->_esi_freq[esdata->_esi_include[i]];
            if(gfreq>0 && efreq>0 )
            {
                if( fabs(gfreq-efreq)>freqthresh )
                {
                    failcount++;
                    tmpstr=esdata->_esi_rs[esdata->_esi_include[i]] + '\t' + esdata->_esi_allele1[esdata->_esi_include[i]] + '\t' + atos(gfreq) + '\t' + atos(efreq) + '\n';
                    fputs(tmpstr.c_str(),tmpfile);
                }
                else
                {
                    pasgid.push_back(etrait->_esi_include[i]);
                    paseid.push_back(esdata->_esi_include[i]);
                }
            }
            else
            {
                pasgid.push_back(etrait->_esi_include[i]);
                paseid.push_back(esdata->_esi_include[i]);
            }
        }
        if(paseid.size()==0)
        {
            printf("Error: No SNP passed Allele Frequency check. please check the frequecny colomn in your GWAS summmay data and your .esi file.\n");
            exit(EXIT_FAILURE);
        }
        etrait->_esi_include = pasgid;
        esdata->_esi_include = paseid;
        double prop=1.0*failcount/snpnum;
        if(prop > percenthresh) {
            success = false;
            printf("%ld SNPs (%0.2f%% > %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        } else  printf("%ld SNPs (%0.2f%% <= %0.2f%%) with allele frequency differences > %0.2f between any pair of the data sets are excluded from the analysis.\n", failcount,100.0*failcount/snpnum,100*percenthresh,freqthresh);
        
        fclose(tmpfile);
        return prop;
    }
    void update_geIndx(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(etrait->_esi_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
        }
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        etrait->_esi_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
    }


    void read_smaslist(vector<string> &smasNames, string eqtlsmaslstName)
    {
        ifstream smas(eqtlsmaslstName.c_str());
        if (!smas) throw ("Error: can not open the file [" + eqtlsmaslstName + "] to read.");
        cout << "Reading eQTL summary file names from [" + eqtlsmaslstName + "]." << endl;
        char buf[MAX_LINE_SIZE];
        map<string, int> probe_map;
        long mapsize=0;
        while (smas.getline(buf, MAX_LINE_SIZE))
        {
            if(buf[0]!='\0')
            {
                vector<string> vs_buf;
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num==0)
                {
                    printf("Blank row found and skipped!\n");
                    continue;
                } else if( col_num >1) {
                    printf("Column number is not correct with this row\n %s\n",buf);
                } else {
                    probe_map.insert(pair<string,int>(vs_buf[0],mapsize));
                    if(mapsize<probe_map.size())
                    {
                        smasNames.push_back(vs_buf[0]);
                        mapsize=probe_map.size();
                    }
                    else
                    {
                        printf("WARNING: duplicate summary file name %s found and skipped.\n",vs_buf[0].c_str());
                    }

                }
                
            }
            
        }
        cout << smasNames.size()<<" eQTL summary file names are included from [" + eqtlsmaslstName + "]." << endl;
        smas.close();
    }
  
    void ld_calc_o2m(VectorXd &ld_v,long targetid, MatrixXd &X, bool centered)
    {
        long size=X.cols();
        long n=X.rows();
        
        VectorXd tmpX(size);
        VectorXd tmpX2(size);
        VectorXd tmpXY(size);
        if(centered)
        {
            #pragma omp parallel for
            for(int i=0;i<size;i++){
                tmpX2[i]=X.col(i).dot(X.col(i));
                tmpXY[i]=X.col(targetid).dot(X.col(i));
            }
            float tmpY2=X.col(targetid).dot(X.col(targetid));
            ld_v=tmpXY.array()/sqrt(tmpX2.array()*tmpY2);
        }
        else
        {
            #pragma omp parallel for
            for(int i=0;i<size;i++){
                tmpX[i]=X.col(i).sum();
                tmpX2[i]=X.col(i).dot(X.col(i));
                tmpXY[i]=X.col(targetid).dot(X.col(i));
            }
            float tmpY=X.col(targetid).sum();;
            float tmpY2=X.col(targetid).dot(X.col(targetid));
            ld_v=(tmpXY*n-tmpX*tmpY).array()/sqrt((tmpX2.array()*n-tmpX.array()*tmpX.array())*(tmpY2*n-tmpY*tmpY));
        }
       
    }
    void init_smr_wk(SMRWK* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear();
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP, int lowerbp,int upperbp,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid=-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                        else smrwk->freq.push_back(esdata->_esi_freq[j]);
                    }
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
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                
                if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid);
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                    else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                    
                }
            }
        }
        return maxid;
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[j]+9>1e-6)
                    {
                        if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                            {
                                smrwk->bxz.push_back(esdata->_bxz[i][j]);
                                smrwk->sexz.push_back(esdata->_sexz[i][j]);
                                smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                smrwk->byz.push_back(gdata->byz[j]);
                                smrwk->seyz.push_back(gdata->seyz[j]);
                                smrwk->pyz.push_back(gdata->pvalue[j]);
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata->_esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                                smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                                smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                                if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                                if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq.push_back(esdata->_esi_freq[j]);
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                        
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            smrwk->bxz.push_back(esdata->_bxz[i][j]);
                            smrwk->sexz.push_back(esdata->_sexz[i][j]);
                            smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            smrwk->byz.push_back(gdata->byz[j]);
                            smrwk->seyz.push_back(gdata->seyz[j]);
                            smrwk->pyz.push_back(gdata->pvalue[j]);
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata->_esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                            if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                            else smrwk->freq.push_back(esdata->_esi_freq[j]);

                        }
                    }
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
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(snpchr==esdata->_epi_chr[i] && abs(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[ge_rowid]+9>1e-6)
                {
                    if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0)
                    {
                        if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                        {
                            smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                            smrwk->sexz.push_back(esdata->_val[se_start+j]);
                            smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                            smrwk->byz.push_back(gdata->byz[ge_rowid]);
                            smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                            smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                            smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                            smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                            if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                            else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                        } else {
                            printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                            double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            double p=pchisq(z*z, 1);
                            string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                            
                            printf("%s\n",tmp.c_str());
                        }
                        
                    } else {
                        smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                        smrwk->sexz.push_back(esdata->_val[se_start+j]);
                        smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                        smrwk->byz.push_back(gdata->byz[ge_rowid]);
                        smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                        smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                        if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                        else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                    }
                  
                }
            }
        }
        
        return maxid;
    }
    long fill_trans_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, vector<int> &topTransRowid, vector<int> &topTransBP,vector<int> &topTransChr,const char* refSNP,int cis_itvl,int trans_itvl, bool heidioffFlag, int tridx)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        if(esdata->_rowid.empty())
        {
            for (int j = topTransRowid[tridx]; j>=0; j--)
            {
                if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==topTransChr[tridx] && ABS(topTransBP[tridx]-snpbp)<=trans_itvl && gdata->seyz[j]+9>1e-6)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        if((refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) || snpbp==topTransBP[tridx]) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                        else smrwk->freq.push_back(esdata->_esi_freq[j]);
                    }
                }
                
            }
            for (int j = topTransRowid[tridx]+1; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==topTransChr[tridx] && ABS(topTransBP[tridx]-snpbp)<=trans_itvl && gdata->seyz[j]+9>1e-6)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        if((refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) || snpbp==topTransBP[tridx]) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                        else smrwk->freq.push_back(esdata->_esi_freq[j]);
                    }
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
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(snpchr==topTransChr[tridx] && abs(topTransBP[tridx]-snpbp)<=trans_itvl && gdata->seyz[ge_rowid]+9>1e-6)
                {
                    //topTransBP[jj]==esdata._esi_bp[esdata._rowid[esdata._cols[i<<1]+topTransRowid[jj]]];
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    if((refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) || snpbp==topTransBP[tridx]) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                    else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                }
            }
        }
        return maxid;
    }
    
    double heidi_test(bInfo* bdata,SMRWK* smrwk, long maxid,double ld_top, double threshold, int m_hetero,long &nsnp )
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        Map<VectorXd> ei_bxz(&smrwk->bxz[0],smrwk->bxz.size());
        Map<VectorXd> ei_sexz(&smrwk->sexz[0],smrwk->sexz.size());
        
        VectorXd zsxz=ei_bxz.array()/ei_sexz.array();
        
        
        make_XMat(bdata,smrwk->curId, _X);
        ld_calc_o2m(ld_v,maxid,_X);
        if(fabs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
        else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
        if(sn_ids.size() < m_hetero) {
            
            printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is less than a threshold (%d).\n", sn_ids.size(), m_hetero);
            return -9;
        }
        
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
            _byz[j]=smrwk->byz[sn_ids[j]];
            _seyz[j]=smrwk->seyz[sn_ids[j]];
            _bxz[j]=smrwk->bxz[sn_ids[j]];
            _sexz[j]=smrwk->sexz[sn_ids[j]];
            _zsxz[j]=zsxz(sn_ids[j]);
            _X_heidi.col(j)=_X.col(sn_ids[j]);
        }
        _X.resize(0,0);
        cor_calc(_LD_heidi, _X_heidi);
        
        _X_heidi.resize(0,0);
        
        nsnp = sn_ids.size();
        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
        
        return pdev;
    }
    void update_snidx(SMRWK* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
    {
        if(sn_ids.size()>max_snp_slct)
        {
            //printf("Top %d SNPs (ordered by eQTL p-value) are used int the %s...\n",max_snp_slct,forwhat.c_str());
            priority_queue< pair<double, int> > q;
            vector<int> sn_slct_ids;
            for(int i=0;i<sn_ids.size();i++) q.push((pair<double, int>(fabs(smrwk->zxz[sn_ids[i]]), sn_ids[i])));
            for(int i=0;i<max_snp_slct;i++){
                sn_slct_ids.push_back(q.top().second);
                q.pop();
            }
            sn_ids.swap(sn_slct_ids);
        }
        
    }
    void extract_smrwk(SMRWK* smrwk,vector<int> &sn_ids,SMRWK* smrwk2)
    {
        smrwk2->cur_chr=smrwk->cur_chr;
        smrwk2->cur_prbidx=smrwk->cur_prbidx;
        for(int i=0;i<sn_ids.size();i++)
        {
            smrwk2->bxz.push_back(smrwk->bxz[sn_ids[i]]);
            smrwk2->sexz.push_back(smrwk->sexz[sn_ids[i]]);
            smrwk2->freq.push_back(smrwk->freq[sn_ids[i]]);
            smrwk2->byz.push_back(smrwk->byz[sn_ids[i]]);
            smrwk2->seyz.push_back(smrwk->seyz[sn_ids[i]]);
            smrwk2->pyz.push_back(smrwk->pyz[sn_ids[i]]);
            smrwk2->zxz.push_back(smrwk->zxz[sn_ids[i]]);
            smrwk2->curId.push_back(smrwk->curId[sn_ids[i]]);
            smrwk2->bpsnp.push_back(smrwk->bpsnp[sn_ids[i]]);
            smrwk2->snpchrom.push_back(smrwk->snpchrom[sn_ids[i]]);
            smrwk2->rs.push_back(smrwk->rs[sn_ids[i]]);
            smrwk2->allele1.push_back(smrwk->allele1[sn_ids[i]]);
            smrwk2->allele2.push_back(smrwk->allele2[sn_ids[i]]);
        }

    }
    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1) {
        //approximate maximum independent set
        //Modified version of rm_cor_indi from grm.cpp
        
        int i = 0, j = 0, i_buf = 0;
        vector<int> rm_ID2;
        
        //float tmpr = 0; //rm_ID1 is the same as indx1 in ld_prune of R Script
        for (i = 0; i < m; i++) {
            for (j = 0; j < i; j++) {
                if (fabs(R(i,j)) > R_cutoff ) {
                    rm_ID1.push_back(j);
                    rm_ID2.push_back(i);
                }
            }
        }
        
        // count the number of appearance of each "position" in the vector, which involves a few steps
        vector<int> rm_uni_ID(rm_ID1);
        rm_uni_ID.insert(rm_uni_ID.end(), rm_ID2.begin(), rm_ID2.end());
        stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
        rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
        map<int, int> rm_uni_ID_count;
        for (i = 0; i < rm_uni_ID.size(); i++) {
            i_buf = count(rm_ID1.begin(), rm_ID1.end(), rm_uni_ID[i]) + count(rm_ID2.begin(), rm_ID2.end(), rm_uni_ID[i]);
            rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
        }
        
        // swapping
        map<int, int>::iterator iter1, iter2;
        for (i = 0; i < rm_ID1.size(); i++) {
            iter1 = rm_uni_ID_count.find(rm_ID1[i]);
            iter2 = rm_uni_ID_count.find(rm_ID2[i]);
            int c1=iter1->second , c2=iter2->second;
            if ( c1<c2 ) {
                i_buf = rm_ID1[i];
                rm_ID1[i] = rm_ID2[i];
                rm_ID2[i] = i_buf;
            }
        }
        stable_sort(rm_ID1.begin(), rm_ID1.end());
        rm_ID1.erase(unique(rm_ID1.begin(), rm_ID1.end()), rm_ID1.end());
    }
    void update_smrwk_x(SMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X)
    {
        vector<double> byz,seyz, bxz,sexz,zxz,freq,pyz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
        MatrixXd _X;
        byz.resize(sn_ids.size());
        seyz.resize(sn_ids.size());
        bxz.resize(sn_ids.size());
        sexz.resize(sn_ids.size());
        zxz.resize(sn_ids.size());
        freq.resize(sn_ids.size());
        pyz.resize(sn_ids.size());
        curId.resize(sn_ids.size());
        bpsnp.resize(sn_ids.size());
        snpchrom.resize(sn_ids.size());
        rs.resize(sn_ids.size());
        allele1.resize(sn_ids.size());
        allele2.resize(sn_ids.size());
        _X.resize(X.rows(), sn_ids.size());
        
        #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            byz[j]=smrwk->byz[sn_ids[j]];
            seyz[j]=smrwk->seyz[sn_ids[j]];
            bxz[j]=smrwk->bxz[sn_ids[j]];
            sexz[j]=smrwk->sexz[sn_ids[j]];
            zxz[j]=smrwk->zxz[sn_ids[j]];
            freq[j]=smrwk->freq[sn_ids[j]];
            pyz[j]=smrwk->pyz[sn_ids[j]];
            curId[j]=smrwk->curId[sn_ids[j]];
            bpsnp[j]=smrwk->bpsnp[sn_ids[j]];
            snpchrom[j]=smrwk->snpchrom[sn_ids[j]];
            rs[j]=smrwk->rs[sn_ids[j]];
            allele1[j]=smrwk->allele1[sn_ids[j]];
            allele2[j]=smrwk->allele2[sn_ids[j]];
            _X.col(j)=X.col(sn_ids[j]);
        }
        smrwk->bxz.swap(bxz);
        smrwk->sexz.swap(sexz);
        smrwk->freq.swap(freq);
        smrwk->byz.swap(byz);
        smrwk->seyz.swap(seyz);
        smrwk->pyz.swap(pyz);
        smrwk->zxz.swap(zxz);
        smrwk->curId.swap(curId);
        smrwk->bpsnp.swap(bpsnp);
        smrwk->snpchrom.swap(snpchrom);
        smrwk->rs.swap(rs);
        smrwk->allele1.swap(allele1);
        smrwk->allele2.swap(allele2);
        X=_X;
    }
    double heidi_test_new(bInfo* bdata,SMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, double theta)
    {
        //the new method would calcualte maxid after each filtering
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        double pthres=pchisq(threshold,1);
        //printf("Filtering SNPs (%ld in total) at eQTL p-value < %e for the HEIDI test.\n",smrwk->zxz.size(), pthres);
        for(int i=0;i<smrwk->zxz.size();i++)
        {
            if(smrwk->zxz[i]*smrwk->zxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {
            
           // printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        //printf("%ld SNPs left after filtering.\n",sn_ids.size());
        update_snidx(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        
        SMRWK smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
        long maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        
        make_XMat(bdata,smrwk_heidi.curId, _X);
        //printf("Removing SNPs with LD r-squared between top-SNP %s > %f or < %f...\n",smrwk_heidi.rs[maxid_heidi].c_str(),ldr2_top,ld_min);
        ld_calc_o2m(ld_v,maxid_heidi,_X);
        
        if(fabs(ldr2_top-1)>1e-6 || ld_min>0) {
            sn_ids.clear();
            for(int i=0;i<smrwk_heidi.zxz.size();i++)
            {
                if(i!= maxid_heidi)
                {
                    double ldr2tmp=ld_v(i)*ld_v(i);
                    if((ldr2tmp<ldr2_top) && (ldr2tmp > ld_min)) sn_ids.push_back(i);
                }
                else{
                    sn_ids.push_back(i);
                }
            }
        }
        //printf("%ld SNPs are removed and %ld SNPs are retained.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
        if(sn_ids.size() < m_hetero) {
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        
        //printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
        int m = (int)smrwk_heidi.bxz.size();
        vector<int> rm_ID1;
        MatrixXd C;
        cor_calc(C, _X);
       
        double ld_top=sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
        //printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
        if(m-rm_ID1.size() < m_hetero) {
            
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
            return -9;
        }
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sn_ids.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        
        update_snidx(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
        if (sn_ids.size() < C.size()) { //Build new matrix
            MatrixXd D(sn_ids.size(),sn_ids.size());
            for (int i = 0 ; i < sn_ids.size() ; i++) {
                for (int j = 0 ; j < sn_ids.size() ; j++) {
                    D(i,j) = C(sn_ids[i],sn_ids[j]);
                }
            }
            C = D;
        }        
        
        VectorXd _byz,_seyz, _bxz,_sexz,_zsxz;
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
        #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            _byz[j]=smrwk_heidi.byz[sn_ids[j]];
            _seyz[j]=smrwk_heidi.seyz[sn_ids[j]];
            _bxz[j]=smrwk_heidi.bxz[sn_ids[j]];
            _sexz[j]=smrwk_heidi.sexz[sn_ids[j]];
            _zsxz[j]=smrwk_heidi.zxz[sn_ids[j]];
        }
        nsnp = sn_ids.size();
        double pdev=-9;
        if(sampleoverlap) pdev=bxy_mltheter_so(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
        else pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp);
        
        //printf("pHeidi is %e with %ld SNPs including in the HEIDI test.\n",pdev,nsnp);
        return pdev;
    }
    int get_num_probes(eqtlInfo* esdata, double p_smr, int cis_itvl)
    {
        int numprb=0;
        if(esdata->_valNum==0)
        {
            for(int ii=0;ii<esdata->_include.size();ii++)
            {
                long i=esdata->_include[ii];
                double pxz_top=1;
             
                int prbbp=esdata->_epi_bp[i];
                int prbchr=esdata->_epi_chr[i];
                int cisR=prbbp+cis_itvl*1000;
                int cisL=(prbbp-cis_itvl*1000)>0?(prbbp-cis_itvl*1000):0;
                for(int jj=0;jj<esdata->_esi_include.size();jj++)
                {
                    long j=esdata->_esi_include[jj];
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==prbchr && snpbp>=cisL && snpbp<=cisR)
                    {
                        double beta=esdata->_bxz[i][j];
                        double se=esdata->_sexz[i][j];
                        if(fabs(se+9)<1e-6) continue;
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=pxz_top) pxz_top=pxz;
                    }
                }
                if(pxz_top<=p_smr) numprb++;
            }
        }
        else
        {
            if(esdata->_val.size()==0)
            {
                printf("Error: No data extracted from the input, please check.\n");
                exit(EXIT_FAILURE);
            }
            
            for(int ii=0;ii<esdata->_include.size();ii++)
            {
                uint64_t proid=esdata->_include[ii];
                uint64_t pos=esdata->_cols[proid<<1];
                uint64_t pos1=esdata->_cols[(proid<<1)+1];
                uint64_t num=pos1-pos;

                double pxz_top=1;
                int prbbp=esdata->_epi_bp[proid];
                int prbchr=esdata->_epi_chr[proid];
                int cisR=prbbp+cis_itvl*1000;
                int cisL=(prbbp-cis_itvl*1000)>0?(prbbp-cis_itvl*1000):0;
                for(int j=0;j<num;j++)
                {
                    int esiid=esdata->_rowid[pos+j];
                    int snpbp=esdata->_esi_bp[esiid];
                    int snpchr=esdata->_esi_chr[esiid];
                    if(snpchr==prbchr && snpbp>=cisL && snpbp<=cisR)
                    {
                        double beta=esdata->_val[pos+j];
                        double se=esdata->_val[pos+j+num];
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=pxz_top) pxz_top=pxz;
                    }
                }
                if(pxz_top<=p_smr) numprb++;
            }
        }
        
        return numprb;
        
    }
    void smr_heidi_func_old(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg)
    {
        
        uint64_t probNum = esdata->_include.size();
        double thresh_heidi= chi_val(1,p_hetero),theta=0;
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        map<string, string>::iterator iter;
        map<string, int>::iterator iter2;
        FILE* smr=NULL;
        long write_count=0, noSNPprb=0, noSNPprbPassthresh=0;
        string outstr="";
        string smrfile="";
        if(outFileName!=NULL) {
            smrfile = string(outFileName)+".smr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
            if(refSNP || targetLstFlg) outstr="probeID\tProbeChr\tGene\tProbe_bp\ttargetSNP\ttargetSNP_chr\ttargetSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            else outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: error in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
        int numsigprobe=get_num_probes(esdata, p_smr, cis_itvl);
        SMRWK smrwk;
        //printf("\nPerforming SMR analysis (SMR and HEIDI tests) for %d probes (with at least a cis-eQTL at p < %0.2e)... \n",numsigprobe, p_smr);
        //printf("For each probe, the analysis will only include SNPs with eQTL p-values < %e,\n",p_hetero);
        //printf("then exclude SNPs with LD r-squared between top-SNP > %0.2f or < %0.2f, and further exclude one of each pair of the remaining SNPs with LD r-squared > %0.2f.\n", ld_top, ld_min,ld_top);
        
       // float progr0=0.0, progr1;
       // progress_print(progr0);
        cis_itvl=cis_itvl*1000;
        double cr=0;
        for(int ii=0;ii<probNum;ii++)
        {
            double desti=1.0*ii/(probNum-1);
            if(desti>=cr)
            {
                //printf("%3.0f%%\r", 100.0*desti);
                //fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            int i=esdata->_include[ii];
            int probebp=esdata->_epi_bp[i];
            int probechr=esdata->_epi_chr[i];
            string probename=esdata->_epi_prbID[i];
            string probegene=esdata->_epi_gene[i];
            char probeorien=esdata->_epi_orien[i];
            string specifiedsnp;
            if(prb_snp.size()>0)
            {
                iter=prb_snp.find(probename);
                if(iter!=prb_snp.end()) {
                    specifiedsnp=iter->second;
                    if(targetLstFlg) {
                        refSNP=specifiedsnp.c_str();
                       // printf("The target SNP %s of probe %s extracted from the target SNP - probe list.\n",refSNP,probename.c_str());
                    } //else printf("The SNP %s of probe %s extracted from the SNP - probe list.\n",specifiedsnp.c_str(),probename.c_str());
                }
                else
                {
                    printf("Can't happen. Please report.\n");
                    continue;
                }
            }
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            smrwk.cur_chr=probechr;
            long maxid =-9;
            if(refSNP || prb_snp.size()>0)
            {
                vector<int> topTransBP;
                vector<int> topTransRowid;
                vector<int> topTransChr;
                if(prb_snp.size()==0) specifiedsnp=string(refSNP);
                iter2=esdata->_snp_name_map.find(specifiedsnp);
                if(iter2!=esdata->_snp_name_map.end())
                {
                    int idx=iter2->second;
                    topTransRowid.push_back(idx);
                    topTransBP.push_back(esdata->_esi_bp[idx]);
                    topTransChr.push_back(esdata->_esi_chr[idx]);
                    maxid =fill_trans_smr_wk(bdata, gdata, esdata, &smrwk, topTransRowid, topTransBP,topTransChr,refSNP, cis_itvl,cis_itvl, heidioffFlag,0);
                } else {
                    //printf("WARNING: can't find target SNP %s for probe %s from .epi file.\n",refSNP, probename.c_str());
                    continue;
                }
                
            } else maxid =fill_smr_wk(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl,heidioffFlag);
            
            
            if(refSNP!=NULL && maxid==-9)
            {
                //printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                continue;
            }
            if (smrwk.bxz.size() == 0) {
                noSNPprb++;
                //printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                continue;
            }
           
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            if(sampleoverlap)
            {
                //printf("Estimating the correlation ...\n");
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
                    //printf("WARNING: less than %d common SNPs obtained from the cis-region of probe %s at a p-value threshold %5.2e.\n ", minCor,probename.c_str(), pmecs);
                    //printf("probe %s is skipped.\n ", probename.c_str());
                    continue;
                }
                else
                {
                    theta=cor(zxz,zyz);
                    //printf("The estimated correlation is %f.\n",theta);
                }
            }
            zsxz=ei_bxz.array()/ei_sexz.array();
            
            if(refSNP==NULL) {
                if(opt) maxid=max_zsmr_id(&smrwk, p_smr);
                else maxid=max_abs_id(zsxz);
            }
            if(maxid==-9) {
                noSNPprbPassthresh++;
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                continue;
            }
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            if(refSNP==NULL && pxz_val>p_smr) {
                noSNPprbPassthresh++;
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                continue;
            } else {
               // printf("Analysing probe %s...\n", probename.c_str());
            }
            
            /*****test***/
           // string failName= probename+".expo";
           // FILE* failfptr=fopen(failName.c_str(),"w");
           // if(failfptr==NULL)
           // {
           //     printf("ERROR: failed in open file %s.\n",failName.c_str()) ;
           //     exit(EXIT_FAILURE);
           // }
           // for(int k=0; k<smrwk.bxz.size();k++)
           // {
           //     string snpstr=smrwk.rs[k] + '\t' + atos(smrwk.bxz[k]) + '\t' + atos(smrwk.sexz[k]) + '\n';
           //     if(fputs_checked(snpstr.c_str(),failfptr))
           //     {
           //         printf("ERROR: in writing file %s .\n", failName.c_str());
           //         exit(EXIT_FAILURE);
           //     }
               
           // }
           // fclose(failfptr);
           // failName= probename+".outo";
           // failfptr=fopen(failName.c_str(),"w");
           // if(failfptr==NULL)
           // {
           //     printf("ERROR: failed in open file %s.\n",failName.c_str()) ;
           //     exit(EXIT_FAILURE);
           // }
           // for(int k=0; k<smrwk.byz.size();k++)
           // {
           //     string snpstr=smrwk.rs[k] + '\t' + atos(smrwk.byz[k]) + '\t' + atos(smrwk.seyz[k]) + '\n';
           //     if(fputs_checked(snpstr.c_str(),failfptr))
           //     {
           //         printf("ERROR: in writing file %s .\n", failName.c_str());
           //         exit(EXIT_FAILURE);
           //     }
           // }
           // fclose(failfptr);
            
            /******/

            double byzt=smrwk.byz[maxid], bxzt=smrwk.bxz[maxid], seyzt=smrwk.seyz[maxid], sexzt=smrwk.sexz[maxid];
            double bxy_val = byzt / bxzt;
            double sexy_val = -9;
            if(sampleoverlap)
            {
                sexy_val = sqrt(seyzt* seyzt/(bxzt*bxzt) + sexzt*sexzt*byzt*byzt/(bxzt*bxzt*bxzt*bxzt) - 2*theta*seyzt*sexzt*byzt/(bxzt*bxzt*bxzt));
            } else {
                sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
            }
            double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
            double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
            
            SMRRLT currlt;
            if(smr)
            {
                outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + smrwk.rs[maxid] + '\t' + atos(smrwk.snpchrom[maxid]) + '\t' + atos(smrwk.bpsnp[maxid]) + '\t' + smrwk.allele1[maxid] + '\t' + smrwk.allele2[maxid] + '\t' + atos(smrwk.freq[maxid]) + '\t';
                outstr += atos(smrwk.byz[maxid]) + '\t' + atos(smrwk.seyz[maxid]) + '\t' + dtos(smrwk.pyz[maxid]) + '\t';
                outstr += atos(smrwk.bxz[maxid]) + '\t' + atos(smrwk.sexz[maxid]) + '\t' + dtos(pxz_val) + '\t';
                outstr += atos(bxy_val) + '\t' + atos(sexy_val) + '\t' + dtos(pxy_val) + '\t';
                
            } else {
                currlt.ProbeID=probename;
                currlt.ProbeChr=probechr;
                currlt.Gene=probegene;
                currlt.Probe_bp=probebp;
                currlt.Orien=probeorien;
                currlt.SNP=smrwk.rs[maxid];
                currlt.SNP_Chr=smrwk.snpchrom[maxid];
                currlt.SNP_bp=smrwk.bpsnp[maxid];
                currlt.A1=smrwk.allele1[maxid];
                currlt.A2=smrwk.allele2[maxid];
                currlt.Freq=smrwk.freq[maxid];
                currlt.b_GWAS=smrwk.byz[maxid];
                currlt.se_GWAS=smrwk.seyz[maxid];
                currlt.p_GWAS=smrwk.pyz[maxid];
                currlt.b_eQTL=smrwk.bxz[maxid];
                currlt.se_eQTL=smrwk.sexz[maxid];
                currlt.p_eQTL=pxz_val;
                currlt.b_SMR=bxy_val;
                currlt.se_SMR=sexy_val;
                currlt.p_SMR=pxy_val;
            }
            
            if(heidioffFlag || pxy_val>threshpsmrest)
            {
              //  printf("INFO: the HEIDI test for probe %s is skipped because HEIDI test is turned off by the --heidi-off option or p_SMR does not pass the %e threshold.\n", probename.c_str(),threshpsmrest);
                if(smr)
                {
                    outstr+= string("NA") + '\t' + string("NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=-9;
                    currlt.nsnp=-9;
                    smrrlts.push_back(currlt);
                }
                
            }
            else
            {
                
                /****Generate the test data***/
                /*
                 MatrixXd _X;
                 make_XMat(bdata,smrwk.curId, _X);
                 string xfilename =string(outFileName)+".mat";
                 FILE* xfile = fopen(xfilename.c_str(), "w");
                 if (!(xfile)) {
                 printf("Open error %s\n", xfilename.c_str());
                 exit(1);
                 }
                 string str=smrwk.rs[0];
                 for( int j=1;j<smrwk.rs.size();j++) str+='\t'+smrwk.rs[j];
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 for(int i=0;i<_X.rows();i++)
                 {
                 string str=atos(_X(i,0));
                 for( int j=1;j<_X.cols();j++) str+='\t'+atos(_X(i,j));
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 }
                 fclose(xfile);
                
                string gwasfname =string(outFileName)+".gwas";
                FILE* gfile = fopen(gwasfname.c_str(), "w");
                if (!(gfile)) {
                    printf("Open error %s\n", gwasfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.byz[j])+'\t'+atos(smrwk.seyz[j])+'\n';
                    fputs_checked(str.c_str(), gfile);
                }
                fclose(gfile);
                
                string eqtlfname =string(outFileName)+".eqtl";
                FILE* efile = fopen(eqtlfname.c_str(), "w");
                if (!(efile)) {
                    printf("Open error %s\n", eqtlfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.bxz[j])+'\t'+atos(smrwk.sexz[j])+'\n';
                    fputs_checked(str.c_str(), efile);
                }
                fclose(efile);
                 */
                /*****/

                long nsnp=-9;
                double pdev=-9;
                if(!heidioffFlag) {
                    if(new_heidi_mth){
                        if(refSNP!=NULL) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid, ld_min,opt_hetero,sampleoverlap, theta);
                        else pdev= heidi_test_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp, ld_min ,opt_hetero,sampleoverlap, theta);
                    } else pdev= heidi_test(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                }
                if(smr)
                {
                    outstr+= (pdev >= 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=pdev;
                    currlt.nsnp=(int)nsnp;
                    smrrlts.push_back(currlt);
                }
            }
        }
        //printf("%ld probes has no SNPs in common.\n", noSNPprb);
        //printf("%ld probes has no SNPs passed the p-value threshold %e for the SMR analysis.\n", noSNPprbPassthresh, p_smr);
        if(smr){
            printf("Results of %ld probes have been saved in file %s.\n",write_count,smrfile.c_str());
            fclose(smr);
        } else {
            //printf("Results of %ld probes have been returned.\n",smrrlts.size());
        }
        
    }

    void smr_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg)
    {
        uint64_t probNum = esdata->_include.size();
        double thresh_heidi= chi_val(1,p_hetero),theta=0;
        // VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        // MatrixXd _X,_LD,_LD_heidi,_X_heidi;    
        FILE* smr=NULL;
        long write_count=0, noSNPprb=0, noSNPprbPassthresh=0;
        string outstr="";
        string smrfile="";
        if(outFileName!=NULL) {
            smrfile = string(outFileName)+".smr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
            if(refSNP || targetLstFlg) outstr="probeID\tProbeChr\tGene\tProbe_bp\ttargetSNP\ttargetSNP_chr\ttargetSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            else outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: error in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
        int numsigprobe=get_num_probes(esdata, p_smr, cis_itvl);    

        vector<SMRRLT> smrrlts_tmp(probNum); vector<int> includ_num(probNum,1);
        cis_itvl=cis_itvl*1000;
        double cr=0;

        #pragma omp parallel for
        for(int ii=0;ii<probNum;ii++)
        {
            SMRWK smrwk; SMRRLT currlt;
            VectorXd zsxz;
            map<string, string>::iterator iter;
            map<string, int>::iterator iter2;
            int i=esdata->_include[ii];
            int probebp=esdata->_epi_bp[i];
            int probechr=esdata->_epi_chr[i];
            string probename=esdata->_epi_prbID[i];
            string probegene=esdata->_epi_gene[i];
            char probeorien=esdata->_epi_orien[i];
            string specifiedsnp;
            if(prb_snp.size()>0)
            {
                iter=prb_snp.find(probename);
                if(iter!=prb_snp.end()) {
                    specifiedsnp=iter->second;
                    if(targetLstFlg) {
                        refSNP=specifiedsnp.c_str();
                       // printf("The target SNP %s of probe %s extracted from the target SNP - probe list.\n",refSNP,probename.c_str());
                    } //else printf("The SNP %s of probe %s extracted from the SNP - probe list.\n",specifiedsnp.c_str(),probename.c_str());
                }
                else
                {
                    includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                    printf("Can't happen. Please report.\n");
                    continue;
                }
            }
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            smrwk.cur_chr=probechr;
            long maxid =-9;
            if(refSNP || prb_snp.size()>0)
            {
                vector<int> topTransBP;
                vector<int> topTransRowid;
                vector<int> topTransChr;
                if(prb_snp.size()==0) specifiedsnp=string(refSNP);
                iter2=esdata->_snp_name_map.find(specifiedsnp);
                if(iter2!=esdata->_snp_name_map.end())
                {
                    int idx=iter2->second;
                    topTransRowid.push_back(idx);
                    topTransBP.push_back(esdata->_esi_bp[idx]);
                    topTransChr.push_back(esdata->_esi_chr[idx]);
                    maxid =fill_trans_smr_wk(bdata, gdata, esdata, &smrwk, topTransRowid, topTransBP,topTransChr,refSNP, cis_itvl,cis_itvl, heidioffFlag,0);
                } else {
                    includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                    //printf("WARNING: can't find target SNP %s for probe %s from .epi file.\n",refSNP, probename.c_str());
                    continue;
                }
                
            } else maxid =fill_smr_wk(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl,heidioffFlag);
        
            
            if(refSNP!=NULL && maxid==-9)
            {
                includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                //printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                continue;
            }
            if (smrwk.bxz.size() == 0) {
                includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                // noSNPprb++;
                //printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                continue;
            }
           
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            if(sampleoverlap)
            {
                //printf("Estimating the correlation ...\n");
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
                    includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                    //printf("WARNING: less than %d common SNPs obtained from the cis-region of probe %s at a p-value threshold %5.2e.\n ", minCor,probename.c_str(), pmecs);
                    //printf("probe %s is skipped.\n ", probename.c_str());
                    continue;
                }
                else
                {
                    theta=cor(zxz,zyz);
                    //printf("The estimated correlation is %f.\n",theta);
                }
            }
            zsxz=ei_bxz.array()/ei_sexz.array();
            
            if(refSNP==NULL) {
                if(opt) maxid=max_zsmr_id(&smrwk, p_smr);
                else maxid=max_abs_id(zsxz);
            }
            if(maxid==-9) {
                includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                //noSNPprbPassthresh++;
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                continue;
            }
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            if(refSNP==NULL && pxz_val>p_smr) {
                includ_num[ii] = 0; //smrrlts_tmp[ii]=currlt;
                //noSNPprbPassthresh++;
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                continue;
            } else {
               // printf("Analysing probe %s...\n", probename.c_str());
            }
            
            /*****test***/
           // string failName= probename+".expo";
           // FILE* failfptr=fopen(failName.c_str(),"w");
           // if(failfptr==NULL)
           // {
           //     printf("ERROR: failed in open file %s.\n",failName.c_str()) ;
           //     exit(EXIT_FAILURE);
           // }
           // for(int k=0; k<smrwk.bxz.size();k++)
           // {
           //     string snpstr=smrwk.rs[k] + '\t' + atos(smrwk.bxz[k]) + '\t' + atos(smrwk.sexz[k]) + '\n';
           //     if(fputs_checked(snpstr.c_str(),failfptr))
           //     {
           //         printf("ERROR: in writing file %s .\n", failName.c_str());
           //         exit(EXIT_FAILURE);
           //     }
               
           // }
           // fclose(failfptr);
           // failName= probename+".outo";
           // failfptr=fopen(failName.c_str(),"w");
           // if(failfptr==NULL)
           // {
           //     printf("ERROR: failed in open file %s.\n",failName.c_str()) ;
           //     exit(EXIT_FAILURE);
           // }
           // for(int k=0; k<smrwk.byz.size();k++)
           // {
           //     string snpstr=smrwk.rs[k] + '\t' + atos(smrwk.byz[k]) + '\t' + atos(smrwk.seyz[k]) + '\n';
           //     if(fputs_checked(snpstr.c_str(),failfptr))
           //     {
           //         printf("ERROR: in writing file %s .\n", failName.c_str());
           //         exit(EXIT_FAILURE);
           //     }
           // }
           // fclose(failfptr);
            
            /******/

            double byzt=smrwk.byz[maxid], bxzt=smrwk.bxz[maxid], seyzt=smrwk.seyz[maxid], sexzt=smrwk.sexz[maxid];
            double bxy_val = byzt / bxzt;
            double sexy_val = -9;
            if(sampleoverlap)
            {
                sexy_val = sqrt(seyzt* seyzt/(bxzt*bxzt) + sexzt*sexzt*byzt*byzt/(bxzt*bxzt*bxzt*bxzt) - 2*theta*seyzt*sexzt*byzt/(bxzt*bxzt*bxzt));
            } else {
                sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
            }
            double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
            double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
            
            if(smr)
            {
                outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + smrwk.rs[maxid] + '\t' + atos(smrwk.snpchrom[maxid]) + '\t' + atos(smrwk.bpsnp[maxid]) + '\t' + smrwk.allele1[maxid] + '\t' + smrwk.allele2[maxid] + '\t' + atos(smrwk.freq[maxid]) + '\t';
                outstr += atos(smrwk.byz[maxid]) + '\t' + atos(smrwk.seyz[maxid]) + '\t' + dtos(smrwk.pyz[maxid]) + '\t';
                outstr += atos(smrwk.bxz[maxid]) + '\t' + atos(smrwk.sexz[maxid]) + '\t' + dtos(pxz_val) + '\t';
                outstr += atos(bxy_val) + '\t' + atos(sexy_val) + '\t' + dtos(pxy_val) + '\t';
                
            } else {
                currlt.ProbeID=probename;
                currlt.ProbeChr=probechr;
                currlt.Gene=probegene;
                currlt.Probe_bp=probebp;
                currlt.Orien=probeorien;
                currlt.SNP=smrwk.rs[maxid];
                currlt.SNP_Chr=smrwk.snpchrom[maxid];
                currlt.SNP_bp=smrwk.bpsnp[maxid];
                currlt.A1=smrwk.allele1[maxid];
                currlt.A2=smrwk.allele2[maxid];
                currlt.Freq=smrwk.freq[maxid];
                currlt.b_GWAS=smrwk.byz[maxid];
                currlt.se_GWAS=smrwk.seyz[maxid];
                currlt.p_GWAS=smrwk.pyz[maxid];
                currlt.b_eQTL=smrwk.bxz[maxid];
                currlt.se_eQTL=smrwk.sexz[maxid];
                currlt.p_eQTL=pxz_val;
                currlt.b_SMR=bxy_val;
                currlt.se_SMR=sexy_val;
                currlt.p_SMR=pxy_val;
            }
            
            if(heidioffFlag || pxy_val>threshpsmrest)
            {
              //  printf("INFO: the HEIDI test for probe %s is skipped because HEIDI test is turned off by the --heidi-off option or p_SMR does not pass the %e threshold.\n", probename.c_str(),threshpsmrest);
                if(smr)
                {
                    outstr+= string("NA") + '\t' + string("NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=-9;
                    currlt.nsnp=-9;
                    smrrlts_tmp[ii] = currlt;
                }
                
            }
            else
            {
                
                /****Generate the test data***/
                /*
                 MatrixXd _X;
                 make_XMat(bdata,smrwk.curId, _X);
                 string xfilename =string(outFileName)+".mat";
                 FILE* xfile = fopen(xfilename.c_str(), "w");
                 if (!(xfile)) {
                 printf("Open error %s\n", xfilename.c_str());
                 exit(1);
                 }
                 string str=smrwk.rs[0];
                 for( int j=1;j<smrwk.rs.size();j++) str+='\t'+smrwk.rs[j];
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 for(int i=0;i<_X.rows();i++)
                 {
                 string str=atos(_X(i,0));
                 for( int j=1;j<_X.cols();j++) str+='\t'+atos(_X(i,j));
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 }
                 fclose(xfile);
                
                string gwasfname =string(outFileName)+".gwas";
                FILE* gfile = fopen(gwasfname.c_str(), "w");
                if (!(gfile)) {
                    printf("Open error %s\n", gwasfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.byz[j])+'\t'+atos(smrwk.seyz[j])+'\n';
                    fputs_checked(str.c_str(), gfile);
                }
                fclose(gfile);
                
                string eqtlfname =string(outFileName)+".eqtl";
                FILE* efile = fopen(eqtlfname.c_str(), "w");
                if (!(efile)) {
                    printf("Open error %s\n", eqtlfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.bxz[j])+'\t'+atos(smrwk.sexz[j])+'\n';
                    fputs_checked(str.c_str(), efile);
                }
                fclose(efile);
                 */
                /*****/

                long nsnp=-9;
                double pdev=-9;
                if(!heidioffFlag) {
                    if(new_heidi_mth){
                        if(refSNP!=NULL) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid, ld_min,opt_hetero,sampleoverlap, theta);
                        else pdev= heidi_test_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp, ld_min ,opt_hetero,sampleoverlap, theta);
                    } else pdev= heidi_test(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                }
                if(smr)
                {
                    outstr+= (pdev >= 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=pdev;
                    currlt.nsnp=(int)nsnp;
                    smrrlts_tmp[ii]=currlt;
                }
            }
        }
        for(int i=0;i<includ_num.size();i++)
        {
            if(includ_num[i] != 0) {
                smrrlts.push_back(smrrlts_tmp[i]);
            }
        }
        //printf("%ld probes has no SNPs in common.\n", noSNPprb);
        //printf("%ld probes has no SNPs passed the p-value threshold %e for the SMR analysis.\n", noSNPprbPassthresh, p_smr);
        if(smr){
            printf("Results of %ld probes have been saved in file %s.\n",write_count,smrfile.c_str());
            fclose(smr);
        } else {
            //printf("Results of %ld probes have been returned.\n",smrrlts.size());
        }
        
    }
    
    void smr(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest, bool new_het_mth,bool opt,char* prbseqregion, double ld_min, bool sampleoverlap, double pmecs, int minCor, char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh)
    {
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetsnpproblstName)
        {
            printf("ERROR: --target-snp and --extract-target-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetsnpproblstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        bool targetLstflg=false;
        map<string, string> prb_snp;
        
        if(!heidioffFlag && bFileName == NULL ) throw("ERROR: please input Plink file for the SMR analysis using the flag --bfile.");
        if(gwasFileName==NULL) throw("ERROR: please input GWAS summary data for the SMR analysis using the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("ERROR: please input eQTL summary data for the SMR analysis using the flag --eqtl-summary.");

        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
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
            if(forcefrqck)
            {
                double prop=freq_check(&bdata, &gdata, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
                
            }
        }else
        {
            allele_check(&gdata, &esdata);
            if(forcefrqck)
            {
                double prop=freq_check(&gdata, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
            }
        }
        update_gwas(&gdata);
       
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        if(prbseqregion!= NULL) {
            read_epistartend(&esdata,prbseqregion);
        }
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(targetsnpproblstName)
        {
            targetLstflg=true;
            extract_targets(&esdata, targetsnpproblstName,  prb_snp);
        }
        if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
       if(esdata._rowid.empty() && esdata._bxz.empty())
       {
           printf("ERROR: no data are included in the analysis.\n");
           exit(EXIT_FAILURE);
       }
        
       vector<SMRRLT> smrrlts;
       smr_heidi_func(smrrlts,  outFileName, &bdata,&gdata,&esdata,  cis_itvl,  heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt, ld_min,opt_hetero, sampleoverlap,pmecs, minCor,prb_snp,targetLstflg);
    }
    double heidi_test_ref_new(bInfo* bdata,SMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, int refid , double ld_min, int opt_hetero , bool sampleoverlap, double theta)
    {
        //refid is the id in smrwk for the target eQTL
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        string refrs=smrwk->rs[refid];
        double pthres=pchisq(threshold,1);
        //printf("Filtering SNPs (%ld in total) at eQTL p-value < %e for the HEIDI test.\n",smrwk->zxz.size(), pthres);
        for(int i=0;i<smrwk->zxz.size();i++)
        {
            if(i==refid || smrwk->zxz[i]*smrwk->zxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {
            
            //printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        //printf("%ld SNPs left after filtering.\n",sn_ids.size());
        update_snidx(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        
        if(find (sn_ids.begin(), sn_ids.end(), refid) == sn_ids.end()) sn_ids.push_back(refid);
        SMRWK smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
        long refid_heidi= find (smrwk_heidi.rs.begin(), smrwk_heidi.rs.end(), refrs) - smrwk_heidi.rs.begin();
        
        make_XMat(bdata,smrwk_heidi.curId, _X);
        //printf("Removing SNPs with LD r-squared between target-SNP %s > %f or < %f...\n",smrwk_heidi.rs[refid_heidi].c_str(),ldr2_top,ld_min);
        ld_calc_o2m(ld_v,refid_heidi,_X);
        
        if(fabs(ldr2_top-1)>1e-6 || ld_min > 0) {
            sn_ids.clear();
            for(int i=0;i<smrwk_heidi.zxz.size();i++)
            {
                if(i!= refid_heidi)
                {
                    double ldtmp=ld_v(i)*ld_v(i);
                    if((ldtmp-ldr2_top)<1e-6 && (ldtmp > ld_min)) sn_ids.push_back(i);
                }
                else{
                    sn_ids.push_back(i);
                }
            }
        }
        //printf("%ld SNPs are removed and %ld SNPs are retained.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
        if(sn_ids.size() < m_hetero) {
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        refid_heidi= find (smrwk_heidi.rs.begin(), smrwk_heidi.rs.end(), refrs) - smrwk_heidi.rs.begin();
        //printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
        int m = (int)smrwk_heidi.bxz.size();
        vector<int> rm_ID1;
        MatrixXd C;
        cor_calc(C, _X);
        double ld_top=sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
        //printf("%ld SNPs are removed and %ld SNPs (including the target SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[refid_heidi].c_str());
        if(m-rm_ID1.size() < m_hetero) {
            
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
            return -9;
        }
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sn_ids.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        
        update_snidx(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
        if(find (sn_ids.begin(), sn_ids.end(), refid_heidi) == sn_ids.end()){
            sn_ids[sn_ids.size()-1]=(int)refid_heidi; //in case of target SNP is not in top 20
        }
        if (sn_ids.size() < C.size()) { //Build new matrix
            MatrixXd D(sn_ids.size(),sn_ids.size());
            for (int i = 0 ; i < sn_ids.size() ; i++) {
                for (int j = 0 ; j < sn_ids.size() ; j++) {
                    D(i,j) = C(sn_ids[i],sn_ids[j]);
                }
            }
            C = D;
        }
        
        VectorXd _byz,_seyz, _bxz,_sexz,_zsxz;
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            _byz[j]=smrwk_heidi.byz[sn_ids[j]];
            _seyz[j]=smrwk_heidi.seyz[sn_ids[j]];
            _bxz[j]=smrwk_heidi.bxz[sn_ids[j]];
            _sexz[j]=smrwk_heidi.sexz[sn_ids[j]];
            _zsxz[j]=smrwk_heidi.zxz[sn_ids[j]];
        }
        
        nsnp = sn_ids.size();
        double pdev=-9;
        if(sampleoverlap) pdev=bxy_mltheter_so(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
        else pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp);
        
        return pdev;
    }
    void smr_heidi_trans_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, int trans_itvl, double p_trans, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, double ld_min,int opt_hetero)
    {
        bool sampleoverlap=false; double theta =0;
        // select trans-eQTL using --peqtl-trans. run SMR analysis using --peqtl-smr.
        uint64_t probNum = esdata->_include.size();
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        FILE* smr=NULL;
        long write_count=0;
        string outstr="";
        string smrfile="";
        if(outFileName!=NULL){
            smrfile = string(outFileName)+".smr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: error in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
        
        
        SMRWK smrwk;
        cout<<endl<<"Performing SMR analysis (SMR and HEIDI tests)..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            
            progr1=1.0*ii/probNum;
            if(progr1-progr0-0.05>1e-6 || ii+1==probNum)
            {
                if(ii+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            int i=esdata->_include[ii];
            int probebp=esdata->_epi_bp[i];
            int probechr=esdata->_epi_chr[i];
            string probename=esdata->_epi_prbID[i];
            string probegene=esdata->_epi_gene[i];
            char probeorien=esdata->_epi_orien[i];
            
            //get top trans-eQTLs
            vector<int> topTransBP;
            vector<int> topTransRowid;
            vector<int> topTransChr;
         
            if(esdata->_rowid.empty())
            {
                for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                    {
                        int snpbp=esdata->_esi_bp[j];
                        int snpchr=esdata->_esi_chr[j];
                        float beta=esdata->_bxz[i][j];
                        float se=esdata->_sexz[i][j];
                        if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>cis_itvl))
                        {
                            float zeqtl=beta/se;
                            float pval=pchisq(zeqtl*zeqtl, 1);
                            if(pval<=p_trans)
                            {
                                topTransBP.push_back(snpbp);
                                topTransRowid.push_back(j);
                                topTransChr.push_back(snpchr);
                            }
                            
                        }
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
                    int snpbp=esdata->_esi_bp[ge_rowid];
                    int snpchr=esdata->_esi_chr[ge_rowid];
                    float beta=esdata->_val[beta_start+j];
                    float se=esdata->_val[se_start+j];
                    
                    if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>cis_itvl))
                    {
                        float zeqtl=beta/se;
                        float pval=pchisq(zeqtl*zeqtl, 1);
                        if(pval<=p_trans)
                        {
                            topTransBP.push_back(snpbp);
                            topTransRowid.push_back(j); // here j is not row id in dense format. it is the value id
                            topTransChr.push_back(snpchr);
                        }
                    }
                }
            }

            printf("%ld trans-eQTLs have been detected of probe %s.\n", topTransBP.size(),probename.c_str());
            for(int jj=0;jj<topTransBP.size();jj++)
            {
                init_smr_wk(&smrwk);
                smrwk.cur_prbidx=i;
                smrwk.cur_chr=probechr;
                
                long maxid =fill_trans_smr_wk(bdata, gdata, esdata, &smrwk, topTransRowid, topTransBP,topTransChr,refSNP, cis_itvl,trans_itvl, heidioffFlag,jj);//return refSNP id or this trans-eQTL id
                if(refSNP!=NULL && maxid==-9)
                {
                    printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                    continue;
                }
                if (smrwk.bxz.size() == 0) {
                    
                    printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                    continue;
                }
                Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
                Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(maxid==-9) {
                    printf("Bug here. please report.\n");
                    exit(EXIT_FAILURE);
                }
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                if(refSNP==NULL && pxz_val>p_smr){
                    printf("WARNING: no trans-eQTL passed a p-value threshold (--peqtl-smr) %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                    continue;
                } else {
                    printf("Analysing probe %s...\n", probename.c_str());
                }
                double bxy_val = smrwk.byz[maxid] / smrwk.bxz[maxid];
                double sexy_val = sqrt((smrwk.seyz[maxid] * smrwk.seyz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] + smrwk.sexz[maxid] * smrwk.sexz[maxid] * smrwk.byz[maxid] * smrwk.byz[maxid]) / (smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                SMRRLT currlt;
                if(smr)
                {
                    outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + smrwk.rs[maxid] + '\t' + atos(smrwk.snpchrom[maxid]) + '\t' + atos(smrwk.bpsnp[maxid]) + '\t' + smrwk.allele1[maxid] + '\t' + smrwk.allele2[maxid] + '\t' + atos(smrwk.freq[maxid]) + '\t';
                    outstr += atos(smrwk.byz[maxid]) + '\t' + atos(smrwk.seyz[maxid]) + '\t' + dtos(smrwk.pyz[maxid]) + '\t';
                    outstr += atos(smrwk.bxz[maxid]) + '\t' + atos(smrwk.sexz[maxid]) + '\t' + dtos(pxz_val) + '\t';
                    outstr += atos(bxy_val) + '\t' + atos(sexy_val) + '\t' + dtos(pxy_val) + '\t';
                    
                } else {
                    currlt.ProbeID=probename;
                    currlt.ProbeChr=probechr;
                    currlt.Gene=probegene;
                    currlt.Probe_bp=probebp;
                    currlt.Orien=probeorien;
                    currlt.SNP=smrwk.rs[maxid];
                    currlt.SNP_Chr=smrwk.snpchrom[maxid];
                    currlt.SNP_bp=smrwk.bpsnp[maxid];
                    currlt.A1=smrwk.allele1[maxid];
                    currlt.A2=smrwk.allele2[maxid];
                    currlt.Freq=smrwk.freq[maxid];
                    currlt.b_GWAS=smrwk.byz[maxid];
                    currlt.se_GWAS=smrwk.seyz[maxid];
                    currlt.p_GWAS=smrwk.pyz[maxid];
                    currlt.b_eQTL=smrwk.bxz[maxid];
                    currlt.se_eQTL=smrwk.sexz[maxid];
                    currlt.p_eQTL=pxz_val;
                    currlt.b_SMR=bxy_val;
                    currlt.se_SMR=sexy_val;
                    currlt.p_SMR=pxy_val;
                }
                if(heidioffFlag || pxy_val>threshpsmrest)
                {
                    printf("INFO: the HEIDI test for probe %s is skipped because HEIDI test is turned off by the --heidi-off option or p_SMR does not pass the %e threshold.\n", probename.c_str(),threshpsmrest);
                    if(smr)
                    {
                        outstr+= string("NA") + '\t' + string("NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        write_count++;
                    } else {
                        currlt.p_HET=-9;
                        currlt.nsnp=-9;
                        smrrlts.push_back(currlt);
                    }
                    
                }
                else
                {
                    
                    long nsnp=-9;
                    double pdev=-9;
                    if(!heidioffFlag) {
                        if(new_heidi_mth) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid ,ld_min,opt_hetero,sampleoverlap, theta);
                        else pdev= heidi_test(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                    }
                    if(smr)
                    {
                        outstr+= (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        write_count++;
                    } else {
                        currlt.p_HET=pdev;
                        currlt.nsnp=(int)nsnp;
                        smrrlts.push_back(currlt);
                    }
                }
            }
        }
        
        if(smr){
            printf("Results of %ld probes have been saved in file %s.\n",write_count,smrfile.c_str());
            fclose(smr);
        } else {
            printf("Results of %ld probes have been returned.\n",smrrlts.size());
        }
        
    }
    // test all the meQTL or eQTL <=threshold.
    void smr_trans(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int cis_itvl,int trans_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest, bool new_het_mth, double p_smr,double ld_min)
    {
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false, 0, prbname);
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
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        
        vector<SMRRLT> smrrlts;
        smr_heidi_trans_func(smrrlts, outFileName, &bdata,&gdata,&esdata, cis_itvl, trans_itvl, p_trans, heidioffFlag, refSNP,p_hetero,ld_top,m_hetero , p_smr,threshpsmrest, new_het_mth,ld_min,opt_hetero);
        
    }
    void slct_trans_per_prb(vector<int> &slct_idx, vector<int> &regionChr, vector<long> &snpNumPerRegion,vector<long> &leftbound, vector<long> &rightbound, probeinfolst* prbifo, vector<info4trans> &snpinfo, long cis_itvl, long trans_itvl,double transThres)
    {
        
        map<string, int> rsa_map;
        map<string, int>::iterator iter;
        long rsNum=0;
        
        vector<bool> extend;
        vector<bool> merge;
        
        //extract info
        string probid=prbifo->probeId;
        long probbp=prbifo->bp;
        long probchr=prbifo->probechr;
        long cisuperBounder= probbp+cis_itvl;
        long cislowerBounder=((probbp-cis_itvl>0)?(probbp-cis_itvl):0);
        
        printf("Extracting trans information of the probe %s...\n", probid.c_str());
        for(int l=0;l<snpinfo.size();l++)
        {
            if(fabs(snpinfo[l].se+9)>1e-6)
            {
                double zsxz=snpinfo[l].beta/snpinfo[l].se;
                double pxz=pchisq(zsxz*zsxz, 1);
                
                if(!(snpinfo[l].snpchr == probchr && snpinfo[l].bp<=cisuperBounder && snpinfo[l].bp>=cislowerBounder) && pxz<=transThres)
                {
                    uint64_t transNum=0;
                    int curChr=snpinfo[l].snpchr;
                    string chckstr=string(snpinfo[l].snprs)+":"+string(snpinfo[l].a1)+":"+string(snpinfo[l].a2);
                    rsa_map.insert(pair<string, int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        transNum++;
                    }
                    
                    long transbp=snpinfo[l].bp;
                    long translowerBounder=((transbp-trans_itvl>0)?(transbp-trans_itvl):0);
                    long transuperBounder=transbp+trans_itvl;
                    bool extended=false;
                    bool merged=false;
                    
                    int startptr=l-1;
                    while(startptr>=0 && curChr == snpinfo[startptr].snpchr && transbp-snpinfo[startptr].bp<=trans_itvl)
                    {
                        if(fabs(snpinfo[startptr].se+9)>1e-6)
                        {
                            translowerBounder=snpinfo[startptr].bp;
                            if(rightbound.size()>0 && translowerBounder<=rightbound[rightbound.size()-1] && regionChr[rightbound.size()-1]==snpinfo[startptr].snpchr) //trans region merges
                            {
                                merged=true;
                                transNum=transNum+snpNumPerRegion[snpNumPerRegion.size()-1];
                                break;
                            }
                            if(snpinfo[startptr].snpchr == probchr && translowerBounder<=cisuperBounder && translowerBounder>=cislowerBounder) // trans touches cis region
                            {
                                break;
                            }
                            string chckstr=string(snpinfo[startptr].snprs)+":"+string(snpinfo[startptr].a1)+":"+string(snpinfo[startptr].a2);
                            
                                rsa_map.insert(pair<string, int>(chckstr,rsNum));
                                if(rsNum<rsa_map.size()){
                                    slct_idx.push_back(startptr);
                                    rsNum=rsa_map.size();
                                    transNum++;
                                }
                            
                        }
                        
                        startptr--;
                    }
                    startptr=l+1;
                    while (startptr<snpinfo.size() && curChr == snpinfo[startptr].snpchr && snpinfo[startptr].bp - transbp <= trans_itvl)
                    {
                        if(fabs(snpinfo[startptr].se+9)>1e-6)
                        {
                            transuperBounder=snpinfo[startptr].bp;
                            if(snpinfo[startptr].snpchr == probchr && transuperBounder>=cislowerBounder && transuperBounder<=cisuperBounder) // trans touches cis region
                            {
                                break;
                                
                            }  else {
                                
                                string chckstr=string(snpinfo[startptr].snprs)+":"+string(snpinfo[startptr].a1)+":"+string(snpinfo[startptr].a2);
                                rsa_map.insert(pair<string, int>(chckstr,rsNum));
                                if(rsNum<rsa_map.size()){
                                    slct_idx.push_back(startptr);
                                    rsNum=rsa_map.size();
                                    transNum++;
                                    double zsxz_tmp=snpinfo[startptr].beta/snpinfo[startptr].se;
                                    double pxz_tmp=pchisq(zsxz_tmp*zsxz_tmp, 1);
                                    if(pxz_tmp<transThres) // trans region extends
                                    {
                                        extended=true;
                                        transbp=snpinfo[startptr].bp;
                                    }
                                    l=startptr;
                                }
                            }
                        }
                        startptr++;
                    }
                    if(merged)
                    {
                        rightbound[rightbound.size()-1]=transuperBounder;
                        snpNumPerRegion[snpNumPerRegion.size()-1]=transNum;
                        extend[snpNumPerRegion.size()-1]=extended;
                        merge[snpNumPerRegion.size()-1]=merged;
                    }
                    else
                    {
                        regionChr.push_back(curChr);
                        rightbound.push_back(transuperBounder);
                        leftbound.push_back(translowerBounder);
                        snpNumPerRegion.push_back(transNum);
                        extend.push_back(extended);
                        merge.push_back(merged);
                    }
                    
                }
            }
        }
        long transnum=snpNumPerRegion.size();
        long transsnpnum=0;
        for(int l=0;l<snpNumPerRegion.size();l++) transsnpnum+=snpNumPerRegion[l];
       
        printf(" %ld SNPs in total %ld trans-region(s) have been extracted.\n",transsnpnum,transnum);
        //log
        if( snpNumPerRegion.size())
        {
            string logstr="{"+probid+","+atos(probchr)+","+atos(probbp)+"}\t";
            if(snpNumPerRegion.size()>0)
            {
                for(int h=0;h<snpNumPerRegion.size();h++)
                {
                    logstr+= "<"+atos(regionChr[h])+","+ atos(leftbound[h])+","+atos(rightbound[h])+","+atos(snpNumPerRegion[h])+">\t";
                }
            }else logstr+="<>\t";
            
            printf("%s.\n",logstr.c_str());

        }
    }
    int max_zsmr_id(SMRWK *smrwk , double p_smr)
    {
        printf("selecting candidiate instrument ...\n");
        vector<int> z_candid;
        for(int i=0;i<smrwk->bxz.size();i++)
        {
            double zxz=smrwk->bxz[i]/smrwk->sexz[i];
            double pxz = pchisq(zxz*zxz, 1);
            if(pxz<=p_smr) z_candid.push_back(i);
        }
        printf("%ld eQTLs found in this region that passed a p value threshold %e.\n", z_candid.size(),p_smr);
        
        int id=-9, idzx=-9;
        double mzxz2=-9,mzxy2=-9;
        for(int i=0;i<z_candid.size();i++)
        {
            int cid=z_candid[i];
            double zxz = smrwk->bxz[cid] / smrwk->sexz[cid];
            double zyz = smrwk->byz[cid] / smrwk->seyz[cid];
            double zxy2=(zxz*zxz*zyz*zyz)/(zxz*zxz+zyz*zyz);
            double zxz2=zxz*zxz;
            if( zxz2 > mzxz2)
            {
                mzxz2=zxz2;
                idzx=cid;
            }
            if( zxy2>mzxy2)
            {
                mzxy2=zxy2;
                id=cid;
            }
        }
        if(id!=-9 && idzx!=-9)
        {
            if(id==idzx)
            {
                printf("SNP %s is the top eQTL and the top SMR as well.\n", smrwk->rs[id].c_str());
            } else {
                printf("SNP %s is the top eQTL and SNP %s is the top SMR. We pick up SNP %s as the instrument.\n",smrwk->rs[idzx].c_str(),smrwk->rs[id].c_str(),smrwk->rs[id].c_str() );
            }

        } else {
            printf("NO eQTLs found in this region that passed a p value threshold %e.\n", p_smr);
        }
        return(id);
    }

    int max_zsmr_id(MTSMRWKEXP *smrwk , double p_smr)
    {
        printf("selecting candidiate instrument ...\n");
        vector<int> z_candid;
        for(int i=0;i<smrwk->bxz[0].size();i++)
        {
            double zxz=smrwk->bxz[0][i]/smrwk->sexz[0][i];
            double pxz = pchisq(zxz*zxz, 1);
            if(pxz<=p_smr) z_candid.push_back(i);
        }
        printf("%ld eQTLs found in this region that passed a p value threshold %e.\n", z_candid.size(),p_smr);
        
        int id=-9, idzx=-9;
        double mzxz2=-9,mzxy2=-9;
        for(int i=0;i<z_candid.size();i++)
        {
            int cid=z_candid[i];
            double zxz = smrwk->bxz[0][cid] / smrwk->sexz[0][cid];
            double zyz = smrwk->byz[cid] / smrwk->seyz[cid];
            double zxy2=(zxz*zxz*zyz*zyz)/(zxz*zxz+zyz*zyz);
            double zxz2=zxz*zxz;
            if( zxz2 > mzxz2)
            {
                mzxz2=zxz2;
                idzx=cid;
            }
            if( zxy2>mzxy2)
            {
                mzxy2=zxy2;
                id=cid;
            }
        }
        if(id!=-9 && idzx!=-9)
        {
            if(id==idzx)
            {
                printf("SNP %s is the top eQTL and the top SMR as well.\n", smrwk->rs[id].c_str());
            } else {
                printf("SNP %s is the top eQTL and SNP %s is the top SMR. We pick up SNP %s as the instrument.\n",smrwk->rs[idzx].c_str(),smrwk->rs[id].c_str(),smrwk->rs[id].c_str() );
            }

        } else {
            printf("NO eQTLs found in this region that passed a p value threshold %e.\n", p_smr);
        }
        return(id);
    }
    
    
    void smr_heidi_trans_region_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, int trans_itvl, double p_trans, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero)
    {
         bool sampleoverlap=false; double theta=0;
        // select trans-eQTL using --peqtl-trans. run SMR analysis using --peqtl-smr.
        uint64_t probNum = esdata->_include.size();
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        FILE* smr=NULL;
        long write_count=0;
        string outstr="";
        string smrfile="";
        if(outFileName!=NULL){
            smrfile = string(outFileName)+".smr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttrans_chr\ttrans_leftBound\ttrans_rightBound\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: error in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
        
        vector<info4trans> snpinfo;
        SMRWK smrwk;
        cout<<endl<<"Performing SMR analysis (SMR and HEIDI tests)..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            
            progr1=1.0*ii/probNum;
            if(progr1-progr0-0.05>1e-6 || ii+1==probNum)
            {
                if(ii+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            int i=esdata->_include[ii];
            
            probeinfolst prbifo;
            prbifo.bp=esdata->_epi_bp[i];
            prbifo.probechr=esdata->_epi_chr[i];
            strcpy2(&prbifo.probeId, esdata->_epi_prbID[i]);
            strcpy2(&prbifo.genename, esdata->_epi_gene[i]);
            prbifo.orien=esdata->_epi_orien[i];
            
            snpinfo.clear();
            
            if(esdata->_rowid.empty())
            {
                for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if (fabs(esdata->_bxz[i][j] + 9) > 1e-6)
                    {
                        if(gdata->seyz[j]+9>1e-6 && (esdata->_esi_chr[j]!=prbifo.probechr || (esdata->_esi_chr[j]==prbifo.probechr && ABS(prbifo.bp-esdata->_esi_bp[j])>cis_itvl)))
                        {
                            info4trans snpinfotmp;
                            strcpy2(&snpinfotmp.snprs, esdata->_esi_rs[j]);
                            snpinfotmp.snpchr=esdata->_esi_chr[j];
                            snpinfotmp.bp=esdata->_esi_bp[j];
                            snpinfotmp.gd=j;
                            strcpy2(&snpinfotmp.a1, esdata->_esi_allele1[j]);
                            strcpy2(&snpinfotmp.a2, esdata->_esi_allele2[j]);
                            if(!heidioffFlag) snpinfotmp.freq=bdata->_mu[bdata->_include[j]] / 2;
                            else snpinfotmp.freq=esdata->_esi_freq[j];
                            snpinfotmp.beta=esdata->_bxz[i][j];
                            snpinfotmp.se=esdata->_sexz[i][j];
                            snpinfotmp.byz=gdata->byz[j];
                            snpinfotmp.seyz=gdata->seyz[j];
                            snpinfotmp.pyz=gdata->pvalue[j];
                            
                            snpinfo.push_back(snpinfotmp);

                        }
                    }
                }
                
            }
            else{
                uint64_t beta_start=esdata->_cols[i<<1];
                uint64_t se_start=esdata->_cols[1+(i<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    //rowid in sparse is not ordered if there are trans region. because when we get trans-eQTL, we would reverse to get the SNPs in trans-region.
                    // it brings me a lot of trouble. I should have sorted the ids after sparse extraction of each probe.
                    int ge_rowid=esdata->_rowid[beta_start+j];
                    if(gdata->seyz[ge_rowid]+9>1e-6 && (esdata->_esi_chr[ge_rowid]!=prbifo.probechr || (esdata->_esi_chr[ge_rowid]==prbifo.probechr && ABS(prbifo.bp-esdata->_esi_bp[ge_rowid])>cis_itvl)))
                    {
                        info4trans snpinfotmp;
                        snpinfotmp.beta=esdata->_val[beta_start+j];
                        snpinfotmp.se=esdata->_val[se_start+j];
                        strcpy2(&snpinfotmp.snprs, esdata->_esi_rs[ge_rowid]);
                        snpinfotmp.snpchr=esdata->_esi_chr[ge_rowid];
                        snpinfotmp.bp=esdata->_esi_bp[ge_rowid];
                        snpinfotmp.gd=ge_rowid; //save snp id of the raw datastruct
                        strcpy2(&snpinfotmp.a1, esdata->_esi_allele1[ge_rowid]);
                        strcpy2(&snpinfotmp.a2, esdata->_esi_allele2[ge_rowid]);
                        if(!heidioffFlag) snpinfotmp.freq=bdata->_mu[bdata->_include[ge_rowid]] / 2;
                        else snpinfotmp.freq=esdata->_esi_freq[ge_rowid];
                        snpinfotmp.byz=gdata->byz[ge_rowid];
                        snpinfotmp.seyz=gdata->seyz[ge_rowid];
                        snpinfotmp.pyz=gdata->pvalue[ge_rowid];
                        
                        snpinfo.push_back(snpinfotmp);
                    }
                }
            }

            info4trans* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(info4trans),comp_i4tran);
            vector<long> snpNumPerRegion;
            vector<long> leftbound;
            vector<long> rightbound;
            vector<int> slct_idx;
            vector<int> regionChr;
            slct_trans_per_prb(slct_idx,regionChr, snpNumPerRegion,leftbound, rightbound,&prbifo, snpinfo,  cis_itvl,  trans_itvl, p_trans);
            // the output would be a little different from the sparse .summary file. because here the SNPs are the common SNPs of the three datasets. furthermore the trans-eqtl is not in common.
            if(snpNumPerRegion.size()==0) {
                 printf("WARNING: no trans-eQTL for probe %s at the p-value threshold %e.\n",prbifo.probeId,p_trans);
            }
            for(int jj=0;jj<snpNumPerRegion.size();jj++)
            {
                init_smr_wk(&smrwk);
                smrwk.cur_prbidx=i;
                smrwk.cur_chr=prbifo.probechr;
                long maxid = -9;
                for(int k=0;k<snpNumPerRegion[jj];k++)
                {
                    smrwk.bxz.push_back(snpinfo[slct_idx[k]].beta);
                    smrwk.sexz.push_back(snpinfo[slct_idx[k]].se);
                    smrwk.zxz.push_back(snpinfo[slct_idx[k]].beta/snpinfo[slct_idx[k]].se);
                    smrwk.byz.push_back(snpinfo[slct_idx[k]].byz);
                    smrwk.seyz.push_back(snpinfo[slct_idx[k]].seyz);
                    smrwk.pyz.push_back(snpinfo[slct_idx[k]].pyz);
                    smrwk.curId.push_back(snpinfo[slct_idx[k]].gd); //save snp id of the raw datastruct
                    smrwk.rs.push_back(snpinfo[slct_idx[k]].snprs);
                    smrwk.snpchrom.push_back(snpinfo[slct_idx[k]].snpchr);
                    smrwk.allele1.push_back(snpinfo[slct_idx[k]].a1);
                    smrwk.allele2.push_back(snpinfo[slct_idx[k]].a2);
                    if((refSNP!=NULL && snpinfo[slct_idx[k]].snprs==string(refSNP)) ) maxid=(smrwk.rs.size()-1);
                    smrwk.bpsnp.push_back(snpinfo[slct_idx[k]].bp);
                    smrwk.freq.push_back(snpinfo[slct_idx[k]].freq);
                    
                }
                
                if(refSNP!=NULL && maxid==-9)
                {
                    printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, prbifo.probeId);
                    continue;
                }
                if (smrwk.bxz.size() == 0) {
                    
                    printf("WARNING: no SNP fetched for probe %s.\n", prbifo.probeId);
                    continue;
                }
                Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
                Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(refSNP==NULL) {
                    if(opt) maxid=max_zsmr_id(&smrwk,p_smr);
                    else maxid=max_abs_id(zsxz);
                }
                if(maxid==-9) continue; //no eQTL in this region passed p_smr. (we seletct trans-region with p_trans which could be bigger than p_smr)
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                if(refSNP==NULL && pxz_val>p_smr){
                    printf("WARNING: current trans-eQTL with p-value %e doen's pass a p-value threshold (--peqtl-smr) %e for the SMR analysis for probe %s.\n", pxz_val, p_smr, prbifo.probeId);
                    continue;
                } else {
                    printf("Analysing probe %s...\n", prbifo.probeId);
                }
                double bxy_val = smrwk.byz[maxid] / smrwk.bxz[maxid];
                double sexy_val = sqrt((smrwk.seyz[maxid] * smrwk.seyz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] + smrwk.sexz[maxid] * smrwk.sexz[maxid] * smrwk.byz[maxid] * smrwk.byz[maxid]) / (smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                SMRRLT currlt;
                if(smr)
                {
                    outstr = atos(prbifo.probeId) + '\t' + atos(prbifo.probechr) + '\t' + atos(prbifo.genename) + '\t' + atos(prbifo.bp) +'\t'+atos(regionChr[jj])+'\t'+atos(leftbound[jj])+'\t'+atos(rightbound[jj]) + '\t' + smrwk.rs[maxid] + '\t' + atos(smrwk.snpchrom[maxid]) + '\t' + atos(smrwk.bpsnp[maxid]) + '\t' + smrwk.allele1[maxid] + '\t' + smrwk.allele2[maxid] + '\t' + atos(smrwk.freq[maxid]) + '\t';
                    outstr += atos(smrwk.byz[maxid]) + '\t' + atos(smrwk.seyz[maxid]) + '\t' + dtos(smrwk.pyz[maxid]) + '\t';
                    outstr += atos(smrwk.bxz[maxid]) + '\t' + atos(smrwk.sexz[maxid]) + '\t' + dtos(pxz_val) + '\t';
                    outstr += atos(bxy_val) + '\t' + atos(sexy_val) + '\t' + dtos(pxy_val) + '\t';
                    
                } else {
                    currlt.ProbeID=prbifo.probeId;
                    currlt.ProbeChr=prbifo.probechr;
                    currlt.Gene=prbifo.genename;
                    currlt.Probe_bp=prbifo.bp;
                    currlt.Orien=prbifo.orien;
                    currlt.SNP=smrwk.rs[maxid];
                    currlt.SNP_Chr=smrwk.snpchrom[maxid];
                    currlt.SNP_bp=smrwk.bpsnp[maxid];
                    currlt.A1=smrwk.allele1[maxid];
                    currlt.A2=smrwk.allele2[maxid];
                    currlt.Freq=smrwk.freq[maxid];
                    currlt.b_GWAS=smrwk.byz[maxid];
                    currlt.se_GWAS=smrwk.seyz[maxid];
                    currlt.p_GWAS=smrwk.pyz[maxid];
                    currlt.b_eQTL=smrwk.bxz[maxid];
                    currlt.se_eQTL=smrwk.sexz[maxid];
                    currlt.p_eQTL=pxz_val;
                    currlt.b_SMR=bxy_val;
                    currlt.se_SMR=sexy_val;
                    currlt.p_SMR=pxy_val;
                }
                if(heidioffFlag || pxy_val>threshpsmrest)
                {
                    printf("INFO: the HEIDI test for probe %s is skipped because HEIDI test is turned off by the --heidi-off option or p_SMR does not pass the %e threshold.\n", prbifo.probeId,threshpsmrest);
                    if(smr)
                    {
                        outstr+= string("NA") + '\t' + string("NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        write_count++;
                    } else {
                        currlt.p_HET=-9;
                        currlt.nsnp=-9;
                        smrrlts.push_back(currlt);
                    }
                    
                }
                else
                {
                    
                    long nsnp=-9;
                    double pdev=-9;
                    if(!heidioffFlag) {
                        if(new_heidi_mth) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid, ld_min,opt_hetero,sampleoverlap, theta);
                        else pdev= heidi_test(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                    }
                    if(smr)
                    {
                        outstr+= (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp+1) : "NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        write_count++;
                    } else {
                        currlt.p_HET=pdev;
                        currlt.nsnp=(int)nsnp;
                        smrrlts.push_back(currlt);
                    }
                }
            }
            free_snplist(snpinfo);
            free2(&prbifo.genename);
            free2(&prbifo.probeId);
        }
        if(smr){
            printf("Results of %ld probes have been saved in file %s.\n",write_count,smrfile.c_str());
            fclose(smr);
        } else {
            printf("Results of %ld probes have been returned.\n",smrrlts.size());
        }
        
    }

    void smr_trans_region(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero ,int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int cis_itvl,int trans_itvl,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest, bool new_het_mth, double p_smr, bool opt, double ld_min,double afthresh,double percenthresh)
    {
        setNbThreads(thread_num);
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        double threshold= chi_val(1,p_hetero);
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false, 0, prbname);
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
            if(forcefrqck)
            {
                double prop=freq_check(&bdata, &gdata, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
                
            }
            
        }else
        {
            allele_check(&gdata, &esdata);
            if(forcefrqck)
            {
                double prop=freq_check(&gdata, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
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
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        vector<SMRRLT> smrrlts;
        smr_heidi_trans_region_func(smrrlts, outFileName, &bdata,&gdata,&esdata, cis_itvl, trans_itvl, p_trans, heidioffFlag, refSNP,p_hetero,ld_top,m_hetero , p_smr,threshpsmrest, new_het_mth,opt,ld_min,opt_hetero);
        
    }
    
   void make_full_besd(char* outFileName, char* eqtlFileName, char* snplstName,char* problstName,bool bFlag,bool make_besd_flag, char* snplst2exclde, char* problst2exclde, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int addn)
    {
        
        eqtlInfo eqtlinfo;
        
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            
           read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.");
        
      
                filter_probe_null(&eqtlinfo); // at the same time, reset the vector _include // checked 20171120, filterring probe here is right.
                //filter_snp_null(&eqtlinfo);
                //update_besd();
                cout<<"\nsaving eQTL data..."<<endl;
                string esdfile = string(outFileName)+".esi";
                ofstream smr(esdfile.c_str());
                if (!smr) throw ("Error: can not open the ESI file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._esi_include.size(); i++) {
                    smr<<eqtlinfo._esi_chr[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_rs[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_gd[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_bp[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_allele1[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_allele2[eqtlinfo._esi_include[i]]<<'\t'<<eqtlinfo._esi_freq[eqtlinfo._esi_include[i]]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._esi_include.size()<<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
                
                esdfile = string(outFileName)+".epi";
                smr.open(esdfile.c_str());
                if (!smr) throw ("Error: can not open the EPI file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._include.size(); i++) {
                    smr<<eqtlinfo._epi_chr[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_prbID[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_gd[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_bp[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_gene[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_orien[eqtlinfo._include[i]]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._include.size()<<" probes have been saved in the file [" + esdfile + "]."<<endl;
            
                //if(make_besd_flag)
               // {
                    esdfile = string(outFileName)+".besd";
                    FILE * smr1;
                    smr1 = fopen (esdfile.c_str(), "wb");
                    if(eqtlinfo._valNum==0)
                    {
                        vector<int> ten_ints(RESERVEDUNITS);
                        ten_ints[0]=DENSE_FILE_TYPE_3;
                        if(addn!=-9)
                        {
                            printf("Adding the sample size  %d to the file %s.\n", addn, esdfile.c_str());
                            ten_ints[1]=addn;
                        }  else {
                            ten_ints[1]=-9;
                        }
                        ten_ints[2]=(int)eqtlinfo._esi_include.size();
                        ten_ints[3]=(int)eqtlinfo._include.size();
                        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                        fwrite(&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                        
                        uint64_t bsize=(eqtlinfo._include.size()*eqtlinfo._snpNum<<1);
                        float* buffer=(float*)malloc (sizeof(float)*bsize);
                        memset(buffer,0,sizeof(float)*bsize);
                        float* ptr=buffer;
                        uint64_t pro_num=eqtlinfo._include.size();
                        uint64_t snp_num=eqtlinfo._snpNum;
                        for(int i=0;i<pro_num;i++)
                        {
                            memcpy(ptr+(i<<1)*snp_num,&eqtlinfo._bxz[eqtlinfo._include[i]][0],sizeof(float)*snp_num);
                            memcpy(ptr+((i<<1)+1)*snp_num,&eqtlinfo._sexz[eqtlinfo._include[i]][0],sizeof(float)*snp_num);
                        }
                        fwrite (buffer,sizeof(float), bsize, smr1);
                        free(buffer);
                    }
                    else
                    {
                        vector<int> ten_ints(RESERVEDUNITS);
                        ten_ints[0]=SPARSE_FILE_TYPE_3;
                        if(addn!=-9)
                        {
                            printf("Adding the sample size  %d to the file %s.\n", addn, esdfile.c_str());
                            ten_ints[1]=addn;
                        }  else {
                            ten_ints[1]=-9;
                        }
                        ten_ints[2]=(int)eqtlinfo._esi_include.size();
                        ten_ints[3]=(int)eqtlinfo._include.size();
                        for(int i=4;i<RESERVEDUNITS;i++) ten_ints[i]=-9;
                        fwrite(&ten_ints[0],sizeof(int), RESERVEDUNITS, smr1);
                        
                        uint64_t colSize=sizeof(uint64_t)*((eqtlinfo._include.size()<<1)+1);
                        uint64_t rowSize=sizeof(uint32_t)*eqtlinfo._valNum;
                        uint64_t valSize=sizeof(float)*eqtlinfo._valNum;
                        uint64_t valNum=eqtlinfo._valNum;
                        uint64_t bufsize=sizeof(uint64_t)+colSize+rowSize+valSize;
                        
                        char* buffer=(char*)malloc (sizeof(char)*bufsize);
                        memset(buffer,0,sizeof(char)*bufsize);
                        char* wptr=buffer;
                        memcpy(wptr,&valNum,sizeof(uint64_t));
                        wptr+=sizeof(uint64_t);
                        uint64_t* uptr=(uint64_t*)wptr; *uptr++=0;
                        for(int i=0;i<eqtlinfo._include.size();i++)
                        {
                            *uptr++=eqtlinfo._cols[(eqtlinfo._include[i]<<1)+1];
                            *uptr++=eqtlinfo._cols[eqtlinfo._include[i]+1<<1];
                        }
                        wptr+=colSize;
                        memcpy(wptr,&eqtlinfo._rowid[0],rowSize);
                        wptr+=rowSize;
                        memcpy(wptr,&eqtlinfo._val[0],valSize);
                        fwrite (buffer,sizeof(char), bufsize, smr1);
                        free(buffer);
                        
                    }
                    fclose (smr1);
                    
                    cout<<"Effect sizes (beta) and SEfor "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
               // }
    
    }
    
 
    int get_besd_format(string besdfName)
    {
        string besdfile = besdfName + ".besd";
        char buf[8];
        ifstream besd(besdfile.c_str(), ios::in|ios::binary);
        if(!besd)
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     besdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        
        besd.read(buf, 4);
        besd.close();
        
        float* flag=(float *)buf;
        return  (int)*flag;

    }
  

    //sort in ascend order
    int comp(const void *a,const void *b){ return (((*(probeinfolst *)a).probechr>(*(probeinfolst *)b).probechr) || ( ((*(probeinfolst *)a).probechr ==(*(probeinfolst *)b).probechr) && ((*(probeinfolst *)a).bp > (*(probeinfolst *)b).bp) ))?1:-1; }
    int comp2(const void *a,const void *b){ return (((*(probeinfolst2 *)a).probechr>(*(probeinfolst2 *)b).probechr) || ( ((*(probeinfolst2 *)a).probechr ==(*(probeinfolst2 *)b).probechr) && ((*(probeinfolst2 *)a).bp > (*(probeinfolst2 *)b).bp) ))?1:-1; }
    int comp_esi(const void *a,const void *b){ return (((*(snpinfolst *)a).snpchr >(*(snpinfolst *)b).snpchr) || ( ((*(snpinfolst *)a).snpchr ==(*(snpinfolst *)b).snpchr) && ((*(snpinfolst *)a).bp > (*(snpinfolst *)b).bp) ))?1:-1; }
    int comp_i4tran(const void *a,const void *b){ return (((*(info4trans *)a).snpchr >(*(info4trans *)b).snpchr) || ( ((*(info4trans *)a).snpchr ==(*(info4trans *)b).snpchr) && ((*(info4trans *)a).bp > (*(info4trans *)b).bp) ))?1:-1; }
    int comp_estn(const void *a,const void *b){ return ((*(snpinfolst *)a).estn >(*(snpinfolst *)b).estn)?1:-1; }
    void free_snplist(vector<snpinfolst> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].snprs) free2(&a[i].snprs);
            if(a[i].a1) free2(&a[i].a1);
            if(a[i].a2) free2(&a[i].a2);
        }
    }
    void free_probelist(vector<probeinfolst> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].probeId) free2(&a[i].probeId);
            if(a[i].genename) free2(&a[i].genename);
            if(a[i].esdpath) free2(&a[i].esdpath);
            if(a[i].bfilepath) free2(&a[i].bfilepath);
        }
    }
    void free_snplist(vector<info4trans> &a){
        for(int i=0;i<a.size();i++)
        {
           if(a[i].snprs) free2(&a[i].snprs);
           if(a[i].a1) free2(&a[i].a1);
           if(a[i].a2) free2(&a[i].a2);
        }
    }
    void free_probelist(vector<probeinfolst2> &a){
        for(int i=0;i<a.size();i++)
        {
            if(a[i].probeId) free2(&a[i].probeId);
            if(a[i].genename) free2(&a[i].genename);
        }
    }
    
     void smr_e2e(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, bool ssmrflag,int expanWind,double ld_top_multi, char* targetsnpproblstName, char* snpproblstName,double afthresh,double percenthresh)
    {
        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetsnpproblstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetsnpproblstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&etrait, problstName);
        if(problst2exclde != NULL) exclude_prob(&etrait, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&etrait, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&etrait, oprobe);
        if(oproblst2exclde != NULL) exclude_prob(&etrait, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&etrait, oprobe2rm);
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); //no need here, already extracted in etrait
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        if(!heidioffFlag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            //if(snplstName != NULL) extract_snp(&bdata, snplstName);
            //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
            allele_check(&bdata, &etrait, &esdata);
            // if no snp left after check
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, &etrait, &esdata);
            }
            if(forcefrqck)
            {
                double prop=freq_check(&bdata, &etrait, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
            }
            
        }else
        {
            allele_check(&etrait, &esdata);
            if(forcefrqck)
            {
                double prop = freq_check( &etrait, &esdata,afthresh,percenthresh);
                if(prop > percenthresh)
                {
                    printf("ERROR: the analysis stopped because more than %0.2f%% of the SNPs were removed by the allele frequency difference check. You can change the proportion threshold by the flag --diff-freq-prop.\n",100*percenthresh);
                    exit(EXIT_FAILURE);
                }
            }
        }
       
        //the etrait is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(targetsnpproblstName)
        {
            targetLstflg=true;
            extract_targets(&esdata, targetsnpproblstName,  prb_snp);
        }
        if(snpproblstName) extract_targets(&esdata, snpproblstName,  prb_snp);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("ERROR: no data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        
      
        long itemcount=0;
        long etraitcount=0;
        string outstr="";
        string smrfile = string(outFileName)+".smr";
        if(ssmrflag) smrfile = string(outFileName)+".msmr";
        FILE* smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("ERROR: open error %s\n", smrfile.c_str());
                exit(1);
            }
        if(ssmrflag)  outstr="Expo_ID\tExpo_Chr\tExpo_Gene\tExpo_bp\tOutco_ID\tOutco_Chr\tOutco_Gene\tOutco_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_Outco\tse_Outco\tp_Outco\tb_Expo\tse_Expo\tp_Expo\tb_SMR\tse_SMR\tp_SMR\tp_SMR_multi\tp_HEIDI\tnsnp_HEIDI\n";
        else outstr="Expo_ID\tExpo_Chr\tExpo_Gene\tExpo_bp\tOutco_ID\tOutco_Chr\tOutco_Gene\tOutco_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_Outco\tse_Outco\tp_Outco\tb_Expo\tse_Expo\tp_Expo\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
       
        vector<int> includebk=esdata._include;
        int outcome_probe_wind=op_wind*1000;
        for( int ii=0;ii<etrait._probNum;ii++)
        {
            
            gwasData gdata;
            gdata.allele_1.resize(etrait._esi_include.size());
            gdata.allele_2.resize(etrait._esi_include.size());
            gdata.byz.resize(etrait._esi_include.size());
            gdata.seyz.resize(etrait._esi_include.size());
            gdata.freq.resize(etrait._esi_include.size());
            gdata.pvalue.resize(etrait._esi_include.size());
            gdata.splSize.resize(etrait._esi_include.size());
            gdata.snpName.resize(etrait._esi_include.size());
            
            for(int j=0;j<etrait._esi_include.size();j++){
                gdata.seyz[j]=-9;
                gdata.pvalue[j]=-9;
            }
            gdata._include.clear();
            
            string traitname=etrait._epi_prbID[ii];
            int traitchr=etrait._epi_chr[ii];
            string traitgene=etrait._epi_gene[ii];
            int traitbp=etrait._epi_bp[ii];
            cout<<"\nPerforming analysis of eTrait [ "+traitname+" ]..."<<endl;
            int count=0;
            if(etrait._rowid.empty())
            {
                for (int j = 0; j<etrait._esi_include.size(); j++)
                {
                    if (fabs(etrait._bxz[ii][etrait._esi_include[j]] + 9) > 1e-6)
                    {
                        gdata.byz[j]=etrait._bxz[ii][etrait._esi_include[j]];
                        gdata.seyz[j]=etrait._sexz[ii][etrait._esi_include[j]];
                        double z=etrait._bxz[ii][etrait._esi_include[j]]/etrait._sexz[ii][etrait._esi_include[j]];
                        double p=pchisq(z*z,1);
                        gdata.pvalue[j]=p;
                        gdata.snpName[j]=etrait._esi_rs[etrait._esi_include[j]];
                        gdata.allele_1[j]=etrait._esi_allele1[etrait._esi_include[j]];
                        gdata.allele_2[j]=etrait._esi_allele2[etrait._esi_include[j]];
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
                    int idx=(int)(find(etrait._esi_include.begin(), etrait._esi_include.end(), ge_rowid)-etrait._esi_include.begin());
                    if(idx<etrait._esi_include.size())
                    {
                        gdata.byz[idx]=etrait._val[beta_start+j];
                        gdata.seyz[idx]=etrait._val[se_start+j];
                        double z=etrait._val[beta_start+j]/etrait._val[se_start+j];
                        double p=pchisq(z*z,1);
                        gdata.pvalue[idx]=p;
                        gdata.snpName[idx]=etrait._esi_rs[ge_rowid];
                        gdata.allele_1[idx]=etrait._esi_allele1[ge_rowid];
                        gdata.allele_2[idx]=etrait._esi_allele2[ge_rowid];
                        gdata._include.push_back(idx);
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
            cout<<gdata.snpNum<<" common SNPs are included from eTrait [ "+traitname+" ] summary."<<endl;
           if(cis2all)
           {
               printf("\n%ld exposure probes are inclued in the analysis with the outcome probe %s.\n", esdata._include.size(),traitname.c_str());
               
           } else
           {
               int lowerbounder=(traitbp-outcome_probe_wind)>0?(traitbp-outcome_probe_wind):0;
               int upperbounder=traitbp+outcome_probe_wind;
               esdata._include.clear();
               for(int j=0;j<includebk.size();j++)
               {
                   int bptmp=esdata._epi_bp[includebk[j]];
                   if(esdata._epi_chr[includebk[j]]==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) esdata._include.push_back(includebk[j]);
               }
               printf("\n%ld exposure probes in the cis region [%d, %d] of outcome probe %s are inclued in the analysis.\n", esdata._include.size(),lowerbounder,upperbounder,traitname.c_str());
               if(esdata._include.size()==0) continue;
           }
            vector<SMRRLT> smrrlts;
            if(ssmrflag) ssmr_heidi_func(smrrlts,  NULL, &bdata,&gdata,&esdata,  cis_itvl,  heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest, ld_min,opt_hetero,expanWind,sampleoverlap,pmecs, minCor,ld_top_multi);
            else smr_heidi_func(smrrlts,  NULL, &bdata,&gdata,&esdata,  cis_itvl,  heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            if(smrrlts.size()>0){
                etraitcount++;
                itemcount+=smrrlts.size();
            }
            for(int j=0;j<smrrlts.size();j++)
            {
                outstr=smrrlts[j].ProbeID+'\t'+atos(smrrlts[j].ProbeChr)+'\t'+smrrlts[j].Gene+'\t'+atos(smrrlts[j].Probe_bp)+'\t'+traitname+'\t'+atos(traitchr)+'\t'+traitgene+'\t'+atos(traitbp)+'\t'+smrrlts[j].SNP+'\t'+atos(smrrlts[j].SNP_Chr)+'\t'+atos(smrrlts[j].SNP_bp)+'\t'+smrrlts[j].A1+'\t'+smrrlts[j].A2+'\t'+atos(smrrlts[j].Freq)+'\t'+atos(smrrlts[j].b_GWAS)+'\t'+atos(smrrlts[j].se_GWAS)+'\t'+dtos(smrrlts[j].p_GWAS)+'\t'+atos(smrrlts[j].b_eQTL)+'\t'+atos(smrrlts[j].se_eQTL)+'\t'+dtos(smrrlts[j].p_eQTL)+'\t'+atos(smrrlts[j].b_SMR)+'\t'+atos(smrrlts[j].se_SMR)+'\t'+dtos(smrrlts[j].p_SMR)+'\t';
                if(ssmrflag) outstr+=dtos(smrrlts[j].p_SSMR)+'\t'+(smrrlts[j].p_HET >= 0 ? dtos(smrrlts[j].p_HET) : "NA") + '\t' + (smrrlts[j].nsnp > 0 ? atos(smrrlts[j].nsnp+1) : "NA") + '\n';
                else outstr+=(smrrlts[j].p_HET >= 0 ? dtos(smrrlts[j].p_HET) : "NA") + '\t' + (smrrlts[j].nsnp > 0 ? atos(smrrlts[j].nsnp+1) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr))
                {
                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(smr);
        if(ssmrflag) printf("muti-SNP SMR and HEIDI analyses completed.\nSMR and heterogeneity analysis results of %ld exposure probes ( %ld outcome probes) have been saved in the file %s.\n",itemcount,etraitcount,smrfile.c_str());
        else printf("SMR and HEIDI analyses completed.\nSMR and heterogeneity analysis results of %ld exposure probes ( %ld outcome probes) have been saved in the file %s.\n",itemcount,etraitcount,smrfile.c_str());
    }

    // OPERA 17.07.2021 
    void multiexposurepi_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag, bool jointsmrflag, int cis_itvl,int piWind, int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included as xQTL data list and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. check input besd file format
        vector<string> besds, multi_bfiles;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }        
        
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // 3. expoNum = besdNum will be used; get prior variance 
        long expoNum; expoNum = besdNum;
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        bool operasmrflag = false;
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        }

        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();        
        // get the priors
        double sigma_def = 0.02;
        vector<string> sigmasplit;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }        
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum) throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        vector<double> sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }

        // 4. define global variables and extract the snp and probe data                
        vector<eqtlInfo> etrait(besdNum), etrait_sig(besdNum);
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        //extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        //extract the probe list
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
        }
        //read the besd
        for(int i=0;i<besdNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the OPERA analysis.\n");
               exit(EXIT_FAILURE);
           }
           // extract probe with a xQTL
           cis_xqtl_probe_include_only(&etrait[i], p_smr, cis_itvl, besds[i]);
        }

        // 5. select the indepedent (no overlap) genomic loci for stage 1 analysis
        // find the molecular trait with minimum num. of probes
        long minprbnum=etrait[0]._include.size();
        int minexpnum=0;
        for(int i=0;i<besdNum;i++) {            
            if(etrait[i]._include.size() < minprbnum) {                
                minprbnum = etrait[i]._include.size();
                minexpnum = i;
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // find the num. of independent loci as the num. of no overlap probes for molecular trait with minimum probes
        int exposure_probe_wind = op_wind*1000;
        int cis_itvl_wind = cis_itvl*1000;
        int indwin = piWind*1000;  // the epi_bp are required to be sorted
        vector<vector<int>> includepi(besdNum);
        includepi[minexpnum].push_back(etrait[minexpnum]._include[0]);
        int tmpbp = etrait[minexpnum]._epi_bp[etrait[minexpnum]._include[0]];
        for(int b=0;b<(etrait[minexpnum]._include.size()-1);b++) {
            int incldidx = etrait[minexpnum]._include[b+1];
            int distbp = abs(etrait[minexpnum]._epi_bp[incldidx] - tmpbp);
            if(distbp > indwin) {
                includepi[minexpnum].push_back(incldidx);
                tmpbp = etrait[minexpnum]._epi_bp[incldidx];
            }
        }
        // find the probes for each molecular trait that are within the window of the target probes
        for( int ii=0;ii<includepi[minexpnum].size();ii++)
        {
            int idxincld=includepi[minexpnum][ii];
            int probechr=etrait[minexpnum]._epi_chr[idxincld];
            int probebp=etrait[minexpnum]._epi_bp[idxincld];
            int lowerbounder=(probebp-indwin/2)>0?(probebp-indwin/2):0;
            int upperbounder=probebp+indwin/2;
           for(int i=0;i<besdNum;i++)
           {
               if(i!=minexpnum) {
                   for(int j=0;j<etrait[i]._include.size();j++)
                   {
                       int idxtmp = etrait[i]._include[j];
                       int bptmp = etrait[i]._epi_bp[idxtmp];
                       if(etrait[i]._epi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                            includepi[i].push_back(idxtmp);
                        }
                   }
               }
               
            }
        }
        // update the _include epi probes list for etrait
        for(int i=0;i<besdNum;i++) {
            getUnique(includepi[i]);
            printf("There are %ld probes included from exposure %ld.\n",includepi[i].size(),i+1);
            etrait[i]._include.clear();
            for(int l=0;l<includepi[i].size();l++) {
                    etrait[i]._include.push_back(includepi[i][l]);
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // select _esi_include SNPs list that are within 2Mb window of the target probes
        vector<vector<int>> slctsnpidx(besdNum);
        for(int i=0;i<besdNum;i++) {
            for( int k=0;k<etrait[i]._include.size();k++)
            {
                int idxtmp=etrait[i]._include[k];
                int probechr=etrait[i]._epi_chr[idxtmp];
                int probebp=etrait[i]._epi_bp[idxtmp];
                int lowerbounder=(probebp-cis_itvl_wind)>0?(probebp-cis_itvl_wind):0;
                int upperbounder=probebp+cis_itvl_wind;
               for(int j=0;j<etrait[i]._esi_include.size();j++)
               {
                   int idxtmp = etrait[i]._esi_include[j];
                   int bptmp = etrait[i]._esi_bp[idxtmp];
                   if(etrait[i]._esi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctsnpidx[i].push_back(idxtmp);}
               }
            }
        }
        // 6. allele checking between data
        // update the esi_include SNPs
        for(int i=0;i<besdNum;i++) {
            getUnique(slctsnpidx[i]);
            printf("There are %ld SNPs included from exposure %ld.\n",slctsnpidx[i].size(),i+1);
            etrait[i]._esi_include.clear();
            for(int m=0;m<slctsnpidx[i].size();m++) {
                    etrait[i]._esi_include.push_back(slctsnpidx[i][m]);
            }
            stable_sort(etrait[i]._esi_include.begin(),etrait[i]._esi_include.end());
        }
        // read the final besd file with updated esi_include and _include
        // for(int i=0;i<besdNum;i++) {
        //    read_besdfile(&etrait[i], string(besds[i])+".besd");
        //    if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
        //    {
        //        printf("ERROR: no data included in the analysis.\n");
        //        exit(EXIT_FAILURE);
        //    }
        // }
        
        // update the xQTL data with _esi_include and _include
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &etrait_sig[i]);
        }
        etrait.clear();

        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            }
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // allele checking between data
        if(!heidioffFlag || jointsmrflag)
        {
            map<string, string> snp_name_per_chr;
            if(bFileName!=NULL) read_famfile(&bdata, string(bFileName)+".fam");
            if(mbFileName!=NULL) read_multi_famfiles(&bdata, multi_bfiles);
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            if(bFileName!=NULL) read_bimfile(&bdata, string(bFileName)+".bim");
            if(mbFileName!=NULL) read_multi_bimfiles(&bdata, multi_bfiles, snp_name_per_chr);
            allele_check_multi_opt(&bdata, etrait_sig, &gdata1);
            //  allele_check_multi(&bdata, etrait_sig, &gdata1);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait_sig, &gdata1);
            }
        } else
        {
            allele_check_multi(etrait_sig,&gdata1);
        }
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata1); ngwas = median(gdata1.splSize);
        }
        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait_sig[i], &esdata[i]);
        }

        // 7. open the output file for pi estimate and sd
        string outstr = "", logpistr = "";
        string pifile = string(outFileName)+".pi";
        FILE* piiter = fopen(pifile.c_str(), "w");
        outstr="Posteriors\t"; logpistr = "Iteration\t";
        for(int i=0; i<combNum;i++)
        {
            outstr+="Pi"+atos(i+1)+"(";
            logpistr+="Pi"+atos(i+1)+"(";
            for(int j=0;j<expoNum;j++)
            {
                outstr+=atos(combins[i][j]);
                if(j<expoNum-1) outstr+=":";
                logpistr+=atos(combins[i][j]);
                if(j<expoNum-1) logpistr+=":";
            }
            if(i<(combNum-1)) outstr+=")\t";
            if(i<(combNum-1)) logpistr+=")\t";
        }
        outstr+=")\n"; logpistr+=")\n";
        if (!(piiter)) {
            printf("ERROR: open error %s\n", pifile.c_str());
            exit(1);
        }
        
        // 8. compute the pairwise SMR effect for all exposure probes
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
            }
        }
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposures with significant instruments (P < 5e-8) are less than the number of specified priors! Please check.\n");
            exit(EXIT_FAILURE);
        }

        // 9. sample a combination of exposure probes from each indepedent loci
        vector<vector<int>> includetmp(expoNum);
        vector<int> NAidx(probNum[minexpnum]);
        for(int ii=0;ii<probNum[minexpnum];ii++)
        {
            int traitchr=smrrlts[minexpnum][ii].ProbeChr;
            int traitbp=smrrlts[minexpnum][ii].Probe_bp;
            int lowerbounder=(traitbp-indwin/2)>0?(traitbp-indwin/2):0;
            int upperbounder=traitbp+indwin/2;
            NAidx[ii] = 0;
            for(int i=0;i<expoNum;i++)
            {
               vector<int> slctprbidx;
               for(int j=0;j<probNum[i];j++)
               {
                   int bptmp=smrrlts[i][j].Probe_bp;
                   if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctprbidx.push_back(j); 
                    }
               }
               if(slctprbidx.size()>0) {
                int randomIndex = rand()%slctprbidx.size();
                includetmp[i].push_back(slctprbidx[randomIndex]);
               } else { 
                    includetmp[i].push_back(-1);
                    NAidx[ii] = 1;
               }
            }
        }
        // remove these loci with missing exposures
        vector<vector<int>> includesmr(expoNum);
        for(int i=0;i<expoNum;i++) {
            for(int l=0;l<includetmp[minexpnum].size();l++) {
                if(NAidx[l]!=1) {
                    includesmr[i].push_back(includetmp[i][l]);
                }
            }
        }

        printf("\nPerforming joint SMR analysis ...\n");
        // joint SMR analysis
        vector<vector<SMRRLT>> smrrlts_joint_all;
        for(int ii=0;ii<includesmr[minexpnum].size();ii++) {
            vector<eqtlInfo> esdatacond(expoNum);
            vector<string> outconamec(besdNum);
            vector<vector<string>> prb_cojolist;
            int findindex = 0;
            for(int t=0; t<expoNum; t++) {
                long idx = includesmr[t][ii];
                outconamec[t] = smrrlts[t][idx].ProbeID;
                // find the target probe esdata
                esdata[t]._include.clear();
                map<string, int>::iterator itt;
                eqtlInfo esdatatmp;
                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                if(itt != esdata[t]._probe_name_map.end()) {
                    esdata[t]._include.push_back(itt->second);
                    e2econvert(&esdata[t], &esdatatmp);
                    esdatacond[t] = esdatatmp;
                    findindex+=1; 
                }
                if(targetcojosnplstName!=NULL) {
                    // find the target probe COJO signals
                    map<string, vector<string>>::iterator prb_pos;
                    prb_pos = prb_cojosnps.find(outconamec[t]);
                    vector<string> navector; navector.push_back("");
                    if(prb_pos!=prb_cojosnps.end()) {
                        prb_cojolist.push_back(prb_pos->second);
                    } else { prb_cojolist.push_back(navector); }
                }
            }
            if(findindex == expoNum) {
                vector<SMRRLT> smrrlts_joint;
                if(!operasmrflag) {
                    if(targetcojosnplstName!=NULL) {
                        multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatacond, ngwas, prb_cojolist, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    } else {multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatacond, ngwas, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);}
                } else {
                    for(int es=0; es<esdatacond.size(); es++) {
                        vector<SMRRLT> smrrlt_esdata;
                        smr_heidi_func(smrrlt_esdata, NULL, &bdata, &gdata1, &esdatacond[es],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
                        smrrlts_joint.push_back(smrrlt_esdata[0]);
                    }
                }          
                if(smrrlts_joint.size() == expoNum) {
                    smrrlts_joint_all.push_back(smrrlts_joint);
                }            
            }
        }

        printf("\nThere are %ld independent loci included in the estimation of the global pi.\n",smrrlts_joint_all.size());

        // 10. start the MCMC sampling with the indepedent loci SMR data
        // MCMC iteration start, define variables;//
        int nloops = 10000, nburn_in = 0.2 * nloops, nsamples = nloops - nburn_in;
        MatrixXd Pr(nloops,combNum);
        VectorXf ngamma(combNum);
        VectorXd alpha(combNum);
        VectorXf Pripi(combNum);
        VectorXd Primean(combNum), Prisd(combNum);
        Stat::Dirichlet Prob;
        printf("\nMCMC sampling starts ...\n");
        // initialize the alpha value as 0.1
        double sumalpha = 0;
        for(int i=0;i<combNum;i++) {
            alpha[i] = 0.1;
            sumalpha += 0.1;
        }
        // output the posterior samples in log file
        printf("\n%s", logpistr.c_str());
        for(int l=0;l<nloops;l++) {
            if(l==0) { //initialize the Pripi as even prior probability
                for(int i=0;i<combNum;i++) {
                    Pripi[i] = (float)1/combNum;
                }
            }
            // set starting values of ngamma as alpha, where ngamma is the sum of PP after sampling pi and alpha is the hyperparamenter
            for(int i=0;i<combNum;i++) {
                ngamma[i] = alpha[i];
            }
            // find the joint-SMR summary data for independent loci exposure combinations
            MatrixXd PP(smrrlts_joint_all.size(),combNum);
            for(int ii=0;ii<smrrlts_joint_all.size();ii++) {
                // Only one combination at each independent locus
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                vector<double> HH(combNum),PO(combNum);
                MatrixXd lh(2,expoNum);                    
                // find the summary-level joint-SMR data
                for(int t=0; t<expoNum; t++) {
                    bxy[t]=smrrlts_joint_all[ii][t].b_SMR;
                    sigma_e[t]=pow(smrrlts_joint_all[ii][t].se_SMR,2);
                    c[t]=1+sigma_e[t]/sigma_b[t];
                }
                // compute the marginal likelood under H0 and H1
                const double PI = 3.141592653589793238463;
                for(int t=0;t<expoNum;t++) {
                    lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                    lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                }
                // compute the joint data likelihood for H_gamma
                for(int i=0;i<combNum;i++) {
                    HH[i]=1.0;
                    for(int t=0;t<expoNum;t++)
                    {
                        HH[i] *= lh(combins[i][t],t);
                    }
                }
                // compute the normalizing constant 
                double POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i]=HH[i]*Pripi[i];
                    POall+=PO[i];
                }
                // compute the posterior probability
                if(POall > 0) {
                    for(int i=0;i<combNum;i++) {                    
                        PP(ii,i)=PO[i]/POall;
                        ngamma[i]+=PP(ii,i);
                    }
                }                                            
            }
            Pripi=Prob.sample(combNum,ngamma);
            logpistr=atos(l);
            for(int c=0;c<combNum;c++) {
                Pr(l,c)=Pripi[c];
                logpistr+='\t'+atos(Pr(l,c));
            }
            logpistr=logpistr+'\n';
            // print the sampled global pi each 100 iteration               
            if(l%100==0) {                                        
                printf("%s",(logpistr).c_str());
            }
        }
        
        // 11. output the posterior mean and SD
        Primean = Pr.bottomRows(nsamples).colwise().mean();
        for(int c=0;c<combNum;c++) {
            ArrayXd vectmp = Pr.bottomRows(nsamples).col(c);
            Prisd[c] = sqrt((vectmp - vectmp.mean()).square().sum()/(vectmp.size()-1));
        }
        outstr=outstr+"Mean"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Primean[c]);}; outstr=outstr+'\n';
        outstr=outstr+"SD"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Prisd[c]);}; outstr=outstr+'\n';
        if(fputs_checked(outstr.c_str(),piiter))
        {
            printf("ERROR: in writing file %s .\n", pifile.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(piiter);

    }

    void multiexposurepi_jointsmr_loci(char* outFileName, char* bFileName, char* mbFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag, bool jointsmrflag, int cis_itvl,int piWind, int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName,char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. checking besd file list
        vector<string> besds, multi_bfiles;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }
        
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // 3. expoNum = besdNum will be used; get prior variance 
        long expoNum; expoNum = besdNum;
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        bool operasmrflag = false;
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        }

        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();        
        // get the priors
        double sigma_def = 0.02;
        vector<string> sigmasplit;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        vector<double> sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }

        // 4. define global variables and extract the snp and probe data                
        vector<eqtlInfo> etrait(besdNum), etrait_sig(besdNum);
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;;
        
        printf("\nReading the exposure summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        //extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        //extract the probe list
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
        }
        //read the besd
        for(int i=0;i<besdNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the analysis.\n");
               exit(EXIT_FAILURE);
           }
           cis_xqtl_probe_include_only(&etrait[i], p_smr, cis_itvl, besds[i]);
        }

        // read the GWAS cojo signals
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
        }

        // 5. select the indepedent (no overlap) genomic loci for stage 1 analysis                
        // find the probes for each exposure that are within the window of the target probes
        int exposure_probe_wind=op_wind*1000;
        int indwin = piWind*1000;  // the epi_bp are required to be sorted
        vector<vector<int>> includepi(besdNum);
        for(int ii=0;ii<ldata._include.size();ii++)
        {
            int probechr=ldata._chr[ii];
            int probebp=ldata._bp[ii];
            int lowerbounder=(probebp-indwin/2)>0?(probebp-indwin/2):0;
            int upperbounder=probebp+indwin/2;
           for(int i=0;i<besdNum;i++)
           {
               for(int j=0;j<etrait[i]._include.size();j++)
               {
                   int idxtmp = etrait[i]._include[j];
                   int bptmp = etrait[i]._epi_bp[idxtmp];
                   if(etrait[i]._epi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        includepi[i].push_back(idxtmp);
                    }
               }
               
            }
        }
        // update the _include epi probes for etrait
        for(int i=0;i<besdNum;i++) {
            getUnique(includepi[i]);
            printf("There are %ld probes included for exposure %ld.\n",includepi[i].size(),i);
            etrait[i]._include.clear();
            for(int l=0;l<includepi[i].size();l++) {
                    etrait[i]._include.push_back(includepi[i][l]);
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // select _esi_include SNPs that are within 2Mb window of the target probes
        vector<vector<int>> slctsnpidx(besdNum);
        for(int i=0;i<besdNum;i++) {
            for( int k=0;k<etrait[i]._include.size();k++)
            {
                int idxtmp=etrait[i]._include[k];
                int probechr=etrait[i]._epi_chr[idxtmp];
                int probebp=etrait[i]._epi_bp[idxtmp];
                int lowerbounder=(probebp-exposure_probe_wind)>0?(probebp-exposure_probe_wind):0;
                int upperbounder=probebp+exposure_probe_wind;
               for(int j=0;j<etrait[i]._esi_include.size();j++)
               {
                   int idxtmp = etrait[i]._esi_include[j];
                   int bptmp = etrait[i]._esi_bp[idxtmp];
                   if(etrait[i]._esi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctsnpidx[i].push_back(idxtmp);}
               }
            }
        }

        // 6. allele checking between data
        // update the esi_include SNPs
        for(int i=0;i<besdNum;i++) {
            getUnique(slctsnpidx[i]);
            printf("There are %ld SNPs included for exposure %ld.\n",slctsnpidx[i].size(),i);
            etrait[i]._esi_include.clear();
            for(int m=0;m<slctsnpidx[i].size();m++) {
                    etrait[i]._esi_include.push_back(slctsnpidx[i][m]);
            }
            stable_sort(etrait[i]._esi_include.begin(),etrait[i]._esi_include.end());
        }
        // read the final besd file with updated esi_include and _include
        // for(int i=0;i<besdNum;i++) {
        //     //cout<<"Reading eQTL summary data..."<<endl;
        //    read_besdfile(&etrait[i], string(besds[i])+".besd");
        //    if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
        //    {
        //        printf("ERROR: no data included in the analysis.\n");
        //        exit(EXIT_FAILURE);
        //    }
        // }

        // extract probes that have at least a xQTL < p_smr
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &etrait_sig[i]);
        }
        etrait.clear();

        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            }
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // allele checking between data
        if(!heidioffFlag || jointsmrflag)
        {
            map<string, string> snp_name_per_chr;
            if(bFileName!=NULL) read_famfile(&bdata, string(bFileName)+".fam");
            if(mbFileName!=NULL) read_multi_famfiles(&bdata, multi_bfiles);
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            if(bFileName!=NULL) read_bimfile(&bdata, string(bFileName)+".bim");
            if(mbFileName!=NULL) read_multi_bimfiles(&bdata, multi_bfiles, snp_name_per_chr);
            allele_check_multi_opt(&bdata, etrait_sig, &gdata1);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait_sig, &gdata1);
            }
        } else
        {
              allele_check_multi(etrait_sig,&gdata1);
        }
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata1); ngwas = median(gdata1.splSize);
        }
        // update the SNPs after allele checking 
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait_sig[i], &esdata[i]);
        }

        // 7. open the output file for pi estimate and sd
        string outstr = "", logpistr = "";
        string pifile = string(outFileName)+".pi";
        FILE* piiter = fopen(pifile.c_str(), "w");
        outstr="Posteriors\t"; logpistr = "Iteration\t";
        for(int i=0; i<combNum;i++)
        {
            outstr+="Pi"+atos(i+1)+"(";
            logpistr+="Pi"+atos(i+1)+"(";
            for(int j=0;j<expoNum;j++)
            {
                outstr+=atos(combins[i][j]);
                if(j<expoNum-1) outstr+=":";
                logpistr+=atos(combins[i][j]);
                if(j<expoNum-1) logpistr+=":";
            }
            if(i<(combNum-1)) outstr+=")\t";
            if(i<(combNum-1)) logpistr+=")\t";
        }
        outstr+=")\n"; logpistr+=")\n";
        if (!(piiter)) {
            printf("ERROR: open error %s\n", pifile.c_str());
            exit(1);
        }

        // 8. compute the pairwise SMR effect for all exposure probes
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
            smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
            }
        }
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposures with significant instruments are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        } 

        // 9. sample a combination of exposure probes from each indepedent loci
        vector<vector<int>> includetmp(expoNum);
        vector<int> NAidx(ldata._include.size());
        for( int ii=0;ii<ldata._include.size();ii++)
        {
            int traitchr=ldata._chr[ii];
            int traitbp=ldata._bp[ii];
            int lowerbounder=(traitbp-indwin/2)>0?(traitbp-indwin/2):0;
            int upperbounder=traitbp+indwin/2;
            NAidx[ii] = 0;
            for(int i=0;i<expoNum;i++)
            {
               vector<int> slctprbidx;
               for(int j=0;j<probNum[i];j++)
               {
                   int bptmp=smrrlts[i][j].Probe_bp;
                   if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctprbidx.push_back(j); 
                    }
               }
               if(slctprbidx.size()>0) {
                int randomIndex = rand()%slctprbidx.size();
                includetmp[i].push_back(slctprbidx[randomIndex]);
               } else { 
                    includetmp[i].push_back(-1);
                    NAidx[ii] = 1;
               }
            }
        }
        // remove these loci with missing exposures
        vector<vector<int>> includesmr(expoNum);
        for(int i=0;i<expoNum;i++) {
            for(int l=0;l<ldata._include.size();l++) {
                if(NAidx[l]!=1) {
                    includesmr[i].push_back(includetmp[i][l]);
                }
            }
        }
        // joint SMR analysis
        vector<vector<SMRRLT>> smrrlts_joint_all;
        for(int ii=0;ii<includesmr[1].size();ii++) {
            vector<eqtlInfo> esdatacond(expoNum);
            vector<string> outconamec(besdNum);
            vector<vector<string>> prb_cojolist;
            int findindex = 0;
            for(int t=0; t<expoNum; t++) {
                long idx = includesmr[t][ii];
                outconamec[t] = smrrlts[t][idx].ProbeID;
                // find the target probe esdata
                esdata[t]._include.clear();
                map<string, int>::iterator itt;
                eqtlInfo esdatatmp;
                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                if(itt != esdata[t]._probe_name_map.end()) {
                    esdata[t]._include.push_back(itt->second);
                    e2econvert(&esdata[t], &esdatatmp);
                    esdatacond[t] = esdatatmp;
                    findindex+=1; 
                }
                if(targetcojosnplstName!=NULL) {
                    // find the target probe COJO signals
                    map<string, vector<string>>::iterator prb_pos;
                    prb_pos = prb_cojosnps.find(outconamec[t]);
                    vector<string> navector; navector.push_back("");
                    if(prb_pos!=prb_cojosnps.end()) {
                        prb_cojolist.push_back(prb_pos->second);
                    } else { prb_cojolist.push_back(navector); }
                }
            }
            if(findindex == expoNum) {
                vector<SMRRLT> smrrlts_joint;
                if(!operasmrflag) { 
                    if(targetcojosnplstName!=NULL) {
                    multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatacond, ngwas, prb_cojolist, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    } else {multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatacond, ngwas, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);}
                } else {
                    for(int es=0; es<esdatacond.size(); es++) {
                        vector<SMRRLT> smrrlt_esdata;
                        smr_heidi_func(smrrlt_esdata, NULL, &bdata, &gdata1, &esdatacond[es],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
                        smrrlts_joint.push_back(smrrlt_esdata[0]);
                    }
                }                                
                if(smrrlts_joint.size() == expoNum) {
                    smrrlts_joint_all.push_back(smrrlts_joint);
                }
            }
        }
        printf("\nThere are %ld independent loci included in the estimation of the global pi.\n",smrrlts_joint_all.size());
        
        // 10. start the MCMC sampling with the indepedent loci SMR data
        // MCMC iteration start, define variables;//
        int nloops = 10000, nburn_in = 0.2 * nloops, nsamples = nloops - nburn_in;
        MatrixXd Pr(nloops,combNum);
        VectorXf ngamma(combNum);
        VectorXd alpha(combNum);
        VectorXf Pripi(combNum);
        VectorXd Primean(combNum), Prisd(combNum);
        Stat::Dirichlet Prob;
        printf("\nMCMC sampling start ...\n");
        //initialize the alpha value as 0.1
        double sumalpha = 0;
        for(int i=0;i<combNum;i++) {
            alpha[i] = 0.1;
            sumalpha += 0.1;
        }
        // output the posterior samples in log file
        printf("\n%s", logpistr.c_str());
        for(int l=0;l<nloops;l++) {
            if(l==0) { //initialize the Pripi as even prior
                for(int i=0;i<combNum;i++) {
                    Pripi[i] = (float)1/combNum;
                }
            }
            // set starting values of ngamma as alpha, where ngamma is the sum of PP after sampling pi and alpha is the hyperparamenter
            for(int i=0;i<combNum;i++) {
                ngamma[i] = alpha[i];
            }
            // find the joint-SMR summary data for independent loci exposure combinations
            MatrixXd PP(smrrlts_joint_all.size(),combNum);
            for(int ii=0;ii<smrrlts_joint_all.size();ii++) {
                // Only one combination at each independent locus
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                vector<double> HH(combNum), PO(combNum);
                MatrixXd lh(2,expoNum);                    
                // find the summary-level joint-SMR data
                for(int t=0; t<expoNum; t++) {
                    bxy[t]=smrrlts_joint_all[ii][t].b_SMR;
                    sigma_e[t]=pow(smrrlts_joint_all[ii][t].se_SMR,2);
                    c[t]=1+sigma_e[t]/sigma_b[t];
                }
                // compute the marginal likelood under H0 and H1
                const double PI = 3.141592653589793238463;
                for(int t=0;t<expoNum;t++) {
                    lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                    lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                }
                // compute the joint likelihood for H_gamma
                for(int i=0;i<combNum;i++) {
                    HH[i]=1.0;
                    for(int t=0;t<expoNum;t++)
                    {
                        HH[i] *= lh(combins[i][t],t);
                    }
                }
                // compute the normalizing constant 
                double POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i]=HH[i]*Pripi[i];
                    POall+=PO[i];
                }
                // compute the posterior probability
                if(POall > 0) {
                    for(int i=0;i<combNum;i++) {
                        PP(ii,i)=PO[i]/POall;
                        ngamma[i]+=PP(ii,i);
                    }
                }                        
            }
            Pripi=Prob.sample(combNum,ngamma);
            logpistr=atos(l);
            for(int c=0;c<combNum;c++) {
                Pr(l,c)=Pripi[c];
                logpistr+='\t'+atos(Pr(l,c));
            }
            logpistr=logpistr+'\n';
            // print the sampled global pi                
            if(l%100==0) {                                        
                printf("%s",(logpistr).c_str());
            }            
        }

        // 11. output the posterior mean and SD
        Primean = Pr.bottomRows(nsamples).colwise().mean();
        for(int c=0;c<combNum;c++) {
            ArrayXd vectmp = Pr.bottomRows(nsamples).col(c);
            Prisd[c] = sqrt((vectmp - vectmp.mean()).square().sum()/(vectmp.size()-1));
        }
        outstr=outstr+"Mean"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Primean[c]);}; outstr=outstr+'\n';
        outstr=outstr+"SD"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Prisd[c]);}; outstr=outstr+'\n';
        if(fputs_checked(outstr.c_str(),piiter))
        {
            printf("ERROR: in writing file %s .\n", pifile.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(piiter);
         
    }

    void multiexposurepi(char* outFileName, char* bFileName, char* mbFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,int piWind, int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. check input besd file format
        vector<string> besds, multi_bfiles;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }        
        
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // 3. expoNum = besdNum will be used; get prior variance 
        long expoNum; expoNum = besdNum;
        printf("There are %ld exposure(s) and 1 outcome are included in OPERA.\n",expoNum);
        
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();        
        // get the priors
        double sigma_def = 0.02;
        vector<string> sigmasplit;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }        
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum) throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        vector<double> sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }

        // 4. define global variables and extract the snp and probe data        
        vector<eqtlInfo> etrait(besdNum), etrait_sig(besdNum);
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if(!heidioffFlag && bFileName == NULL && mbFileName == NULL) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        //extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        //extract the probe list
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
        }
        //read the besd
        for(int i=0;i<besdNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the analysis.\n");
               exit(EXIT_FAILURE);
           }
           cis_xqtl_probe_include_only(&etrait[i], p_smr, cis_itvl, besds[i]);
        }
        
        // 5. select the indepedent (no overlap) genomic loci for stage 1 analysis
        // find the molecular trait with minimum num. of probes
        long minprbnum=etrait[0]._include.size();
        int minexpnum=0;
        for(int i=0;i<besdNum;i++) {
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
            if(etrait[i]._include.size() < minprbnum) {
                minprbnum = etrait[i]._include.size();
                minexpnum = i;
            }
        }
        // find the num. of independent loci as the num. of no overlap probes for the smallest exposure
        int exposure_probe_wind=op_wind*1000;
        int indwin = piWind*1000;  // the epi_bp are required to be sorted
        vector<vector<int>> includepi(besdNum);
        includepi[minexpnum].push_back(etrait[minexpnum]._include[0]);
        int tmpbp = etrait[minexpnum]._epi_bp[etrait[minexpnum]._include[0]];
        for(int b=0;b<(etrait[minexpnum]._include.size()-1);b++) {
            int incldidx = etrait[minexpnum]._include[b+1];
            int distbp = abs(etrait[minexpnum]._epi_bp[incldidx] - tmpbp);
            if(distbp > indwin) {
                includepi[minexpnum].push_back(incldidx);
                tmpbp = etrait[minexpnum]._epi_bp[incldidx];
            }
        }
        // find the probes for each exposure that are within the window of the target probes
        for( int ii=0;ii<includepi[minexpnum].size();ii++)
        {
            int idxincld=includepi[minexpnum][ii];
            int probechr=etrait[minexpnum]._epi_chr[idxincld];
            int probebp=etrait[minexpnum]._epi_bp[idxincld];
            int lowerbounder=(probebp-indwin/2)>0?(probebp-indwin/2):0;
            int upperbounder=probebp+indwin/2;
           for(int i=0;i<besdNum;i++)
           {
               if(i!=minexpnum) {
                   for(int j=0;j<etrait[i]._include.size();j++)
                   {
                       int idxtmp = etrait[i]._include[j];
                       int bptmp = etrait[i]._epi_bp[idxtmp];
                       if(etrait[i]._epi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                            includepi[i].push_back(idxtmp);
                        }
                   }
               }
               
            }
        }
        // update the _include epi probes list for etrait
        for(int i=0;i<besdNum;i++) {
            getUnique(includepi[i]);
            printf("There are %ld probes included for exposure %ld.\n",includepi[i].size(),i+1);
            etrait[i]._include.clear();
            for(int l=0;l<includepi[i].size();l++) {
                    etrait[i]._include.push_back(includepi[i][l]);
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // select _esi_include SNPs that are within 2Mb window of the target probes
        vector<vector<int>> slctsnpidx(besdNum);
        for(int i=0;i<besdNum;i++) {
            for( int k=0;k<etrait[i]._include.size();k++)
            {
                int idxtmp=etrait[i]._include[k];
                int probechr=etrait[i]._epi_chr[idxtmp];
                int probebp=etrait[i]._epi_bp[idxtmp];
                int lowerbounder=(probebp-exposure_probe_wind)>0?(probebp-exposure_probe_wind):0;
                int upperbounder=probebp+exposure_probe_wind;
               for(int j=0;j<etrait[i]._esi_include.size();j++)
               {
                   int idxtmp = etrait[i]._esi_include[j];
                   int bptmp = etrait[i]._esi_bp[idxtmp];
                   if(etrait[i]._esi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctsnpidx[i].push_back(idxtmp);}
               }
            }
        }

        // 6. allele checking between data
        // update the _esi_include SNPs
        for(int i=0;i<besdNum;i++) {
            getUnique(slctsnpidx[i]);
            printf("There are %ld SNPs included for exposure %ld.\n",slctsnpidx[i].size(),i+1);
            etrait[i]._esi_include.clear();
            for(int m=0;m<slctsnpidx[i].size();m++) {
                    etrait[i]._esi_include.push_back(slctsnpidx[i][m]);
            }
            stable_sort(etrait[i]._esi_include.begin(),etrait[i]._esi_include.end());
        }
        // read the final besd file with updated esi_include and _include
        // for(int i=0;i<besdNum;i++) {
        //    read_besdfile(&etrait[i], string(besds[i])+".besd");
        //    if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
        //    {
        //        printf("ERROR: no data included in the analysis.\n");
        //        exit(EXIT_FAILURE);
        //    }
        // }

        // extract probes that have at least a xQTL < p_smr
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &etrait_sig[i]);
        }
        etrait.clear();

        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            }
        }        
        // allele checking between data
        if(!heidioffFlag)
        {
            map<string, string> snp_name_per_chr;
            if(bFileName!=NULL) read_famfile(&bdata, string(bFileName)+".fam");
            if(mbFileName!=NULL) read_multi_famfiles(&bdata, multi_bfiles);
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            if(bFileName!=NULL) read_bimfile(&bdata, string(bFileName)+".bim");
            if(mbFileName!=NULL) read_multi_bimfiles(&bdata, multi_bfiles, snp_name_per_chr);
            allele_check_multi_opt(&bdata, etrait_sig, &gdata1);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait_sig, &gdata1);
            }
        } else
        {
            allele_check_multi(etrait_sig,&gdata1);
        }
        if(gwasFileName!=NULL)  update_gwas(&gdata1);
        // update the SNPs after allele checking 
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait_sig[i], &esdata[i]);
        }

        // 7. open the output file for pi estimate and sd
        string outstr = "", logpistr = "";
        string pifile = string(outFileName)+".pi";
        FILE* piiter = fopen(pifile.c_str(), "w");
        outstr="Posteriors\t"; logpistr = "Iteration\t";
        for(int i=0; i<combNum;i++)
        {
            outstr+="Pi"+atos(i+1)+"(";
            logpistr+="Pi"+atos(i+1)+"(";
            for(int j=0;j<expoNum;j++)
            {
                outstr+=atos(combins[i][j]);
                if(j<expoNum-1) outstr+=":";
                logpistr+=atos(combins[i][j]);
                if(j<expoNum-1) logpistr+=":";
            }
            if(i<(combNum-1)) outstr+=")\t";
            if(i<(combNum-1)) logpistr+=")\t";
        }
        outstr+=")\n"; logpistr+=")\n";
        if (!(piiter)) {
            printf("ERROR: open error %s\n", pifile.c_str());
            exit(1);
        }

        // 8. compute the pairwise SMR effect for all exposure probes
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
            }
        }
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposures with significant instruments are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        } 

        // 9. sample a combination of exposure probes from each indepedent loci
        vector<vector<int>> includetmp(expoNum);
        vector<int> NAidx(probNum[minexpnum]);
        for( int ii=0;ii<probNum[minexpnum];ii++)
        {
            int traitchr=smrrlts[minexpnum][ii].ProbeChr;
            int traitbp=smrrlts[minexpnum][ii].Probe_bp;
            int lowerbounder=(traitbp-indwin/2)>0?(traitbp-indwin/2):0;
            int upperbounder=traitbp+indwin/2;
            NAidx[ii] = 0;
            for(int i=0;i<expoNum;i++)
            {
               vector<int> slctprbidx;
               for(int j=0;j<probNum[i];j++)
               {
                   int bptmp=smrrlts[i][j].Probe_bp;
                   if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctprbidx.push_back(j); 
                    }
               }
               if(slctprbidx.size()>0) {
                int randomIndex = rand()%slctprbidx.size();
                includetmp[i].push_back(slctprbidx[randomIndex]);
               } else { 
                    includetmp[i].push_back(-1);
                    NAidx[ii] = 1;
               }
            }
        }
        // remove these loci with missing exposures
        vector<vector<int>> includesmr(expoNum);
        for(int i=0;i<expoNum;i++) {
            for(int l=0;l<includetmp[minexpnum].size();l++) {
                if(NAidx[l]!=1) {
                    includesmr[i].push_back(includetmp[i][l]);
                }
            }
        }            
        printf("\nThere are %ld independent loci included in the estimation of the global pi.\n",includesmr[0].size());
        
        // 10. start the MCMC sampling with the indepedent loci SMR data
        // MCMC iteration start, define variables;//
        int nloops = 10000, nburn_in = 0.2 * nloops, nsamples = nloops - nburn_in;
        MatrixXd Pr(nloops,combNum);
        VectorXf ngamma(combNum);
        VectorXd alpha(combNum);
        VectorXf Pripi(combNum);
        VectorXd Primean(combNum), Prisd(combNum);
        Stat::Dirichlet Prob;
        printf("\nMCMC sampling start ...\n");
        //initialize the alpha value as 0.1
        double sumalpha = 0;
        for(int i=0;i<combNum;i++) {
            alpha[i] = 0.1;
            sumalpha += 0.1;
        }
        // output the posterior samples in log file
        printf("\n%s", logpistr.c_str());
        for(int l=0;l<nloops;l++) {
            if(l==0) { //initialize the Pripi as even prior
                for(int i=0;i<combNum;i++) {
                    Pripi[i] = (float)1/combNum;
                }
            }
            // set starting values of ngamma as alpha, where ngamma is the sum of PP after sampling pi and alpha is the hyperparamenter
            for(int i=0;i<combNum;i++) {
                ngamma[i] = alpha[i];
            }
            // find the SMR summary data for independent loci exposure combinations
            MatrixXd PP(includesmr[minexpnum].size(),combNum);
            for(int ii=0;ii<includesmr[minexpnum].size();ii++) {
                // Only one combination at each independent locus
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                vector<double> HH(combNum),PO(combNum);
                MatrixXd lh(2,expoNum);                    
                // get the summary-level SMR data
                for(int t=0; t<expoNum; t++) {
                    long idx = includesmr[t][ii];
                    bxy[t]=smrrlts[t][idx].b_SMR;
                    sigma_e[t]=pow(smrrlts[t][idx].se_SMR,2);
                    c[t]=1+sigma_e[t]/sigma_b[t];
                }
                // compute the marginal likelood under H0 and H1
                const double PI = 3.141592653589793238463;
                for(int t=0;t<expoNum;t++) {
                    lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                    lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                }
                // compute the joint data likelihood for H_gamma
                for(int i=0;i<combNum;i++) {
                    HH[i]=1.0;
                    for(int t=0;t<expoNum;t++)
                    {
                        HH[i] *= lh(combins[i][t],t);
                    }
                }
                // compute the normalizing constant
                double POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i]=HH[i]*Pripi[i];
                    POall+=PO[i];
                }
                // compute the posterior probability
                if(POall > 0) {
                    for(int i=0;i<combNum;i++) {
                        PP(ii,i)=PO[i]/POall;
                        ngamma[i]+=PP(ii,i);
                    }
                }                            
            }
            Pripi=Prob.sample(combNum,ngamma);
            logpistr=atos(l);
            for(int c=0;c<combNum;c++) {
                Pr(l,c)=Pripi[c];
                logpistr+='\t'+atos(Pr(l,c));
            }
            logpistr=logpistr+'\n';
            // output the sampled global pi each 100 iteration               
            if(l%100==0) {                                        
                printf("%s",(logpistr).c_str());
            }
            
        }

        // 11. output the posterior mean and SD
        Primean = Pr.bottomRows(nsamples).colwise().mean();
        for(int c=0;c<combNum;c++) {
            ArrayXd vectmp = Pr.bottomRows(nsamples).col(c);
            Prisd[c] = sqrt((vectmp - vectmp.mean()).square().sum()/(vectmp.size()-1));
        }
        outstr=outstr+"Mean"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Primean[c]);}; outstr=outstr+'\n';
        outstr=outstr+"SD"; for(int c=0;c<combNum;c++) { outstr+='\t'+atos(Prisd[c]);}; outstr=outstr+'\n';
        if(fputs_checked(outstr.c_str(),piiter))
        {
            printf("ERROR: in writing file %s .\n", pifile.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(piiter);

    }

    // not updated anymore
    void multiexposurepi_condsmr(char* outFileName, char* bFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag, bool condismrflag, int cis_itvl,int piWind, int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        // checking besd file list
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"eQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
        //printf("Checking the BESD format...\n");
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // this expoNum is the total number of exposures included
        long expoNum;
        expoNum = besdNum;
        printf("There are %ld exposures and 1 outcome are included in OPERA.\n",expoNum);
        if(expoNum<2) throw("Error: The program can not perform the conditional SMR analysis because there is only one exposure included. Please remove the flag --conditional-smr to analyze single exposure.");
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++){
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // get the priors
        vector<string> sigmasplit;
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        vector<double> sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        //define variables
        bool targetLstflg=false;
        map<string, string> prb_snp;
        vector<eqtlInfo> etrait(besdNum);
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        bool heidiFlag=false;
        
        printf("Reading the exposure summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL) || (condismrflag && bFileName == NULL)) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        //if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        //extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }

        //extract the probe list
        for(int i=0;i<besdNum;i++) {
            //cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
        }
        /////////////////////////////////////////////////
        /////////select the no overlap exposures/////////
        /////////////////////////////////////////////////
        // find the exposure with smallest num. of probes
        long minprbnum=etrait[0]._include.size();
        int minexpnum=0;
        for(int i=0;i<besdNum;i++) {
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
            if(etrait[i]._include.size() < minprbnum) {
                minprbnum = etrait[i]._include.size();
                minexpnum = i;
            }
        }
        // find the num. of independent loci as the num. of no overlap probes for the smallest exposure
        int exposure_probe_wind=op_wind*1000;
        int indwin = piWind*1000;  // the epi_bp are required to be sorted
        vector<vector<int>> includepi(besdNum);
        includepi[minexpnum].push_back(etrait[minexpnum]._include[0]);
        int tmpbp = etrait[minexpnum]._epi_bp[etrait[minexpnum]._include[0]];
        for(int b=0;b<(etrait[minexpnum]._include.size()-1);b++) {
            int incldidx = etrait[minexpnum]._include[b+1];
            int distbp = abs(etrait[minexpnum]._epi_bp[incldidx] - tmpbp);
            if(distbp > indwin) {
                includepi[minexpnum].push_back(incldidx);
                tmpbp = etrait[minexpnum]._epi_bp[incldidx];
            }
        }
        // find the probes for each exposure that are within the window of the target probes
        for( int ii=0;ii<includepi[minexpnum].size();ii++)
        {
            int idxincld=includepi[minexpnum][ii];
            int probechr=etrait[minexpnum]._epi_chr[idxincld];
            int probebp=etrait[minexpnum]._epi_bp[idxincld];
            int lowerbounder=(probebp-indwin/2)>0?(probebp-indwin/2):0;
            int upperbounder=probebp+indwin/2;
           for(int i=0;i<besdNum;i++)
           {
               if(i!=minexpnum) {
                   for(int j=0;j<etrait[i]._include.size();j++)
                   {
                       int idxtmp = etrait[i]._include[j];
                       int bptmp = etrait[i]._epi_bp[idxtmp];
                       if(etrait[i]._epi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                            includepi[i].push_back(idxtmp);
                        }
                   }
               }
               
            }
        }
        // update the _include epi probes for etrait
        for(int i=0;i<besdNum;i++) {
            getUnique(includepi[i]);
            printf("There are %ld probes included for exposure %ld.\n",includepi[i].size(),i);
            etrait[i]._include.clear();
            for(int l=0;l<includepi[i].size();l++) {
                    etrait[i]._include.push_back(includepi[i][l]);
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // select esi_include SNPs that are within 2Mb window of the target probes
        vector<vector<int>> slctsnpidx(besdNum);
        for(int i=0;i<besdNum;i++) {
            for( int k=0;k<etrait[i]._include.size();k++)
            {
                int idxtmp=etrait[i]._include[k];
                int probechr=etrait[i]._epi_chr[idxtmp];
                int probebp=etrait[i]._epi_bp[idxtmp];
                int lowerbounder=(probebp-exposure_probe_wind)>0?(probebp-exposure_probe_wind):0;
                int upperbounder=probebp+exposure_probe_wind;
               for(int j=0;j<etrait[i]._esi_include.size();j++)
               {
                   int idxtmp = etrait[i]._esi_include[j];
                   int bptmp = etrait[i]._esi_bp[idxtmp];
                   if(etrait[i]._esi_chr[idxtmp]==probechr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctsnpidx[i].push_back(idxtmp);}
               }
            }
        }
        // update the esi_include SNPs
        for(int i=0;i<besdNum;i++) {
            getUnique(slctsnpidx[i]);
            printf("There are %ld SNPs included for exposure %ld.\n",slctsnpidx[i].size(),i);
            etrait[i]._esi_include.clear();
            for(int m=0;m<slctsnpidx[i].size();m++) {
                    etrait[i]._esi_include.push_back(slctsnpidx[i][m]);
            }
            stable_sort(etrait[i]._esi_include.begin(),etrait[i]._esi_include.end());
        }
        // read the final besd file with updated esi_include and _include
        for(int i=0;i<besdNum;i++) {
            //cout<<"Reading eQTL summary data..."<<endl;
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the analysis.\n");
               exit(EXIT_FAILURE);
           }
        }
        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            }
        }        
        // allele checking between data
        if(!heidioffFlag || condismrflag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            allele_check_multi(&bdata, etrait, &gdata1);
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        } else
        {
              allele_check_multi(etrait,&gdata1);
        }
        if(gwasFileName!=NULL)  update_gwas(&gdata1);
        // update the SNPs after allele checking 
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        // open output file for pi iterations
        string outstr="";
        string pifile = string(outFileName)+".pi";
        FILE* piiter = fopen(pifile.c_str(), "w");
            outstr="Iteration\t";
            for(int i=0; i<combNum;i++)
            {
                outstr+="Pi"+atos(i+1)+"(";
                for(int j=0;j<expoNum;j++)
                {
                    outstr+=atos(combins[i][j]);
                    if(j<expoNum-1) outstr+=":";
                }
                if(i<(combNum-1)) outstr+=")\t";
            }
            outstr+=")\n";
            if (!(piiter)) {
                printf("ERROR: open error %s\n", pifile.c_str());
                exit(1);
            }
            if(fputs_checked(outstr.c_str(),piiter))
            {
                printf("ERROR: in writing file %s .\n", pifile.c_str());
                exit(EXIT_FAILURE);
            }
        // caculate the SMR effect size for all exposure probes
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
            smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
            } else {
                //printf("WARNING: No common SNP passed the p-value threshold %e for the SMR analysis of complex trait.\n", p_smr);
                continue;
            }
        }
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposures with significant instruments are less than the number of specified priors.\n");
            // exit(EXIT_FAILURE);
        } else {
            // sample the probe from each exposure at each independent loci
            vector<vector<int>> includetmp(expoNum);
            vector<int> NAidx(probNum[minexpnum]);
            for( int ii=0;ii<probNum[minexpnum];ii++)
            {
                int traitchr=smrrlts[minexpnum][ii].ProbeChr;
                int traitbp=smrrlts[minexpnum][ii].Probe_bp;
                int lowerbounder=(traitbp-indwin/2)>0?(traitbp-indwin/2):0;
                int upperbounder=traitbp+indwin/2;
                NAidx[ii] = 0;
                for(int i=0;i<expoNum;i++)
                {
                   vector<int> slctprbidx;
                   for(int j=0;j<probNum[i];j++)
                   {
                       int bptmp=smrrlts[i][j].Probe_bp;
                       if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                            slctprbidx.push_back(j); 
                        }
                   }
                   if(slctprbidx.size()>0) {
                    int randomIndex = rand()%slctprbidx.size();
                    includetmp[i].push_back(slctprbidx[randomIndex]);
                   } else { 
                        includetmp[i].push_back(-1);
                        NAidx[ii] = 1;
                   }
                }
            }
            // remove these loci with missing exposures
            vector<vector<int>> includesmr(expoNum);
            for(int i=0;i<expoNum;i++) {
                for(int l=0;l<includetmp[minexpnum].size();l++) {
                    if(NAidx[l]!=1) {
                        includesmr[i].push_back(includetmp[i][l]);
                    }
                }
            }
            // conditional SMR analysis
            vector<vector<SMRRLT>> smrrlts_condi_all;
            for(int ii=0;ii<includesmr[minexpnum].size();ii++) {
                vector<eqtlInfo> esdatacond(expoNum);
                vector<string> outconamec(besdNum);
                int findindex = 0;
                for(int t=0; t<expoNum; t++) {
                    long idx = includesmr[t][ii];
                    outconamec[t] = smrrlts[t][idx].ProbeID;
                    // find the target probe esdata
                    esdata[t]._include.clear();
                    map<string, int>::iterator itt;
                    eqtlInfo esdatatmp;
                    itt = esdata[t]._probe_name_map.find(outconamec[t]);
                    if(itt != esdata[t]._probe_name_map.end()) {
                        esdata[t]._include.push_back(itt->second);
                        e2econvert(&esdata[t], &esdatatmp);
                        esdatacond[t] = esdatatmp;
                        findindex+=1; 
                    }                    
                }
                if(findindex == expoNum) {
                    vector<SMRRLT> smrrlts_condi;
                    multi_cond_smr_func(smrrlts_condi, NULL, &bdata, &gdata1, esdatacond, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    if(smrrlts_condi.size() == expoNum) {
                        smrrlts_condi_all.push_back(smrrlts_condi);
                    }
                }
            }
            printf("\nThere are %ld independent loci included in the estimation of the global pi.\n",smrrlts_condi_all.size());
            ////////////////////////////////////////////
            // MCMC iteration start, define variables;//
            ////////////////////////////////////////////
            int nloops = 1000, nburn_in = 0.2 * nloops;
            MatrixXd Pr(nloops,combNum);
            VectorXf ngamma(combNum);
            VectorXd alpha(combNum);
            VectorXf Pripi(combNum);
            Stat::Dirichlet Prob;
            printf("MCMC sampling start ...\n");
            //initialize the alpha value as 1
            double sumalpha = 0;
            for(int i=0;i<combNum;i++) {
                //alpha[i] = (rand()%10 + 1);
                alpha[i] = 0.1;
                sumalpha += 0.1;
            }
            for(int l=0;l<nloops;l++) {
                if(l==0) { //initialize the Pripi as even prior
                    for(int i=0;i<combNum;i++) {
                        Pripi[i] = (float)1/combNum;
                    }
                }
                // if(l >= 1 && l <= nburn_in) {
                //      alpha = Pr.topRows(l).colwise().mean() * sumalpha;
                // }
                // if(l >= nburn_in) {
                //    alpha = Pr.topRows(nburn_in).colwise().mean() * sumalpha;
                // }
                // set ngamma as alpha as starting values, where ngamma is the sum of PP after sampling pi and alpha is the hyperparamenter
                for(int i=0;i<combNum;i++) {
                    ngamma[i] = alpha[i];
                }
                // find the SMR summary data for independent loci exposure combinations
                MatrixXd PP(smrrlts_condi_all.size(),combNum);
                for(int ii=0;ii<smrrlts_condi_all.size();ii++) {
                    // Only one combination at each independent locus
                    vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                    vector<double> HH(combNum),PO(combNum);
                    MatrixXd lh(2,expoNum);                    
                    // get the summary-level SMR data
                    string sampled_prb_name;
                    for(int t=0; t<expoNum; t++) {
                        //long idx = includesmr[t][ii];
                        bxy[t]=smrrlts_condi_all[ii][t].b_SMR;
                        sigma_e[t]=pow(smrrlts_condi_all[ii][t].se_SMR,2);
                        c[t]=1+sigma_e[t]/sigma_b[t];
                        // debug code
                      // sampled_prb_name += smrrlts[t][idx].ProbeID + '\t';
                    }
                        // compute the marginal likelood under H0 and H1
                        const double PI = 3.141592653589793238463;
                        for(int t=0;t<expoNum;t++) {
                            lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                            lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                        }
                        // compute the joint likelihood for H_gamma
                        for(int i=0;i<combNum;i++) {
                            HH[i]=1.0;
                            for(int t=0;t<expoNum;t++)
                            {
                                HH[i] *= lh(combins[i][t],t);
                            }
                        }
                        // compute the posterior odd
                        double POall = 0;
                        for(int i=0;i<combNum;i++) {
                            PO[i]=HH[i]*Pripi[i];
                            POall+=PO[i];
                        }
                        // compute the posterior probability
                     // if(POall>0){
                        for(int i=0;i<combNum;i++) {
                            PP(ii,i)=PO[i]/POall;
                            ngamma[i]+=PP(ii,i);
                        }
                    // }
                            
                }
                Pripi=Prob.sample(combNum,ngamma);
                // output the estimated global pi
                outstr=atos(l);
                for(int c=0;c<combNum;c++) {
                    Pr(l,c)=Pripi[c];
                    outstr+='\t'+atos(Pr(l,c));
                }
                outstr=outstr+'\n';
                if(fputs_checked(outstr.c_str(),piiter))
                {
                  printf("ERROR: in writing file %s .\n", pifile.c_str());
                  exit(EXIT_FAILURE);
                }
                
            }
        }         
        fclose(piiter);
    }
    
    // not updated anymore
    void multiexposuresmr_loci(char* outFileName, char* bFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP,char* refSNP, bool heidioffFlag,int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        // checking besd file list
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"eQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
        //printf("Checking the BESD format...\n");
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // this expoNum is the total number of exposures 
        long expoNum;
        expoNum = besdNum;
        printf("There are %ld exposures and 1 outcome are included in OPERA.\n",expoNum);
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++){
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // get the priors
        vector<string> priorsplit, sigmasplit;
        split_string(priorstr,priorsplit);
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of possible configurations.");
        vector<double> prior, sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        for(int t=0; t<combNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        vector<eqtlInfo> etrait(besdNum);        //vector<eqtlInfo*> etrait;
        vector<eqtlInfo> esdata(besdNum);
        //eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata1;
        bool heidiFlag=false;
        
        printf("Reading the exposure summary data file ...\n");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        //if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        printf("Reading the outcome summary data file ...\n");
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }

        //the etrait is not updated, so from now on _esi_include should be used always.
        for(int i=0;i<besdNum;i++) {
            //cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            //cout<<etrait[i]._rowid<<endl;
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }
        // Yang Wu read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // read the GWAS cojo signals
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
        }
        // allele checking between data
        if(!heidioffFlag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            allele_check_multi(&bdata, etrait, &gdata1);
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        }else
        {
            allele_check_multi(etrait,&gdata1);
        }
        if(gwasFileName!=NULL)  update_gwas(&gdata1);
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        // open multiple outcome files to write
        long itemcountsmr=0,itercountmlt=0;
        string outstr="";
        // GWAS trait
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL) {
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
        outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        string smrfile2 = string(outFileName)+".multismr";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
            if (!(smr2)) {
                printf("ERROR: open error %s\n", smrfile2.c_str());
                exit(1);
            }
            //outstr="Chr\tExpo_ID\tExpo_Gene\tExpo_bp\tOutco1_ID\tOutco1_Gene\tOutco1_bp\tOutco2_ID\tOutco2_Gene\tOutco2_bp\tOutco3_ID\tPP1(0:0:0)\tPP2(0:0:1)\tPP3(0:1:0)\tPP4(0:1:1)\tPP5(1:0:0)\tPP6(1:0:1)\tPP7(1:1:0)\tPP8(1:1:1)\n";
            outstr="Chr\t";
            int j = 0;
            for(int i=0; i<besdNum;i++)
            {
                j = i+1;
                outstr+="Expo"+atos(j)+"_ID"+'\t'+"Expo"+atos(j)+"_Gene"+'\t'+"Expo"+atos(j)+"_bp"+'\t';
            }
            //if(gwasFileName!=NULL) outstr+="Outco"+atos(j+1)+"_ID"+'\t';
        
            for(int i=0; i<combins.size();i++)
            {
                outstr+="PP"+atos(i+1)+"(";
                for(int j=0;j<expoNum;j++)
                {
                    outstr+=atos(combins[i][j]);
                    if(j<expoNum-1) outstr+=":";
                }
                if(i<(combins.size()-1)) outstr+=")\t";
            }
            outstr+=")\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr2))
            {
                printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                exit(EXIT_FAILURE);
            }           
        // compute SMR effect size for all exposure probes and output
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0){
            	smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                // for(int j=0;j<smrrltstmp.size();j++)
                // {
                //     outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                //     if(fputs_checked(outstr.c_str(),smr0))
                //     {
                //         printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                //         exit(EXIT_FAILURE);
                //     }
                // }                
            } else {
                continue;
            }
        }        
        // for( int k=0;k<probNum.size();k++)
        // {   
        //     itemcount = itemcount + probNum[k]; 
        // }
        // printf("SMR analysis results of %ld exposure probes have been saved in the file %s .\n",itemcount,smrfile0.c_str());
        // fclose(smr0);
        printf("\nPerforming multiple exposures OPERA analysis (including multi-exposure HEIDI tests) ... \n");
        if(probNum.size()!=expoNum){
            throw("ERROR: The number of exposures with significant instruments are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        }
        // loop with GWAS COJO loci
        double cr=0;
        for(int ii=0;ii<ldata._include.size();ii++)
        {   
            double desti=1.0*ii/(ldata._include.size());
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            vector<SMRRLT> smrrltsbf;
            vector<long> probNumbf;
            vector<eqtlInfo> esdatabf(expoNum);
            
            int locichr=ldata._chr[ii];
            int locibp=ldata._bp[ii];

            int lowerbounder=(locibp-exposure_probe_wind)>0?(locibp-exposure_probe_wind):0;
            int upperbounder=locibp+exposure_probe_wind;

            for(int i=0;i<besdNum;i++)
            {
                int countNum = 0;
                for(int j=0;j<probNum[i];j++)
                {
                   int bptmp=smrrlts[i][j].Probe_bp;
                   if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        smrrltsbf.push_back(smrrlts[i][j]); 
                        countNum = countNum + 1; itemcountsmr = itemcountsmr + 1;
                        // output the pairwise SMR results
                        outstr=smrrlts[i][j].ProbeID+'\t'+atos(smrrlts[i][j].ProbeChr)+'\t'+smrrlts[i][j].Gene+'\t'+atos(smrrlts[i][j].Probe_bp)+'\t'+smrrlts[i][j].SNP+'\t'+atos(smrrlts[i][j].SNP_Chr)+'\t'+atos(smrrlts[i][j].SNP_bp)+'\t'+smrrlts[i][j].A1+'\t'+smrrlts[i][j].A2+'\t'+atos(smrrlts[i][j].Freq)+'\t'+atos(smrrlts[i][j].b_GWAS)+'\t'+atos(smrrlts[i][j].se_GWAS)+'\t'+dtos(smrrlts[i][j].p_GWAS)+'\t'+atos(smrrlts[i][j].b_eQTL)+'\t'+atos(smrrlts[i][j].se_eQTL)+'\t'+dtos(smrrlts[i][j].p_eQTL)+'\t'+atos(smrrlts[i][j].b_SMR)+'\t'+atos(smrrlts[i][j].se_SMR)+'\t'+dtos(smrrlts[i][j].p_SMR)+'\t'+(smrrlts[i][j].p_HET >= 0 ? dtos(smrrlts[i][j].p_HET) : "NA") + '\t' + (smrrlts[i][j].nsnp > 0 ? atos(smrrlts[i][j].nsnp+1) : "NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr0))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                probNumbf.push_back(countNum);
                //printf("\n%ld probes from outcome%d in the cis region [%d, %d] of exposure probe %s are inclued in the analysis.\n", etrait[i]._include.size(),i+1,lowerbounder,upperbounder,traitname.c_str());
            }

            int expocout=0; 
            for(int i=0;i<probNumbf.size();i++){
                if(probNumbf[i]>0) expocout+=1;
            }

            if(expocout==expoNum) {
                // illustrate all the combinations
                vector<vector<int>> indexall;
                for(int i=0;i<probNumbf.size();i++){
                    vector<int> index(probNumbf[i]);
                    std::iota(index.begin(),index.end(),0);
                    indexall.push_back(index);
                }
                vector<vector<int>> combines;
                permute_vector(indexall, combines);
                //if (combines.size()==0) printf("\nWARNING: Less than %ld outcomes included in the analysis. \nThe multi-SMR and multi-HEIDI test will skip for exposure probe %s\n",expoNum,traitname.c_str());
                printf("\nThere are %ld possible combinations to test in a %ldKb window at %ldbp on chromosome %ld in the joint SMR model\n",combines.size(),op_wind,locibp,locichr);
                for(int i=0; i<combines.size();i++)
                {
                    vector<float> bxy(expoNum), sigma_e(expoNum),c(expoNum);
                    vector<string> outconamec(besdNum), outcogenec(besdNum); 
                    vector<long> outcobpc(besdNum);
                    vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum);
                    vector<gwasData> esdatain(expoNum);
                    MatrixXd lh(2,expoNum);
                    // get the probe and gene information for output
                    outstr=atos(locichr)+'\t';
                    long postmp = 0;                    
                    for(int t=0; t<besdNum; t++)
                    {   
                        long idxtmp = combines[i][t]+postmp;
                        outconamec[t] = smrrltsbf[idxtmp].ProbeID;
                        outcogenec[t] = smrrltsbf[idxtmp].Gene;
                        outcobpc[t] = smrrltsbf[idxtmp].Probe_bp;
                        postmp = postmp+probNumbf[t];
                        if(! heidioffFlag){
                            esdata[t]._include.clear();
                            map<string, int>::iterator itt;
                            eqtlInfo esdatatmp;
                            itt = esdata[t]._probe_name_map.find(outconamec[t]);
                            if(itt != esdata[t]._probe_name_map.end()){
                                esdata[t]._include.push_back(itt->second);
                                e2econvert(&esdata[t], &esdatatmp);
                                esdatabf[t] = esdatatmp;
                            }
                        }
                        outstr+=outconamec[t]+'\t'+outcogenec[t]+'\t'+atos(outcobpc[t])+'\t';
                    }
                    // Bayesian SMR test here
                    long pos=0;
                    for(int t=0; t<expoNum; t++)
                    {
                        long idx = combines[i][t]+pos;
                        bxy[t]=smrrltsbf[idx].b_SMR;
                        sigma_e[t]=pow(smrrltsbf[idx].se_SMR,2);
                        c[t]=1+sigma_e[t]/sigma_b[t];
                        pos=pos+probNumbf[t];
                    }
                    // multi-HEIDI test
                    vector<SMRRLT> smrrltsheidi;
                    if(! heidioffFlag){
                    	multi_heidi_func(smrrltsheidi, NULL, &bdata, &gdata1, esdatabf, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    }
                    // get the H0 and H1 prior pi and likelihood
                    const double PI = 3.141592653589793238463;
                    for(int t=0;t<expoNum;t++) {
                        lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                        lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                    }
                    // caculate the posterier probablity
                    for(int i=0;i<combNum;i++){
                        HH[i]=1.0;
                        for(int t=0;t<expoNum;t++)
                        {
                            HH[i] *= lh(combins[i][t],t);
                        }
                    }
                    float POall=0;
                    for(int i=0;i<combNum;i++){
                        PO[i]=HH[i]*prior[i];
                        POall+=PO[i];
                    }
                    for(int i=0;i<combNum;i++){
                        PP[i]=PO[i]/POall;
                        outstr=outstr+atos(PP[i])+'\t';
                    }
                    bool sigflag = false;
                    for(int i=1;i<combNum;i++) {
                        if(PP[i]>=thresh_PP) sigflag = true;
                    }
                    if(! heidioffFlag){
                    	outstr=outstr+(smrrltsheidi[0].p_HET >= 0 ? dtos(smrrltsheidi[0].p_HET) : "NA") + '\t' + (smrrltsheidi[0].nsnp > 0 ? atos(smrrltsheidi[0].nsnp+1) : "NA") + '\n';
                    } else{
                    	outstr=outstr + "NA" + '\t' + "NA" + '\n';
                    }
                    if(sigflag) {
                        itercountmlt+=1;
                        if(fputs_checked(outstr.c_str(),smr2))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                
                }

            } else {
                continue;
            }
            
        }
        fclose(smr0);            
        fclose(smr2);
        printf("\nPairwise SMR and HEIDI analyses for %ld exposure probes have been saved in the file %s .\n",itemcountsmr,smrfile0.c_str());
        printf("\nOPERA analyses for %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s.\n",expoNum,itercountmlt,smrfile2.c_str());
    }

    // not updated anymore
    void multiexposuresmr(char* outFileName, char* bFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        // checking besd file list
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"eQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
        //printf("Checking the BESD format...\n");
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // this outcoNum is the total number of exposures 
        long expoNum;
        expoNum = besdNum;
        printf("There are %ld exposures and 1 outcome are included in OPERA.\n",expoNum);
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++){
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // get the priors
        vector<string> priorsplit, sigmasplit;
        split_string(priorstr,priorsplit);
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of possible configurations.");
        vector<double> prior, sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        for(int t=0; t<combNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        vector<eqtlInfo> etrait(besdNum);        //vector<eqtlInfo*> etrait;
        vector<eqtlInfo> esdata(besdNum);
        //eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata1;
        bool heidiFlag=false;
        
        printf("Reading the exposure summary data file ...\n");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        //if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        printf("Reading the outcome summary data file ...\n");
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }

        //the etrait is not updated, so from now on _esi_include should be used always.
        for(int i=0;i<besdNum;i++) {
            //cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            //cout<<etrait[i]._rowid<<endl;
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }
        // Yang Wu read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // allele checking between data
        if(!heidioffFlag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            allele_check_multi(&bdata, etrait, &gdata1);
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        }else
        {
            allele_check_multi(etrait,&gdata1);
        }
        if(gwasFileName!=NULL)  update_gwas(&gdata1);
        //eqtlInfo esdatatmp;
        //e2econvert(&etrait[2], &esdatatmp);
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        //read_besdfile(&esdata, string(eqtlFileName)+".besd");
        // if(esdata._rowid.empty() && esdata._bxz.empty())
        // {
        //     printf("No data included in the analysis.\n");
        //     exit(EXIT_FAILURE);
        // }
        // open multiple outcome files to write
        long itemcount=0,itercountmlt=0;
        string outstr="";
        // GWAS trait
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL) {
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
        outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        string smrfile2 = string(outFileName)+".multismr";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
            if (!(smr2)) {
                printf("ERROR: open error %s\n", smrfile2.c_str());
                exit(1);
            }
            //outstr="Chr\tExpo_ID\tExpo_Gene\tExpo_bp\tOutco1_ID\tOutco1_Gene\tOutco1_bp\tOutco2_ID\tOutco2_Gene\tOutco2_bp\tOutco3_ID\tPP1(0:0:0)\tPP2(0:0:1)\tPP3(0:1:0)\tPP4(0:1:1)\tPP5(1:0:0)\tPP6(1:0:1)\tPP7(1:1:0)\tPP8(1:1:1)\n";
            outstr="Chr\t";
            int j = 0;
            for(int i=0; i<besdNum;i++)
            {
                j = i+1;
                outstr+="Expo"+atos(j)+"_ID"+'\t'+"Expo"+atos(j)+"_Gene"+'\t'+"Expo"+atos(j)+"_bp"+'\t';
            }
            //if(gwasFileName!=NULL) outstr+="Outco"+atos(j+1)+"_ID"+'\t';
        
            for(int i=0; i<combins.size();i++)
            {
                outstr+="PP"+atos(i+1)+"(";
                for(int j=0;j<expoNum;j++)
                {
                    outstr+=atos(combins[i][j]);
                    if(j<expoNum-1) outstr+=":";
                }
                if(i<(combins.size()-1)) outstr+=")\t";
            }
            outstr+=")\tp_HEIDI\tnsnp_HEIDI\n";
            //outstr+=")\n";
            if(fputs_checked(outstr.c_str(),smr2))
            {
                printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                exit(EXIT_FAILURE);
            }   
        // end of opening files
        int exposure_probe_wind=op_wind*1000;
        // loop with exposure probes
        //printf("\nPerforming multi-exposure-SMR analysis (multi-SMR and multi-HEIDI tests) for %d exposure probes... \n",esdata[0]._probNum);
        // caculate the SMR effect size for all exposure probes
        vector<vector<SMRRLT>> smrrlts;
        //vector<eqtlInfo> esdatall;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<besdNum;i++)
        {
            //eqtlInfo esdata;
            // e2econvert(&etrait[i],&esdata);
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0){
            smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                // esdatall.push_back(esdata[i]);
                for(int j=0;j<smrrltstmp.size();j++)
                {
                    outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr0))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                        exit(EXIT_FAILURE);
                    }
                }                
            } else {
                //printf("WARNING: No common SNP passed the p-value threshold %e for the SMR analysis of complex trait.\n", p_smr);
                continue;
            }
        }        

        for( int k=0;k<probNum.size();k++)
        {   
            itemcount = itemcount + probNum[k]; 
        }
        if(probNum.size()!=expoNum){
                throw("ERROR: The number of exposures with significant instruments are less than the number of specified priors.\n");
                // exit(EXIT_FAILURE);
        } else {
            double cr=0;
            for( int ii=0;ii<probNum[0];ii++)
            {   
                double desti=1.0*ii/(probNum[0]);
                if(desti>=cr)
                {
                    printf("%3.0f%%\r", 100.0*desti);
                    fflush(stdout);
                    if(cr==0) cr+=0.05;
                    else if(cr==0.05) cr+=0.2;
                    else if(cr==0.25) cr+=0.5;
                    else cr+=0.25;
                }
                vector<SMRRLT> smrrltsbf;
                vector<long> probNumbf;
                smrrltsbf.push_back(smrrlts[0][ii]); 
                probNumbf.push_back(1);
                // //gwasData gdata;
                // etrait[1]._include.clear();
                // etrait[1]._include.push_back(ii);
                string traitname=smrrlts[0][ii].ProbeID;
                int traitchr=smrrlts[0][ii].ProbeChr;
                string traitgene=smrrlts[0][ii].Gene;
                int traitbp=smrrlts[0][ii].Probe_bp;
                vector<eqtlInfo> esdatabf(expoNum);
                
                if(! heidioffFlag){
                    map<string, int>::iterator itt;
                    eqtlInfo esdatatmp;
                    esdata[0]._include.clear();
                    itt = esdata[0]._probe_name_map.find(traitname);
                    if(itt != esdata[0]._probe_name_map.end()){
                        esdata[0]._include.push_back(itt->second);
                        e2econvert(&esdata[0], &esdatatmp);
                        esdatabf[0]=esdatatmp;
                    } else {
                        continue;
                    }
                }
                //cout<<"\nPerforming analysis of exposure [ "+traitname+" ]..."<<endl;
        
               if(cis2all)
               {
                   //printf("\n%ld outcome probes are inclued in the analysis with the exposure probe %s.\n", esdata._include.size(),traitname.c_str());
                   
               } else
               {
                   int lowerbounder=(traitbp-exposure_probe_wind)>0?(traitbp-exposure_probe_wind):0;
                   int upperbounder=traitbp+exposure_probe_wind;

                   for(int i=1;i<besdNum;i++)
                   {
                       int countNum = 0;
                       //esdata[i]._include.clear();
                       for(int j=0;j<probNum[i];j++)
                       {
                           int bptmp=smrrlts[i][j].Probe_bp;
                           if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                                smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                                //esdatabf.push_back(esdatall[i][j]);
                            }
                       }
                       probNumbf.push_back(countNum);
                       //printf("\n%ld probes from outcome%d in the cis region [%d, %d] of exposure probe %s are inclued in the analysis.\n", etrait[i]._include.size(),i+1,lowerbounder,upperbounder,traitname.c_str());
                    }
                }

                int outcomcout=0; 
                for(int i=0;i<probNumbf.size();i++){
                    if(probNumbf[i]>0) outcomcout+=1;
                }

                if(outcomcout==expoNum) {
                    // illustrate all the combinations
                    vector<vector<int>> indexall;
                    for(int i=0;i<probNumbf.size();i++){
                        vector<int> index(probNumbf[i]);
                        std::iota(index.begin(),index.end(),0);
                        indexall.push_back(index);
                    }
                    vector<vector<int>> combines;
                    permute_vector(indexall, combines);
                    // Bayesian factor caculation
                    //printf("\nPerforming Multi-SMR tests.....\n");
                    //if (combines.size()==0) printf("\nWARNING: Less than %ld outcomes included in the analysis. \nThe multi-SMR and multi-HEIDI test will skip for exposure probe %s\n",expoNum,traitname.c_str());
                    for(int i=0; i<combines.size();i++)
                    {
                        vector<float> bxy(expoNum), sigma_e(expoNum),c(expoNum);
                        vector<string> outconamec(besdNum), outcogenec(besdNum); vector<long> outcobpc(besdNum);
                        vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum);
                        vector<gwasData> esdatain(expoNum);
                        MatrixXd lh(2,expoNum);
                        // get the probe and gene information for output
                        outstr=atos(traitchr)+'\t'+traitname+'\t'+traitgene+'\t'+atos(traitbp)+'\t';
                        long postmp=1;
                        
                        for(int t=1; t<besdNum; t++)
                        {   
                            long idxtmp = combines[i][t]+postmp;
                            outconamec[t] = smrrltsbf[idxtmp].ProbeID;
                            outcogenec[t] = smrrltsbf[idxtmp].Gene;
                            outcobpc[t] = smrrltsbf[idxtmp].Probe_bp;
                            postmp = postmp+probNumbf[t];
                            if(! heidioffFlag){
                                esdata[t]._include.clear();
                                map<string, int>::iterator itt;
                                eqtlInfo esdatatmp;
                                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                                if(itt != esdata[t]._probe_name_map.end()){
                                    esdata[t]._include.push_back(itt->second);
                                    e2econvert(&esdata[t], &esdatatmp);
                                    esdatabf[t] = esdatatmp;
                                }
                            }
                            outstr+=outconamec[t]+'\t'+outcogenec[t]+'\t'+atos(outcobpc[t])+'\t';
                        }
                        //if(gwasFileName!=NULL) outstr+="Trait\t";
                        //get the bxy, sigma_b and sigma_e from SMR
                        long pos=0;
                        for(int t=0; t<expoNum; t++)
                        {
                            long idx = combines[i][t]+pos;
                            bxy[t]=smrrltsbf[idx].b_SMR;
                            sigma_e[t]=pow(smrrltsbf[idx].se_SMR,2);
                            c[t]=1+sigma_e[t]/sigma_b[t];
                            //esdatain[t]=esdatabf[idx];
                            pos=pos+probNumbf[t];
                        }
                        // multi-HEIDI test
                        vector<SMRRLT> smrrltsheidi;
                        if(! heidioffFlag){
                        multi_heidi_func(smrrltsheidi, NULL, &bdata, &gdata1, esdatabf, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                        }
                        // get the H0 and H1 prior pi and likelihood
                        const double PI = 3.141592653589793238463;
                        for(int t=0;t<expoNum;t++) {
                            lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                            lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                        }
                        // caculate the posterier probablity
                        for(int i=0;i<combNum;i++){
                            HH[i]=1.0;
                            for(int t=0;t<expoNum;t++)
                            {
                                HH[i] *= lh(combins[i][t],t);
                            }
                        }
                        float POall=0;
                        for(int i=0;i<combNum;i++){
                            PO[i]=HH[i]*prior[i];
                            POall+=PO[i];
                        }
                        for(int i=0;i<combNum;i++){
                            PP[i]=PO[i]/POall;
                            outstr=outstr+atos(PP[i])+'\t';
                        }
                        bool sigflag = false;
                        for(int i=1;i<combNum;i++) {
                            if(PP[i]>=0.8) sigflag = true;
                        }
                        if(! heidioffFlag){
                        outstr=outstr+(smrrltsheidi[0].p_HET >= 0 ? dtos(smrrltsheidi[0].p_HET) : "NA") + '\t' + (smrrltsheidi[0].nsnp > 0 ? atos(smrrltsheidi[0].nsnp+1) : "NA") + '\n';
                        } else{
                        outstr=outstr + "NA" + '\t' + "NA" + '\n';
                        }
                        if(sigflag) {
                            itercountmlt+=1;
                            if(fputs_checked(outstr.c_str(),smr2))
                            {
                                printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                                exit(EXIT_FAILURE);
                            }
                        }
                    
                    }

                } else {
                    continue;
                }
                
            }
        }         
        fclose(smr0);
        fclose(smr2);
        //printf("\nSMR and HEIDI analyses for molecular traits completed.\nSMR and heterogeneity analysis results of %ld outcome probes (%ld exposure probe) have been saved in the file %s.\n",itemcount,etraitcount,smrfile1.c_str());
        printf("\nOPERA analyses for %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations (%ld exposure probes) have been saved in the file %s.\n",expoNum,itercountmlt,itemcount,smrfile2.c_str());
    }

    void multiexposure_jointsmr_loci(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP, double thresh_smr,double thresh_heidi,char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. check input besd file format
        vector<string> besds, multi_bfiles;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }
        
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // 3. expoNum = besdNum will be used; get prior variance and PIP header
        long expoNum;  expoNum = besdNum;
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure file included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        } 
            
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // list marginal PIP combins
        vector<vector<int>> combmarg, comborg(combins.size()), idxmarg(combins.size());
        comborg[0].push_back(0); idxmarg[0].push_back(0);
        for(int i=1;i<combNum;i++) {
            for(int j=0;j<expoNum;j++) {
                if(combins[i][j] == 1) comborg[i].push_back(j+1);
            }
        }
        // reorder the output PIP
        for(int i=1;i<=expoNum;i++) {
            vector<vector<int>> combtmp;
            for(int j=0;j<combNum;j++) {
                if(comborg[j].size()==i) combtmp.push_back(comborg[j]);
            }
            sort(combtmp.begin(),combtmp.end());
            for(int k=0;k<combtmp.size();k++) {
                combmarg.push_back(combtmp[k]);
            }
        }
        // find the PIP correspoding configuration index
        for(int i=1;i<combmarg.size();i++) {
            for(int c=0;c<combNum;c++) {
                int sumcount=0;
                for(int j=0;j<combmarg[i].size();j++) {
                    int tmpidx = combmarg[i][j] - 1;
                    sumcount += combins[c][tmpidx];
                }
                if(sumcount==combmarg[i].size()) idxmarg[i].push_back(c);
            }
        }
        // get the priors
        vector<string> priorsplit, sigmasplit;
        double sigma_def = 0.02;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }
        if(piFileName!=NULL) {
            read_pifile(piFileName, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the total number of possible configurations.");
        vector<double> prior, sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        for(int t=0; t<combNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }

        // 4. define global variables and extract the snp and probe data 
        vector<eqtlInfo> etrait(besdNum); 
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        // extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        // extract the probe list for exposures
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }

        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // read the GWAS cojo signals
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
            if(prbchr!=0) extract_ldata_by_chr(&ldata,prbchr);
        }

        // 5. allele checking between data
        if(!heidioffFlag || jointsmrflag)
        {
            map<string, string> snp_name_per_chr;
            if(bFileName!=NULL) read_famfile(&bdata, string(bFileName)+".fam");
            if(mbFileName!=NULL) read_multi_famfiles(&bdata, multi_bfiles);
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            if(bFileName!=NULL) read_bimfile(&bdata, string(bFileName)+".bim");
            if(mbFileName!=NULL) read_multi_bimfiles(&bdata, multi_bfiles, snp_name_per_chr);
            allele_check_multi_opt(&bdata, etrait, &gdata1);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        } else
        {
            allele_check_multi(etrait,&gdata1);
        }
        // update the SNPs after allele checking
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata1); ngwas = median(gdata1.splSize);
        }
        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }

        // 6. open .smr and .ppa for writing output
        long itemcountsmr=0,itercountmlt=0,itercounttest=0;
        string outstr="";
        // header for .smr
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL) {
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        // header for .ppa
        string smrfile2 = string(outFileName)+".ppa";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
        if (!(smr2)) {
            printf("ERROR: open error %s\n", smrfile2.c_str());
            exit(1);
        }
        outstr="Chr\t";
        int j = 0;
        for(int i=0; i<besdNum;i++)
        {
            j = i+1;
            outstr+="Expo"+atos(j)+"_ID"+'\t'+"Expo"+atos(j)+"_bp"+'\t';
        }
        for(int i=0; i<combmarg.size();i++)
        {
            outstr+="PPA(";
            for(int j=0;j<combmarg[i].size();j++)
            {
                outstr+=atos(combmarg[i][j]);
                if(j<combmarg[i].size()-1) outstr+=",";
            }
            outstr+=")\t";
        }
        for(int i=1; i<combmarg.size();i++)
        {
            outstr+="p_HEIDI(";
            for(int j=0;j<combmarg[i].size();j++)
            {
                outstr+=atos(combmarg[i][j]);
                if(j<combmarg[i].size()-1) outstr+=",";
            }
            if(i<(combmarg.size()-1)) outstr+=")\t";
        }
        outstr+=")\n";
        if(fputs_checked(outstr.c_str(),smr2))
        {
            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
            exit(EXIT_FAILURE);
        }

        // 7. compute the pairwise SMR effect for all exposure probes
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        map<string, double> hdirlts;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                for(int j=0;j<smrrltstmp.size();j++)
                {
                    hdirlts.insert(pair<string, double>(smrrltstmp[j].ProbeID,smrrltstmp[j].p_HET));
                }
            }
        }
        printf("\nPerforming multi-exposure OPERA analysis (including multi-exposure HEIDI tests) ... \n");        

        // 8. loop with GWAS COJO loci; test all possible combinations at each loci
        double cr=0;
        for(int ii=0;ii<ldata._include.size();ii++)
        {
            double desti=1.0*ii/(ldata._include.size());
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            vector<SMRRLT> smrrltsbf;
            vector<long> probNumbf(besdNum,0); 
            vector<long> expoNumbf; //for missing exposures

            int locichr=ldata._chr[ldata._include[ii]];
            int locibp=ldata._bp[ldata._include[ii]];
            int lowerbounder=(locibp-exposure_probe_wind)>0?(locibp-exposure_probe_wind):0;
            int upperbounder=locibp+exposure_probe_wind;
            
            // find all the probes across exposures within the window
            for(int i=0;i<probNum.size();i++)
            {
                int countNum = 0;
                for(int j=0;j<probNum[i];j++)
                {
                    // distance between probe and GWAS loci
                    int bptmp=smrrlts[i][j].Probe_bp;
                    if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                    
                    // select significant probes based on SMR only
                    // if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder && smrrlts[i][j].p_SMR<=thresh_smr) {

                    // select significant probes based on SMR + HEIDI
                    //if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder && smrrlts[i][j].p_SMR<=0.05 && smrrlts[i][j].p_HET>=0.01) {
                    
                    // distance between top-SNP and GWAS loci
                    // int bptmp=smrrlts[i][j].SNP_bp;
                    // if(smrrlts[i][j].SNP_Chr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        // if(smrrlts[i][j].p_SMR<=thresh_smr) {
                        if(smrrlts[i][j].p_SMR<=thresh_smr && smrrlts[i][j].p_HET>=thresh_heidi) {
                            smrrltsbf.push_back(smrrlts[i][j]);
                            countNum = countNum + 1;
                        } 
                        itemcountsmr = itemcountsmr + 1;
                        // output the pairwise SMR results
                        outstr=smrrlts[i][j].ProbeID+'\t'+atos(smrrlts[i][j].ProbeChr)+'\t'+smrrlts[i][j].Gene+'\t'+atos(smrrlts[i][j].Probe_bp)+'\t'+smrrlts[i][j].SNP+'\t'+atos(smrrlts[i][j].SNP_Chr)+'\t'+atos(smrrlts[i][j].SNP_bp)+'\t'+smrrlts[i][j].A1+'\t'+smrrlts[i][j].A2+'\t'+atos(smrrlts[i][j].Freq)+'\t'+atos(smrrlts[i][j].b_GWAS)+'\t'+atos(smrrlts[i][j].se_GWAS)+'\t'+dtos(smrrlts[i][j].p_GWAS)+'\t'+atos(smrrlts[i][j].b_eQTL)+'\t'+atos(smrrlts[i][j].se_eQTL)+'\t'+dtos(smrrlts[i][j].p_eQTL)+'\t'+atos(smrrlts[i][j].b_SMR)+'\t'+atos(smrrlts[i][j].se_SMR)+'\t'+dtos(smrrlts[i][j].p_SMR)+'\t'+(smrrlts[i][j].p_HET >= 0 ? dtos(smrrlts[i][j].p_HET) : "NA") + '\t' + (smrrlts[i][j].nsnp > 0 ? atos(smrrlts[i][j].nsnp+1) : "NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr0))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                probNumbf[i]=countNum;
            }
            // update esdata._esi_include/gdata1._include/bdata._include with only SNPs in the cis-window; 
            lowerbounder=(locibp-cis_itvl*1000)>0?(locibp-cis_itvl*1000):0;
            upperbounder=locibp+cis_itvl*1000;
            for(int i=0;i<expoNum;i++) esdata[i]._esi_include.clear();                 
            gdata1._include.clear(); bdata._include.clear();
            for(int j=0;j<esdata[0]._snpNum;j++) {
               int bptmp=esdata[0]._esi_bp[j];
               if(esdata[0]._esi_chr[j]==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                   for(int i=0;i<expoNum;i++) {
                        esdata[i]._esi_include.push_back(j); 
                   } 
                   gdata1._include.push_back(j);
                   bdata._include.push_back(j);
               }
            }
            // skip the GWAS loci without an exposure probe with significant instrument
            int expocout = 0, missNum = 0; 
            vector<int> missexpoNum(expoNum);
            for(int i=0;i<probNumbf.size();i++) {
                if(probNumbf[i]>0) { expocout+=1; expoNumbf.push_back(i);
                } else { missNum += 1;}
                missexpoNum[i] = missNum;
            }
            // if(expocout == 0) { continue; }
            
            // illustrate all the combinations
            vector<vector<int>> indexall;
            for(int i=0;i<probNumbf.size();i++) {
                if(probNumbf[i] > 0) {
                    vector<int> index(probNumbf[i]);
                    std::iota(index.begin(),index.end(),0);
                    indexall.push_back(index);
                } else {
                    vector<int> index;
                    index.push_back(0);
                    indexall.push_back(index);
                }
            }
            vector<vector<int>> combines;
            permute_vector(indexall, combines);
            printf("\nThere are %ld possible combinations to test in a %ldKb window at %ldbp on chromosome %ld in the joint SMR model\n",combines.size(),op_wind,locibp,locichr);

            // loop with all the possible combinations
            vector<long> idxcomb_smrrltsbf(expoNumbf.size()), idxcomb_smrrltsbf_last(expoNumbf.size());
            vector<eqtlInfo> esdatabf(expoNumbf.size());
            double crcomb=0;
            for(int i=0; i<combines.size(); i++)
            {
                double desticomb=1.0*i/(combines.size());
                if(desticomb>=crcomb)
                {
                    printf("%3.0f%%\r", 100.0*desticomb);
                    fflush(stdout);
                    if(crcomb==0) crcomb+=0.05;
                    else if(crcomb==0.05) crcomb+=0.2;
                    else if(crcomb==0.25) crcomb+=0.5;
                    else crcomb+=0.25;
                }
                vector<vector<string>> prb_cojolist;
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                vector<string> outconamec(besdNum), outcogenec(besdNum); vector<long> outcobpc(besdNum);
                vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum),PIP(combNum);
                MatrixXd lh(2,expoNum);
                // find the probe in smrrltsbf and esdata
                long postmp = 0; idxcomb_smrrltsbf.clear(), idxcomb_smrrltsbf.resize(expoNumbf.size());
                if(operasmrflag) jointsmrflag = false;
                for(int t=0; t<besdNum; t++)
                {   
                    if(probNumbf[t] > 0) {
                        int t_new = t - missexpoNum[t];
                        long idxtmp = combines[i][t]+postmp;
                        idxcomb_smrrltsbf[t_new] = idxtmp;
                        outconamec[t] = smrrltsbf[idxtmp].ProbeID;
                        outcogenec[t] = smrrltsbf[idxtmp].Gene;
                        outcobpc[t] = smrrltsbf[idxtmp].Probe_bp;
                        postmp = postmp + probNumbf[t];

                        if(!heidioffFlag || jointsmrflag) {
                            if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new] || i==0) {
                                esdata[t]._include.clear();
                                map<string, int>::iterator itt;
                                eqtlInfo esdatatmp;
                                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                                if(itt != esdata[t]._probe_name_map.end()) {
                                    esdata[t]._include.push_back(itt->second);
                                    e2econvert(&esdata[t], &esdatatmp);
                                    esdatabf[t_new] = esdatatmp;
                                }
                            }                            
                            if(targetcojosnplstName!=NULL) {
                                // find the target probe COJO signals
                                map<string, vector<string>>::iterator prb_pos;
                                prb_pos = prb_cojosnps.find(outconamec[t]);
                                vector<string> navector; navector.push_back("");
                                if(prb_pos!=prb_cojosnps.end()) {
                                    prb_cojolist.push_back(prb_pos->second);
                                } else { prb_cojolist.push_back(navector); }
                            }
                        }
                    
                    } else {
                        outconamec[t] = "NA"; outcogenec[t] = "NA"; outcobpc[t] = 0;
                    }                                            
                }
                idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;

                if(esdatabf.size()!=0 && esdatabf.size() != expocout) { continue; } 

                // perform the SMR or joint-SMR analysis, including unbalanced exposures
                vector<SMRRLT> smrrlts_joint;
                if(!operasmrflag) { // compute the joint SMR effect
                    if(targetcojosnplstName!=NULL) {
                     multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, prb_cojolist, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    } else { multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor); }
                } else { // compute the marginal SMR effect
                    for(int s=0; s<idxcomb_smrrltsbf.size(); s++) {
                        smrrlts_joint.push_back(smrrltsbf[idxcomb_smrrltsbf[s]]);
                    }
                }
                // skip no joint SMR effect due to no common SNPs
                if(smrrlts_joint.size() != expocout) { continue; }
                // output the probe information;
                outstr=atos(locichr)+'\t';
                for(int t=0; t<besdNum; t++) {
                    outstr+=outconamec[t]+'\t'+atos(outcobpc[t])+'\t';
                }
                //get the bxy, sigma_b and sigma_e from joint-SMR
                int k_joint = 0;
                for(int t=0; t<expoNum; t++)
                {
                    if(probNumbf[t] > 0) {
                        bxy[t] = smrrlts_joint[k_joint].b_SMR;
                        sigma_e[t] = pow(smrrlts_joint[k_joint].se_SMR,2);
                        c[t] = 1+sigma_e[t]/sigma_b[t];
                        k_joint = k_joint + 1;
                    } else {
                        bxy[t] = 0; sigma_e[t] = 0; c[t] = 0;
                    }
                }
                // get the H0 and H1 prior pi and likelihood
                const double PI = 3.141592653589793238463;
                for(int t=0;t<expoNum;t++) {
                    if(probNumbf[t] > 0) {
                        lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                        lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                    } else {
                        lh(0,t) = 1; lh(1,t) = 0;
                    }
                }
                // caculate the posterior probablity
                for(int i=0;i<combNum;i++){
                    HH[i]=1.0;
                    for(int t=0;t<expoNum;t++)
                    {
                        HH[i] *= lh(combins[i][t],t);
                    }
                }
                float POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i] = HH[i]*prior[i];
                    POall+=PO[i];
                }
                for(int i=0;i<combNum;i++) {
                    PP[i] = PO[i]/POall;
                }
                for(int i=0;i<combmarg.size();i++) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        PIP[i] += PP[idxmarg[i][j]];
                    }
                }
                for(int i=0;i<PIP.size();i++) {
                    outstr = outstr + atos(PIP[i])+'\t';
                }
                bool sigflag = false; vector<int> sigcomb;  
                for(int i=1;i<combmarg.size();i++) {
                    if(PIP[i]>=thresh_PP) { sigflag = true; sigcomb.push_back(i); }
                }

                // multi-exposure HEIDI test
                vector<vector<SMRRLT>> smrrltsheidi(combmarg.size());
                if(! heidioffFlag && sigflag) {
                    for(int h=0;h<sigcomb.size();h++) {
                        vector<eqtlInfo> esdataheidi;
                        int tmpidx = sigcomb[h];
                        string prbname;
                        for(int c=0;c<combmarg[tmpidx].size();c++) {
                            int combidx = combmarg[tmpidx][c]-1;
                            for(int k=0;k<expoNumbf.size();k++) {
                                if(expoNumbf[k] == combidx) {
                                    esdataheidi.push_back(esdatabf[k]);
                                    prbname.append(esdatabf[k]._epi_prbID[0]);
                                }
                            }                                                            
                        }
                        map<string, double>::iterator itmp;
                        itmp = hdirlts.find(prbname);
                        if(itmp != hdirlts.end()) {
                            SMRRLT currlt;
                            currlt.p_HET = itmp->second;
                            currlt.ProbeID = prbname;
                            smrrltsheidi[tmpidx].push_back(currlt);
                        } else {
                            multi_heidi_func(smrrltsheidi[tmpidx], NULL, &bdata, &gdata1, esdataheidi, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                            hdirlts.insert(pair<string, double>(prbname,smrrltsheidi[tmpidx][0].p_HET));    
                        }                        
                    }
                    // output heidi pvalue
                    for(int i=1;i<combmarg.size();i++) {
                        if(smrrltsheidi[i].size()>0) {
                            if(i<(combmarg.size()-1)) {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\t' ;
                            } else {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\n';}
                        } else {
                            if(i<(combmarg.size()-1)) {outstr=outstr + "NA" + '\t';
                            } else {outstr=outstr + "NA" + '\n';}
                        }
                    }                
                } else {
                     for(int i=1;i<combmarg.size();i++) {
                         if(i<(combmarg.size()-1)) { outstr=outstr + "NA" + '\t';
                         } else { outstr=outstr + "NA" + '\n'; }
                    }
                }
                if(sigflag) {
                    itercountmlt+=1;
                    if(fputs_checked(outstr.c_str(),smr2))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                itercounttest+=1;
            }
            
        }        
        fclose(smr0);            
        fclose(smr2);
        printf("\nPairwise SMR and HEIDI analyses for %ld exposure probes have been saved in the file %s.\n",itemcountsmr,smrfile0.c_str());
        printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s.\n",itercounttest,expoNum,itercountmlt,smrfile2.c_str());
    }

    void multiexposure_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP, double thresh_smr,double thresh_heidi, char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. check input besd file format
        vector<string> besds, multi_bfiles;;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());
        
        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }

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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // 3. expoNum = besdNum will be used; get prior variance and PIP header
        long expoNum; expoNum = besdNum;
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        }

        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // list marginal PIP combins
        vector<vector<int>> combmarg, comborg(combins.size()), idxmarg(combins.size());
        comborg[0].push_back(0); idxmarg[0].push_back(0);
        for(int i=1;i<combNum;i++) {
            for(int j=0;j<expoNum;j++) {
                if(combins[i][j] == 1) comborg[i].push_back(j+1);
            }
        }
        // reorder the output PIP
        for(int i=1;i<=expoNum;i++) {
            vector<vector<int>> combtmp;
            for(int j=0;j<combNum;j++) {
                if(comborg[j].size()==i) combtmp.push_back(comborg[j]);
            }
            sort(combtmp.begin(),combtmp.end());
            for(int k=0;k<combtmp.size();k++) {
                combmarg.push_back(combtmp[k]);
            }
        }
        // find the PIP correspoding configuration index
        for(int i=1;i<combmarg.size();i++) {
            for(int c=0;c<combNum;c++) {
                int sumcount=0;
                for(int j=0;j<combmarg[i].size();j++) {
                    int tmpidx = combmarg[i][j] - 1;
                    sumcount += combins[c][tmpidx];
                }
                if(sumcount==combmarg[i].size()) idxmarg[i].push_back(c);
            }
        }
        // get the priors (variance and pi)
        vector<string> priorsplit, sigmasplit;
        double sigma_def = 0.02;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }
        if(piFileName!=NULL) {
            read_pifile(piFileName, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }        
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of test configurations.");
        vector<double> prior, sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        for(int t=0; t<combNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }

        // 4. define global variables and extract the snp and probe data         
        vector<eqtlInfo> etrait(besdNum); 
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        // extract the SNP list for exposures
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        // extract the probe list for exposures
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }
        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }

        // 5. allele checking between data
        if(!heidioffFlag || jointsmrflag)
        {
            map<string, string> snp_name_per_chr;
            if(bFileName!=NULL) read_famfile(&bdata, string(bFileName)+".fam");
            if(mbFileName!=NULL) read_multi_famfiles(&bdata, multi_bfiles);
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            if(bFileName!=NULL) read_bimfile(&bdata, string(bFileName)+".bim");
            if(mbFileName!=NULL) read_multi_bimfiles(&bdata, multi_bfiles, snp_name_per_chr);
            allele_check_multi_opt(&bdata, etrait, &gdata1);
            //            allele_check_multi(&bdata, etrait, &gdata1);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        } else
        {
            allele_check_multi(etrait,&gdata1);
        }
        // update the SNPs after allele checking
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata1); ngwas = median(gdata1.splSize);
        }

        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }

        // 6. open .smr and .ppa for writing output
        long itemcount=0,itercountmlt=0,itercounttest=0;
        string outstr="";
        // header for .smr
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL) {
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        // header for .ppa
        string smrfile2 = string(outFileName)+".ppa";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
        if (!(smr2)) {
            printf("ERROR: open error %s\n", smrfile2.c_str());
            exit(1);
        }
        outstr="Chr\t";
        int j = 0;
        for(int i=0; i<besdNum;i++)
        {
            j = i+1;
            outstr+="Expo"+atos(j)+"_ID"+'\t'+"Expo"+atos(j)+"_bp"+'\t';
        }
        for(int i=0; i<combmarg.size();i++)
        {
            outstr+="PPA(";
            for(int j=0;j<combmarg[i].size();j++)
            {
                outstr+=atos(combmarg[i][j]);
                if(j<combmarg[i].size()-1) outstr+=",";
            }
            outstr+=")\t";
        }
        for(int i=1; i<combmarg.size();i++)
        {
            outstr+="p_HEIDI(";
            for(int j=0;j<combmarg[i].size();j++)
            {
                outstr+=atos(combmarg[i][j]);
                if(j<combmarg[i].size()-1) outstr+=",";
            }
            if(i<(combmarg.size()-1)) outstr+=")\t";
        }
        outstr+=")\n";
        if(fputs_checked(outstr.c_str(),smr2))
        {
            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
            exit(EXIT_FAILURE);
        }

        // 7. compute the pairwise SMR effect for all exposure probes
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        map<string, double> hdirlts;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                for(int j=0;j<smrrltstmp.size();j++)
                {
                    outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr0))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                        exit(EXIT_FAILURE);
                    }
                    hdirlts.insert(pair<string, double>(smrrltstmp[j].ProbeID,smrrltstmp[j].p_HET));
                }                
            } 
        }
        for(int k=0;k<probNum.size();k++)
        {   
            itemcount = itemcount + probNum[k]; 
        }
        printf("SMR analysis results of %ld exposure probes have been saved in the file %s .\n",itemcount,smrfile0.c_str());
        fclose(smr0);

        printf("\nPerforming multi-exposure OPERA analysis (including multi-exposure HEIDI tests) ... \n");
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposure probes with significant instrument are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        }

        // 8. loop with each SMR < 0.05 loci and test all possible combinations at each probe loci;
        int process = 0;
        map<string, long> combname_set;
        for(int tt=0;tt<expoNum;tt++) {            
            double cr=0;
            
            for(int ii=0;ii<probNum[tt];ii++) {
                process = process + 1;
                double desti=1.0*process/itemcount;
                if(desti>=cr)
                {
                    printf("%3.0f%%\r", 100.0*desti);
                    fflush(stdout);
                    if(cr==0) cr+=0.05;
                    else if(cr==0.05) cr+=0.2;
                    else if(cr==0.25) cr+=0.5;
                    else cr+=0.25;
                }
                vector<SMRRLT> smrrltsbf;
                vector<long> probNumbf(expoNum);
                vector<long> expoNumbf; //for missing exposures
                // skip exposures P_SMR > 0.05
                // if(smrrlts[tt][ii].p_SMR > thresh_smr) { continue; }
                if(smrrlts[tt][ii].p_SMR > thresh_smr || smrrlts[tt][ii].p_HET < thresh_heidi) { continue; }
                
                int traitchr=smrrlts[tt][ii].ProbeChr;
                int traitbp=smrrlts[tt][ii].Probe_bp;
                int lowerbounder=(traitbp-exposure_probe_wind)>0?(traitbp-exposure_probe_wind):0;
                int upperbounder=traitbp+exposure_probe_wind;
                // find other exposure probes in the cis-window of target exposure
                for(int i=0;i<expoNum;i++)
                {
                    int countNum = 0;
                    if(i != tt) {
                        for(int j=0;j<probNum[i];j++)
                        {
                           int bptmp=smrrlts[i][j].Probe_bp;
                           // select probes with SMR pvalue < 0.05 & HEIDI pvalue > 1e-5 for OPERA analysis
                           // if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder && smrrlts[i][j].p_SMR<=thresh_smr) {
                           if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder && smrrlts[i][j].p_SMR<=thresh_smr && smrrlts[i][j].p_HET>=thresh_heidi) {
                                smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                           }
                        }                        
                    } 
                    if(i == tt) {                        
                        smrrltsbf.push_back(smrrlts[tt][ii]); countNum = countNum + 1;
                    }
                    probNumbf[i] = countNum;                                                
                }
                // update esdata._esi_include/gdata1._include/bdata._include with only SNPs in the cis-window; 
                // use the exposure tt bp as gold standard; gdata1, bdata and esdata included SNPs are the same;
                lowerbounder=(traitbp-cis_itvl*1000)>0?(traitbp-cis_itvl*1000):0;
                upperbounder=traitbp+cis_itvl*1000;
                for(int i=0;i<expoNum;i++) esdata[i]._esi_include.clear();                 
                gdata1._include.clear(); bdata._include.clear();
                for(int j=0;j<esdata[tt]._snpNum;j++) {
                   int bptmp=esdata[tt]._esi_bp[j];
                   if(esdata[tt]._esi_chr[j]==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                       for(int i=0;i<expoNum;i++) {
                            esdata[i]._esi_include.push_back(j); 
                       } 
                       gdata1._include.push_back(j);
                       bdata._include.push_back(j);
                   }
                }
                // skip the exposure probe where no nearby probes from other exposure with SMR effect
                int expocout = 0, missNum = 0;
                vector<int> missexpoNum(expoNum);
                for(int i=0;i<expoNum;i++) {
                    if(probNumbf[i]>0) { expocout+=1; expoNumbf.push_back(i);
                    } else { missNum += 1;}
                    missexpoNum[i] = missNum;
                }
                // if(expocout == 0) { continue; } // not possible

                // illustrate all the combinations
                vector<vector<int>> indexall;
                for(int i=0;i<probNumbf.size();i++) {
                    if(probNumbf[i] > 0) {
                        vector<int> index(probNumbf[i]);
                        std::iota(index.begin(),index.end(),0);
                        indexall.push_back(index);
                    } else {
                        vector<int> index;
                        index.push_back(0);
                        indexall.push_back(index);
                    }
                }
                vector<vector<int>> combines;
                permute_vector(indexall, combines);

                 // loop with all the possible combinations
                vector<long> idxcomb_smrrltsbf(expoNumbf.size()), idxcomb_smrrltsbf_last(expoNumbf.size());
                vector<eqtlInfo> esdatabf(expoNumbf.size());
                for(int cc=0; cc<combines.size(); cc++)
                {
                    vector<vector<string>> prb_cojolist;
                    vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);                    
                    vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum),PIP(combNum);
                    vector<string> outconamec(besdNum), outcogenec(besdNum); vector<long> outcobpc(besdNum);
                    MatrixXd lh(2,expoNum);                    
                    // get the target combination name
                    long postmp = 0; string combname;
                    for(int t=0; t<besdNum; t++)
                    {   
                        if(probNumbf[t] > 0) {
                            long idxtmp = combines[cc][t] + postmp;
                            outconamec[t] = smrrltsbf[idxtmp].ProbeID;
                            outcogenec[t] = smrrltsbf[idxtmp].Gene;
                            outcobpc[t] = smrrltsbf[idxtmp].Probe_bp;
                            postmp = postmp + probNumbf[t];
                        } else {
                            outconamec[t] = "NA"; outcogenec[t] = "NA"; outcobpc[t] = 0;
                        }
                        combname.append(outconamec[t]);                        
                    }
                    // skip the tested combination
                    map<string, long>::iterator comb_pos;
                    comb_pos = combname_set.find(combname);
                    if(comb_pos != combname_set.end()) { continue; } 

                    combname_set.insert(pair<string, long>(combname, itercounttest)); 
                    // find the probe in smrrltsbf and esdata                    
                    postmp = 0; idxcomb_smrrltsbf.clear(); idxcomb_smrrltsbf.resize(expoNumbf.size());
                    if(operasmrflag) jointsmrflag = false;
                    for(int t=0; t<besdNum; t++)
                    {   
                        if(probNumbf[t] > 0) {
                            int t_new = t - missexpoNum[t];
                            long idxtmp = combines[cc][t] + postmp;
                            idxcomb_smrrltsbf[t_new] = idxtmp;
                            postmp = postmp + probNumbf[t];
                            if(!heidioffFlag || jointsmrflag) {
                                if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new] || cc==0) {
                                    esdata[t]._include.clear();
                                    map<string, int>::iterator itt;
                                    eqtlInfo esdatatmp;
                                    itt = esdata[t]._probe_name_map.find(outconamec[t]);
                                    if(itt != esdata[t]._probe_name_map.end()) {
                                        esdata[t]._include.push_back(itt->second);
                                        e2econvert(&esdata[t], &esdatatmp);                                
                                        esdatabf[t_new] = esdatatmp;
                                    }
                                }                            
                                if(targetcojosnplstName!=NULL) {
                                    // find the target probe COJO signals
                                    map<string, vector<string>>::iterator prb_pos;
                                    prb_pos = prb_cojosnps.find(outconamec[t]);
                                    vector<string> navector; navector.push_back("");
                                    if(prb_pos!=prb_cojosnps.end()) {
                                        prb_cojolist.push_back(prb_pos->second);
                                    } else { prb_cojolist.push_back(navector); }
                                }
                            }                        
                        }                                                
                    }
                    idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;
                    
                    if(esdatabf.size()!=0 && esdatabf.size() != expocout) { continue; } 

                    // perform joint-SMR analysis or extract the SMR effect                
                    vector<SMRRLT> smrrlts_joint;
                    if(!operasmrflag) { // compute the joint SMR effect
                        if(targetcojosnplstName!=NULL) {
                            multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, prb_cojolist, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                        } else { multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor); }
                    } else { // compute the marginal SMR effect
                        for(int s=0; s<idxcomb_smrrltsbf.size(); s++) {
                            smrrlts_joint.push_back(smrrltsbf[idxcomb_smrrltsbf[s]]);
                        }
                    }
                    // skip no joint SMR effect due to no common SNPs
                    if(smrrlts_joint.size() != expocout) { continue; }
                    // output the probe information;
                    outstr=atos(traitchr)+'\t';
                    for(int t=0; t<besdNum; t++) {
                        outstr+=outconamec[t]+'\t'+atos(outcobpc[t])+'\t';
                    }
                    //get the bxy, sigma_b and sigma_e from joint-SMR
                    int k_joint = 0;
                    for(int t=0; t<expoNum; t++)
                    {
                        if(probNumbf[t] > 0) {
                            bxy[t] = smrrlts_joint[k_joint].b_SMR;
                            sigma_e[t] = pow(smrrlts_joint[k_joint].se_SMR,2);
                            c[t] = 1+sigma_e[t]/sigma_b[t];
                            k_joint = k_joint + 1;
                        } else {
                            bxy[t] = 0; sigma_e[t] = 0; c[t] = 0;
                        }
                    }
                    // get the H0 and H1 prior pi and likelihood
                    const double PI = 3.141592653589793238463;
                    for(int t=0;t<expoNum;t++) {
                        if(probNumbf[t] > 0) {
                            lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                            lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                        } else {
                            lh(0,t) = 1; lh(1,t) = 0;
                        }
                    }
                    // caculate the posterier probablity
                    for(int i=0;i<combNum;i++) {
                        HH[i]=1.0;
                        for(int t=0;t<expoNum;t++)
                        {
                            HH[i] *= lh(combins[i][t],t);
                        }
                    }
                    float POall = 0;
                    for(int i=0;i<combNum;i++) {
                        PO[i] = HH[i]*prior[i];
                        POall+=PO[i];
                    }
                    for(int i=0;i<combNum;i++) {
                        PP[i] = PO[i]/POall;
                    }
                    for(int i=0;i<combmarg.size();i++) {
                        for(int j=0;j<idxmarg[i].size();j++) {
                            PIP[i] += PP[idxmarg[i][j]];
                        }
                    }
                    for(int i=0;i<PIP.size();i++) {
                        outstr = outstr + atos(PIP[i])+'\t';
                    }
                    bool sigflag = false; vector<int> sigcomb;  
                    for(int i=1;i<combmarg.size();i++) {
                        if(PIP[i]>=thresh_PP) { sigflag = true; sigcomb.push_back(i); }
                    }
                    
                    // multi-exposure HEIDI test
                    vector<vector<SMRRLT>> smrrltsheidi(combmarg.size());
                    if(! heidioffFlag && sigflag) {
                        for(int h=0;h<sigcomb.size();h++) {
                            vector<eqtlInfo> esdataheidi;
                            int tmpidx = sigcomb[h];
                            string prbname;
                            for(int c=0;c<combmarg[tmpidx].size();c++) {
                                int combidx = combmarg[tmpidx][c]-1;
                                for(int k=0;k<expoNumbf.size();k++) {
                                    if(expoNumbf[k] == combidx) {
                                        esdataheidi.push_back(esdatabf[k]);
                                        prbname.append(esdatabf[k]._epi_prbID[0]);
                                    }
                                }                                                            
                            }
                            map<string, double>::iterator itmp;
                            itmp = hdirlts.find(prbname);
                            if(itmp != hdirlts.end()) {
                                SMRRLT currlt;
                                currlt.p_HET = itmp->second;
                                currlt.ProbeID = prbname;
                                smrrltsheidi[tmpidx].push_back(currlt);
                            } else {
                                multi_heidi_func(smrrltsheidi[tmpidx], NULL, &bdata, &gdata1, esdataheidi, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                                hdirlts.insert(pair<string, double>(prbname,smrrltsheidi[tmpidx][0].p_HET));
                            }
                        }
                        // output heidi pvalue
                        for(int i=1;i<combmarg.size();i++) {
                            if(smrrltsheidi[i].size()>0) {
                                if(i<(combmarg.size()-1)) {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\t' ;
                                } else {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\n';}
                            } else {
                                if(i<(combmarg.size()-1)) {outstr=outstr + "NA" + '\t';
                                } else {outstr=outstr + "NA" + '\n';}
                            }
                        }                
                    } else {
                         for(int i=1;i<combmarg.size();i++) {
                             if(i<(combmarg.size()-1)) { outstr=outstr + "NA" + '\t';
                             } else { outstr=outstr + "NA" + '\n'; }
                        }
                    }
                    if(sigflag) {
                        itercountmlt+=1;
                        if(fputs_checked(outstr.c_str(),smr2))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }   
                    itercounttest+=1;

                }
                
            }
        }
        
        fclose(smr2);
        printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s.\n",itercounttest,expoNum,itercountmlt,smrfile2.c_str());
    }

    void multioutcomesmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* eqtlFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        //here eqtlFileName is the outcome and eqtlFileName2 is the exposure. in the main we pass the outcome (eqtlFileName2) to eqtlFileName and the exposure (eqtlFileName) to eqtlFileName2
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        // checking besd file list
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"eQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
        //printf("Checking the BESD format...\n");
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        long outcoNum;
        if(gwasFileName!=NULL)
        {
            outcoNum = besdNum +1;
        } else {outcoNum = besdNum;}
        printf("There are %ld outcomes and 1 exposure are included in multiple-SMR and multiple-HEIDI test.\n",outcoNum);
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<outcoNum;i++){
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        
        // get the priors
        vector<string> priorsplit, sigmasplit;
        split_string(priorstr,priorsplit);
        split_string(sigmastr,sigmasplit);
        if(priorsplit.size()!=sigmasplit.size() || priorsplit.size()!=outcoNum || sigmasplit.size()!=outcoNum)
            throw("Error: The number of input prior is not consistent with the number of input outcomes.");
        vector<double> prior, sigma_b;
        for(int t=0; t<outcoNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }
        
        bool targetLstflg=false;
        map<string, string> prb_snp;
        vector<eqtlInfo> etrait(besdNum);        //vector<eqtlInfo*> etrait;
        eqtlInfo esdata;
        bInfo bdata;
        gwasData gdata1;
        bool heidiFlag=false;
        
        printf("Reading the exposure summary data file ...\n");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        printf("Reading the outcome summary data file ...\n");
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }

        //the etrait is not updated, so from now on _esi_include should be used always.
        for(int i=0;i<besdNum;i++) {
            //cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(oproblstName != NULL ) extract_prob(&etrait[i], oproblstName);
            else if(oprobe != NULL) extract_eqtl_single_probe(&etrait[i], oprobe);
            if(oproblst2exclde != NULL) exclude_prob(&etrait[i], oproblst2exclde);
            else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], oprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            //cout<<etrait[i]._rowid<<endl;
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }
        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // allele checking between data
        if(!heidioffFlag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            if(gwasFileName!=NULL){
                allele_check_multi(&bdata, etrait, &gdata1, &esdata);
            }else { allele_check_multi(&bdata, etrait, &esdata);}
            // if no snp left after check
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                if(gwasFileName!=NULL)
                {
                    update_geIndx(&bdata, etrait, &gdata1, &esdata);
                } else {update_geIndx(&bdata, etrait, &esdata);}
            }
        }else
        {
            if(gwasFileName!=NULL)
            {
                allele_check_multi(etrait,&gdata1,&esdata);
            } else {allele_check_multi(etrait,&esdata);}
            
        }
        if(gwasFileName!=NULL)  update_gwas(&gdata1);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included in the analysis.\n");
            exit(EXIT_FAILURE);
        }
        // open multiple outcome files to write
        long itemcount=0,itercountmlt=0,itercounttest=0;
        long etraitcount=0;
        string outstr="";
        // GWAS trait
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL){
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        // molecular trait
        string smrfile1 = string(outFileName)+".msmr";
        FILE* smr1 = fopen(smrfile1.c_str(), "w");
            if (!(smr1)) {
                printf("ERROR: open error %s\n", smrfile1.c_str());
                exit(1);
            }
            outstr="Expo_ID\tExpo_Chr\tExpo_Gene\tExpo_bp\tOutco_ID\tOutco_Chr\tOutco_Gene\tOutco_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_Outco\tse_Outco\tp_Outco\tb_Expo\tse_Expo\tp_Expo\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr1))
            {
                printf("ERROR: in writing file %s .\n", smrfile1.c_str());
                exit(EXIT_FAILURE);
            }
        // posterior probability
        string smrfile2 = string(outFileName)+".multismr";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
            if (!(smr2)) {
                printf("ERROR: open error %s\n", smrfile2.c_str());
                exit(1);
            }
            //outstr="Chr\tExpo_ID\tExpo_Gene\tExpo_bp\tOutco1_ID\tOutco1_Gene\tOutco1_bp\tOutco2_ID\tOutco2_Gene\tOutco2_bp\tOutco3_ID\tPP1(0:0:0)\tPP2(0:0:1)\tPP3(0:1:0)\tPP4(0:1:1)\tPP5(1:0:0)\tPP6(1:0:1)\tPP7(1:1:0)\tPP8(1:1:1)\n";
            outstr="Chr\tExpo_ID\tExpo_Gene\tExpo_bp\t";
            int j = 0;
            for(int i=0; i<besdNum;i++)
            {
                j = i+1;
                outstr+="Outco"+atos(j)+"_ID"+'\t'+"Outco"+atos(j)+"_Gene"+'\t'+"Outco"+atos(j)+"_bp"+'\t';
            }
            if(gwasFileName!=NULL) outstr+="Outco"+atos(j+1)+"_ID"+'\t';
        
            for(int i=0; i<combins.size();i++)
            {
                outstr+="PP"+atos(i+1)+"(";
                for(int j=0;j<outcoNum;j++)
                {
                    outstr+=atos(combins[i][j]);
                    if(j<outcoNum-1) outstr+=":";
                }
                if(i<(combins.size()-1)) outstr+=")\t";
            }
            outstr+=")\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr2))
            {
                printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                exit(EXIT_FAILURE);
            }   
        // end of opening files
        int exposure_probe_wind=op_wind*1000;
        vector<vector<int>> includebk(besdNum);
        for(int t=0;t<besdNum;t++)
        {
            includebk[t]=etrait[t]._include;
        }
        // loop with exposure probes
        printf("\nPerforming multi-outcome-SMR analysis (multi-SMR and multi-HEIDI tests) for %d exposure probes... \n",esdata._probNum);
        double cr=0;
        for( int ii=0;ii<esdata._probNum;ii++)
        {   
            double desti=1.0*ii/(esdata._probNum-1);
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //gwasData gdata;
            esdata._include.clear();
            esdata._include.push_back(ii);
            string traitname=esdata._epi_prbID[ii];
            int traitchr=esdata._epi_chr[ii];
            string traitgene=esdata._epi_gene[ii];
            int traitbp=esdata._epi_bp[ii];
            //cout<<"\nPerforming analysis of exposure [ "+traitname+" ]..."<<endl;
    
           if(cis2all)
           {
               printf("\n%ld outcome probes are inclued in the analysis with the exposure probe %s.\n", esdata._include.size(),traitname.c_str());
               
           } else
           {
               int lowerbounder=(traitbp-exposure_probe_wind)>0?(traitbp-exposure_probe_wind):0;
               int upperbounder=traitbp+exposure_probe_wind;

               for(int i=0;i<besdNum;i++)
               {
                   etrait[i]._include.clear();
                   for(int j=0;j<includebk[i].size();j++)
                   {
                       int bptmp=etrait[i]._epi_bp[includebk[i][j]];
                       if(etrait[i]._epi_chr[includebk[i][j]]==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) etrait[i]._include.push_back(includebk[i][j]);
                   }
                   //printf("\n%ld probes from outcome%d in the cis region [%d, %d] of exposure probe %s are inclued in the analysis.\n", etrait[i]._include.size(),i+1,lowerbounder,upperbounder,traitname.c_str());
                }
            }
            vector<long> probNum;
            vector<vector<SMRRLT>> smrrlts;
            vector<gwasData> gdatall;
            vector<string> outconame;
            vector<int> outcochr;
            vector<string> outcogene;
            vector<int> outcobp;
            // molecular traits SMR
            for(int i=0;i<besdNum;i++)
            {
                long count=0, cisNum = etrait[i]._include.size();
                for(int p=0;p<cisNum;p++)
                {
                gwasData gdata;
                vector<SMRRLT> smrrlts1;
                int j = etrait[i]._include[p];
                e2gconvert(&etrait[i], &gdata, j);
                if(gdata.snpName.size()< m_hetero)
                {
                    //cout<<"\n"<<gdata.snpNum<<" common SNPs (less than parameter m_hetero: "+atos(m_hetero)+" ) are included from eTrait [ "+etrait[i]._epi_prbID[j]+" ] summary."<<endl;
                    continue;
                }
                gdata.snpNum=gdata._include.size();
                //cout<<"\n"<<gdata.snpNum<<" common SNPs are included from eTrait [ "+etrait[i]._epi_prbID[j]+" ] summary."<<endl;
                smr_heidi_func(smrrlts1, NULL, &bdata, &gdata, &esdata, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);                
                if(smrrlts1.size()>0)
                    {
                        count+=1;
                        itemcount+=1;
                        if(!heidioffFlag) {gdatall.push_back(gdata);}
                        smrrlts.push_back(smrrlts1);
                        outcochr.push_back(etrait[i]._epi_chr[j]);
                        outconame.push_back(etrait[i]._epi_prbID[j]);
                        outcogene.push_back(etrait[i]._epi_gene[j]);
                        outcobp.push_back(etrait[i]._epi_bp[j]);
                    }
                }
                probNum.push_back(count);
            }
            // GWAS SMR
            vector<SMRRLT> smrrltsgwas;
            if(gwasFileName!=NULL)
            {
            smr_heidi_func(smrrltsgwas, NULL, &bdata,&gdata1,&esdata, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            if(smrrltsgwas.size()>0)
                {
                    
                    probNum.push_back(1);
                    smrrlts.push_back(smrrltsgwas);
                    if(!heidioffFlag) {gdatall.push_back(gdata1);}
                } else {
                    //printf("WARNING: No common SNP passed the p-value threshold %e for the SMR analysis of complex trait.\n", p_smr);
                    continue;
                }
            }
            if(smrrlts.size()>0) etraitcount+=1;

            int outcomcout=0; 
            for(int i=0;i<probNum.size();i++){
                if(probNum[i]>0) outcomcout+=1;
            }

            if(outcomcout==outcoNum) {
                // illustrate all the combinations
                vector<vector<int>> indexall;
                for(int i=0;i<probNum.size();i++){
                    vector<int> index(probNum[i]);
                    std::iota(index.begin(),index.end(),0);
                    indexall.push_back(index);
                }
                vector<vector<int>> combines;
                permute_vector(indexall, combines);
                // Bayesian factor caculation
                //printf("\nPerforming Multi-SMR tests.....\n");
                //if (combines.size()==0) printf("\nWARNING: Less than %ld outcomes included in the analysis. \nThe multi-SMR and multi-HEIDI test will skip for exposure probe %s\n",outcoNum,traitname.c_str());
                long combNum=round(pow(2,outcoNum));
                for(int i=0; i<combines.size();i++)
                {
                    vector<float> bxy(outcoNum), sigma_e(outcoNum),c(outcoNum);
                    vector<string> outconamec(besdNum), outcogenec(besdNum); vector<long> outcobpc(besdNum);
                    vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum);
                    MatrixXd pie(2,outcoNum), lh(2,outcoNum);

                    vector<gwasData> gdatain(outcoNum);

                    // get the probe and gene information for output
                    outstr=atos(traitchr)+'\t'+traitname+'\t'+traitgene+'\t'+atos(traitbp)+'\t';
                    long postmp=0;
                    for(int t=0; t<besdNum; t++)
                    {
                        long idxtmp = combines[i][t]+postmp;
                        outconamec[t] = outconame[idxtmp];
                        outcogenec[t] = outcogene[idxtmp];
                        outcobpc[t] = outcobp[idxtmp];
                        postmp=postmp+probNum[t];
                        outstr+=outconamec[t]+'\t'+outcogenec[t]+'\t'+atos(outcobpc[t])+'\t';
                    }
                    if(gwasFileName!=NULL) outstr+="Trait\t";
                    //get the bxy, sigma_b and sigma_e from SMR
                    long pos=0;
                    for(int t=0; t<outcoNum; t++)
                    {
                        long idx = combines[i][t]+pos;
                        bxy[t]=smrrlts[idx][0].b_SMR;
                        sigma_e[t]=pow(smrrlts[idx][0].se_SMR,2);
                        c[t]=1+sigma_e[t]/sigma_b[t];
                        pos=pos+probNum[t];
                        if(!heidioffFlag) {gdatain[t]=gdatall[idx];}

                    }
                    // multi-HEIDI test
                    vector<SMRRLT> smrrltsheidi;
                    if(!heidioffFlag){
                    multi_heidi_func(smrrltsheidi, NULL, &bdata, gdatain, &esdata, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    }
                    // get the H0 and H1 prior pi and likelihood
                    const double PI = 3.141592653589793238463;
                    for(int t=0;t<outcoNum;t++) {
                        pie(0,t)=1-prior[t]; pie(1,t)=prior[t];
                        lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                        lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                    }
                    // caculate the posterier probablity
                    for(int i=0;i<combNum;i++){
                        Pr[i]=1.0;  HH[i]=1.0;
                        for(int t=0;t<outcoNum;t++)
                        {
                            Pr[i] *= pie(combins[i][t],t);
                            HH[i] *= lh(combins[i][t],t);
                        }
                    }
                    float POall=0;
                    for(int i=0;i<combNum;i++){
                        PO[i]=HH[i]*Pr[i];
                        POall+=PO[i];
                    }
                    for(int i=0;i<combNum;i++){
                        PP[i]=PO[i]/POall;
                        outstr=outstr+atos(PP[i])+'\t';
                    }
                    if(!heidioffFlag){
                    outstr=outstr+(smrrltsheidi[0].p_HET >= 0 ? dtos(smrrltsheidi[0].p_HET) : "NA") + '\t' + (smrrltsheidi[0].nsnp > 0 ? atos(smrrltsheidi[0].nsnp+1) : "NA") + '\n';
                    } else {
                        outstr=outstr + "NA" + '\t' + "NA" + '\n';
                    }
                    // if(PP[combNum-1]>=0.6){
                    if(1-PP[0]>=0.5){
                        itercountmlt+=1;
                        if(fputs_checked(outstr.c_str(),smr2))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                    itercounttest+=1;                    
                }

            } else {
                continue;
            }
            
            // output the GWAS SMR results if included
            if(gwasFileName!=NULL){
                outstr=smrrltsgwas[0].ProbeID+'\t'+atos(smrrltsgwas[0].ProbeChr)+'\t'+smrrltsgwas[0].Gene+'\t'+atos(smrrltsgwas[0].Probe_bp)+'\t'+smrrltsgwas[0].SNP+'\t'+atos(smrrltsgwas[0].SNP_Chr)+'\t'+atos(smrrltsgwas[0].SNP_bp)+'\t'+smrrltsgwas[0].A1+'\t'+smrrltsgwas[0].A2+'\t'+atos(smrrltsgwas[0].Freq)+'\t'+atos(smrrltsgwas[0].b_GWAS)+'\t'+atos(smrrltsgwas[0].se_GWAS)+'\t'+dtos(smrrltsgwas[0].p_GWAS)+'\t'+atos(smrrltsgwas[0].b_eQTL)+'\t'+atos(smrrltsgwas[0].se_eQTL)+'\t'+dtos(smrrltsgwas[0].p_eQTL)+'\t'+atos(smrrltsgwas[0].b_SMR)+'\t'+atos(smrrltsgwas[0].se_SMR)+'\t'+dtos(smrrltsgwas[0].p_SMR)+'\t'+(smrrltsgwas[0].p_HET >= 0 ? dtos(smrrltsgwas[0].p_HET) : "NA") + '\t' + (smrrltsgwas[0].nsnp > 0 ? atos(smrrltsgwas[0].nsnp+1) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr0))
                {
                    printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                    exit(EXIT_FAILURE);
                }
                
            }
            // output molecular trait SMR results
            long mNum=smrrlts.size();
            if (gwasFileName!=NULL) mNum=mNum-1;
            for(int j=0;j<mNum;j++)
            {
                outstr=smrrlts[j][0].ProbeID+'\t'+atos(smrrlts[j][0].ProbeChr)+'\t'+smrrlts[j][0].Gene+'\t'+atos(smrrlts[j][0].Probe_bp)+'\t'+outconame[j]+'\t'+atos(outcochr[j])+'\t'+outcogene[j]+'\t'+atos(outcobp[j])+'\t'+smrrlts[j][0].SNP+'\t'+atos(smrrlts[j][0].SNP_Chr)+'\t'+atos(smrrlts[j][0].SNP_bp)+'\t'+smrrlts[j][0].A1+'\t'+smrrlts[j][0].A2+'\t'+atos(smrrlts[j][0].Freq)+'\t'+atos(smrrlts[j][0].b_GWAS)+'\t'+atos(smrrlts[j][0].se_GWAS)+'\t'+dtos(smrrlts[j][0].p_GWAS)+'\t'+atos(smrrlts[j][0].b_eQTL)+'\t'+atos(smrrlts[j][0].se_eQTL)+'\t'+dtos(smrrlts[j][0].p_eQTL)+'\t'+atos(smrrlts[j][0].b_SMR)+'\t'+atos(smrrlts[j][0].se_SMR)+'\t'+dtos(smrrlts[j][0].p_SMR)+'\t'+(smrrlts[j][0].p_HET >= 0 ? dtos(smrrlts[j][0].p_HET) : "NA") + '\t' + (smrrlts[j][0].nsnp > 0 ? atos(smrrlts[j][0].nsnp+1) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr1))
                {
                    printf("ERROR: in writing file %s .\n", smrfile1.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }
        if(gwasFileName!=NULL)
        {
            fclose(smr0);
            printf("\n\nSMR and HEIDI analyses for complex traits completed.\nSMR and heterogeneity analysis results of %ld exposure probes have been saved in the file %s.\n",etraitcount,smrfile0.c_str());
        }
        fclose(smr1);
        fclose(smr2);
        printf("\nSMR and HEIDI analyses for molecular traits completed.\nSMR and heterogeneity analysis results of %ld outcome probes (%ld exposure probe) have been saved in the file %s.\n",itemcount,etraitcount,smrfile1.c_str());
        printf("\nMulti-SMR and Multi-HEIDI analyses for %ld between 1 exposure and %ld outcomes completed.\nPosterior probability and HEIDI results of %ld combinations (%ld outcome probes and %ld exposure probes) have been saved in the file %s.\n",itercounttest,outcoNum,itercountmlt,itemcount,etraitcount,smrfile2.c_str());
    }

    // jointSMR balanced, no updated anymore
    void multiexposure_jointsmr_balanced(char* outFileName, char* bFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag, bool jointsmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
    {   
        // 1. check flags; eqtlsmaslstName is the included exposure probes and gwasFileName will be the outcome 
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(ld_min>ld_top) {
            printf("ERROR: --ld-min %f is larger than --ld-top %f.\n",ld_min,ld_top);
            exit(EXIT_FAILURE);
        }
        if(refSNP && targetcojosnplstName)
        {
            printf("ERROR: --target-snp and --target-snp-probe-list both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }
        if(snpproblstName && targetcojosnplstName)
        {
            printf("ERROR: --extract-target-snp-probe and --extract-snp-probe both found in your command. please disable one.\n");
            exit(EXIT_FAILURE);
        }

        // 2. check input besd file format
        vector<string> besds;
        vector<smr_snpinfo> snpinfo;
        vector<smr_probeinfo> probeinfo;
        vector<uint64_t> nprb,nsnp;
        vector< vector<int> > lookup;
        
        read_msglist(eqtlsmaslstName, besds,"eQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld eQTL summary file names are included.\n",besds.size());
        
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
        long besdNum=besds.size();
        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besdNum;i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening file %s .\n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
        // this expoNum is the total number of exposures 
        long expoNum;
        expoNum = besdNum;
        printf("There are %ld exposures and 1 outcome are included in OPERA.\n",expoNum);
        if(expoNum<2) throw("Error: The program can not perform the joint SMR analysis because there is only one exposure included. Please remove the flag --joint-smr to analyze single exposure.");
        // illustrate all the combinations
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++){
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        vector<vector<int>> combins;
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // get the priors
        vector<string> priorsplit, sigmasplit;
        split_string(priorstr,priorsplit);
        split_string(sigmastr,sigmasplit);
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of test configurations.");
        vector<double> prior, sigma_b;
        for(int t=0; t<expoNum; t++)
        {
            sigma_b.push_back(atof(sigmasplit[t].c_str()));
            if(sigma_b[t]<0 || sigma_b[t]>1) throw("Error: --prior-sigma. Prior variance values should be betweeen 0 and 1.");
        }
        for(int t=0; t<combNum; t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }
        bool targetLstflg=false;
        map<string, string> prb_snp;
        vector<eqtlInfo> etrait(besdNum); 
        vector<eqtlInfo> esdata(besdNum);
        bInfo bdata;
        gwasData gdata1;
        bool heidiFlag=false;
        
        printf("Reading the exposure summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL) || (jointsmrflag && bFileName == NULL)) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        //if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        printf("Reading the exposure summary data file ...\n");
        for(int i=0;i<besdNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }

        //the etrait is not updated, so from now on _esi_include should be used always.
        for(int i=0;i<besdNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");
            if(problstName != NULL) extract_prob(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);
            read_besdfile(&etrait[i], string(besds[i])+".besd");
            if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
            {
                printf("ERROR: no data included in the analysis.\n");
                exit(EXIT_FAILURE);
            }
        }
        // Yang Wu read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata1, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata1, snplstName);
                update_gwas(&gdata1);
            } 
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // read the GWAS cojo signals
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
        }
        // allele checking between data
        if(!heidioffFlag || jointsmrflag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            allele_check_multi(&bdata, etrait, &gdata1);
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, etrait, &gdata1);
            }
        } else
        {
            allele_check_multi(etrait,&gdata1);
        }
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata1); ngwas = median(gdata1.splSize);
        }
        for(int i=0;i<besdNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        // open .smr and .multismr for writing output
        long itemcount=0,itercountmlt=0;
        string outstr="";
        // header for .smr
        string smrfile0 = string(outFileName)+".smr";
        FILE* smr0;
        if(gwasFileName!=NULL) {
            smr0 = fopen(smrfile0.c_str(), "w");
            if (!(smr0)) {
                printf("ERROR: open error %s\n", smrfile0.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n";
            if(fputs_checked(outstr.c_str(),smr0))
            {
                printf("ERROR: error in writing file %s .\n", smrfile0.c_str());
                exit(EXIT_FAILURE);
            }
        }
        // header for .ppa
        string smrfile2 = string(outFileName)+".ppa";
        FILE* smr2 = fopen(smrfile2.c_str(), "w");
        if (!(smr2)) {
            printf("ERROR: open error %s\n", smrfile2.c_str());
            exit(1);
        }
        outstr="Chr\t";
        int j = 0;
        for(int i=0; i<besdNum;i++)
        {
            j = i+1;
            outstr+="Expo"+atos(j)+"_ID"+'\t'+"Expo"+atos(j)+"_Gene"+'\t'+"Expo"+atos(j)+"_bp"+'\t';
        }
        for(int i=0; i<combins.size();i++)
        {
            outstr+="PP"+atos(i+1)+"(";
            for(int j=0;j<expoNum;j++)
            {
                outstr+=atos(combins[i][j]);
                if(j<expoNum-1) outstr+=":";
            }
            if(i<(combins.size()-1)) outstr+=")\t";
        }
        outstr+=")\tp_HEIDI\tnsnp_HEIDI\n";
        if(fputs_checked(outstr.c_str(),smr2))
        {
            printf("ERROR: in writing file %s .\n", smrfile2.c_str());
            exit(EXIT_FAILURE);
        }
        // compute SMR effect size for all exposure probes and output
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        for(int i=0;i<besdNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func(smrrltstmp, NULL, &bdata,&gdata1,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                for(int j=0;j<smrrltstmp.size();j++)
                {
                    outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr0))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile0.c_str());
                        exit(EXIT_FAILURE);
                    }
                }                
            } else {
                continue;
            }
        }
        for(int k=0;k<probNum.size();k++)
        {   
            itemcount = itemcount + probNum[k]; 
        }
        printf("SMR analysis results of %ld exposure probes have been saved in the file %s .\n",itemcount,smrfile0.c_str());
        fclose(smr0);
        printf("\nPerforming multiple exposures OPERA analysis (including multi-exposure HEIDI tests) ... \n");
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposure probes with significant instruments are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        }        
        // loop with GWAS COJO loci
        double cr=0;
        for(int ii=0;ii<ldata._include.size();ii++)
        {
            double desti=1.0*ii/(ldata._include.size());
            if(desti>=cr)
            {
                printf("%3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            vector<SMRRLT> smrrltsbf;
            vector<long> probNumbf;

            int locichr=ldata._chr[ii];
            int locibp=ldata._bp[ii];

            int lowerbounder=(locibp-exposure_probe_wind)>0?(locibp-exposure_probe_wind):0;
            int upperbounder=locibp+exposure_probe_wind;
            
            for(int i=0;i<besdNum;i++)
            {
                int countNum = 0;
                for(int j=0;j<probNum[i];j++)
                {
                   int bptmp=smrrlts[i][j].Probe_bp;
                   if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                    }
                }
                probNumbf.push_back(countNum);
                //printf("\n%ld probes from outcome%d in the cis region [%d, %d] of exposure probe %s are inclued in the analysis.\n", etrait[i]._include.size(),i+1,lowerbounder,upperbounder,traitname.c_str());
            }

            int expocout = 0; 
            for(int i=0;i<probNumbf.size();i++) {
                if(probNumbf[i]>0) expocout+=1;
            }

            if(expocout==expoNum) {
                // illustrate all the combinations
                vector<vector<int>> indexall;
                for(int i=0;i<probNumbf.size();i++) {
                    vector<int> index(probNumbf[i]);
                    std::iota(index.begin(),index.end(),0);
                    indexall.push_back(index);
                }
                vector<vector<int>> combines;
                permute_vector(indexall, combines);
                //if (combines.size()==0) printf("\nWARNING: Less than %ld outcomes included in the analysis. \nThe multi-SMR and multi-HEIDI test will skip for exposure probe %s\n",expoNum,traitname.c_str());
                printf("\nThere are %ld possible combinations to test in a %ldKb window at %ldbp on chromosome %ld in the joint SMR model\n",combines.size(),op_wind,locibp,locichr);
                for(int i=0; i<combines.size(); i++)
                {                    
                    vector<eqtlInfo> esdatabf(expoNum);
                    vector<vector<string>> prb_cojolist(expoNum);
                    vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);
                    vector<string> outconamec(besdNum), outcogenec(besdNum);
                    vector<long> outcobpc(besdNum);
                    vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum);
                    vector<gwasData> esdatain(expoNum);
                    MatrixXd lh(2,expoNum);
                    // get the probe and gene information for output
                    outstr=atos(locichr)+'\t';
                    long postmp = 0;
                    for(int t=0; t<besdNum; t++)
                    {   
                        long idxtmp = combines[i][t]+postmp;
                        outconamec[t] = smrrltsbf[idxtmp].ProbeID;
                        outcogenec[t] = smrrltsbf[idxtmp].Gene;
                        outcobpc[t] = smrrltsbf[idxtmp].Probe_bp;
                        postmp = postmp + probNumbf[t];
                        if(!heidioffFlag || jointsmrflag) {
                            esdata[t]._include.clear();
                            map<string, int>::iterator itt;
                            eqtlInfo esdatatmp;
                            itt = esdata[t]._probe_name_map.find(outconamec[t]);
                            if(itt != esdata[t]._probe_name_map.end()) {
                                esdata[t]._include.push_back(itt->second);
                                e2econvert(&esdata[t], &esdatatmp);
                                esdatabf[t] = esdatatmp;
                            }
                            if(targetcojosnplstName!=NULL) {
                                // find the target probe COJO signals
                                map<string, vector<string>>::iterator prb_pos;
                                prb_pos = prb_cojosnps.find(outconamec[t]);
                                vector<string> navector; navector.push_back("");
                                if(prb_pos!=prb_cojosnps.end()) {
                                    prb_cojolist[t]=prb_pos->second;
                                } else { prb_cojolist[t]=navector; }
                            }
                        }
                        outstr+=outconamec[t]+'\t'+outcogenec[t]+'\t'+atos(outcobpc[t])+'\t';
                    }
                    vector<SMRRLT> smrrlts_joint;
                    if(esdatabf.size() == expoNum) {
                        if(targetcojosnplstName!=NULL) {
                            multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, prb_cojolist, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                        } else { multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata1, esdatabf, ngwas, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor); }
                    }
                    if(smrrlts_joint.size() == expoNum) {
                        //get the bxy, sigma_b and sigma_e from joint-SMR
                        for(int t=0; t<expoNum; t++)
                        {
                            bxy[t] = smrrlts_joint[t].b_SMR;
                            sigma_e[t] = pow(smrrlts_joint[t].se_SMR,2);
                            c[t] = 1+sigma_e[t]/sigma_b[t];
                        }
                        // multi-HEIDI test
                        vector<SMRRLT> smrrltsheidi;
                        if(! heidioffFlag){
                            multi_heidi_func(smrrltsheidi, NULL, &bdata, &gdata1, esdatabf, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                        }
                        // get the H0 and H1 prior pi and likelihood
                        const double PI = 3.141592653589793238463;
                        for(int t=0;t<expoNum;t++) {
                            lh(0,t)=pow(2*PI,-0.5)*pow(sigma_e[t],-0.5)*exp(-1*pow(bxy[t],2)/(2*sigma_e[t]));
                            lh(1,t)=pow(2*PI,-0.5)*pow(c[t]*sigma_b[t],-0.5)*exp((1/c[t]-1)*pow(bxy[t],2)/(2*sigma_e[t]));
                        }
                        // caculate the posterier probablity
                        for(int i=0;i<combNum;i++){
                            HH[i]=1.0;
                            for(int t=0;t<expoNum;t++)
                            {
                                HH[i] *= lh(combins[i][t],t);
                            }
                        }
                        float POall = 0;
                        for(int i=0;i<combNum;i++) {
                            PO[i] = HH[i]*prior[i];
                            POall+=PO[i];
                        }
                        for(int i=0;i<combNum;i++) {
                            PP[i] = PO[i]/POall;
                            outstr = outstr + atos(PP[i])+'\t';
                        }
                        bool sigflag = false;
                        for(int i=1;i<combNum;i++) {
                            if(PP[i]>=0.8) sigflag = true;
                        }
                        if(! heidioffFlag) {
                            outstr=outstr+(smrrltsheidi[0].p_HET >= 0 ? dtos(smrrltsheidi[0].p_HET) : "NA") + '\t' + (smrrltsheidi[0].nsnp > 0 ? atos(smrrltsheidi[0].nsnp+1) : "NA") + '\n';
                        } else {
                            outstr=outstr + "NA" + '\t' + "NA" + '\n';
                        }
                        if(sigflag) {
                            itercountmlt+=1;
                            if(fputs_checked(outstr.c_str(),smr2))
                            {
                                printf("ERROR: in writing file %s .\n", smrfile2.c_str());
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                
                }

            } else {
                continue;
            }
            
        }         
        fclose(smr2);
        //printf("\nSMR and HEIDI analyses for molecular traits completed.\nSMR and heterogeneity analysis results of %ld outcome probes (%ld exposure probe) have been saved in the file %s.\n",itemcount,etraitcount,smrfile1.c_str());
        printf("\nOPERA analyses for %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s .\n",expoNum,itercountmlt,smrfile2.c_str());
    }

    void multi_joint_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        MatrixXd _X;
        MTSMRWKEXP smrwk;
        int expoNum = esdata.size();
        uint64_t probNum = esdata[0]._include.size();
        smrwk.bxz.resize(expoNum); smrwk.sexz.resize(expoNum);
        smrwk.zxz.resize(expoNum); smrwk.freq.resize(expoNum);
        // compute the average GWAS sample size
        // double ngwas = 0.0;
        // ngwas = median(gdata->splSize);
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            int i=esdata[0]._include[ii];
            int probechr=esdata[0]._epi_chr[i];
            string probename=esdata[0]._epi_prbID[i];
            string probegene=esdata[0]._epi_gene[i];
            
            init_smr_wk_mlt(&smrwk);
            smrwk.cur_prbidx = i;
            smrwk.cur_chr = probechr;        
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0) {
                SMRRLT currlt;
                smrrlts.push_back(currlt);
                continue;
            }

            vector<int> max_ids, max_ids_ref;
            vector<string> maxSNPs;
            for(int t=0; t<smrwk.zxz.size();t++)
            {
                long maxid=max_abs_id(smrwk.zxz[t]);
                max_ids_ref.push_back(maxid);
                max_ids.push_back(maxid);
                maxSNPs.push_back(smrwk.rs[maxid]);
            }

            // joint bzy only include top SNPs for each exposure
            MTSMRWKEXP smrwk_joint;
            extract_smrwk(&smrwk,max_ids,&smrwk_joint);
            bool minus_2p = true;
            make_XMat(bdata,smrwk_joint.curId, _X, minus_2p);
            // get the summary data
            vector<VectorXd> _bxz(expoNum),_sexz(expoNum),_zsxz(expoNum);
            VectorXd  _byz(expoNum),_seyz(expoNum), _byz_adj(expoNum),_seyz_adj(expoNum);
            for(int t=0; t<expoNum; t++) {
                _bxz[t].resize(expoNum);
                _sexz[t].resize(expoNum);
                _zsxz[t].resize(expoNum);
            }
            for(int j=0;j<max_ids.size();j++)
            {
                for(int t=0;t<smrwk_joint.bxz.size();t++)
                {
                _bxz[t][j]=smrwk_joint.bxz[t][j];
                _sexz[t][j]=smrwk_joint.sexz[t][j];
                _zsxz[t][j]=smrwk_joint.zxz[t][j];
                }
                _byz[j]=smrwk_joint.byz[j];
                _seyz[j]=smrwk_joint.seyz[j];
            }
            
            CommFunc::getUnique(max_ids);
            if(max_ids.size() > 1) {
                vector<int> for_idx, back_idx;
                CommFunc::match_only(max_ids,max_ids_ref,for_idx);
                CommFunc::match_only(max_ids_ref,max_ids,back_idx);

                VectorXd _byz_joint(max_ids.size()), _seyz_joint(max_ids.size()),_byz_tmp(max_ids.size()), _seyz_tmp(max_ids.size());
                // run joint analysis on _bzy
                VectorXd DJ(max_ids.size());
                MatrixXd _Xtmp(_X.rows(),max_ids.size()), D(max_ids.size(),max_ids.size()), V(max_ids.size(),max_ids.size()), Vinv(max_ids.size(),max_ids.size());
                int nld = _X.rows();
                // compute variance y
                double ypy=0.0, sigma_c=0.0;
                for(int s=0; s<max_ids.size(); s++) {
                    _byz_tmp[s] =_byz[for_idx[s]];
                    _seyz_tmp[s] = _seyz[for_idx[s]];
                    _Xtmp.col(s) = _X.col(for_idx[s]);
                    double xpx = _X.col(for_idx[s]).transpose() * _X.col(for_idx[s]);
                    DJ[s] = xpx * ngwas/(double)nld;
                    ypy += DJ[s] * pow(_seyz_tmp[s],2) * (ngwas-1) + DJ[s] * pow(_byz_tmp[s],2);
                }
                ypy/=double(max_ids.size());
                V = (_Xtmp.transpose() * _Xtmp) * (ngwas/(double)nld);
                D = V.diagonal().asDiagonal();
                // do eigen decomposition
                SelfAdjointEigenSolver<MatrixXd> eigensolver(V);
                MatrixXd evec = eigensolver.eigenvectors();
                VectorXd eval = eigensolver.eigenvalues();
                vector<int> eigen_index;
                float eigen_thresh = 1e-8;
                for(int i=0; i<eval.size(); i++) {
                    if(eval[i] > eigen_thresh) eigen_index.push_back(i);
                }
                MatrixXd evec_trunc(evec.rows(),eigen_index.size());
                VectorXd eval_trunc(eigen_index.size());
                for(int s=0; s<eigen_index.size(); s++) {
                    evec_trunc.col(s) = evec.col(eigen_index[s]);
                    eval_trunc(s) = eval(eigen_index[s]);
                }
                Vinv = evec_trunc * eval_trunc.asDiagonal().inverse() * evec_trunc.transpose();
                _byz_joint = Vinv * D * _byz_tmp;
                sigma_c = (ypy - _byz_joint.transpose() * D * _byz_tmp)/(ngwas - max_ids.size());
                if(sigma_c > 0) {
                    _seyz_joint = Vinv.diagonal().cwiseSqrt() * sqrt(sigma_c);
                } else {
                    _seyz_joint = _seyz_tmp;
                }
                
                for(int b=0; b<back_idx.size(); b++) {
                    _byz_adj[b] = _byz_joint[back_idx[b]];
                    _seyz_adj[b] = _seyz_joint[back_idx[b]];
                }

            } else {
                _byz_adj = _byz;
                _seyz_adj = _seyz;
            }
            // compute the joint bxy
            for(int t=0; t<expoNum; t++)
            {
                SMRRLT currlt;
                // compute and output conditional SMR
                double byzt=_byz_adj[t], seyzt=_seyz_adj[t];
                double bxzt=_bxz[t][t], sexzt=_sexz[t][t];
                double bxy_val = byzt / bxzt;
                double sexy_val = -9;
                sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                currlt.b_SMR=bxy_val;
                currlt.se_SMR=sexy_val;
                currlt.p_SMR=pxy_val;
                smrrlts.push_back(currlt);
            }
        }
    }

    void multi_joint_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, vector<vector<string>> &prbcojolist, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        MatrixXd _X;
        MTSMRWKEXP smrwk;
        int expoNum = esdata.size();
        uint64_t probNum = esdata[0]._include.size();
        smrwk.bxz.resize(expoNum); smrwk.sexz.resize(expoNum);
        smrwk.zxz.resize(expoNum); smrwk.freq.resize(expoNum);
        // compute the average GWAS sample size
        // double ngwas = 0.0;
        // ngwas = median(gdata->splSize);
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            int i=esdata[0]._include[0];
            int probechr=esdata[0]._epi_chr[i];            
            init_smr_wk_mlt(&smrwk);
            smrwk.cur_prbidx = i;
            smrwk.cur_chr = probechr;        
            // get all the summary data
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0 || smrwk.zxz.size()!=expoNum) {
                SMRRLT currlt;
                smrrlts.push_back(currlt);
                continue;
            }
            // compute the mean ypy
            vector<double> ypy(smrwk.byz.size());
            double ypymedian = 0.0;
            for(int s=0; s<smrwk.byz.size(); s++) {
                double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*ngwas;
                ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (ngwas-1) + pow(smrwk.byz[s],2));
            }
            ypymedian = median(ypy);
            // find the top and COJO SNPs index in smrwk
            vector<int> max_ids, cojo_ids;
            vector<vector<int>> all_ids(expoNum);
            for(int t=0; t<expoNum;t++)
            {
                long maxid=max_abs_id(smrwk.zxz[t]);
                max_ids.push_back(maxid);
                all_ids[t].push_back(maxid);
                cojo_ids.push_back(maxid);
                for(int j=0; j<prbcojolist[t].size(); j++) {
                    int cojoidx = find(smrwk.rs.begin(),smrwk.rs.end(),prbcojolist[t][j])-smrwk.rs.begin();
                    if(cojoidx!=smrwk.rs.size()) {
                        all_ids[t].push_back(cojoidx);
                        cojo_ids.push_back(cojoidx);
                    }
                }
                CommFunc::getUnique(all_ids[t]);
            }
            CommFunc::getUnique(cojo_ids);
            // only include COJO SNPs for summary data
            MTSMRWKEXP smrwk_joint;
            extract_smrwk(&smrwk,cojo_ids,&smrwk_joint);
            bool minus_2p = true;
            make_XMat(bdata,smrwk_joint.curId, _X, minus_2p);
            // find the top and conditional SNP positions in cojo_ids
            vector<int> max_pos(expoNum);
            vector<vector<int>> cond_pos(expoNum);
            for(int t=0; t<expoNum; t++)
            {
                max_pos[t] = find(cojo_ids.begin(),cojo_ids.end(),max_ids[t])-cojo_ids.begin();                
                for(int k=0; k<cojo_ids.size(); k++) cond_pos[t].push_back(k);
                for(int m=0; m<all_ids[t].size(); m++) {
                    long tmp_ids = find(cojo_ids.begin(),cojo_ids.end(),all_ids[t][m])-cojo_ids.begin();
                    long cond_ids = find(cond_pos[t].begin(),cond_pos[t].end(),tmp_ids)-cond_pos[t].begin();
                    cond_pos[t].erase(cond_pos[t].begin()+cond_ids);
                }
            }
            // check the top and COJO SNP list
            vector<string> maxSNPs;
            vector<vector<string>> condSNPs(expoNum);
            for(int t=0; t<expoNum; t++)
            {
                maxSNPs.push_back(smrwk_joint.rs[max_pos[t]]);
                for(int m=0; m<cond_pos[t].size(); m++) {
                    condSNPs[t].push_back(smrwk_joint.rs[cond_pos[t][m]]);
                }
            }
            // get the summary data
            vector<VectorXd> _bxz(expoNum),_sexz(expoNum),_zsxz(expoNum),_byz_adj(expoNum),_seyz_adj(expoNum);
            VectorXd  _byz(cojo_ids.size()),_seyz(cojo_ids.size());            
            for(int t=0; t<expoNum; t++) {
                _bxz[t].resize(cojo_ids.size());
                _sexz[t].resize(cojo_ids.size());
                _zsxz[t].resize(cojo_ids.size());
                _byz_adj[t].resize(cojo_ids.size());
                _seyz_adj[t].resize(cojo_ids.size());
            }
            for(int j=0;j<cojo_ids.size();j++)
            {
                for(int t=0;t<expoNum;t++)
                {
                _bxz[t][j]=smrwk_joint.bxz[t][j];
                _sexz[t][j]=smrwk_joint.sexz[t][j];
                _zsxz[t][j]=smrwk_joint.zxz[t][j];
                }
                _byz[j]=smrwk_joint.byz[j];
                _seyz[j]=smrwk_joint.seyz[j];
            }
            // run joint GWAS effect            
            for(int t=0; t<expoNum; t++)
            {
                run_joint_effect_func(max_pos[t],cond_pos[t],_byz,_seyz,_byz_adj[t],_seyz_adj[t],_X,ngwas,ypymedian);
            }
            // compute the joint bxy
            for(int t=0; t<expoNum; t++)
            {
                SMRRLT currlt;
                // compute and output joint SMR effect
                double byzt=_byz_adj[t][max_pos[t]], seyzt=_seyz_adj[t][max_pos[t]];
                double bxzt=_bxz[t][max_pos[t]], sexzt=_sexz[t][max_pos[t]];
                double bxy_val = byzt / bxzt;
                double sexy_val = -9;
                sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                currlt.ProbeID = esdata[t]._epi_prbID[0];
                currlt.ProbeChr = esdata[t]._epi_chr[0];
                currlt.SNP = maxSNPs[t];
                currlt.b_SMR = bxy_val;
                currlt.se_SMR = sexy_val;
                currlt.p_SMR = pxy_val;
                smrrlts.push_back(currlt);
            }
        }
    }

    void run_joint_effect_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy)
    {
        int nld = X.rows(), jointsize_tmp = condpos.size() + 1;
        MatrixXd Xtmp(nld,jointsize_tmp);
        Xtmp.col(0) = X.col(maxpos);
        for(int s=1; s<jointsize_tmp; s++) {
            Xtmp.col(s) = X.col(condpos[s-1]);
        }
        // LD pruning with top > 0.9
        float ldr2_thresh = 0.9;
        vector<int> sn_ids, snp_ids;
        VectorXd ld_v;
        ld_calc_o2m(ld_v,0,Xtmp);
        for(int i=0;i<ld_v.size();i++)
        {
           double ldr2tmp=ld_v(i)*ld_v(i);
           if(ldr2tmp < ldr2_thresh) sn_ids.push_back(i);
        }
        // initialize the adjusted value as marginal effect
        byz_adj = byz; seyz_adj = seyz;
        if(sn_ids.size() > 0) {
            MatrixXd X_tmp(Xtmp.rows(),sn_ids.size());
            for(int s=0; s<sn_ids.size(); s++) {
                X_tmp.col(s) = Xtmp.col(sn_ids[s]);
            }
            // LD pruning with pairwise LD > 0.9
            MatrixXd C;
            vector<int> rm_ID1;
            cor_calc(C, X_tmp);
            int qi = 0; int m = sn_ids.size();
            if (ldr2_thresh < 1) rm_cor_sbat(C, sqrt(ldr2_thresh), m, rm_ID1);
            for (int i=0 ; i<m ; i++) {
                if (rm_ID1.size() == 0) snp_ids.push_back(sn_ids[i]);
                else {
                    if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                    else snp_ids.push_back(sn_ids[i]);
                }
            }
            // extract snp_ids for joint analysis
            vector<int> jointpos;
            jointpos.push_back(maxpos);
            if(snp_ids.size() > 0) {
                for(int s=0; s<snp_ids.size(); s++) {
                    jointpos.push_back(condpos[snp_ids[s]-1]);
                }
                int jointsize = jointpos.size();
                VectorXd byz_joint(jointsize), seyz_joint(jointsize), byz_tmp(jointsize), seyz_tmp(jointsize);
                // run joint analysis on bzy
                MatrixXd X_joint(X.rows(),jointsize), D(jointsize,jointsize), V(jointsize,jointsize), Vinv(jointsize,jointsize);
                double sigma_c=0.0;
                for(int s=0; s<jointsize; s++) {
                    byz_tmp[s] = byz[jointpos[s]];
                    seyz_tmp[s] = seyz[jointpos[s]];
                    X_joint.col(s) = X.col(jointpos[s]);
                }
                V = (X_joint.transpose() * X_joint) * (ngwas/(double)nld);
                D = V.diagonal().asDiagonal();
                // do eigen decomposition
                SelfAdjointEigenSolver<MatrixXd> eigensolver(V);
                MatrixXd evec = eigensolver.eigenvectors();
                VectorXd eval = eigensolver.eigenvalues();
                vector<int> eigen_index;
                float eigen_thresh = 1e-8;
                for(int i=0; i<eval.size(); i++) {
                    if(eval[i] > eigen_thresh) eigen_index.push_back(i);
                }
                MatrixXd evec_trunc(evec.rows(),eigen_index.size());
                VectorXd eval_trunc(eigen_index.size());
                for(int s=0; s<eigen_index.size(); s++) {
                    evec_trunc.col(s) = evec.col(eigen_index[s]);
                    eval_trunc(s) = eval(eigen_index[s]);
                }
                Vinv = evec_trunc * eval_trunc.asDiagonal().inverse() * evec_trunc.transpose();
                // Vinv = V.inverse();
                // compute joint effect and se
                byz_joint = Vinv * D * byz_tmp;
                sigma_c = (ypy - byz_joint.transpose() * D * byz_tmp)/(ngwas - jointsize);
                if(sigma_c > 0) {
                    seyz_joint = Vinv.diagonal().cwiseSqrt() * sqrt(sigma_c);
                } else {
                    seyz_joint = seyz_tmp;
                }
                // push back the joint effect
                for(int b=0; b<jointsize; b++) {
                    byz_adj[jointpos[b]] = byz_joint[b];
                    seyz_adj[jointpos[b]] = seyz_joint[b];
                }

            } 
        } 
    }

    void multi_cond_smr_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        MatrixXd _X;
        MTSMRWKEXP smrwk;
        int expoNum = esdata.size();
        uint64_t probNum = esdata[0]._include.size();
        smrwk.bxz.resize(expoNum); smrwk.sexz.resize(expoNum);
        smrwk.zxz.resize(expoNum); smrwk.freq.resize(expoNum);
        // compute the average GWAS sample size
        double ngwas;
        for(int i=0; i<gdata->splSize.size(); i++) ngwas+=gdata->splSize[i];
        ngwas/=double(gdata->splSize.size());
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            int i=esdata[0]._include[ii];
            int probechr=esdata[0]._epi_chr[i];
            string probename=esdata[0]._epi_prbID[i];
            string probegene=esdata[0]._epi_gene[i];
            
            init_smr_wk_mlt(&smrwk);
            smrwk.cur_prbidx = i;
            smrwk.cur_chr = probechr;        
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0) {
                SMRRLT currlt;
                smrrlts.push_back(currlt);
                continue;
            }

            vector<int> max_ids;
            vector<string> maxSNPs;
            for(int t=0; t<smrwk.zxz.size();t++)
            {
                long maxid=max_abs_id(smrwk.zxz[t]);
                max_ids.push_back(maxid);
                maxSNPs.push_back(smrwk.rs[maxid]);
            }
            // conditional SMR for each exposure, only include top SNPs
            MTSMRWKEXP smrwk_condi;
            extract_smrwk(&smrwk,max_ids,&smrwk_condi);
            bool minus_2p = true;
            make_XMat(bdata,smrwk_condi.curId, _X,minus_2p);
            // get the summary data
            vector<VectorXd> _bxz(expoNum),_sexz(expoNum),_zsxz(expoNum);
            VectorXd  _byz(expoNum),_seyz(expoNum), _byz_adj(expoNum),_seyz_adj(expoNum);
            for(int t=0; t<expoNum; t++) {
                _bxz[t].resize(expoNum);
                _sexz[t].resize(expoNum);
                _zsxz[t].resize(expoNum);
            }
            for(int j=0;j<max_ids.size();j++)
            {
                for(int t=0;t<smrwk_condi.bxz.size();t++)
                {
                _bxz[t][j]=smrwk_condi.bxz[t][j];
                _sexz[t][j]=smrwk_condi.sexz[t][j];
                _zsxz[t][j]=smrwk_condi.zxz[t][j];
                }
                _byz[j]=smrwk_condi.byz[j];
                _seyz[j]=smrwk_condi.seyz[j];
            }
            // run conditional analysis on _bzy
            float ldr2_thresh = 0.9;
            for(int t=0; t<expoNum; t++)
            {
                SMRRLT currlt;
                // remove snps in LD with top > 0.9
                vector<int> sn_ids, snp_ids;
                VectorXd ld_v;
                ld_calc_o2m(ld_v,t,_X);
                for(int i=0;i<ld_v.size();i++)
                {
                   double ldr2tmp=ld_v(i)*ld_v(i);
                   if(ldr2tmp < ldr2_thresh) sn_ids.push_back(i);
                }
                if(sn_ids.size() > 0) {
                    MatrixXd _X_tmp(_X.rows(),sn_ids.size());
                    for(int s=0; s<sn_ids.size(); s++) {
                        _X_tmp.col(s) = _X.col(sn_ids[s]);
                    }
                    // remove snps in pairwise LD > 0.9
                    MatrixXd C;
                    vector<int> rm_ID1;
                    cor_calc(C, _X_tmp);
                    int qi = 0; int m = sn_ids.size();
                    if (ldr2_thresh < 1) rm_cor_sbat(C, sqrt(ldr2_thresh), m, rm_ID1);
                    for (int i=0 ; i<m ; i++) {
                        if (rm_ID1.size() == 0) snp_ids.push_back(sn_ids[i]);
                        else {
                            if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                            else snp_ids.push_back(sn_ids[i]);
                        }
                    }
                    // extract snp_ids for conditional analysis
                    int nld = _X.rows();
                    if(snp_ids.size() > 0) {
                    MatrixXd _X_condi(nld,snp_ids.size());
                    for(int s=0; s<snp_ids.size(); s++) {
                        _X_condi.col(s) = _X.col(snp_ids[s]);
                    }
                    VectorXd _byz_condi(snp_ids.size()), _seyz_condi(snp_ids.size()), DJ(snp_ids.size());
                    MatrixXd D22(snp_ids.size(),snp_ids.size());
                    double ypy=0.0, bvar_condi=0.0, sigma_c=0.0;
                    for(int s=0; s<snp_ids.size(); s++) {
                        _byz_condi[s] = _byz[snp_ids[s]];
                        _seyz_condi[s] = _seyz[snp_ids[s]];
                        double xpx = _X_condi.col(s).transpose() * _X_condi.col(s);
                        DJ[s] = xpx * ngwas/(double)nld;
                        ypy += DJ[s] * pow(_seyz_condi[s],2) * (ngwas-1) + DJ[s] * pow(_byz_condi[s],2);
                    }
                    D22 = DJ.asDiagonal();
                    ypy/=double(snp_ids.size());
                    double tmp = 1/(_X.col(t).transpose() * _X.col(t)) * ((double)nld/ngwas);
                    _byz_adj[t] = _byz[t] -  1/(_X.col(t).transpose() * _X.col(t)) * (_X.col(t).transpose() * _X_condi) * (_X_condi.transpose() * _X_condi).inverse() * (_byz_condi.cwiseProduct((_X_condi.transpose() * _X_condi).diagonal()));
                    bvar_condi = _byz_condi.transpose() * D22 * (_X_condi.transpose() * _X_condi).inverse() * (nld/ngwas) * D22 * _byz_condi;
                    sigma_c = ypy - bvar_condi - _byz_adj[t] * (1/tmp) * _byz[t];
                        if(sigma_c > 0){
                            _seyz_adj[t] = sqrt((sigma_c * tmp)/(ngwas - 1 - _byz_condi.size()));
                        } else{
                            _seyz_adj[t] = _seyz[t];
                        }
                    }   else {
                        _byz_adj[t] = _byz[t];
                        _seyz_adj[t] = _seyz[t];
                    }
                }   else {
                    _byz_adj[t] = _byz[t];
                    _seyz_adj[t] = _seyz[t];
                }
                // compute and output conditional SMR
                double byzt=_byz_adj[t], seyzt=_seyz_adj[t];
                double bxzt=_bxz[t][t], sexzt=_sexz[t][t];
                double bxy_val = byzt / bxzt;
                double sexy_val = -9;
                sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                currlt.b_SMR=bxy_val;
                currlt.se_SMR=sexy_val;
                currlt.p_SMR=pxy_val;
                smrrlts.push_back(currlt);
            }
        }
    }
    
    void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
    {
           long numRows = matrix.rows();
           long numCols = matrix.cols()-1;

           if( colToRemove < numCols )
               matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
           matrix.conservativeResize(numRows,numCols);
    }

    void e2gconvert(eqtlInfo* etrait, gwasData* gdata, int &ii)
    {
        gdata->allele_1.resize(etrait->_esi_include.size());
        gdata->allele_2.resize(etrait->_esi_include.size());
        gdata->byz.resize(etrait->_esi_include.size());
        gdata->seyz.resize(etrait->_esi_include.size());
        gdata->freq.resize(etrait->_esi_include.size());
        gdata->pvalue.resize(etrait->_esi_include.size());
        gdata->splSize.resize(etrait->_esi_include.size());
        gdata->snpName.resize(etrait->_esi_include.size());
        
        for(int j=0;j<etrait->_esi_include.size();j++){
            gdata->seyz[j]=-9;
            gdata->pvalue[j]=-9;
        }
        gdata->_include.clear();
        
        int count=0;
        if(etrait->_rowid.empty())
        {
            for (int j = 0; j<etrait->_esi_include.size(); j++)
            {
                if (fabs(etrait->_bxz[ii][etrait->_esi_include[j]] + 9) > 1e-6)
                {
                    gdata->byz[j]=etrait->_bxz[ii][etrait->_esi_include[j]];
                    gdata->seyz[j]=etrait->_sexz[ii][etrait->_esi_include[j]];
                    double z=etrait->_bxz[ii][etrait->_esi_include[j]]/etrait->_sexz[ii][etrait->_esi_include[j]];
                    double p=pchisq(z*z,1);
                    gdata->pvalue[j]=p;
                    gdata->snpName[j]=etrait->_esi_rs[etrait->_esi_include[j]];
                    gdata->allele_1[j]=etrait->_esi_allele1[etrait->_esi_include[j]];
                    gdata->allele_2[j]=etrait->_esi_allele2[etrait->_esi_include[j]];
                    gdata->_include.push_back(etrait->_esi_include[j]); // row id selected
                    count++;
                }
            }
        }
        else
        {
            uint64_t beta_start=etrait->_cols[ii<<1];
            uint64_t se_start=etrait->_cols[1+(ii<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=etrait->_rowid[beta_start+j];
                int idx=(int)(find(etrait->_esi_include.begin(), etrait->_esi_include.end(), ge_rowid)-etrait->_esi_include.begin());
                //if(idx!=etrait->_esi_include.size()) cout<<idx<<endl;
                if(idx<etrait->_esi_include.size())
                {
                    gdata->byz[idx]=etrait->_val[beta_start+j];
                    gdata->seyz[idx]=etrait->_val[se_start+j];
                    double z=etrait->_val[beta_start+j]/etrait->_val[se_start+j];
                    double p=pchisq(z*z,1);
                    gdata->pvalue[idx]=p;
                    gdata->snpName[idx]=etrait->_esi_rs[ge_rowid];
                    gdata->allele_1[idx]=etrait->_esi_allele1[ge_rowid];
                    gdata->allele_2[idx]=etrait->_esi_allele2[ge_rowid];
                    gdata->_include.push_back(idx);
                    count++;
                }
                
            }
        }
    }
    void e2econvert(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        esdata->_esi_rs.resize(etrait->_esi_include.size());
        esdata->_esi_allele1.resize(etrait->_esi_include.size());
        esdata->_esi_allele2.resize(etrait->_esi_include.size());
        esdata->_esi_gd.resize(etrait->_esi_include.size());
        esdata->_esi_bp.resize(etrait->_esi_include.size());
        esdata->_esi_chr.resize(etrait->_esi_include.size());
        esdata->_esi_freq.resize(etrait->_esi_include.size());
        esdata->_esi_include.resize(etrait->_esi_include.size());
        esdata->_snp_name_map.clear();
        map<int,int> id_map;
        for (int j = 0; j<etrait->_esi_include.size(); j++)
        {
            string rs=etrait->_esi_rs[etrait->_esi_include[j]];
            esdata->_esi_rs[j]=rs;
            esdata->_esi_allele1[j]=etrait->_esi_allele1[etrait->_esi_include[j]];
            esdata->_esi_allele2[j]=etrait->_esi_allele2[etrait->_esi_include[j]];
            esdata->_esi_gd[j]=etrait->_esi_gd[etrait->_esi_include[j]];
            esdata->_esi_bp[j]=etrait->_esi_bp[etrait->_esi_include[j]];
            esdata->_esi_chr[j]=etrait->_esi_chr[etrait->_esi_include[j]];
            esdata->_esi_freq[j]=etrait->_esi_freq[etrait->_esi_include[j]];
            esdata->_esi_include[j]=j;
            esdata->_snp_name_map.insert(pair<string,int>(rs,j));
            id_map.insert(pair<int,int>(etrait->_esi_include[j],j));
        }
        esdata->_snpNum=etrait->_esi_include.size();
        
        // update the epi include information.
        esdata->_epi_bp.resize(etrait->_include.size());
        esdata->_epi_gd.resize(etrait->_include.size());
        esdata->_epi_chr.resize(etrait->_include.size());
        esdata->_epi_gene.resize(etrait->_include.size());
        esdata->_epi_orien.resize(etrait->_include.size());
        esdata->_epi_prbID.resize(etrait->_include.size());
        esdata->_include.resize(etrait->_include.size());
        esdata->_probe_name_map.clear();
        map<int,int> probe_map;
        for (int j = 0; j<etrait->_include.size(); j++)
        {
            string probe=etrait->_epi_prbID[etrait->_include[j]];
            esdata->_epi_prbID[j]=probe;
            esdata->_epi_bp[j] =  etrait->_epi_bp[etrait->_include[j]];
            esdata->_epi_gd[j] = etrait->_epi_gd[etrait->_include[j]];
            esdata->_epi_chr[j] = etrait->_epi_chr[etrait->_include[j]];
            esdata->_epi_gene[j] = etrait->_epi_gene[etrait->_include[j]];
            esdata->_epi_orien[j] = etrait->_epi_orien[etrait->_include[j]];
            esdata->_include[j] = j;
            esdata->_probe_name_map.insert(pair<string,int>(probe,j));
            probe_map.insert(pair<int,int>(etrait->_include[j],j));
        }
        esdata->_probNum = etrait->_include.size();       
        
        map<int, int>::iterator iter;
        if(etrait->_rowid.empty())
        {
            esdata->_bxz.resize(etrait->_include.size());
            esdata->_sexz.resize(etrait->_include.size());
            for(int i=0;i<etrait->_include.size();i++)
            {
                esdata->_bxz[i].resize(etrait->_esi_include.size());
                esdata->_sexz[i].resize(etrait->_esi_include.size());
            }
            
            for (int ii = 0; ii<etrait->_include.size(); ii++)
            {
                for (int j = 0; j<etrait->_esi_include.size(); j++)
                {
                    // Here the values in etrait->_include should be 0,1,2,....
                    esdata->_bxz[ii][j]=etrait->_bxz[ii][etrait->_esi_include[j]];
                    esdata->_sexz[ii][j]=etrait->_sexz[ii][etrait->_esi_include[j]];
                    
                }
            }
        }
        else
        {
            // Here the values in etrait->_include should be 0,1,2,....
            esdata->_cols.resize((etrait->_include.size()<<1)+1);
            esdata->_cols[0]=0;
            for (int i = 0; i<etrait->_include.size(); i++)
            {
                long ii = etrait->_include[i];
                uint64_t beta_start=etrait->_cols[ii<<1];
                uint64_t se_start=etrait->_cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                int real_num=0;
                for(int j=0;j<numsnps<<1;j++)
                {
                    int ge_rowid=etrait->_rowid[beta_start+j];
                    iter=id_map.find(ge_rowid);
                    if(iter!=id_map.end())
                    {
                        int sid=iter->second;
                        esdata->_rowid.push_back(sid);
                        esdata->_val.push_back(etrait->_val[beta_start+j]);
                        real_num++;
                    }
                }
                esdata->_cols[(i<<1)+1]=(real_num>>1)+esdata->_cols[i<<1];
                esdata->_cols[i+1<<1]=real_num+esdata->_cols[i<<1];
            }
            esdata->_valNum=esdata->_val.size();
        }
    }
    void e2econvert_old2(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        esdata->_esi_rs.resize(etrait->_esi_include.size());
        esdata->_esi_allele1.resize(etrait->_esi_include.size());
        esdata->_esi_allele2.resize(etrait->_esi_include.size());
        esdata->_esi_gd.resize(etrait->_esi_include.size());
        esdata->_esi_bp.resize(etrait->_esi_include.size());
        esdata->_esi_chr.resize(etrait->_esi_include.size());
        esdata->_esi_freq.resize(etrait->_esi_include.size());
        esdata->_esi_include.resize(etrait->_esi_include.size());
        esdata->_snp_name_map.clear();
        map<int,int> id_map;
        for (int j = 0; j<etrait->_esi_include.size(); j++)
        {
            string rs=etrait->_esi_rs[etrait->_esi_include[j]];
            esdata->_esi_rs[j]=rs;
            esdata->_esi_allele1[j]=etrait->_esi_allele1[etrait->_esi_include[j]];
            esdata->_esi_allele2[j]=etrait->_esi_allele2[etrait->_esi_include[j]];
            esdata->_esi_gd[j]=etrait->_esi_gd[etrait->_esi_include[j]];
            esdata->_esi_bp[j]=etrait->_esi_bp[etrait->_esi_include[j]];
            esdata->_esi_chr[j]=etrait->_esi_chr[etrait->_esi_include[j]];
            esdata->_esi_freq[j]=etrait->_esi_freq[etrait->_esi_include[j]];
            esdata->_esi_include[j]=j;
            esdata->_snp_name_map.insert(pair<string,int>(rs,j));
            id_map.insert(pair<int,int>(etrait->_esi_include[j],j));
        }
        esdata->_snpNum=etrait->_esi_include.size();
        
        // update the epi include information.
        esdata->_epi_bp.resize(etrait->_include.size());
        esdata->_epi_gd.resize(etrait->_include.size());
        esdata->_epi_chr.resize(etrait->_include.size());
        esdata->_epi_gene.resize(etrait->_include.size());
        esdata->_epi_orien.resize(etrait->_include.size());
        esdata->_epi_prbID.resize(etrait->_include.size());
        esdata->_include.resize(etrait->_include.size());
        esdata->_probe_name_map.clear();
        map<int,int> probe_map;
        for (int j = 0; j<etrait->_include.size(); j++)
        {
            string probe=etrait->_epi_prbID[etrait->_include[j]];
            esdata->_epi_prbID[j]=probe;
            esdata->_epi_bp[j] =  etrait->_epi_bp[etrait->_include[j]];
            esdata->_epi_gd[j] = etrait->_epi_gd[etrait->_include[j]];
            esdata->_epi_chr[j] = etrait->_epi_chr[etrait->_include[j]];
            esdata->_epi_gene[j] = etrait->_epi_gene[etrait->_include[j]];
            esdata->_epi_orien[j] = etrait->_epi_orien[etrait->_include[j]];
            esdata->_include[j] = j;
            esdata->_probe_name_map.insert(pair<string,int>(probe,j));
            probe_map.insert(pair<int,int>(etrait->_include[j],j));
        }
        esdata->_probNum = etrait->_include.size();       
        
        map<int, int>::iterator iter;
        if(etrait->_rowid.empty())
        {
            esdata->_bxz.resize(etrait->_include.size());
            esdata->_sexz.resize(etrait->_include.size());
            for(int i=0;i<etrait->_include.size();i++)
            {
                esdata->_bxz[i].resize(etrait->_esi_include.size());
                esdata->_sexz[i].resize(etrait->_esi_include.size());
            }
            
            for (int j = 0; j<etrait->_esi_include.size(); j++)
            {
                for (int ii = 0; ii<etrait->_include.size(); ii++)
                {
                    // Here the values in etrait->_include should be 0,1,2,....
                    esdata->_bxz[ii][j]=etrait->_bxz[ii][etrait->_esi_include[j]];
                    esdata->_sexz[ii][j]=etrait->_sexz[ii][etrait->_esi_include[j]];
                    
                }
            }
        }
        else
        {
            // Here the values in etrait->_include should be 0,1,2,....
            esdata->_cols.resize((etrait->_include.size()<<1)+1);
            esdata->_cols[0]=0;
            for (int i = 0; i<etrait->_include.size(); i++)
            {
                long ii = etrait->_include[i];
                uint64_t beta_start=etrait->_cols[ii<<1];
                uint64_t se_start=etrait->_cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                int real_num=0;
                for(int j=0;j<numsnps<<1;j++)
                {
                    int ge_rowid=etrait->_rowid[beta_start+j];
                    iter=id_map.find(ge_rowid);
                    if(iter!=id_map.end())
                    {
                        int sid=iter->second;
                        esdata->_rowid.push_back(sid);
                        esdata->_val.push_back(etrait->_val[beta_start+j]);
                        real_num++;
                    }
                }
                esdata->_cols[(i<<1)+1]=(real_num>>1)+esdata->_cols[i<<1];
                esdata->_cols[i+1<<1]=real_num+esdata->_cols[i<<1];
            }
            esdata->_valNum=esdata->_val.size();
        }
    }
    void e2econvert_old(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        esdata->_esi_rs.resize(etrait->_esi_include.size());
        esdata->_esi_allele1.resize(etrait->_esi_include.size());
        esdata->_esi_allele2.resize(etrait->_esi_include.size());
        esdata->_esi_gd.resize(etrait->_esi_include.size());
        esdata->_esi_bp.resize(etrait->_esi_include.size());
        esdata->_esi_chr.resize(etrait->_esi_include.size());
        esdata->_esi_freq.resize(etrait->_esi_include.size());
        esdata->_esi_include.resize(etrait->_esi_include.size());
        esdata->_snp_name_map.clear();
        map<int,int> id_map;
        for (int j = 0; j<etrait->_esi_include.size(); j++)
        {
            string rs=etrait->_esi_rs[etrait->_esi_include[j]];
            esdata->_esi_rs[j]=rs;
            esdata->_esi_allele1[j]=etrait->_esi_allele1[etrait->_esi_include[j]];
            esdata->_esi_allele2[j]=etrait->_esi_allele2[etrait->_esi_include[j]];
            esdata->_esi_gd[j]=etrait->_esi_gd[etrait->_esi_include[j]];
            esdata->_esi_bp[j]=etrait->_esi_bp[etrait->_esi_include[j]];
            esdata->_esi_chr[j]=etrait->_esi_chr[etrait->_esi_include[j]];
            esdata->_esi_freq[j]=etrait->_esi_freq[etrait->_esi_include[j]];
            esdata->_esi_include[j]=j;
            esdata->_snp_name_map.insert(pair<string,int>(rs,j));
            id_map.insert(pair<int,int>(etrait->_esi_include[j],j));
        }
        esdata->_snpNum=etrait->_esi_include.size();
        
        // Also do not forget to copy epi information from etrait to esdata.
        esdata->_epi_bp =  etrait->_epi_bp;
        esdata->_epi_gd = etrait->_epi_gd;
        esdata->_epi_chr = etrait->_epi_chr;
        esdata->_epi_gene = etrait->_epi_gene;
        esdata->_epi_orien = etrait->_epi_orien;
        esdata->_epi_prbID = etrait->_epi_prbID;
        esdata->_include = etrait->_include;
        esdata->_probe_name_map = etrait->_probe_name_map;
        esdata->_probNum = etrait->_probNum;
        
        map<int, int>::iterator iter;
        if(etrait->_rowid.empty())
        {
            esdata->_bxz.resize(etrait->_include.size());
            esdata->_sexz.resize(etrait->_include.size());
            for(int i=0;i<etrait->_include.size();i++)
            {
                esdata->_bxz[i].resize(etrait->_esi_include.size());
                esdata->_sexz[i].resize(etrait->_esi_include.size());
            }
            
            for (int j = 0; j<etrait->_esi_include.size(); j++)
            {
                for (int ii = 0; ii<etrait->_include.size(); ii++)
                {
                    // Here the values in etrait->_include should be 0,1,2,....
                    esdata->_bxz[ii][j]=etrait->_bxz[ii][etrait->_esi_include[j]];
                    esdata->_sexz[ii][j]=etrait->_sexz[ii][etrait->_esi_include[j]];
                    
                }
            }
        }
        else
        {
            // Here the values in etrait->_include should be 0,1,2,....
            esdata->_cols.resize((etrait->_include.size()<<1)+1);
            esdata->_cols[0]=0;
            for (int ii = 0; ii<etrait->_include.size(); ii++)
            {
                uint64_t beta_start=etrait->_cols[ii<<1];
                uint64_t se_start=etrait->_cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                int real_num=0;
                for(int j=0;j<numsnps<<1;j++)
                {
                    int ge_rowid=etrait->_rowid[beta_start+j];
                    iter=id_map.find(ge_rowid);
                    if(iter!=id_map.end())
                    {
                        int sid=iter->second;
                        esdata->_rowid.push_back(sid);
                        esdata->_val.push_back(etrait->_val[beta_start+j]);
                        real_num++;
                    }
                }
                esdata->_cols[(ii<<1)+1]=(real_num>>1)+esdata->_cols[ii<<1];
                esdata->_cols[ii+1<<1]=real_num+esdata->_cols[ii<<1];
            }
            esdata->_valNum=esdata->_val.size();
        }
    }
    void cis_xqtl_probe_include_only(eqtlInfo* eqtlinfo, double p_smr, int cis_itvl, string eqtlFileName)
    {
        vector<int> newIcld;
        if(eqtlinfo->_valNum==0)
        {
            for(uint32_t i=0;i<eqtlinfo->_probNum;i++)
            {
                double beta_top=0;
                double se_top=1;
                double pxz_top=1;
                int esi_id_top=-1;
                int epi_id_top=-1;
                int prbbp=eqtlinfo->_epi_bp[i];
                int prbchr=eqtlinfo->_epi_chr[i];
                int cisR=prbbp+cis_itvl*1000;
                int cisL=(prbbp-cis_itvl*1000)>0?(prbbp-cis_itvl*1000):0;
                
                for(uint32_t j=0;j<eqtlinfo->_snpNum;j++)
                {
                    int snpbp=eqtlinfo->_esi_bp[j];
                    int snpchr=eqtlinfo->_esi_chr[j];
                    if(snpchr==prbchr && snpbp>=cisL && snpbp<=cisR)
                    {
                        double beta=eqtlinfo->_bxz[i][j];
                        double se=eqtlinfo->_sexz[i][j];
                        if(ABS(se+9)<1e-6) continue;
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=pxz_top)
                        {
                            esi_id_top=j;
                            epi_id_top=i;
                            beta_top=beta;
                            se_top=se;
                            pxz_top=pxz;
                        }
                        
                    }
                }
                
                if(pxz_top<=p_smr)
                {
                    newIcld.push_back(eqtlinfo->_include[i]);
                }
            }
        }
        else
        {
            if(eqtlinfo->_val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.\n");
            }
            
            for(uint32_t i=0;i<eqtlinfo->_probNum;i++)
            {
                uint64_t proid=eqtlinfo->_include[i];
                uint64_t pos=eqtlinfo->_cols[proid<<1];
                uint64_t pos1=eqtlinfo->_cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                double beta_top=0;
                double se_top=1;
                double pxz_top=1;
                int esi_id_top=-1;
                int epi_id_top=-1;
                int prbbp=eqtlinfo->_epi_bp[proid];
                int prbchr=eqtlinfo->_epi_chr[proid];
                int cisR=prbbp+cis_itvl*1000;
                int cisL=(prbbp-cis_itvl*1000)>0?(prbbp-cis_itvl*1000):0;
                for(int j=0;j<num;j++)
                {
                    int esiid=eqtlinfo->_rowid[pos+j];
                    int snpbp=eqtlinfo->_esi_bp[esiid];
                    int snpchr=eqtlinfo->_esi_chr[esiid];
                    if(snpchr==prbchr && snpbp>=cisL && snpbp<=cisR)
                    {
                        double beta=eqtlinfo->_val[pos+j];
                        double se=eqtlinfo->_val[pos+j+num];
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=pxz_top)
                        {
                            esi_id_top=eqtlinfo->_rowid[pos+j];
                            epi_id_top=i;
                            beta_top=beta;
                            se_top=se;
                            pxz_top=pxz;
                        }
                    }
                    
                }
                if(pxz_top<=p_smr)
                {
                    newIcld.push_back(eqtlinfo->_include[i]);
                }

            }
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout<<eqtlinfo->_include.size()<<" probes with at least a cis-xQTL at p-value threshold of "<<p_smr<<" are extracted from ["+eqtlFileName+"]."<<endl;    
    }
    // multiple outcomes HEIDI test
    void multi_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, vector<gwasData> &gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        
        uint64_t probNum = esdata->_include.size();
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        MTSMRWK smrwk;
        long outcoNum=gdata.size();
        smrwk.byz.resize(outcoNum);smrwk.seyz.resize(outcoNum);smrwk.pyz.resize(outcoNum);
        //cout<<endl<<"Performing Multi-HEIDI tests..... "<<endl;
        float progr0=0.0 , progr1;
        //progress_print(progr0);
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            
            progr1=1.0*ii/probNum;
            if(progr1-progr0-0.05>1e-6 || ii+1==probNum)
            {
                if(ii+1==probNum) progr1=1.0;
                //progress_print(progr1);
                progr0=progr1;
            }
            int i=esdata->_include[ii];
            int probechr=esdata->_epi_chr[i];
            string probename=esdata->_epi_prbID[i];
            string probegene=esdata->_epi_gene[i];
            
            init_smr_wk_mlt(&smrwk);
            smrwk.cur_prbidx=i;
            smrwk.cur_chr=probechr;
            long maxid =fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, heidioffFlag);
            SMRRLT currlt;
             if(refSNP!=NULL && maxid==-9)
             {
                 //printf("WARNING: can't find target SNP %s for probe %s.\n",refSNP, probename.c_str());
                 currlt.p_HET=-9;
                 currlt.nsnp=-9;
                 smrrlts.push_back(currlt);
                 continue;
             }
             if (smrwk.bxz.size() == 0) {
                
                 //printf("WARNING: no SNP fetched for probe %s.\n", probename.c_str());
                 currlt.p_HET=-9;
                 currlt.nsnp=-9;
                 smrrlts.push_back(currlt);
                 continue;
             }
           
             Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
             Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
             vector<double> theta(outcoNum);
             if(sampleoverlap)
             {
                 printf("Estimating the correlation ...\n");
                 double z2mecs=qchisq(pmecs,1);
                 double zmecs=sqrt(z2mecs);
                 for(int t=0;t<outcoNum;t++)
                 {
                     vector<double> zxz, zyz;
                     for(int k=0;k<smrwk.bxz.size();k++)
                     {
                         
                         double z1=smrwk.bxz[k]/smrwk.sexz[k];
                         double z2=smrwk.byz[t][k]/smrwk.seyz[t][k];
                         if(abs(z1)<zmecs && abs(z2)<zmecs )
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
                         theta[t]=cor(zxz,zyz);
                         printf("The estimated correlation between exposure and outcome%d is %f.\n",t+1,theta[t]);
                     }
                 }
             }
             zsxz=ei_bxz.array()/ei_sexz.array();
            if(refSNP==NULL) {
             //if(opt) maxid=max_zsmr_id(&smrwk, p_smr);
             //else maxid=max_abs_id(zsxz);
             maxid=max_abs_id(zsxz);
            }
            if(maxid==-9) {
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
                 continue;
            }
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);

            if(refSNP==NULL && pxz_val>p_smr){
                //printf("WARNING: no SNP passed the p-value threshold %e for the SMR analysis for probe %s.\n", p_smr, probename.c_str());
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
            continue;
            } else {
            //printf("Analysing probe %s...\n", probename.c_str());
            }

            if(heidioffFlag)
            {
                //printf("INFO: the HEIDI test for probe %s is skipped because HEIDI test is turned off by the --heidi-off option or p_SMR does not pass the %e threshold.\n", probename.c_str(),threshpsmrest);
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
                
            }
            else
            {
                long nsnp=-9;
                double pdev=-9;
                if(!heidioffFlag) {
                    if(new_heidi_mth){
                        //if(refSNP!=NULL) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid, ld_min,opt_hetero,sampleoverlap, theta);
                        pdev= multi_heidi_test_new(bdata, &smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp, ld_min ,opt_hetero,sampleoverlap, theta);
                    } //else pdev= heidi_test(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                }
                    currlt.p_HET=pdev;
                    currlt.nsnp=(int)nsnp;
                    smrrlts.push_back(currlt);
             }
        }
        
    }
    // multiple exposures HEIDI test
    void multi_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        
        uint64_t probNum = esdata[0]._include.size();
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        MTSMRWKEXP smrwk;
        long expoNum = esdata.size();
        smrwk.bxz.resize(expoNum);smrwk.sexz.resize(expoNum);smrwk.zxz.resize(expoNum);smrwk.freq.resize(expoNum);
        float progr0=0.0 , progr1;
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            
            progr1=1.0*ii/probNum;
            if(progr1-progr0-0.05>1e-6 || ii+1==probNum)
            {
                if(ii+1==probNum) progr1=1.0;
                progr0=progr1;
            }
            int i=esdata[0]._include[ii];
            int probechr=esdata[0]._epi_chr[i];
            string probename=esdata[0]._epi_prbID[i];
            string probegene=esdata[0]._epi_gene[i];
            
            init_smr_wk_mlt(&smrwk);
            smrwk.cur_prbidx=i;
            smrwk.cur_chr=probechr;
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, heidioffFlag);
            
            SMRRLT currlt;
            if(refSNP!=NULL && maxid==-9)
            {
              currlt.p_HET=-9;
              currlt.nsnp=-9;
              smrrlts.push_back(currlt);
              continue;
            }
            if (smrwk.bxz[0].size() == 0) {
              currlt.p_HET=-9;
              currlt.nsnp=-9;
              smrrlts.push_back(currlt);
              continue;
            }
        
            Map<VectorXd> ei_bxz(&smrwk.bxz[0][0],smrwk.bxz[0].size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0][0],smrwk.sexz[0].size());
            vector<double> theta(expoNum);
            if(sampleoverlap)
            {
                printf("Estimating the correlation ...\n");
                double z2mecs=qchisq(pmecs,1);
                double zmecs=sqrt(z2mecs);
                for(int t=0;t<expoNum;t++)
                {
                    vector<double> zxz, zyz;
                    for(int k=0;k<smrwk.bxz.size();k++)
                    {                     
                        double z1=smrwk.bxz[t][k]/smrwk.sexz[t][k];
                        double z2=smrwk.byz[k]/smrwk.seyz[k];
                        if(abs(z1)<zmecs && abs(z2)<zmecs )
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
                        theta[t]=cor(zxz,zyz);
                        printf("The estimated correlation between exposure and outcome%d is %f.\n",t+1,theta[t]);
                    }
                }
            }
            zsxz=ei_bxz.array()/ei_sexz.array();
            if(refSNP==NULL) {
                if(opt) maxid=max_zsmr_id(&smrwk, p_smr);
                else maxid=max_abs_id(zsxz);
                maxid=max_abs_id(zsxz);
            }
            if(maxid==-9) {
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
                continue;
            }
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);

            if(refSNP==NULL && pxz_val>p_smr){
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
                continue;
            }

            if(heidioffFlag) {
                currlt.p_HET=-9;
                currlt.nsnp=-9;
                smrrlts.push_back(currlt);
            } else {
                long nsnp=-9;
                double pdev=-9;
                if(!heidioffFlag) {
                    if(new_heidi_mth){
                        //if(refSNP!=NULL) pdev= heidi_test_ref_new(bdata,&smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp,(int)maxid, ld_min,opt_hetero,sampleoverlap, theta);
                        pdev = multi_heidi_test_new_v2(bdata, &smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp, ld_min ,opt_hetero,sampleoverlap, theta);
                    }
                }
                currlt.p_HET=pdev;
                currlt.nsnp=(int)nsnp;
                smrrlts.push_back(currlt);
            }
        }
        
    }
    // multi-outcome HEIDI
    double multi_heidi_test_new(bInfo* bdata,MTSMRWK* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, vector<double> theta)
    {
        //the new method would calcualte maxid after each filtering
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        double pthres=pchisq(threshold,1);
        //printf("Filtering SNPs (%ld in total) at eQTL p-value < %e for the HEIDI test.\n",smrwk->zxz.size(), pthres);
        for(int i=0;i<smrwk->zxz.size();i++)
        {
            if(smrwk->zxz[i]*smrwk->zxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {
            
            //printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        //printf("%ld SNPs left after filtering.\n",sn_ids.size());
        update_snidx(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        MTSMRWK smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
        long maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        make_XMat(bdata,smrwk_heidi.curId, _X);
        //printf("Removing SNPs with LD r-squared between top-SNP %s > %f or < %f...\n",smrwk_heidi.rs[maxid_heidi].c_str(),ldr2_top,ld_min);
        ld_calc_o2m(ld_v,maxid_heidi,_X);
        if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
            sn_ids.clear();
            for(int i=0;i<smrwk_heidi.zxz.size();i++)
            {
                if(i!= maxid_heidi)
                {
                    double ldr2tmp=ld_v(i)*ld_v(i);
                    if((ldr2tmp<ldr2_top) && (ldr2tmp > ld_min)) sn_ids.push_back(i);
                }
                else{
                    sn_ids.push_back(i);
                }
            }
        }
        //printf("%ld SNPs are removed and %ld SNPs are retained.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
        if(sn_ids.size() < m_hetero) {
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        //printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
        int m = (int)smrwk_heidi.bxz.size();
        vector<int> rm_ID1;
        MatrixXd C;
        cor_calc(C, _X);
        double ld_top=sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
        //printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
        if(m-rm_ID1.size() < m_hetero) {
            
            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
            return -9;
        }
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sn_ids.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        update_snidx(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
        if (sn_ids.size() < C.size()) { //Build new matrix
            MatrixXd D(sn_ids.size(),sn_ids.size());
            for (int i = 0 ; i < sn_ids.size() ; i++) {
                for (int j = 0 ; j < sn_ids.size() ; j++) {
                    D(i,j) = C(sn_ids[i],sn_ids[j]);
                }
            }
            C = D;
        }
        VectorXd _bxz,_sexz,_zsxz;
        vector<VectorXd> _byz,_seyz;
        _byz.resize(smrwk_heidi.byz.size());
        _seyz.resize(smrwk_heidi.seyz.size());
        for(int t=0;t<smrwk_heidi.byz.size();t++)
        {
        _byz[t].resize(sn_ids.size());
        _seyz[t].resize(sn_ids.size());
        }
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
        //#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            for(int t=0;t<smrwk_heidi.byz.size();t++)
            {
            _byz[t][j]=smrwk_heidi.byz[t][sn_ids[j]];
            _seyz[t][j]=smrwk_heidi.seyz[t][sn_ids[j]];
            }
            _bxz[j]=smrwk_heidi.bxz[sn_ids[j]];
            _sexz[j]=smrwk_heidi.sexz[sn_ids[j]];
            _zsxz[j]=smrwk_heidi.zxz[sn_ids[j]];
            //cout<<_bxz[j]<<endl;
        }
        nsnp = sn_ids.size();
        double pdev=-9;
        //if(sampleoverlap) pdev=bxy_mltheter_so(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
        //else pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp);
        pdev=bxy_mltheter(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
        //printf("pHeidi is %e with %ld SNPs including in the HEIDI test.\n",pdev,nsnp);
        return pdev;
    }
    // multi-exposure HEIDI
    double multi_heidi_test_new(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp ,double ld_min,int opt_hetero, bool sampleoverlap, vector<double> theta)
    {
        //the new method would calcualte maxid after each filtering
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        //printf("Filtering SNPs (%ld in total) at eQTL p-value < %e for the HEIDI test.\n",smrwk->zxz.size(), pthres);
        vector<string> selectSNPs;
        vector<vector<string>> selectall(smrwk->zxz.size());
        for(int t=0; t<smrwk->zxz.size();t++)
        {
            selectSNPs.clear();
            multi_heidi_pruning(bdata, smrwk, t, ldr2_top, threshold,  m_hetero, ld_min, opt_hetero, selectSNPs);
            selectall[t]=selectSNPs;
        }
        
        vector<string> maxSNPs;
        for(int t=0; t<smrwk->zxz.size();t++)
        {
            long maxid=max_abs_id(smrwk->zxz[t]);
            maxSNPs.push_back(smrwk->rs[maxid]);
        }
        
        vector<int> cmmnID;
        vector<string> cmmnSNPs;
        vector<string> snptmp0;
        snptmp0.resize(selectall[0].size());
        for(int i=0;i<selectall[0].size();i++) snptmp0[i]=selectall[0][i];
        if(smrwk->zxz.size()>1) {
            for(int t=1; t<smrwk->zxz.size();t++)
            {
                cmmnID.clear();
                cmmnSNPs.clear();
                vector<string> snptmp(selectall[t].size());
                for(int i=0;i<selectall[t].size();i++)
                snptmp[i]=selectall[t][i];
                match_only(snptmp,snptmp0,cmmnID);
                if(cmmnID.empty())
                {
                    break;
                } else {
                cmmnSNPs.resize(cmmnID.size());
                for(int i=0;i<cmmnID.size();i++)
                    cmmnSNPs[i]=snptmp0[cmmnID[i]];
                snptmp0.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                    snptmp0[i]=cmmnSNPs[i];
                }
            }
        } else {
            cmmnSNPs.resize(snptmp0.size());
            for(int i=0;i<snptmp0.size();i++)
                cmmnSNPs[i]=snptmp0[i];
        }
        
        if(cmmnSNPs.size() < m_hetero) {
            
            //printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            return -9;
        }
        for(int t=0; t<maxSNPs.size();t++) cmmnSNPs.push_back(maxSNPs[t]);
        
        match_only(cmmnSNPs,smrwk->rs,sn_ids);
        CommFunc::getUnique(sn_ids);
//        vector<int> tmp_ids;
//        match_only(maxSNPs,cmmnSNPs,tmp_ids);
        
//        for(int i=0;i<smrwk->zxz[0].size();i++)
//        {
//            if(smrwk->zxz[0][i]*smrwk->zxz[0][i]-threshold>1e-6) sn_ids.push_back(i);
//        }
     
        //printf("%ld SNPs left after filtering.\n",sn_ids.size());
//        update_snidx(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
//        long maxid_heidi=max_abs_id(smrwk_heidi.zxz[0]);
        make_XMat(bdata,smrwk_heidi.curId, _X);
//        //printf("Removing SNPs with LD r-squared between top-SNP %s > %f or < %f...\n",smrwk_heidi.rs[maxid_heidi].c_str(),ldr2_top,ld_min);
//        ld_calc_o2m(ld_v,maxid_heidi,_X);
//
//        if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
//            sn_ids.clear();
//            for(int i=0;i<smrwk_heidi.zxz[0].size();i++)
//            {
//                if(i!= maxid_heidi)
//                {
//                    double ldr2tmp=ld_v(i)*ld_v(i);
//                    if((ldr2tmp<ldr2_top) && (ldr2tmp > ld_min)) sn_ids.push_back(i);
//                }
//                else{
//                    sn_ids.push_back(i);
//                }
//            }
//        }
//        //printf("%ld SNPs are removed and %ld SNPs are retained.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
//        if(sn_ids.size() < m_hetero) {
//            //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
//            return -9;
//        }
       for(int i=0;i<sn_ids.size();i++) sn_ids[i]=i;
       update_smrwk_x(&smrwk_heidi,sn_ids,_X);
//       maxid_heidi=max_abs_id(smrwk_heidi.zxz[0]);
       //printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
//       int m = (int)smrwk_heidi.bxz[0].size();
//       vector<int> rm_ID1;
       MatrixXd C;
       cor_calc(C, _X);
//       double ld_top=sqrt(ldr2_top);
//       if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
//       //printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
//       if(m-rm_ID1.size() < m_hetero) {
//
//           //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
//           return -9;
//       }
//       //Create new index
//       sn_ids.clear();
//       int qi=0;
//       for (int i=0 ; i<m ; i++) {
//           if (rm_ID1.size() == 0) sn_ids.push_back(i);
//           else {
//               if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
//               else sn_ids.push_back(i);
//           }
//       }
//       update_snidx(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
       
        if (sn_ids.size() < C.size()) { //Build new matrix
           MatrixXd D(sn_ids.size(),sn_ids.size());
           for (int i = 0 ; i < sn_ids.size() ; i++) {
               for (int j = 0 ; j < sn_ids.size() ; j++) {
                   D(i,j) = C(sn_ids[i],sn_ids[j]);
               }
           }
           C = D;
        }
        vector<VectorXd> _bxz,_sexz,_zsxz;
        VectorXd  _byz,_seyz;
        _bxz.resize(smrwk_heidi.bxz.size());
        _sexz.resize(smrwk_heidi.sexz.size());
        _zsxz.resize(smrwk_heidi.zxz.size());
        for(int t=0;t<smrwk_heidi.bxz.size();t++)
        {
            _bxz[t].resize(sn_ids.size());
            _sexz[t].resize(sn_ids.size());
            _zsxz[t].resize(sn_ids.size());
        }
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
       #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            for(int t=0;t<smrwk_heidi.bxz.size();t++)
            {
            _bxz[t][j]=smrwk_heidi.bxz[t][sn_ids[j]];
            _sexz[t][j]=smrwk_heidi.sexz[t][sn_ids[j]];
            _zsxz[t][j]=smrwk_heidi.zxz[t][sn_ids[j]];
            }
            _byz[j]=smrwk_heidi.byz[sn_ids[j]];
            _seyz[j]=smrwk_heidi.seyz[sn_ids[j]];
            //cout<<_bxz[j]<<endl;
        }
       nsnp = sn_ids.size();
       double pdev=-9;
       //if(sampleoverlap) pdev=bxy_mltheter_so(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
       //else pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp);
       pdev=bxy_mltheter(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
       //printf("pHeidi is %e with %ld SNPs including in the HEIDI test.\n",pdev,nsnp);
       return pdev;
    }

    // multi-exposure HEIDI with new pruning
    double multi_heidi_test_new_v2(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, double ld_min, int opt_hetero, bool sampleoverlap, vector<double> theta)
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;        
        vector<string> slctsnps;

        // topSNP and pairwise LD pruning
        multi_heidi_pruning_v2(bdata, smrwk, ldr2_top, threshold,  m_hetero, ld_min, opt_hetero, slctsnps);
        if(slctsnps.size() < m_hetero) {
            return -9;
        }
        match_only(slctsnps,smrwk->rs,sn_ids);
        CommFunc::getUnique(sn_ids);

        MTSMRWKEXP smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
        make_XMat(bdata,smrwk_heidi.curId, _X);
        for(int i=0;i<sn_ids.size();i++) sn_ids[i]=i;
        //update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        MatrixXd C;
        cor_calc(C, _X);
       
       // if (sn_ids.size() < C.size()) { //Build new matrix
       //    MatrixXd D(sn_ids.size(),sn_ids.size());
       //    for (int i = 0 ; i < sn_ids.size() ; i++) {
       //        for (int j = 0 ; j < sn_ids.size() ; j++) {
       //            D(i,j) = C(sn_ids[i],sn_ids[j]);
       //        }
       //    }
       //    C = D;
       // }

        vector<VectorXd> _bxz,_sexz,_zsxz;
        VectorXd  _byz,_seyz;
        _bxz.resize(smrwk_heidi.bxz.size());
        _sexz.resize(smrwk_heidi.sexz.size());
        _zsxz.resize(smrwk_heidi.zxz.size());
        for(int t=0;t<smrwk_heidi.bxz.size();t++)
        {
            _bxz[t].resize(sn_ids.size());
            _sexz[t].resize(sn_ids.size());
            _zsxz[t].resize(sn_ids.size());
        }
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());

        #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            for(int t=0;t<smrwk_heidi.bxz.size();t++)
            {
            _bxz[t][j]=smrwk_heidi.bxz[t][sn_ids[j]];
            _sexz[t][j]=smrwk_heidi.sexz[t][sn_ids[j]];
            _zsxz[t][j]=smrwk_heidi.zxz[t][sn_ids[j]];
            //cout<<_bxz[t][j]<<endl;
            }
            _byz[j]=smrwk_heidi.byz[sn_ids[j]];
            _seyz[j]=smrwk_heidi.seyz[sn_ids[j]];
        }
        nsnp = sn_ids.size();
        double pdev=-9;
        pdev = bxy_mltheter(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
        return pdev;
    }

    void multi_heidi_pruning_v2(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero ,double ld_min,int opt_hetero,vector<string> &slctsnps)
    {
        int expoNum = smrwk->zxz.size();
        MatrixXd _X;
        vector<int> sn_ids, rmsn_ids;
        MTSMRWKEXP* smrwk_tmp = smrwk;
        
        // 1. HEIDI instrument at p-value < 1e-3;
        for(int i=0;i<smrwk_tmp->zxz[0].size();i++)
        {
            int inldNum = 0;
            for(int t=0; t<expoNum; t++) {
                if(smrwk_tmp->zxz[t][i]*smrwk_tmp->zxz[t][i]-threshold>1e-6) {
                    inldNum = inldNum + 1;
                }
            }
            if(inldNum == expoNum) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {
            slctsnps.push_back("NA");
            return;
        }
        // insert the maxids
        for(int t=0; t<expoNum; t++) {
            long maxid = max_abs_id(smrwk_tmp->zxz[t]);
            sn_ids.push_back(maxid);
        }
        CommFunc::getUnique(sn_ids);

        // 2. LD pruning with the top SNPs;
        // update_snidx(smrwk_tmp,t,sn_ids,MAX_NUM_LD,"LD pruning");
        sort(sn_ids.begin(),sn_ids.end());
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk(smrwk_tmp,sn_ids,&smrwk_heidi);
        for(int i=0;i<smrwk_heidi.byz.size();i++) sn_ids[i] = i;
        for(int t=0; t<expoNum; t++) {
            VectorXd ld_v; MatrixXd _Xtmp;
            long maxid = max_abs_id(smrwk_heidi.zxz[t]);
            // LD r-squared between top-SNP and others;
            make_XMat(bdata,smrwk_heidi.curId, _Xtmp);
            ld_calc_o2m(ld_v,maxid,_Xtmp);
            
            if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
                for(int i=0;i<smrwk_heidi.zxz[t].size();i++)
                {
                    if(i!= maxid)
                    {
                        double ldr2tmp=ld_v(i)*ld_v(i);
                        if((ldr2tmp > ldr2_top) || (ldr2tmp < ld_min)) rmsn_ids.push_back(i);
                    }
                }
            }
        }
        CommFunc::getUnique(rmsn_ids);
        for(int i=0;i<rmsn_ids.size();i++) sn_ids.erase(std::remove(sn_ids.begin(), sn_ids.end(), rmsn_ids[i]), sn_ids.end());

        if(sn_ids.size() < m_hetero) {
            slctsnps.push_back("NA");
            return;
        }

        // insert the maxids
        for(int t=0; t<expoNum; t++) {
            long maxid = max_abs_id(smrwk_heidi.zxz[t]);
            sn_ids.push_back(maxid);
        }
        CommFunc::getUnique(sn_ids);

        // 3. pairwise LD pruning
        make_XMat(bdata,smrwk_heidi.curId, _X);
        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        int m = (int)smrwk_heidi.byz.size();
        vector<int> rm_ID;
        MatrixXd C;
        cor_calc(C, _X);
        double ld_top = sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID);

        if(m-rm_ID.size() < m_hetero) {
           slctsnps.push_back("NA");
           return;
        }
        
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0; i<m; i++) {
            if (rm_ID.size() == 0) sn_ids.push_back(i);
            else {
                if (qi<rm_ID.size() && rm_ID[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        // insert the maxids
        for(int t=0; t<expoNum; t++) {
            long maxid = max_abs_id(smrwk_heidi.zxz[t]);
            sn_ids.push_back(maxid);
        }
        CommFunc::getUnique(sn_ids);
        
        // return select SNPs
        slctsnps.resize(sn_ids.size());
        for(int i=0;i<sn_ids.size();i++)
        {
            slctsnps[i]=smrwk_heidi.rs[sn_ids[i]];
        }
    }

    void multi_heidi_pruning(bInfo* bdata,MTSMRWKEXP* smrwk,int t,double ldr2_top, double threshold, int m_hetero ,double ld_min,int opt_hetero,vector<string> &selectSNPs)
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        MTSMRWKEXP* smrwk_tmp = smrwk;
        //printf("Filtering SNPs (%ld in total) at eQTL p-value < %e for the HEIDI test.\n",smrwk->zxz.size(), pthres);
        for(int i=0;i<smrwk_tmp->zxz[t].size();i++)
        {
            if(smrwk_tmp->zxz[t][i]*smrwk_tmp->zxz[t][i]-threshold>1e-6) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {           
            //printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            selectSNPs.push_back("NA");
            return;
        }
        //printf("%ld SNPs left after filtering.\n",sn_ids.size());
        update_snidx(smrwk_tmp,t,sn_ids,MAX_NUM_LD,"LD pruning");
        sort(sn_ids.begin(),sn_ids.end());
        
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk(smrwk_tmp,sn_ids,&smrwk_heidi);
        long maxid_heidi = max_abs_id(smrwk_heidi.zxz[t]);
        make_XMat(bdata,smrwk_heidi.curId, _X);
        //printf("Removing SNPs with LD r-squared between top-SNP %s > %f or < %f...\n",smrwk_heidi.rs[maxid_heidi].c_str(),ldr2_top,ld_min);
        ld_calc_o2m(ld_v,maxid_heidi,_X);
        
        if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
            sn_ids.clear();
            for(int i=0;i<smrwk_heidi.zxz[t].size();i++)
            {
                if(i!= maxid_heidi)
                {
                    double ldr2tmp=ld_v(i)*ld_v(i);
                    if((ldr2tmp<ldr2_top) && (ldr2tmp > ld_min)) sn_ids.push_back(i);
                }
                else{
                    sn_ids.push_back(i);
                }
            }
        }

        if(sn_ids.size() < m_hetero) {           
            //printf("INFO: the HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
            selectSNPs.push_back("NA");
            return;
        }

        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        //maxid_heidi=max_abs_id(smrwk_heidi.zxz[t]);
        //printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
        int m = (int)smrwk_heidi.bxz[t].size();
        vector<int> rm_ID1;
        MatrixXd C;
        cor_calc(C, _X);
        double ld_top=sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);

        if(m-rm_ID1.size() < m_hetero) {
           //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
           selectSNPs.push_back("NA");
           return;
        }
        //printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sn_ids.push_back(i);
            else {
                if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        
//        update_snidx(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
        
        selectSNPs.resize(sn_ids.size());
        for(int i=0;i<sn_ids.size();i++)
        {   
            selectSNPs[i]=smrwk_heidi.rs[sn_ids[i]];
        }
    
    }

    void est_cov_dxy(MatrixXd &covdxy,VectorXd &vardev,int maxid, vector<VectorXd> &_byz, VectorXd &_bxz, vector<VectorXd> &_seyz,VectorXd &_sexz, MatrixXd &_LD_heidi, vector<double> theta)
    {
        long nsnp =_bxz.size();
        long outcoNum = _byz.size();
        vardev.resize((nsnp-1)*outcoNum);
        if(nsnp>1)
        {
            VectorXd ivect = VectorXd::Ones(nsnp);
            MatrixXd bexpo=_bxz*_bxz.transpose();
            MatrixXd seexpo=_sexz*_sexz.transpose();
            VectorXd _bxz2=_bxz.array()*_bxz.array();
            for(int k=0;k<outcoNum;k++)
            {
                for(int j=0;j<outcoNum;j++)
                {
                    MatrixXd covbxy(nsnp,nsnp);
                    MatrixXd vdev((nsnp-1),(nsnp-1));
                    VectorXd tmp3;
                    MatrixXd cov1 = _LD_heidi.array()*seexpo.array()*(_byz[k]*_byz[j].transpose()).array()/(bexpo.array()*bexpo.array());
                    MatrixXd cov2 = _LD_heidi.array()*(_seyz[k]*_seyz[j].transpose()).array()/bexpo.array();
                    if(j==k) {covbxy = cov1 + cov2;
                    } else {covbxy = cov1;}
                    double tmp1=covbxy(maxid,maxid);
                    tmp3.resize(nsnp-1);
                    for(int i=0; i<maxid; i++) tmp3[i]=covbxy(maxid,i);
                    for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=covbxy(maxid,i);
                    // vdev as tmp2
                    vdev.block(0, 0, maxid, maxid) = covbxy.block(0, 0, maxid, maxid);
                    vdev.block(0, maxid, maxid, nsnp - maxid - 1) = covbxy.block(0, maxid + 1, maxid, nsnp - maxid - 1);
                    vdev.block(maxid, 0, nsnp - maxid - 1, maxid) = covbxy.block(maxid + 1, 0, nsnp - maxid - 1, maxid);
                    vdev.block(maxid, maxid, nsnp - maxid - 1, nsnp - maxid - 1) = covbxy.block(maxid + 1, maxid + 1, nsnp - maxid - 1, nsnp - maxid - 1);
                    
                    // get vdev
                    VectorXd v1 = VectorXd::Zero(nsnp - 1);
                    v1 = v1.array() + 1.0;
                    
                    vdev = tmp1 + vdev.array() - (v1*tmp3.transpose()).array() - (tmp3*v1.transpose()).array();
                    for (int i = 0; i<nsnp - 1; i++)  vdev(i,i) += 1e-8; // in R code
                    //tmp3 as vardev
                    for(int i=0; i<maxid; i++) tmp3[i]=tmp1+covbxy(i,i)-2*tmp3[i]+ 1e-8;
                    for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy(i,i)-2*tmp3[i-1]+ 1e-8;
                    if(j==k)
                    {
                        for(int m=0;m<(nsnp-1);m++) vardev[m+j*(nsnp-1)]=tmp3[m];
                    }
                    long rowid1=(k)*(nsnp-1);
                    long colid1=(j)*(nsnp-1);
                    covdxy.block(rowid1,colid1,nsnp-1,nsnp-1)=vdev.block(0,0,nsnp-1,nsnp-1);
                }
            }
        }
    }
     void est_cov_dxy(MatrixXd &covdxy,VectorXd &vardev,vector<int> maxid, VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz,vector<VectorXd> &_sexz, MatrixXd &_LD_heidi, vector<double> theta)
    {
        long nsnp =_byz.size();
        long outcoNum = _bxz.size();
        vardev.resize((nsnp-1)*outcoNum);
        if(nsnp>1)
        {
            // VectorXd ivect = VectorXd::Ones(nsnp);
            MatrixXd boutco=_byz*_byz.transpose();
            MatrixXd seoutco=_seyz*_seyz.transpose();
            for(int k=0;k<outcoNum;k++)
            {
                for(int j=0;j<outcoNum;j++)
                {
                    MatrixXd covbxy(nsnp,nsnp);
                    MatrixXd vdev((nsnp-1),(nsnp-1));
                    VectorXd tmp3, tmp4;
                    MatrixXd cov1 = _LD_heidi.array()*(_sexz[k]*_sexz[j].transpose()).array()*boutco.array()/((_bxz[k]*_bxz[j].transpose()).array()*(_bxz[k]*_bxz[j].transpose()).array()).array();
                    MatrixXd cov2 = _LD_heidi.array()*seoutco.array()/(_bxz[k]*_bxz[j].transpose()).array();
                    if(j==k) {covbxy = cov1 + cov2;
                    } else {covbxy = cov2;}
                    
                    double tmp1=covbxy(maxid[k],maxid[j]);
                    tmp3.resize(nsnp-1);
                    // maxid[k] row
                    for(int i=0; i<maxid[j]; i++) tmp3[i]=covbxy(maxid[k],i);
                    for(int i=maxid[j]+1; i<nsnp; i++) tmp3[i-1]=covbxy(maxid[k],i);
                    tmp4.resize(nsnp-1);
                    // maxid[j] column
                    for(int i=0; i<maxid[k]; i++) tmp4[i]=covbxy(i,maxid[j]);
                    for(int i=maxid[k]+1; i<nsnp; i++) tmp4[i-1]=covbxy(i,maxid[j]);
                    // vdev as tmp2
                    vdev.block(0, 0, maxid[k], maxid[j]) = covbxy.block(0, 0, maxid[k], maxid[j]);
                    vdev.block(0, maxid[j], maxid[k], nsnp - maxid[j] - 1) = covbxy.block(0, maxid[j] + 1, maxid[k], nsnp - maxid[j] - 1);
                    vdev.block(maxid[k], 0, nsnp - maxid[k] - 1, maxid[j]) = covbxy.block(maxid[k] + 1, 0, nsnp - maxid[k] - 1, maxid[j]);
                    vdev.block(maxid[k], maxid[j], nsnp - maxid[k] - 1, nsnp - maxid[j] - 1) = covbxy.block(maxid[k] + 1, maxid[j] + 1, nsnp - maxid[k] - 1, nsnp - maxid[j] - 1);
                    // get vdev
                    VectorXd v1 = VectorXd::Zero(nsnp - 1);
                    v1 = v1.array() + 1.0;
                    
                    vdev = tmp1 + vdev.array() - (v1 * tmp3.transpose()).array() - (tmp4 * v1.transpose()).array();
                    
                    for (int i = 0; i<nsnp - 1; i++)  vdev(i,i) += 1e-8; // in R code
                    //tmp3 as vardev
                    for(int i=0; i<maxid[j]; i++) tmp3[i]=tmp1+covbxy(i,i)-2*tmp3[i]+ 1e-8;
                    for(int i=maxid[j]+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy(i,i)-2*tmp3[i-1]+ 1e-8;
                    if(j==k)
                    {
                        for(int m=0;m<(nsnp-1);m++) vardev[m+j*(nsnp-1)]=tmp3[m];
                    }
                    long rowid1=(k)*(nsnp-1);
                    long colid1=(j)*(nsnp-1);
                    covdxy.block(rowid1,colid1,nsnp-1,nsnp-1)=vdev.block(0,0,nsnp-1,nsnp-1);
                }
            }
        }
    }
    
    float bxy_mltheter(vector<VectorXd> &_byz, VectorXd &_bxz, vector<VectorXd> &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num, vector<double> theta)
    {
        vector<VectorXd> _bxy, _sexy;
        VectorXd dev,chisqdev;
        long nsnp=*snp_num;
        int maxid;
        float pdev=-1.0;
        long outcoNum=_byz.size();
        _bxy.resize(outcoNum),_sexy.resize(outcoNum);
        MatrixXd covbxy((nsnp-1)*outcoNum,(nsnp-1)*outcoNum);
        MatrixXd covdxy((nsnp-1)*outcoNum,(nsnp-1)*outcoNum);
        VectorXd vardev;
        VectorXd bxz2;
        dev.resize((nsnp-1)*outcoNum);
        chisqdev.resize((nsnp-1)*outcoNum);
        if(nsnp>1)
        {
            bxz2=_bxz.array()*_bxz.array();
            for(int t=0; t<outcoNum;t++){
                _bxy[t]=_byz[t].array()/_bxz.array();
                _sexy[t]=((_seyz[t].array()*_seyz[t].array()*bxz2.array()+_sexz.array()*_sexz.array()*_byz[t].array()*_byz[t].array())/(bxz2.array()*bxz2.array())).sqrt();
            }
            maxid=max_abs_id(_zsxz);
            
            for(int t=0; t<outcoNum;t++){
            for(int j=0;j<maxid;j++) dev[t*(nsnp-1)+j]=_bxy[t][maxid]-_bxy[t][j];
            for(int j=maxid+1;j<nsnp;j++) dev[t*(nsnp-1)+j-1]=_bxy[t][maxid]-_bxy[t][j];
            }
            est_cov_dxy(covdxy,vardev,maxid, _byz,_bxz,_seyz,_sexz,_LD_heidi,theta);

            //dev as chisq_dev
            for(int i=0;i<(nsnp-1)*outcoNum;i++) chisqdev[i]=dev[i]*dev[i]/vardev[i];
           
            double sumChisq_dev=0.0; double sumdev=0.0, sumvar=0.0;
            for(int i=0;i<(nsnp-1)*outcoNum;i++)
            {
                //cout<<chisqdev[i]<<endl;
                sumChisq_dev+=chisqdev[i];
                sumdev+=dev[i];
                sumvar+=vardev[i];
            }
            
            //using covbxy to store corr_dev
            covbxy.resize((nsnp - 1)*outcoNum, (nsnp - 1)*outcoNum);
            covbxy = covdxy.array() / sqrt((covdxy.diagonal()*covdxy.diagonal().transpose()).array());
            covbxy.triangularView<Lower>() = covbxy.transpose(); // they make the matrix symmetrical in the Rscript
           
            // using Eigen Library
            SelfAdjointEigenSolver<MatrixXd> es(covbxy);
            VectorXd lambda;
            lambda=es.eigenvalues();
            
            /*
             EigenSolver<MatrixXd> es(A);
             cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
             cout<< es.eigenvalues().transpose()<<endl;
             */
            
            /*
             MatrixXd D = es.pseudoEigenvalueMatrix();
             cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
             */
            
            pdev= pchisqsum(sumChisq_dev,lambda);
            //*snp_num=lambda.size();
            *snp_num = nsnp - 1;
            
        }else *snp_num=-9;      
        
        return(pdev);
    }
    float bxy_mltheter(VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz, vector<VectorXd> &_sexz, vector<VectorXd> &_zsxz, MatrixXd &_LD_heidi, long* snp_num, vector<double> theta)
    {
        vector<VectorXd> _bxy, _sexy, bxz2;
        VectorXd dev,chisqdev;
        long nsnp=*snp_num;
        float pdev=-1.0;
        long expoNum=_bxz.size();
        vector<int> maxid(expoNum);
        _bxy.resize(expoNum),_sexy.resize(expoNum),bxz2.resize(expoNum);
        MatrixXd covbxy((nsnp-1)*expoNum,(nsnp-1)*expoNum);
        MatrixXd covdxy((nsnp-1)*expoNum,(nsnp-1)*expoNum);
        VectorXd vardev;
        dev.resize((nsnp-1)*expoNum);
        chisqdev.resize((nsnp-1)*expoNum);
        if(nsnp>1)
        {
            for(int t=0; t<expoNum; t++){
                maxid[t]=max_abs_id(_zsxz[t]);
                bxz2[t]=_bxz[t].array()*_bxz[t].array();
            }
            for(int t=0; t<expoNum; t++){
                _bxy[t]=_byz.array()/_bxz[t].array();
                _sexy[t]=((_seyz.array()*_seyz.array()*bxz2[t].array()+_sexz[t].array()*_sexz[t].array()*_byz.array()*_byz.array())/(bxz2[t].array()*bxz2[t].array())).sqrt();
            }
            for(int t=0; t<expoNum; t++){
                for(int j=0;j<maxid[t];j++) dev[t*(nsnp-1)+j]=_bxy[t][maxid[t]]-_bxy[t][j];
                for(int j=maxid[t]+1;j<nsnp;j++) dev[t*(nsnp-1)+j-1]=_bxy[t][maxid[t]]-_bxy[t][j];
            }
            
            est_cov_dxy(covdxy,vardev,maxid,_byz,_bxz,_seyz,_sexz,_LD_heidi,theta);
            
            //dev as chisq_dev
            for(int i=0;i<(nsnp-1)*expoNum;i++) chisqdev[i]=dev[i]*dev[i]/vardev[i];
            
            double sumChisq_dev=0.0; double sumdev=0.0, sumvar=0.0;
            for(int i=0;i<(nsnp-1)*expoNum;i++)
            {
                //cout<<chisqdev[i]<<endl;
                sumChisq_dev+=chisqdev[i];
                sumdev+=dev[i];
                sumvar+=vardev[i];
            }
            
            //using covbxy to store corr_dev
            covbxy.resize((nsnp - 1)*expoNum, (nsnp - 1)*expoNum);
            covbxy = covdxy.array() / sqrt((covdxy.diagonal()*covdxy.diagonal().transpose()).array());
            covbxy.triangularView<Lower>() = covbxy.transpose(); // they make the matrix symmetrical in the Rscript
           
            // using Eigen Library
            SelfAdjointEigenSolver<MatrixXd> eigensolver(covbxy);
            VectorXd lambda;
            lambda = eigensolver.eigenvalues();
            
            /*
             EigenSolver<MatrixXd> es(A);
             cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
             cout<< es.eigenvalues().transpose()<<endl;
             */
            
            /*
             MatrixXd D = es.pseudoEigenvalueMatrix();
             cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
             */
            
            pdev= pchisqsum(sumChisq_dev,lambda);
            //*snp_num=lambda.size();
            *snp_num = nsnp - 1;
            
        }else *snp_num=-9;      
        
        return(pdev);
    }
    void allele_compare(string a1, string a2, string s1,string s2, int &Id, int &flip)
    {
        if(a1 == s1 && a2 == s2)
        {
            Id = 1; flip = 0 ;
        }
        else if(a1 == s2 && a2 == s1)
        {
            Id = 1; flip = 1 ;
        }
        else 
        {
            Id = 0; flip = 0 ;  
        }
    }
    void allele_check_multi(vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId1;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among GWAS summary, exposure summary and outcome summary data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> gsnp(gdata->_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum ){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                #pragma omp parallel for
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with multiple dataset //
        if(esdata->_esi_include.size()<esdata->_snpNum || gdata->_include.size()< gdata->snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[gdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(gsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
            //cout<<slctSNPs.size()<<endl;
        }else
        {
            match_only(gdata->snpName, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];  
            //cout<<slctSNPs.size()<<endl;
        }
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        //alleles check
        match(slctSNPs, gdata->snpName, gdId1);
        //cout<<gdId1.size()<<endl;
        //match(slctSNPs, bdata->_snp_name, bdId);
        //match(slctSNPs, etrait->_esi_rs, gdId);
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        gdata->_include.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();    
        }
        //etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {   
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string ta1, ta2, ea1, ea2, ga1, ga2;
            ga1 = gdata->allele_1[gdId1[i]];
            ga2 = gdata->allele_2[gdId1[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            //ta1 = etrait->_esi_allele1[gdId[i]];
            //ta2 = etrait->_esi_allele2[gdId[i]];
            int Id3,flip3;
            //allele_compare(ea1,ea2,a1,a2,Id1,flip1);
            //allele_compare(ea1,ea2,ta1,ta2,Id2,flip2);
            allele_compare(ea1,ea2,ga1,ga2,Id3,flip3);
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=1,Idsrt=1; //start with all Id = 1;
            for(int t = 0; t < etraitNum; t++)
            {   
                //gdId.clear();
                //match(slctSNPs, etrait[t]._esi_rs, gdId);
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ea1,ea2,ta1,ta2,Id[t],flip[t]);
                Idall=Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            #pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                esdata->_esi_include.push_back(edId[i]);
                gdata->_include.push_back(gdId1[i]);
                if(flip3==1) {gdata->byz[gdId1[i]]=-gdata->byz[gdId1[i]];}
            }
        }
        //cout<<etrait[t]->_esi_include.size()<<endl;
        //cout<<bdata->_include.size()<<endl;
        //cout<<gdata->_include.size()<<endl;
        logstr=itos(gdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
    }
    void allele_check_multi(vector<eqtlInfo> &etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among exposure summary and outcome summary data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> essnp(esdata->_esi_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum ){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with multiple dataset //
        if(esdata->_esi_include.size()<esdata->_snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(etsnp, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }else
        {
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];  
        }
    
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();    
        }
        esdata->_esi_include.clear();
        // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {   
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string ta1, ta2, ea1, ea2;
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=0,Idsrt=0; //start with all Id = 0;
            for(int t = 0; t < etraitNum; t++)
            {   
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ea1,ea2,ta1,ta2,Id[t],flip[t]);
                Idall=Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            #pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                esdata->_esi_include.push_back(edId[i]);
            }
        }
        logstr=itos(esdata->_esi_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
    }
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId1;
        //vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among GWAS summary, exposure summary, outcome summary and LD reference data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> gsnp(gdata->_include.size());
        vector<string> bsnp(bdata->_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum ){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with multiple dataset //
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum || gdata->_include.size()< gdata->snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[gdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(bsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            //edId.clear();
            match_only(cmmnSNPs,gsnp,edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gsnp[edId[i]];
            //edId.clear();
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
            //cout<<slctSNPs.size()<<endl;
        }else
        {
            match_only(bdata->_snp_name, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gdata->snpName[edId[i]];
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];  
            //cout<<slctSNPs.size()<<endl;
        }
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        //alleles check
        match(slctSNPs, gdata->snpName, gdId1);
        //cout<<gdId1.size()<<endl;
        match(slctSNPs, bdata->_snp_name, bdId);
        //match(slctSNPs, etrait->_esi_rs, gdId);
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();    
        }
        //etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string a1, a2, ta1, ta2, ea1, ea2, ga1, ga2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId1[i]];
            ga2 = gdata->allele_2[gdId1[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            //ta1 = etrait->_esi_allele1[gdId[i]];
            //ta2 = etrait->_esi_allele2[gdId[i]];
            int Id1,Id3,flip1,flip3;
            allele_compare(ea1,ea2,a1,a2,Id1,flip1);
            //allele_compare(ea1,ea2,ta1,ta2,Id2,flip2);
            allele_compare(ea1,ea2,ga1,ga2,Id3,flip3);
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=2,Idsrt=2; //start with all Id = 1;
            for(int t = 0; t < etraitNum; t++)
            {   
                //gdId.clear();
                //match(slctSNPs, etrait[t]._esi_rs, gdId);
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ea1,ea2,ta1,ta2,Id[t],flip[t]);
                Idall=Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            #pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                bdata->_include.push_back(bdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                gdata->_include.push_back(gdId1[i]);
                if(flip1==1){
                    //cout<<"test1"<<endl;
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                }
                if(flip3==1) {gdata->byz[gdId1[i]]=-gdata->byz[gdId1[i]];}
            }
        }
        //cout<<etrait[t]->_esi_include.size()<<endl;
        //cout<<bdata->_include.size()<<endl;
        //cout<<gdata->_include.size()<<endl;
        logstr=atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId1;
        //vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among GWAS summary, exposure summary and LD reference data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> gsnp(gdata->_include.size());
        vector<string> bsnp(bdata->_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with multiple dataset //
        if(bdata->_include.size()< bdata->_snp_num || gdata->_include.size()< gdata->snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[gdata->_include[i]];
            match_only(bsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between plink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            //edId.clear();
            match_only(cmmnSNPs,gsnp,edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gsnp[edId[i]];
            //edId.clear();
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=cmmnSNPs[i];
            //cout<<slctSNPs.size()<<endl;
        } else
        {
            match_only(bdata->_snp_name, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between plink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gdata->snpName[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=cmmnSNPs[i];  
            //cout<<slctSNPs.size()<<endl;
        }
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        //alleles check
        match(slctSNPs, gdata->snpName, gdId1);
        //cout<<gdId1.size()<<endl;
        match(slctSNPs, bdata->_snp_name, bdId);
        //match(slctSNPs, etrait->_esi_rs, gdId);
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();    
        }
        //etrait->_esi_include.clear();
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string a1, a2, ta1, ta2, ea1, ea2, ga1, ga2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId1[i]];
            ga2 = gdata->allele_2[gdId1[i]];
            // ea1 = esdata->_esi_allele1[edId[i]];
            // ea2 = esdata->_esi_allele2[edId[i]];
            //ta1 = etrait->_esi_allele1[gdId[i]];
            //ta2 = etrait->_esi_allele2[gdId[i]];
            int Id1,flip1;
            allele_compare(ga1,ga2,a1,a2,Id1,flip1);
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=1,Idsrt=1; //start with all Id = 1;
            for(int t = 0; t < etraitNum; t++)
            {   
                //gdId.clear();
                //match(slctSNPs, etrait[t]._esi_rs, gdId);
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ga1,ga2,ta1,ta2,Id[t],flip[t]);
                Idall=Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            #pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                bdata->_include.push_back(bdId[i]);
                gdata->_include.push_back(gdId1[i]);
                if(flip1==1){
                    //cout<<"test1"<<endl;
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                }
                //if(flip3==1) {gdata->byz[gdId1[i]]=-gdata->byz[gdId1[i]];}
            }
        }
        //cout<<etrait[t]->_esi_include.size()<<endl;
        //cout<<bdata->_include.size()<<endl;
        //cout<<gdata->_include.size()<<endl;
        logstr=atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    void allele_check_multi_opt_old(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        string logstr="Checking consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the xQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        int etraitNum = etrait.size(), cmpnum = etraitNum + 1;
        map<string, int>::iterator iter, iter1, iter2;
        vector<int> bin,gin,eintmp(etraitNum);
        vector<vector<int>> ein(etraitNum);
        // use 1st etrait as reference; [0]:etrait[0]vsGWAS; [1:t]:etrait[0]:etrait[1:t]; [t+1]:etrait[0]:LDref; 
        vector<int> incld(cmpnum), flip(cmpnum);
        double disp=0;
        for(int i=0;i<gdata->_include.size();i++) //no map in gwas data
        {
            progress(i, disp, (int)gdata->_include.size());      
            int rid = gdata->_include[i];
            string grs = gdata->snpName[rid];
            string ba1, ba2, ga1, ga2, ea1, ea2, ea1tmp, ea2tmp;
            ga1 = gdata->allele_1[rid];
            ga2 = gdata->allele_2[rid];
            iter1 = bdata-> _snp_name_map.find(grs);
            if (iter1 != bdata->_snp_name_map.end()) {
                eintmp.clear(); eintmp.resize(etraitNum);
                incld.clear(); incld.resize(cmpnum);
                flip.clear(); flip.resize(cmpnum);
                for(int t=0; t<etraitNum; t++) {
                    iter2 = etrait[t]. _snp_name_map.find(grs);
                    if(t == 0 && iter2 != etrait[t]._snp_name_map.end()) {
                        ea1 = etrait[t]._esi_allele1[iter2->second];
                        ea2 = etrait[t]._esi_allele2[iter2->second];
                        allele_compare(ea1,ea2,ga1,ga2,incld[t],flip[t]);
                        ba1 = bdata->_allele1[iter1->second];
                        ba2 = bdata->_allele2[iter1->second];
                        allele_compare(ea1,ea2,ba1,ba2,incld[etraitNum],flip[etraitNum]);
                        eintmp[t] = iter2->second;
                    }
                    if (t > 0 && iter2 != etrait[t]._snp_name_map.end()) {
                        ea1tmp = etrait[t]._esi_allele1[iter2->second];
                        ea2tmp = etrait[t]._esi_allele2[iter2->second];
                        allele_compare(ea1,ea2,ea1tmp,ea2tmp,incld[t],flip[t]);
                        eintmp[t] = iter2->second;
                    }
                }
                int incld_sum = accumulate(incld.begin(), incld.end(), 0);
                if(incld_sum == cmpnum) {
                    // GWAS
                    if(flip[0]==1) {
                        gdata->byz[rid]=-1.0*gdata->byz[rid];
                    }
                    gin.push_back(rid);                
                    // xQTL[1:t]
                    for(int t=0; t<etraitNum; t++) {
                        ein[t].push_back(eintmp[t]);
                        if(flip[t]==1 && t > 0) {
                            if(etrait[t]._val.size()>0)
                            {
                                
                                int count=0;
                                for(int j=0;j<etrait[t]._rowid.size();j++)
                                {
                                    if(etrait[t]._rowid[j]==eintmp[t])
                                    {
                                        count++;
                                        if(count & 1)
                                            etrait[t]._val[j]=-1*etrait[t]._val[j];
                                    }
                                }
                                
                            }
                            else
                            {
                                #pragma omp parallel for private(i)
                                for(int j=0;j<etrait[t]._include.size();j++)   
                                    if( etrait[t]._bxz[j][eintmp[t]]+9 > 1e-6 ) etrait[t]._bxz[j][eintmp[t]]=-1*etrait[t]._bxz[j][eintmp[t]];
                            }
                        }
                    }
                    // LD reference
                    if(flip[etraitNum]==1) {
                        string tmpch=bdata->_ref_A[iter1->second];
                        bdata->_ref_A[iter1->second]=bdata->_other_A[iter1->second];
                        bdata->_other_A[iter1->second]=tmpch;
                    }
                    bin.push_back(iter1->second);
                }

            }
        }

        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        for(int t=0; t<etraitNum; t++) {
            etrait[t]._esi_include.swap(ein[t]);
        }
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check_multi_opt_old2(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        string logstr="Checking consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the xQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        int etraitNum = etrait.size(), cmpnum = etraitNum + 1;
        map<string, int>::iterator iter, iter1, iter2;
        vector<int> bin,gin,eintmp(etraitNum);
        vector<vector<int>> ein(etraitNum);
        // use 1st etrait as reference; [0]:etrait[0]vsGWAS; [1:t]:etrait[0]:etrait[1:t]; [t+1]:etrait[0]:LDref;         
        double disp=0;
        for(int i=0;i<gdata->_include.size();i++) //no map in gwas data
        {
            progress(i, disp, (int)gdata->_include.size());      
            int rid = gdata->_include[i];
            string grs = gdata->snpName[rid];
            string ba1, ba2, ga1, ga2, ea1, ea2, ea1tmp, ea2tmp;
            ga1 = gdata->allele_1[rid];
            ga2 = gdata->allele_2[rid];
            iter1 = bdata-> _snp_name_map.find(grs);
            if (iter1 != bdata->_snp_name_map.end()) {                
                vector<int> rsincld;
                eintmp.clear(); eintmp.resize(etraitNum);
                for(int t=0; t<etraitNum; t++) {
                    iter2 = etrait[t]. _snp_name_map.find(grs);
                    if(iter2 != etrait[t]._snp_name_map.end()) {
                        eintmp[t] = iter2->second;
                        rsincld.push_back(1);
                    }
                }
                if(rsincld.size() == etraitNum) {
                    vector<int> incld(cmpnum), flip(cmpnum);
                    for(int t=0; t<etraitNum; t++) {
                        if(t == 0) {
                            ea1 = etrait[t]._esi_allele1[eintmp[t]];
                            ea2 = etrait[t]._esi_allele2[eintmp[t]];
                            allele_compare(ea1,ea2,ga1,ga2,incld[t],flip[t]);
                            ba1 = bdata->_allele1[iter1->second];
                            ba2 = bdata->_allele2[iter1->second];
                            allele_compare(ea1,ea2,ba1,ba2,incld[etraitNum],flip[etraitNum]);
                        }
                        if (t > 0) {
                            ea1tmp = etrait[t]._esi_allele1[eintmp[t]];
                            ea2tmp = etrait[t]._esi_allele2[eintmp[t]];
                            allele_compare(ea1,ea2,ea1tmp,ea2tmp,incld[t],flip[t]);
                        }
                    }
                    int incld_sum = accumulate(incld.begin(), incld.end(), 0);
                    if(incld_sum == cmpnum) {
                        // GWAS
                        if(flip[0]==1) {
                            gdata->byz[rid]=-1.0*gdata->byz[rid];
                        }
                        gin.push_back(rid);                
                        // xQTL[1:t]
                        for(int t=0; t<etraitNum; t++) {
                            ein[t].push_back(eintmp[t]);
                            if(flip[t]==1 && t > 0) {
                                if(etrait[t]._val.size()>0)
                                {
                                    
                                    int count=0;
                                    for(int j=0;j<etrait[t]._rowid.size();j++)
                                    {
                                        if(etrait[t]._rowid[j]==eintmp[t])
                                        {
                                            count++;
                                            if(count & 1)
                                                etrait[t]._val[j]=-1*etrait[t]._val[j];
                                        }
                                    }
                                    
                                }
                                else
                                {
                                    #pragma omp parallel for private(i)
                                    for(int j=0;j<etrait[t]._include.size();j++)   
                                        if( etrait[t]._bxz[j][eintmp[t]]+9 > 1e-6 ) etrait[t]._bxz[j][eintmp[t]]=-1*etrait[t]._bxz[j][eintmp[t]];
                                }
                            }
                        }
                        // LD reference
                        if(flip[etraitNum]==1) {
                            string tmpch=bdata->_ref_A[iter1->second];
                            bdata->_ref_A[iter1->second]=bdata->_other_A[iter1->second];
                            bdata->_other_A[iter1->second]=tmpch;
                        }
                        bin.push_back(iter1->second);
                    }

                }       

            }
        }

        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        for(int t=0; t<etraitNum; t++) {
            etrait[t]._esi_include.swap(ein[t]);
        }
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check_multi_opt_test(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        string logstr="Checking consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the xQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        int etraitNum = etrait.size(), cmpnum = etraitNum + 1;
        map<string, int>::iterator iter, iter0, iter1, iter2;
        vector<int> bin,gin,eintmp(etraitNum);
        vector<vector<int>> ein(etraitNum);
        vector<string> gsnp(gdata->_include.size()),bsnp(bdata->_include.size()), esnp, cmmsnp, slctsnp;
        // extract SNPs in common between dataset and save in slctsnp;
        #pragma omp parallel for
        for(int i=0; i<gdata->_include.size(); i++)
            gsnp[i]=gdata->snpName[gdata->_include[i]];
        #pragma omp parallel for
        for(int i=0; i<bdata->_include.size(); i++)
            bsnp[i]=bdata->_snp_name[bdata->_include[i]];
        sort(gsnp.begin(), gsnp.end()); sort(bsnp.begin(), bsnp.end());
        set_intersection(gsnp.begin(), gsnp.end(), bsnp.begin(), bsnp.end(), back_inserter(slctsnp));
        for(int t=0; t<etraitNum; t++) {
            esnp.clear(); esnp.resize(etrait[t]._esi_include.size());
            #pragma omp parallel for    
            for(int i=0; i<etrait[t]._esi_include.size(); i++)
                esnp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
            cmmsnp.clear();
            sort(slctsnp.begin(), slctsnp.end()); sort(esnp.begin(), esnp.end());
            set_intersection(slctsnp.begin(), slctsnp.end(), esnp.begin(), esnp.end(), back_inserter(cmmsnp));
            slctsnp.clear(); slctsnp.resize(cmmsnp.size());
            #pragma omp parallel for
            for(int i=0; i<cmmsnp.size(); i++)
                slctsnp[i]=cmmsnp[i];        
        }
        // allele compare between SNPs in common; use [0:t-1]:etrait[0:t-1]vsGWAS; [t]:LDrefvsGWAS;         
        double disp=0;
        for(int i=0; i<slctsnp.size(); i++) //select SNPs only
        {
            progress(i, disp, (int)slctsnp.size());
            string rs = slctsnp[i];
            iter0 = gdata->_snp_name_map.find(rs);
            iter1 = bdata-> _snp_name_map.find(rs);
            if(iter0 != gdata->_snp_name_map.end() && iter1 != bdata-> _snp_name_map.end()) {
                string ga1, ga2, ba1, ba2, ea1, ea2;
                ga1 = gdata->allele_1[iter0->second];
                ga2 = gdata->allele_2[iter0->second];
                int rsincld = 0;
                eintmp.clear(); eintmp.resize(etraitNum);
                for(int t=0; t<etraitNum; t++) {
                    iter2 = etrait[t]. _snp_name_map.find(rs);
                    if(iter2 != etrait[t]._snp_name_map.end()) {
                        eintmp[t] = iter2->second;
                        rsincld = rsincld + 1;
                    }
                }
                if(rsincld == etraitNum) {
                    int incld_sum = 0;
                    vector<int> incld(cmpnum), flip(cmpnum);                
                    for(int t=0; t<etraitNum; t++) {
                        ea1 = etrait[t]._esi_allele1[eintmp[t]];
                        ea2 = etrait[t]._esi_allele2[eintmp[t]];                        
                        allele_compare(ga1,ga2,ea1,ea2,incld[t],flip[t]);
                        incld_sum = incld_sum + incld[t];
                    }
                    ba1 = bdata->_allele1[iter1->second];
                    ba2 = bdata->_allele2[iter1->second];
                    allele_compare(ga1,ga2,ba1,ba2,incld[etraitNum],flip[etraitNum]);
                    incld_sum = incld_sum + incld[etraitNum];
                    // int incld_sum = accumulate(incld.begin(), incld.end(), 0);
                    if(incld_sum == cmpnum) {
                        // GWAS
                        gin.push_back(iter0->second);                
                        // xQTL[0:t-1]
                        for(int t=0; t<etraitNum; t++) {                        
                            if(flip[t]==1) {
                                if(etrait[t]._val.size()>0)
                                {
                                    
                                    int count=0;
                                    for(int j=0;j<etrait[t]._rowid.size();j++)
                                    {
                                        if(etrait[t]._rowid[j]==eintmp[t])
                                        {
                                            count++;
                                            if(count & 1)
                                                etrait[t]._val[j]=-1*etrait[t]._val[j];
                                        }
                                    }
                                    
                                }
                                else
                                {
                                    //#pragma omp parallel for private(i)
                                    for(int j=0;j<etrait[t]._include.size();j++)   
                                        if( etrait[t]._bxz[j][eintmp[t]]+9 > 1e-6 ) etrait[t]._bxz[j][eintmp[t]]=-1*etrait[t]._bxz[j][eintmp[t]];
                                }
                            }
                            ein[t].push_back(eintmp[t]);
                        }
                        // LD reference
                        if(flip[etraitNum]==1) {
                            string tmpch=bdata->_ref_A[iter1->second];
                            bdata->_ref_A[iter1->second]=bdata->_other_A[iter1->second];
                            bdata->_other_A[iter1->second]=tmpch;
                        }
                        bin.push_back(iter1->second);
                    }

                }

            }

        }

        stable_sort(bin.begin(),bin.end());
        stable_sort(gin.begin(),gin.end());
        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        for(int t=0; t<etraitNum; t++) {
            stable_sort(ein[t].begin(),ein[t].end());
            etrait[t]._esi_include.swap(ein[t]);
        }
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check_multi_opt_old3(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        string logstr="Checking consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the xQTL summary data and the LD reference data).\n ";
        cout<<logstr<<endl;
        int etraitNum = etrait.size(), cmpnum = etraitNum + 1;
        map<string, int>::iterator iter, iter1, iter2;
        vector<int> bin,gin,eintmp(etraitNum);
        vector<vector<int>> ein(etraitNum);
        // use [0:t-1]:etrait[0:t-1]vsGWAS; [t]:LDrefvsGWAS;         
        double disp=0;
        for(int i=0;i<gdata->_include.size();i++) // gwas data
        {
            progress(i, disp, (int)gdata->_include.size());      
            int rid = gdata->_include[i];
            string grs = gdata->snpName[rid];
            string ba1, ba2, ga1, ga2, ea1, ea2;
            ga1 = gdata->allele_1[rid];
            ga2 = gdata->allele_2[rid];
            iter1 = bdata-> _snp_name_map.find(grs);
            if (iter1 != bdata->_snp_name_map.end()) {                
                int rsincld = 0;
                eintmp.clear(); eintmp.resize(etraitNum);
                for(int t=0; t<etraitNum; t++) {
                    iter2 = etrait[t]. _snp_name_map.find(grs);
                    if(iter2 != etrait[t]._snp_name_map.end()) {
                        eintmp[t] = iter2->second;
                        rsincld = rsincld + 1;
                    }
                }
                if(rsincld == etraitNum) {
                    int incld_sum = 0;
                    vector<int> incld(cmpnum), flip(cmpnum);                
                    for(int t=0; t<etraitNum; t++) {
                        ea1 = etrait[t]._esi_allele1[eintmp[t]];
                        ea2 = etrait[t]._esi_allele2[eintmp[t]];                        
                        allele_compare(ga1,ga2,ea1,ea2,incld[t],flip[t]);
                        incld_sum = incld_sum + incld[t];
                    }
                    ba1 = bdata->_allele1[iter1->second];
                    ba2 = bdata->_allele2[iter1->second];
                    allele_compare(ga1,ga2,ba1,ba2,incld[etraitNum],flip[etraitNum]);
                    incld_sum = incld_sum + incld[etraitNum];
                    // int incld_sum = accumulate(incld.begin(), incld.end(), 0);
                    if(incld_sum == cmpnum) {
                        gin.push_back(rid);                
                        // xQTL[0:t-1]
                        for(int t=0; t<etraitNum; t++) {                            
                            if(flip[t]==1) {
                                if(etrait[t]._val.size()>0)
                                {
                                    
                                    int count=0;
                                    for(int j=0;j<etrait[t]._rowid.size();j++)
                                    {
                                        if(etrait[t]._rowid[j]==eintmp[t])
                                        {
                                            count++;
                                            if(count & 1)
                                                etrait[t]._val[j]=-1*etrait[t]._val[j];
                                        }
                                    }
                                    
                                }
                                else
                                {
                                    #pragma omp parallel for private(i)
                                    for(int j=0;j<etrait[t]._include.size();j++)   
                                        if( etrait[t]._bxz[j][eintmp[t]]+9 > 1e-6 ) etrait[t]._bxz[j][eintmp[t]]=-1*etrait[t]._bxz[j][eintmp[t]];
                                }
                            }
                            ein[t].push_back(eintmp[t]);
                        }
                        // LD reference
                        if(flip[etraitNum]==1) {
                            string tmpch=bdata->_ref_A[iter1->second];
                            bdata->_ref_A[iter1->second]=bdata->_other_A[iter1->second];
                            bdata->_other_A[iter1->second]=tmpch;
                        }
                        bin.push_back(iter1->second);
                    }

                }       

            }
        }

        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        for(int t=0; t<etraitNum; t++) {
            etrait[t]._esi_include.swap(ein[t]);
        }
        printf("%ld SNPs are included after allele checking.\n ",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check_multi_opt(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        string logstr="\nChecking consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the xQTL summary data and the LD reference data).\nNote: running parallel jobs with multiple threads by '--thread-num' can accelarate this process.";
        cout<<logstr<<endl;
        int etraitNum = etrait.size(), cmpnum = etraitNum + 1;
        vector<int> bin,gin; vector<vector<int>> ein(etraitNum);
        vector<string> gsnp(gdata->_include.size()),bsnp(bdata->_include.size()), esnp, cmmsnp, slctsnp;
         // find the exposure Num with largest sites
        int t_max = 0; long num_tmp = etrait[0]._include.size();
        for(int t=0; t<etraitNum; t++) {
            if(etrait[t]._include.size() > num_tmp) {
                num_tmp = etrait[t]._include.size();
                t_max = t;
            }
        }
        // step 1: extract SNPs in common between dataset and save in slctsnp;
        #pragma omp parallel for
        for(int i=0; i<gdata->_include.size(); i++)
            gsnp[i]=gdata->snpName[gdata->_include[i]];
        #pragma omp parallel for
        for(int i=0; i<bdata->_include.size(); i++)
            bsnp[i]=bdata->_snp_name[bdata->_include[i]];
        sort(gsnp.begin(), gsnp.end()); sort(bsnp.begin(), bsnp.end());
        set_intersection(gsnp.begin(), gsnp.end(), bsnp.begin(), bsnp.end(), back_inserter(slctsnp));
        for(int t=0; t<etraitNum; t++) {
            esnp.clear(); esnp.resize(etrait[t]._esi_include.size());
            #pragma omp parallel for    
            for(int i=0; i<etrait[t]._esi_include.size(); i++)
                esnp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
            cmmsnp.clear();
            sort(slctsnp.begin(), slctsnp.end()); sort(esnp.begin(), esnp.end());
            set_intersection(slctsnp.begin(), slctsnp.end(), esnp.begin(), esnp.end(), back_inserter(cmmsnp));
            slctsnp.clear(); slctsnp.resize(cmmsnp.size());
            #pragma omp parallel for
            for(int i=0; i<cmmsnp.size(); i++)
                slctsnp[i]=cmmsnp[i];        
        }
        // reorder the SNP rsIDs based on bp positions using bfile
        #pragma omp parallel for
        for(int i=0; i<bdata->_include.size(); i++)
            bsnp[i]=bdata->_snp_name[bdata->_include[i]];
        vector<int> bidtmp;
        match_only(slctsnp, bsnp, bidtmp);        
        if(bidtmp.empty()) throw("Error: no common SNPs found between data! Please check your LD reference data.");
        sort(bidtmp.begin(), bidtmp.end());
        slctsnp.clear(); slctsnp.resize(bidtmp.size());
        for(int i=0; i<bidtmp.size(); i++) {
            slctsnp[i] = bsnp[bidtmp[i]];
        }        
        // step 2: allele compare between SNPs in common; use [0:t-1]:etrait[0:t-1]vsGWAS; [t]:LDrefvsGWAS;         
        vector<int> incld_slct_tmp(slctsnp.size()), bintmp(slctsnp.size()), gintmp(slctsnp.size());
        vector<vector<int>> flip_slct(cmpnum), flip_slct_tmp(cmpnum , vector<int>(slctsnp.size())), eintmp(etraitNum, vector<int>(slctsnp.size()));
        double disp=0;
        #pragma omp parallel for
        for(int i=0; i<slctsnp.size(); i++) //select SNPs only
        {
            //progress(i, disp, (int)slctsnp.size());
            map<string, int>::iterator iter0, iter1, iter2;
            string rs = slctsnp[i];
            vector<int> einsgl(etraitNum);
            iter0 = gdata->_snp_name_map.find(rs);
            iter1 = bdata-> _snp_name_map.find(rs);
            if(iter0 != gdata->_snp_name_map.end() && iter1 != bdata-> _snp_name_map.end()) {
                string ga1, ga2, ba1, ba2, ea1, ea2, ea1tmp, ea2tmp;
                ga1 = gdata->allele_1[iter0->second];
                ga2 = gdata->allele_2[iter0->second];
                for(int t=0; t<etraitNum; t++) {
                    iter2 = etrait[t]. _snp_name_map.find(rs);
                    if(iter2 != etrait[t]._snp_name_map.end()) {
                        einsgl[t] = iter2->second;
                        eintmp[t][i] = einsgl[t];
                    }
                }
                int incld_sum = 0;
                vector<int> incld(cmpnum), flip(cmpnum);
                ea1 = etrait[t_max]._esi_allele1[einsgl[t_max]];
                ea2 = etrait[t_max]._esi_allele2[einsgl[t_max]];                        
                allele_compare(ea1,ea2,ga1,ga2,incld[t_max],flip[t_max]);  
                incld_sum = incld_sum + incld[t_max];
                flip_slct_tmp[t_max][i] = flip[t_max];              
                for(int t=0; t<etraitNum; t++) {
                    if(t != t_max) {
                        ea1tmp = etrait[t]._esi_allele1[einsgl[t]];
                        ea2tmp = etrait[t]._esi_allele2[einsgl[t]];                        
                        allele_compare(ea1,ea2,ea1tmp,ea2tmp,incld[t],flip[t]);
                        incld_sum = incld_sum + incld[t];
                        flip_slct_tmp[t][i] = flip[t];
                    }
                }
                ba1 = bdata->_allele1[iter1->second];
                ba2 = bdata->_allele2[iter1->second];
                allele_compare(ea1,ea2,ba1,ba2,incld[etraitNum],flip[etraitNum]);
                incld_sum = incld_sum + incld[etraitNum];
                flip_slct_tmp[etraitNum][i] = flip[etraitNum];
                gintmp[i] = iter0->second;
                bintmp[i] = iter1->second;
                incld_slct_tmp[i] = incld_sum;
            }
        }
        // step 3: remove select SNPs with different alleles
        for(int i=0; i<incld_slct_tmp.size(); i++) // SNPs with matched alleles between data
        {
            if(incld_slct_tmp[i] == cmpnum) {
                gin.push_back(gintmp[i]);
                bin.push_back(bintmp[i]);
                for(int t=0; t<etraitNum; t++) {
                    ein[t].push_back(eintmp[t][i]);
                    flip_slct[t].push_back(flip_slct_tmp[t][i]);
                }
                flip_slct[etraitNum].push_back(flip_slct_tmp[etraitNum][i]);
            }
        }
        // step 4: switch the effect size direction for xQTL[-t_max], GWAS and LD reference
        #pragma omp parallel for
        for(int i=0; i<bin.size(); i++) //SNPs with matched rsID and alleles only
        {
            // xQTL[0:t-1]
            for(int t=0; t<etraitNum; t++) {                        
                if(flip_slct[t][i]==1 && t != t_max) {
                    if(etrait[t]._val.size()>0)
                    {
                        
                        int count=0;
                        for(int j=0;j<etrait[t]._rowid.size();j++)
                        {
                            if(etrait[t]._rowid[j]==ein[t][i])
                            {
                                count++;
                                if(count & 1)
                                    etrait[t]._val[j]=-1*etrait[t]._val[j];
                            }
                        }

                    }
                    else
                    {
                        for(int j=0;j<etrait[t]._include.size();j++)   
                            if( etrait[t]._bxz[j][ein[t][i]]+9 > 1e-6 ) etrait[t]._bxz[j][ein[t][i]]=-1*etrait[t]._bxz[j][ein[t][i]];
                    }
                    // switch the effect allele
                    string tmpch=etrait[t]._esi_allele1[ein[t][i]];
                    etrait[t]._esi_allele1[ein[t][i]]=etrait[t]._esi_allele2[ein[t][i]];
                    etrait[t]._esi_allele2[ein[t][i]]=tmpch;
                }
            }
            // GWAS [t_max]
            if(flip_slct[t_max][i]==1) {
                gdata->byz[gin[i]]=-1.0*gdata->byz[gin[i]];
                // switch the effect allele
                string tmpch=gdata->allele_1[gin[i]];
                gdata->allele_1[gin[i]]=gdata->allele_2[gin[i]];
                gdata->allele_2[gin[i]]=tmpch;
            }
            // LD reference etraitNum
            if(flip_slct[etraitNum][i]==1) {
                string tmpch=bdata->_ref_A[bin[i]];
                bdata->_ref_A[bin[i]]=bdata->_other_A[bin[i]];
                bdata->_other_A[bin[i]]=tmpch;
            }
        }

        // step 5: update _include with common SNPs index
        bdata->_include.swap(bin);
        gdata->_include.swap(gin);
        for(int t=0; t<etraitNum; t++) {
            etrait[t]._esi_include.swap(ein[t]);
        }
        printf("%ld SNPs are included after allele checking.\n",bdata->_include.size());
     
        // only update _snp_name_map which would be used in read bed file.
        map<string, int>::iterator iter;
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void allele_check_multi(vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId1;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among GWAS summary, exposure summary and LD reference data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> gsnp(gdata->_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum ){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                #pragma omp parallel for
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                #pragma omp parallel for
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with GWAS dataset //
        if(gdata->_include.size()< gdata->snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<gdata->_include.size();i++)
                gsnp[i]=gdata->snpName[gdata->_include[i]];
            match_only(gsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between GWAS data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=cmmnSNPs[i];
            //cout<<slctSNPs.size()<<endl;
        }else
        {
            match_only(etsnp, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs found between eTratis and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=gdata->snpName[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=cmmnSNPs[i];  
            //cout<<slctSNPs.size()<<endl;
        }
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        //alleles check
        match(slctSNPs, gdata->snpName, gdId1);
        //cout<<gdId1.size()<<endl;
        //match(slctSNPs, etrait->_esi_rs, gdId);
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        gdata->_include.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();
        }
        //etrait->_esi_include.clear();
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string ta1, ta2, ga1, ga2;
            ga1 = gdata->allele_1[gdId1[i]];
            ga2 = gdata->allele_2[gdId1[i]];
            //int Id1,flip1;
            //allele_compare(ga1,ga2,a1,a2,Id1,flip1);
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=0,Idsrt=0; //start with all Id = 1;
            for(int t = 0; t < etraitNum; t++)
            {   
                //gdId.clear();
                //match(slctSNPs, etrait[t]._esi_rs, gdId);
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ga1,ga2,ta1,ta2,Id[t],flip[t]);
                Idall = Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            //#pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                gdata->_include.push_back(gdId1[i]);
                // if(flip1==1){
                //     //cout<<"test1"<<endl;
                //     string tmpch=bdata->_ref_A[bdId[i]];
                //     bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                //     bdata->_other_A[bdId[i]]=tmpch;
                // }
                // if(flip3==1) {gdata->byz[gdId1[i]]=-gdata->byz[gdId1[i]];}
            }
        }
        //cout<<etrait[t]->_esi_include.size()<<endl;
        //cout<<bdata->_include.size()<<endl;
        //cout<<gdata->_include.size()<<endl;
        logstr=atos(gdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        // map<string, int> id_map_buf(bdata->_snp_name_map);
        // for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        // map<string, int>::iterator iter;
        // for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    void allele_check_multi(bInfo* bdata, vector<eqtlInfo> &etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Checking the consistency of SNP alleles among exposure summary, outcome summary and LD reference data. ";
        cout<<logstr<<endl;
        long etraitNum=etrait.size();
        vector<string> bsnp(bdata->_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        vector<string> etsnp(etrait[0]._esi_include.size());
        // matching with the multiple etrait first //
        if(etrait[0]._esi_include.size()<etrait[0]._snpNum ){
            #pragma omp parallel for    
            for(int i=0;i<etrait[0]._esi_include.size();i++)
            etsnp[i]=etrait[0]._esi_rs[etrait[0]._esi_include[i]];
        } else {
            etsnp=etrait[0]._esi_rs;
        }       
        for (int t = 1; t<etraitNum; t++)
        {   
            if(etrait[t]._esi_include.size()<etrait[t]._snpNum)
            {   
                vector<string> etsnptmp(etrait[t]._esi_include.size());
                #pragma omp parallel for    
                for(int i=0;i<etrait[t]._esi_include.size();i++)
                etsnptmp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
                match_only(etsnptmp,etsnp,edId);
                if(edId.empty()) throw("Error: no common SNPs between eTraits data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            } else {
                match_only(etrait[t]._esi_rs, etsnp, edId);
                if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
                cmmnSNPs.resize(edId.size());
                #pragma omp parallel for
                for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
                etsnp.resize(cmmnSNPs.size());
                for(int i=0;i<cmmnSNPs.size();i++)
                etsnp[i]=cmmnSNPs[i];
            }
        }
        // matching with multiple dataset //
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum)
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(bsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }else
        {
            match_only(bdata->_snp_name, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];  
        }
        match(slctSNPs, bdata->_snp_name, bdId);
        long snpcmn=slctSNPs.size();
        vector<int> gdId(etraitNum*snpcmn);
        for(int t = 0; t < etraitNum; t++){
            vector<int> gdIdtmp;
            match(slctSNPs, etrait[t]._esi_rs, gdIdtmp);
            for(int i=0;i<slctSNPs.size();i++)
            gdId[t*snpcmn+i]=gdIdtmp[i];
        }
        cmmnSNPs.clear();
        bdata->_include.clear();
        for(int t = 0; t < etraitNum; t++){
            etrait[t]._esi_include.clear();    
        }
        esdata->_esi_include.clear();
        // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
        double cr=0;
        for (int i = 0; i<edId.size(); i++)
        {   
            double desti=1.0*i/(edId.size()-1);
            if(desti>=cr)
            {
                printf("Checking...  %3.0f%%\r", 100.0*desti);
                fflush(stdout);
                if(cr==0) cr+=0.05;
                else if(cr==0.05) cr+=0.2;
                else if(cr==0.25) cr+=0.5;
                else cr+=0.25;
            }
            //if(i%500==0) printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            //fflush(stdout);
            string a1, a2, ta1, ta2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            int Id1,flip1;
            allele_compare(ea1,ea2,a1,a2,Id1,flip1);
            vector<int> Id(etraitNum), flip(etraitNum);
            int Idall=1,Idsrt=1; //start with all Id = 1;
            for(int t = 0; t < etraitNum; t++)
            {   
                ta1 = etrait[t]._esi_allele1[gdId[t*snpcmn+i]];
                ta2 = etrait[t]._esi_allele2[gdId[t*snpcmn+i]];
                allele_compare(ea1,ea2,ta1,ta2,Id[t],flip[t]);
                Idall=Idall+Id[t];
            }

            for(int t = 0; t < etraitNum; t++)
            {      
                if(Idall==(etraitNum+Idsrt))
                {
                    etrait[t]._esi_include.push_back(gdId[t*snpcmn+i]);
                    if(flip[t]==1)
                    {
                        if(etrait[t]._val.size()>0)
                        {
                            
                            int count=0;
                            for(int j=0;j<etrait[t]._rowid.size();j++)
                            {
                                if(etrait[t]._rowid[j]==gdId[t*snpcmn+i])
                                {
                                    count++;
                                    if(count & 1)
                                        etrait[t]._val[j]=-1*etrait[t]._val[j];
                                }
                            }
                            
                        }
                        else
                        {
                            #pragma omp parallel for private(i)
                            for(int j=0;j<etrait[t]._include.size();j++)   
                                if( etrait[t]._bxz[j][gdId[t*snpcmn+i]]+9 > 1e-6 ) etrait[t]._bxz[j][gdId[t*snpcmn+i]]=-1*etrait[t]._bxz[j][gdId[t*snpcmn+i]];
                        }
                    }
                }
            }
            if (Idall==(etraitNum+Idsrt)) {
                bdata->_include.push_back(bdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                if(flip1==1){
                    //cout<<"test1"<<endl;
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                }
            }
        }
        logstr=atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    void update_geIndx(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
    {
        long etraitNum=etrait.size();
        vector<int> gtmpIdx;
        vector<int> etmpIdx;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            gtmpIdx.push_back(gdata->_include[bdata->_include[i]]);
        }
        gdata->_include.clear();
        gdata->_include = gtmpIdx;
        
        for(int t=0;t<etraitNum;t++)
        {
            for (int i = 0; i < bdata->_include.size(); i++)
            {
                etmpIdx.push_back(etrait[t]._esi_include[bdata->_include[i]]);
            }
            etrait[t]._esi_include.clear();
            etrait[t]._esi_include = etmpIdx;
            etmpIdx.clear();
        }
    }
    void update_geIndx(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata)
    {
        long etraitNum=etrait.size();
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        vector<int> tmpIdx3;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            //tmpIdx1.push_back(etrait->_esi_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
            tmpIdx3.push_back(gdata->_include[bdata->_include[i]]);
        }

        //etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        gdata->_include.clear();
        //etrait->_esi_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
        gdata->_include = tmpIdx3;
        for(int t=0;t<etraitNum;t++)
        {
            for (int i = 0; i < bdata->_include.size(); i++)
            {
                tmpIdx1.push_back(etrait[t]._esi_include[bdata->_include[i]]);
            }
            etrait[t]._esi_include.clear();
            etrait[t]._esi_include = tmpIdx1;
            tmpIdx1.clear();
        }
    }
    void update_geIndx(bInfo* bdata, vector<eqtlInfo> &etrait, eqtlInfo* esdata)
    {
        long etraitNum=etrait.size();
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
        }

        esdata->_esi_include.clear();
        esdata->_esi_include = tmpIdx2;
        for(int t=0;t<etraitNum;t++)
        {
            for (int i = 0; i < bdata->_include.size(); i++)
            {
                tmpIdx1.push_back(etrait[t]._esi_include[bdata->_include[i]]);
            }
            etrait[t]._esi_include.clear();
            etrait[t]._esi_include = tmpIdx1;
            tmpIdx1.clear();
        }
    }
    void init_smr_wk_mlt(MTSMRWK* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear();
    }
    void init_smr_wk_mlt(MTSMRWKEXP* smrwk)
    {
        long outcoNum = smrwk->bxz.size();
        for(int t=0;t<outcoNum;t++)
        {
            smrwk->bxz[t].clear();
            smrwk->sexz[t].clear();
            smrwk->freq[t].clear();
            smrwk->zxz[t].clear();
        }
        smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear();
    }
    long fill_smr_wk_mlt(bInfo* bdata,vector<gwasData> &gdata,eqtlInfo* esdata,MTSMRWK* smrwk, const char* refSNP, int lowerbp,int upperbp,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid=-9;
        long outcoNum=gdata.size();
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (abs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        for( int t=0;t<outcoNum;t++)
                        {
                        smrwk->byz[t].push_back(gdata[t].byz[j]);
                        smrwk->seyz[t].push_back(gdata[t].seyz[j]);
                        smrwk->pyz[t].push_back(gdata[t].pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        }
                        if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                        else smrwk->freq.push_back(esdata->_esi_freq[j]);
                    }
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
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                
                if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    for( int t=0;t<outcoNum;t++)
                    {
                    smrwk->byz[t].push_back(gdata[t].byz[ge_rowid]);
                    smrwk->seyz[t].push_back(gdata[t].seyz[ge_rowid]);
                    smrwk->pyz[t].push_back(gdata[t].pvalue[ge_rowid]);
                    }
                    smrwk->curId.push_back(ge_rowid);
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                    else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                    
                }
            }
        }
        return maxid;
    }
    long fill_smr_wk_mlt(bInfo* bdata,vector<gwasData> &gdata,eqtlInfo* esdata,MTSMRWK* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        long outcoNum=gdata.size();
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (abs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    int sumtrue=0;
                    for( int t=0;t<outcoNum;t++)
                    { if(gdata[t].seyz[j]+9>1e-6) sumtrue+=1;}
                    if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl && sumtrue==outcoNum)
                    {
                        if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                            {
                                smrwk->bxz.push_back(esdata->_bxz[i][j]);
                                smrwk->sexz.push_back(esdata->_sexz[i][j]);
                                smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                for( int t=0;t<outcoNum;t++)
                                {
                                smrwk->byz[t].push_back(gdata[t].byz[j]);
                                smrwk->seyz[t].push_back(gdata[t].seyz[j]);
                                smrwk->pyz[t].push_back(gdata[t].pvalue[j]);
                                }
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata->_esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                                smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                                smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                                if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                                if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq.push_back(esdata->_esi_freq[j]);
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                        
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            smrwk->bxz.push_back(esdata->_bxz[i][j]);
                            smrwk->sexz.push_back(esdata->_sexz[i][j]);
                            smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            for( int t=0;t<outcoNum;t++)
                            {
                            smrwk->byz[t].push_back(gdata[t].byz[j]);
                            smrwk->seyz[t].push_back(gdata[t].seyz[j]);
                            smrwk->pyz[t].push_back(gdata[t].pvalue[j]);
                            }
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata->_esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                            if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                            else smrwk->freq.push_back(esdata->_esi_freq[j]);

                        }
                    }
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
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                int sumtmp=0;
                for( int t=0;t<outcoNum;t++)
                {
                    if(gdata[t].seyz[ge_rowid]+9>1e-6) sumtmp+=1;
                }
                if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl && sumtmp==outcoNum)
                {
                    if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0)
                    {
                        if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                        {
                            smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                            smrwk->sexz.push_back(esdata->_val[se_start+j]);
                            smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                            for( int t=0;t<outcoNum;t++)
                            {
                            smrwk->byz[t].push_back(gdata[t].byz[ge_rowid]);
                            smrwk->seyz[t].push_back(gdata[t].seyz[ge_rowid]);
                            smrwk->pyz[t].push_back(gdata[t].pvalue[ge_rowid]);
                            }
                            smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                            smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                            if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                            else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                        } else {
                            printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                            double z=(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                            double p=pchisq(z*z, 1);
                            string tmp=atos(esdata->_esi_rs[j])+"\t"+ atos(esdata->_esi_chr[j])+"\t"+ atos(esdata->_esi_bp[j])+"\t"+ atos(esdata->_esi_allele1[j])+"\t"+ atos(esdata->_esi_allele2[j])+"\t"+ atos(esdata->_esi_freq[j])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][j])+"\t"+ atos(esdata->_sexz[i][j])+"\t"+ dtos(p)+"\n";
                            
                            printf("%s\n",tmp.c_str());
                        }
                        
                    } else {
                        smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                        smrwk->sexz.push_back(esdata->_val[se_start+j]);
                        smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                        smrwk->byz.resize(outcoNum);smrwk->seyz.resize(outcoNum);smrwk->pyz.resize(outcoNum);
                        for( int t=0;t<outcoNum;t++)
                        {
                        smrwk->byz[t].push_back(gdata[t].byz[ge_rowid]);
                        smrwk->seyz[t].push_back(gdata[t].seyz[ge_rowid]);
                        smrwk->pyz[t].push_back(gdata[t].pvalue[ge_rowid]);
                        }
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                        if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                        else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                    }
                  
                }
            }
        }
        
        return maxid;
    }

    // need to update the function to extract the workspace for cis variants
    long fill_smr_wk_mlt_old(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        long outcoNum=esdata.size();

        if(esdata[0]._rowid.empty())
        {
            for (int j = 0; j<esdata[0]._esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata[0]._bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata[0]._esi_bp[j];
                    int snpchr=esdata[0]._esi_chr[j];
                    if(snpchr==esdata[0]._epi_chr[i] && ABS(esdata[0]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[j]+9>1e-6)
                    {
                        if(esdata[0]._epi_start.size()>0 && esdata[0]._epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata[0]._epi_end[i]==-9 || (snpbp>esdata[0]._epi_end[i] && snpbp<esdata[0]._epi_start[i]))
                            {
                                for( int t=0; t<outcoNum; t++)
                                {
                                    smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                    smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                    smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                                }
                                smrwk->byz.push_back(gdata->byz[j]);
                                smrwk->seyz.push_back(gdata->seyz[j]);
                                smrwk->pyz.push_back(gdata->pvalue[j]);                                    
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata[0]._esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata[0]._esi_chr[j]);
                                smrwk->allele1.push_back(esdata[0]._esi_allele1[j]);
                                smrwk->allele2.push_back(esdata[0]._esi_allele2[j]);
                                if(refSNP!=NULL && esdata[0]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata[0]._esi_bp[j]);
                                
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata[0]._bxz[i][j]/esdata[0]._sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata[0]._esi_rs[j])+"\t"+ atos(esdata[0]._esi_chr[j])+"\t"+ atos(esdata[0]._esi_bp[j])+"\t"+ atos(esdata[0]._esi_allele1[j])+"\t"+ atos(esdata[0]._esi_allele2[j])+"\t"+ atos(esdata[0]._esi_freq[j])+"\t"+ atos(esdata[0]._epi_prbID[i])+"\t"+ atos(esdata[0]._epi_chr[i])+"\t"+ atos(esdata[0]._epi_bp[i])+"\t" + atos(esdata[0]._epi_gene[i])+"\t"+ atos(esdata[0]._epi_orien[i])+"\t"+ atos(esdata[0]._bxz[i][j])+"\t"+ atos(esdata[0]._sexz[i][j])+"\t"+ dtos(p)+"\n";
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            for( int t=0; t<outcoNum; t++)
                            {
                                smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                            }
                            smrwk->byz.push_back(gdata->byz[j]);
                            smrwk->seyz.push_back(gdata->seyz[j]);
                            smrwk->pyz.push_back(gdata->pvalue[j]);                            
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata[0]._esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata[0]._esi_chr[j]);
                            smrwk->allele1.push_back(esdata[0]._esi_allele1[j]);
                            smrwk->allele2.push_back(esdata[0]._esi_allele2[j]);
                            if(refSNP!=NULL && esdata[0]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata[0]._esi_bp[j]);
                            
                        }
                    }
                }
            }                     
        }
        
        else{
                vector<uint32_t> ge_rowid_common;
                vector<uint32_t> ge_rowid_first;
                vector<int> common_idx;
                for( int t=0;t<outcoNum;t++)
                {
                    uint64_t beta_start=esdata[t]._cols[i<<1];
                    uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    vector<uint32_t> ge_rowid_tmp;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata[t]._rowid[beta_start+j];
                        int snpbp=esdata[t]._esi_bp[ge_rowid];
                        int snpchr=esdata[t]._esi_chr[ge_rowid];
                        if(snpchr==esdata[t]._epi_chr[i] && ABS(esdata[t]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[ge_rowid]+9>1e-6)
                        {
                            if(t==0) {ge_rowid_first.push_back(ge_rowid);}
                            ge_rowid_tmp.push_back(ge_rowid);
                        }
                    }
                    common_idx.clear();
                    CommFunc::match_only(ge_rowid_tmp,ge_rowid_first,common_idx);
                    //if(common_idx.empty()) throw("Error: no common SNPs between multiple traits.");
                    ge_rowid_common.resize(common_idx.size());
                    #pragma omp parallel for
                    for(int i=0;i<common_idx.size();i++)
                    ge_rowid_common[i]=ge_rowid_first[common_idx[i]];
                    ge_rowid_first.resize(ge_rowid_common.size());
                    for(int i=0;i<ge_rowid_common.size();i++)
                    ge_rowid_first[i]=ge_rowid_common[i];
                }
            
                uint64_t numsnps=ge_rowid_common.size();
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=ge_rowid_common[j];
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);                            
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->rs.push_back(esdata[0]._esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata[0]._esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata[0]._esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata[0]._esi_allele2[ge_rowid]);
                    if(refSNP!=NULL && esdata[0]._esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata[0]._esi_bp[ge_rowid]);
                }
            
                vector<int> val_idx;
                for( int t=0;t<outcoNum;t++)
                {
                    uint64_t beta_start=esdata[t]._cols[i<<1];
                    uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                    uint64_t numsnps=ge_rowid_common.size();
                    val_idx.clear();
                    match_only(ge_rowid_common,esdata[t]._rowid,val_idx);
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=ge_rowid_common[j];
                        int snpbp=esdata[t]._esi_bp[ge_rowid];

                            if(esdata[t]._epi_start.size()>0 && esdata[t]._epi_end.size()>0)
                            {
                                if(esdata[t]._epi_end[i]==-9 || (snpbp>esdata[t]._epi_end[i] && snpbp<esdata[t]._epi_start[i]))
                                {
                                    smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                    smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                    smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);
                                    
                                } else {
                                    printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                    double z=(esdata[t]._bxz[i][val_idx[j]]/esdata[t]._sexz[i][val_idx[j]]);
                                    double p=pchisq(z*z, 1);
                                    // string tmp=atos(esdata[t]._esi_rs[j])+"\t"+ atos(esdata[t]._esi_chr[j])+"\t"+ atos(esdata[t]._esi_bp[j])+"\t"+ atos(esdata[t]._esi_allele1[j])+"\t"+ atos(esdata[t]._esi_allele2[j])+"\t"+ atos(esdata[t]._esi_freq[j])+"\t"+ atos(esdata[t]._epi_prbID[i])+"\t"+ atos(esdata[t]._epi_chr[i])+"\t"+ atos(esdata[t]._epi_bp[i])+"\t" + atos(esdata[t]._epi_gene[i])+"\t"+ atos(esdata[t]._epi_orien[i])+"\t"+ atos(esdata[t]._bxz[i][j])+"\t"+ atos(esdata[t]._sexz[i][j])+"\t"+ dtos(p)+"\n";                            
                                    // printf("%s\n",tmp.c_str());
                                }
                                
                            } else {
                                    smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                    smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                    smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);                                
                                }                 
                        // }
                    }
                }                
            }       
        return maxid;
    }
    // updated; need to update the function to extract the workspace for cis variants
    long fill_smr_wk_mlt_old2(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        long outcoNum=esdata.size();

        if(esdata[0]._rowid.empty())
        {
            for (int j = 0; j<esdata[0]._esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata[0]._bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata[0]._esi_bp[j];
                    int snpchr=esdata[0]._esi_chr[j];
                    if(snpchr==esdata[0]._epi_chr[i] && ABS(esdata[0]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[j]]+9>1e-6)
                    {
                        if(esdata[0]._epi_start.size()>0 && esdata[0]._epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata[0]._epi_end[i]==-9 || (snpbp>esdata[0]._epi_end[i] && snpbp<esdata[0]._epi_start[i]))
                            {
                                for( int t=0; t<outcoNum; t++)
                                {
                                    smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                    smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                    smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                                }
                                smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                                smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                                smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);                                                                        
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata[0]._esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata[0]._esi_chr[j]);
                                smrwk->allele1.push_back(esdata[0]._esi_allele1[j]);
                                smrwk->allele2.push_back(esdata[0]._esi_allele2[j]);
                                if(refSNP!=NULL && esdata[0]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata[0]._esi_bp[j]);
                                
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata[0]._bxz[i][j]/esdata[0]._sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata[0]._esi_rs[j])+"\t"+ atos(esdata[0]._esi_chr[j])+"\t"+ atos(esdata[0]._esi_bp[j])+"\t"+ atos(esdata[0]._esi_allele1[j])+"\t"+ atos(esdata[0]._esi_allele2[j])+"\t"+ atos(esdata[0]._esi_freq[j])+"\t"+ atos(esdata[0]._epi_prbID[i])+"\t"+ atos(esdata[0]._epi_chr[i])+"\t"+ atos(esdata[0]._epi_bp[i])+"\t" + atos(esdata[0]._epi_gene[i])+"\t"+ atos(esdata[0]._epi_orien[i])+"\t"+ atos(esdata[0]._bxz[i][j])+"\t"+ atos(esdata[0]._sexz[i][j])+"\t"+ dtos(p)+"\n";
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            for( int t=0; t<outcoNum; t++)
                            {
                                smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                            }
                            smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);                             
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata[0]._esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata[0]._esi_chr[j]);
                            smrwk->allele1.push_back(esdata[0]._esi_allele1[j]);
                            smrwk->allele2.push_back(esdata[0]._esi_allele2[j]);
                            if(refSNP!=NULL && esdata[0]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata[0]._esi_bp[j]);
                            
                        }
                    }
                }
            }                     
        }
        
        else{
                vector<uint32_t> ge_rowid_common;
                vector<uint32_t> ge_rowid_first;
                vector<int> common_idx;
                for( int t=0;t<outcoNum;t++)
                {
                    uint64_t beta_start=esdata[t]._cols[i<<1];
                    uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    vector<uint32_t> ge_rowid_tmp;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata[t]._rowid[beta_start+j];
                        int snpbp=esdata[t]._esi_bp[ge_rowid];
                        int snpchr=esdata[t]._esi_chr[ge_rowid];
                        if(snpchr==esdata[t]._epi_chr[i] && ABS(esdata[t]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[ge_rowid]]+9>1e-6)
                        {
                            if(t==0) {ge_rowid_first.push_back(ge_rowid);}
                            ge_rowid_tmp.push_back(ge_rowid);
                        }
                    }
                    common_idx.clear();
                    CommFunc::match_only(ge_rowid_tmp,ge_rowid_first,common_idx);
                    //if(common_idx.empty()) throw("Error: no common SNPs between multiple traits.");
                    ge_rowid_common.resize(common_idx.size());
                    #pragma omp parallel for
                    for(int i=0;i<common_idx.size();i++)
                    ge_rowid_common[i]=ge_rowid_first[common_idx[i]];
                    ge_rowid_first.resize(ge_rowid_common.size());
                    for(int i=0;i<ge_rowid_common.size();i++)
                    ge_rowid_first[i]=ge_rowid_common[i];
                }
                if(ge_rowid_common.size()!=0) {
                    uint64_t numsnps=ge_rowid_common.size();
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=ge_rowid_common[j];
                        smrwk->byz.push_back(gdata->byz[gdata->_include[ge_rowid]]);
                        smrwk->seyz.push_back(gdata->seyz[gdata->_include[ge_rowid]]);
                        smrwk->pyz.push_back(gdata->pvalue[gdata->_include[ge_rowid]]);                           
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata[0]._esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata[0]._esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata[0]._esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata[0]._esi_allele2[ge_rowid]);
                        if(refSNP!=NULL && esdata[0]._esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata[0]._esi_bp[ge_rowid]);
                    }
                
                    vector<int> val_idx;
                    for( int t=0;t<outcoNum;t++)
                    {
                        uint64_t beta_start=esdata[t]._cols[i<<1];
                        uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                        uint64_t numsnps=ge_rowid_common.size();
                        val_idx.clear();
                        match_only(ge_rowid_common,esdata[t]._rowid,val_idx);
                        for(int j=0;j<numsnps;j++)
                        {
                            int ge_rowid=ge_rowid_common[j];
                            int snpbp=esdata[t]._esi_bp[ge_rowid];

                                if(esdata[t]._epi_start.size()>0 && esdata[t]._epi_end.size()>0)
                                {
                                    if(esdata[t]._epi_end[i]==-9 || (snpbp>esdata[t]._epi_end[i] && snpbp<esdata[t]._epi_start[i]))
                                    {
                                        smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                        smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                        smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                        if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                        else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);
                                        
                                    } else {
                                        printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                        double z=(esdata[t]._bxz[i][val_idx[j]]/esdata[t]._sexz[i][val_idx[j]]);
                                        double p=pchisq(z*z, 1);
                                        // string tmp=atos(esdata[t]._esi_rs[j])+"\t"+ atos(esdata[t]._esi_chr[j])+"\t"+ atos(esdata[t]._esi_bp[j])+"\t"+ atos(esdata[t]._esi_allele1[j])+"\t"+ atos(esdata[t]._esi_allele2[j])+"\t"+ atos(esdata[t]._esi_freq[j])+"\t"+ atos(esdata[t]._epi_prbID[i])+"\t"+ atos(esdata[t]._epi_chr[i])+"\t"+ atos(esdata[t]._epi_bp[i])+"\t" + atos(esdata[t]._epi_gene[i])+"\t"+ atos(esdata[t]._epi_orien[i])+"\t"+ atos(esdata[t]._bxz[i][j])+"\t"+ atos(esdata[t]._sexz[i][j])+"\t"+ dtos(p)+"\n";                            
                                        // printf("%s\n",tmp.c_str());
                                    }
                                    
                                } else {
                                        smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                        smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                        smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                        if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                        else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);                                
                                    }                 
                            // }
                        }
                    } 

                }                              
            }       
        return maxid;
    }

    // updated; need to update the function to extract the workspace for cis variants
    long fill_smr_wk_mlt(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        long outcoNum=esdata.size();
        // find the exposure Num with smallest number of cis-SNPs
        int t_min = 0; 

        if(esdata[0]._rowid.empty())
        {
            long num_tmp = esdata[0]._bxz[i].size();
            for(int t=0; t<outcoNum; t++) {
                if(esdata[t]._bxz[i].size() < num_tmp) {
                    num_tmp = esdata[t]._bxz[i].size();
                    t_min = t;
                }
            }
            for (int j = 0; j<esdata[t_min]._esi_include.size() ; j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata[t_min]._bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata[t_min]._esi_bp[j];
                    int snpchr=esdata[t_min]._esi_chr[j];
                    if(snpchr==esdata[t_min]._epi_chr[i] && ABS(esdata[t_min]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[j]]+9>1e-6)
                    {
                        if(esdata[t_min]._epi_start.size()>0 && esdata[t_min]._epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata[t_min]._epi_end[i]==-9 || (snpbp>esdata[t_min]._epi_end[i] && snpbp<esdata[t_min]._epi_start[i]))
                            {
                                for( int t=0; t<outcoNum; t++)
                                {
                                    smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                    smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                    smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                                }
                                smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                                smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                                smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);                                    
                                smrwk->curId.push_back(j);
                                smrwk->rs.push_back(esdata[t_min]._esi_rs[j]);
                                smrwk->snpchrom.push_back(esdata[t_min]._esi_chr[j]);
                                smrwk->allele1.push_back(esdata[t_min]._esi_allele1[j]);
                                smrwk->allele2.push_back(esdata[t_min]._esi_allele2[j]);
                                if(refSNP!=NULL && esdata[t_min]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata[t_min]._esi_bp[j]);
                                
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata[t_min]._bxz[i][j]/esdata[t_min]._sexz[i][j]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata[t_min]._esi_rs[j])+"\t"+ atos(esdata[t_min]._esi_chr[j])+"\t"+ atos(esdata[t_min]._esi_bp[j])+"\t"+ atos(esdata[t_min]._esi_allele1[j])+"\t"+ atos(esdata[t_min]._esi_allele2[j])+"\t"+ atos(esdata[t_min]._esi_freq[j])+"\t"+ atos(esdata[t_min]._epi_prbID[i])+"\t"+ atos(esdata[t_min]._epi_chr[i])+"\t"+ atos(esdata[t_min]._epi_bp[i])+"\t" + atos(esdata[t_min]._epi_gene[i])+"\t"+ atos(esdata[t_min]._epi_orien[i])+"\t"+ atos(esdata[t_min]._bxz[i][j])+"\t"+ atos(esdata[t_min]._sexz[i][j])+"\t"+ dtos(p)+"\n";
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            for( int t=0; t<outcoNum; t++)
                            {
                                smrwk->bxz[t].push_back(esdata[t]._bxz[i][j]);
                                smrwk->sexz[t].push_back(esdata[t]._sexz[i][j]);
                                smrwk->zxz[t].push_back(esdata[t]._bxz[i][j]/esdata[t]._sexz[i][j]);
                                if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq[t].push_back(esdata[t]._esi_freq[j]);
                            }
                            smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);                            
                            smrwk->curId.push_back(j);
                            smrwk->rs.push_back(esdata[t_min]._esi_rs[j]);
                            smrwk->snpchrom.push_back(esdata[t_min]._esi_chr[j]);
                            smrwk->allele1.push_back(esdata[t_min]._esi_allele1[j]);
                            smrwk->allele2.push_back(esdata[t_min]._esi_allele2[j]);
                            if(refSNP!=NULL && esdata[t_min]._esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata[t_min]._esi_bp[j]);
                            
                        }
                    }
                }
            }                     
        }
        
        else{
                long num_tmp = esdata[0]._valNum;
                for(int t=0; t<outcoNum; t++) {
                    if(esdata[t]._valNum < num_tmp) {
                        num_tmp = esdata[t]._valNum;
                        t_min = t;
                    }
                }
                vector<uint32_t> ge_rowid_common;
                vector<uint32_t> ge_rowid_first;
                vector<int> common_idx;
                for( int t=0;t<outcoNum;t++)
                {
                    uint64_t beta_start=esdata[t]._cols[i<<1];
                    uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    vector<uint32_t> ge_rowid_tmp;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata[t]._rowid[beta_start+j];
                        int snpbp=esdata[t]._esi_bp[ge_rowid];
                        int snpchr=esdata[t]._esi_chr[ge_rowid];
                        if(snpchr==esdata[t]._epi_chr[i] && ABS(esdata[t]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[ge_rowid]]+9>1e-6)
                        {
                            if(t==0) {ge_rowid_first.push_back(ge_rowid);}
                            ge_rowid_tmp.push_back(ge_rowid);
                        }
                    }
                    common_idx.clear();
                    CommFunc::match_only(ge_rowid_tmp,ge_rowid_first,common_idx);
                    //if(common_idx.empty()) throw("Error: no common SNPs between multiple traits.");
                    ge_rowid_common.resize(common_idx.size());
                    #pragma omp parallel for
                    for(int i=0;i<common_idx.size();i++)
                    ge_rowid_common[i]=ge_rowid_first[common_idx[i]];
                    ge_rowid_first.resize(ge_rowid_common.size());
                    for(int i=0;i<ge_rowid_common.size();i++)
                    ge_rowid_first[i]=ge_rowid_common[i];
                }
                if(ge_rowid_common.size()!=0) {
                    uint64_t numsnps=ge_rowid_common.size();
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=ge_rowid_common[j];
                        smrwk->byz.push_back(gdata->byz[gdata->_include[ge_rowid]]);
                        smrwk->seyz.push_back(gdata->seyz[gdata->_include[ge_rowid]]);
                        smrwk->pyz.push_back(gdata->pvalue[gdata->_include[ge_rowid]]);                            
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata[t_min]._esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata[t_min]._esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata[t_min]._esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata[t_min]._esi_allele2[ge_rowid]);
                        if(refSNP!=NULL && esdata[t_min]._esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata[t_min]._esi_bp[ge_rowid]);
                    }
                
                    vector<int> val_idx;
                    for( int t=0;t<outcoNum;t++)
                    {
                        uint64_t beta_start=esdata[t]._cols[i<<1];
                        uint64_t se_start=esdata[t]._cols[1+(i<<1)];
                        uint64_t numsnps=ge_rowid_common.size();
                        val_idx.clear();
                        match_only(ge_rowid_common,esdata[t]._rowid,val_idx);
                        for(int j=0;j<numsnps;j++)
                        {
                            int ge_rowid=ge_rowid_common[j];
                            int snpbp=esdata[t]._esi_bp[ge_rowid];

                                if(esdata[t]._epi_start.size()>0 && esdata[t]._epi_end.size()>0)
                                {
                                    if(esdata[t]._epi_end[i]==-9 || (snpbp>esdata[t]._epi_end[i] && snpbp<esdata[t]._epi_start[i]))
                                    {
                                        smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                        smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                        smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                        if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                        else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);
                                        
                                    } else {
                                        printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                        double z=(esdata[t]._bxz[i][val_idx[j]]/esdata[t]._sexz[i][val_idx[j]]);
                                        double p=pchisq(z*z, 1);
                                        // string tmp=atos(esdata[t]._esi_rs[j])+"\t"+ atos(esdata[t]._esi_chr[j])+"\t"+ atos(esdata[t]._esi_bp[j])+"\t"+ atos(esdata[t]._esi_allele1[j])+"\t"+ atos(esdata[t]._esi_allele2[j])+"\t"+ atos(esdata[t]._esi_freq[j])+"\t"+ atos(esdata[t]._epi_prbID[i])+"\t"+ atos(esdata[t]._epi_chr[i])+"\t"+ atos(esdata[t]._epi_bp[i])+"\t" + atos(esdata[t]._epi_gene[i])+"\t"+ atos(esdata[t]._epi_orien[i])+"\t"+ atos(esdata[t]._bxz[i][j])+"\t"+ atos(esdata[t]._sexz[i][j])+"\t"+ dtos(p)+"\n";                            
                                        // printf("%s\n",tmp.c_str());
                                    }
                                    
                                } else {
                                        smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                        smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                        smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                        if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                                        else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);                                
                                    }                 
                            // }
                        }
                    }

                }                                          
            }       
        return maxid;
    }    
    
    void update_snidx(MTSMRWK* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
    {
        if(sn_ids.size()>max_snp_slct)
        {
            //printf("Top %d SNPs (ordered by eQTL p-value) are used int the %s...\n",max_snp_slct,forwhat.c_str());
            priority_queue< pair<double, int> > q;
            vector<int> sn_slct_ids;
            for(int i=0;i<sn_ids.size();i++) q.push((pair<double, int>(abs(smrwk->zxz[sn_ids[i]]), sn_ids[i])));
            for(int i=0;i<max_snp_slct;i++){
                sn_slct_ids.push_back(q.top().second);
                q.pop();
            }
            sn_ids.swap(sn_slct_ids);
        }
        
    }
    void update_snidx(MTSMRWKEXP* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
    {
        if(sn_ids.size()>max_snp_slct)
        {
            //printf("Top %d SNPs (ordered by eQTL p-value) are used int the %s...\n",max_snp_slct,forwhat.c_str());
            priority_queue< pair<double, int> > q;
            vector<int> sn_slct_ids;
            for(int i=0;i<sn_ids.size();i++) q.push((pair<double, int>(abs(smrwk->zxz[0][sn_ids[i]]), sn_ids[i])));
            for(int i=0;i<max_snp_slct;i++){
                sn_slct_ids.push_back(q.top().second);
                q.pop();
            }
            sn_ids.swap(sn_slct_ids);
        }        
    }
    void update_snidx(MTSMRWKEXP* smrwk,int t, vector<int> &sn_ids,int max_snp_slct, string forwhat)
    {
        if(sn_ids.size()>max_snp_slct)
        {
            //printf("Top %d SNPs (ordered by eQTL p-value) are used int the %s...\n",max_snp_slct,forwhat.c_str());
            priority_queue< pair<double, int> > q;
            vector<int> sn_slct_ids;
            for(int i=0;i<sn_ids.size();i++) q.push((pair<double, int>(abs(smrwk->zxz[t][sn_ids[i]]), sn_ids[i])));
            for(int i=0;i<max_snp_slct;i++){
                sn_slct_ids.push_back(q.top().second);
                q.pop();
            }
            sn_ids.swap(sn_slct_ids);
        }        
    }
    void extract_smrwk(MTSMRWK* smrwk,vector<int> &sn_ids,MTSMRWK* smrwk2)
    {
        long outcoNum=smrwk->byz.size();
        smrwk2->cur_chr=smrwk->cur_chr;
        smrwk2->cur_prbidx=smrwk->cur_prbidx;
        smrwk2->byz.resize(outcoNum);
        smrwk2->seyz.resize(outcoNum);
        smrwk2->pyz.resize(outcoNum);
        for(int i=0;i<sn_ids.size();i++)
        {
            smrwk2->bxz.push_back(smrwk->bxz[sn_ids[i]]);
            smrwk2->sexz.push_back(smrwk->sexz[sn_ids[i]]);
            smrwk2->freq.push_back(smrwk->freq[sn_ids[i]]);
            for(int t=0;t<outcoNum;t++)
            {
            smrwk2->byz[t].push_back(smrwk->byz[t][sn_ids[i]]);
            smrwk2->seyz[t].push_back(smrwk->seyz[t][sn_ids[i]]);
            smrwk2->pyz[t].push_back(smrwk->pyz[t][sn_ids[i]]);
            }
            smrwk2->zxz.push_back(smrwk->zxz[sn_ids[i]]);
            smrwk2->curId.push_back(smrwk->curId[sn_ids[i]]);
            smrwk2->bpsnp.push_back(smrwk->bpsnp[sn_ids[i]]);
            smrwk2->snpchrom.push_back(smrwk->snpchrom[sn_ids[i]]);
            smrwk2->rs.push_back(smrwk->rs[sn_ids[i]]);
            smrwk2->allele1.push_back(smrwk->allele1[sn_ids[i]]);
            smrwk2->allele2.push_back(smrwk->allele2[sn_ids[i]]);
        }
    }
    void extract_smrwk(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MTSMRWKEXP* smrwk2)
    {
        long outcoNum=smrwk->bxz.size();
        smrwk2->cur_chr=smrwk->cur_chr;
        smrwk2->cur_prbidx=smrwk->cur_prbidx;
        smrwk2->bxz.resize(outcoNum);
        smrwk2->sexz.resize(outcoNum);
        smrwk2->zxz.resize(outcoNum);
        smrwk2->freq.resize(outcoNum);
        for(int i=0;i<sn_ids.size();i++)
        {
            smrwk2->byz.push_back(smrwk->byz[sn_ids[i]]);
            smrwk2->seyz.push_back(smrwk->seyz[sn_ids[i]]);
            smrwk2->pyz.push_back(smrwk->pyz[sn_ids[i]]);
            for(int t=0;t<outcoNum;t++)
            {
                smrwk2->bxz[t].push_back(smrwk->bxz[t][sn_ids[i]]);
                smrwk2->sexz[t].push_back(smrwk->sexz[t][sn_ids[i]]);
                smrwk2->freq[t].push_back(smrwk->freq[t][sn_ids[i]]);
                smrwk2->zxz[t].push_back(smrwk->zxz[t][sn_ids[i]]);
            }            
            smrwk2->curId.push_back(smrwk->curId[sn_ids[i]]);
            smrwk2->bpsnp.push_back(smrwk->bpsnp[sn_ids[i]]);
            smrwk2->snpchrom.push_back(smrwk->snpchrom[sn_ids[i]]);
            smrwk2->rs.push_back(smrwk->rs[sn_ids[i]]);
            smrwk2->allele1.push_back(smrwk->allele1[sn_ids[i]]);
            smrwk2->allele2.push_back(smrwk->allele2[sn_ids[i]]);
        }
    }
    void update_smrwk_x(MTSMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X)
    {
        vector<vector<double>> byz,seyz,pyz;
        vector<double> bxz,sexz,zxz,freq;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
        MatrixXd _X;
        byz.resize(smrwk->byz.size()); //byz.resize(sn_ids.size());
        seyz.resize(smrwk->seyz.size()); //seyz.resize(sn_ids.size());
        pyz.resize(smrwk->pyz.size()); //pyz.resize(sn_ids.size());
        bxz.resize(sn_ids.size());
        sexz.resize(sn_ids.size());
        zxz.resize(sn_ids.size());
        freq.resize(sn_ids.size());
        curId.resize(sn_ids.size());
        bpsnp.resize(sn_ids.size());
        snpchrom.resize(sn_ids.size());
        rs.resize(sn_ids.size());
        allele1.resize(sn_ids.size());
        allele2.resize(sn_ids.size());
        _X.resize(X.rows(), sn_ids.size());
        
        //#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            for(int t=0;t<smrwk->byz.size();t++){
                byz[t].push_back(smrwk->byz[t][sn_ids[j]]);
                seyz[t].push_back(smrwk->seyz[t][sn_ids[j]]);
                pyz[t].push_back(smrwk->pyz[t][sn_ids[j]]);
            }
            bxz[j]=smrwk->bxz[sn_ids[j]];
            sexz[j]=smrwk->sexz[sn_ids[j]];
            zxz[j]=smrwk->zxz[sn_ids[j]];
            freq[j]=smrwk->freq[sn_ids[j]];
            curId[j]=smrwk->curId[sn_ids[j]];
            bpsnp[j]=smrwk->bpsnp[sn_ids[j]];
            snpchrom[j]=smrwk->snpchrom[sn_ids[j]];
            rs[j]=smrwk->rs[sn_ids[j]];
            allele1[j]=smrwk->allele1[sn_ids[j]];
            allele2[j]=smrwk->allele2[sn_ids[j]];
            _X.col(j)=X.col(sn_ids[j]);
        }
        smrwk->bxz.swap(bxz);
        smrwk->sexz.swap(sexz);
        smrwk->freq.swap(freq);
        for(int t=0;t<smrwk->byz.size();t++){
        smrwk->byz[t].swap(byz[t]);
        smrwk->seyz[t].swap(seyz[t]);
        smrwk->pyz[t].swap(pyz[t]);
        }
        smrwk->zxz.swap(zxz);
        smrwk->curId.swap(curId);
        smrwk->bpsnp.swap(bpsnp);
        smrwk->snpchrom.swap(snpchrom);
        smrwk->rs.swap(rs);
        smrwk->allele1.swap(allele1);
        smrwk->allele2.swap(allele2);
        X=_X;
    }
    void update_smrwk_x(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MatrixXd &X)
    {
        vector<double> byz,seyz,pyz;
        vector<vector<double>> bxz,sexz,zxz,freq;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
        MatrixXd _X;
        
        bxz.resize(smrwk->bxz.size());
        sexz.resize(smrwk->sexz.size());
        zxz.resize(smrwk->zxz.size());
        freq.resize(smrwk->freq.size());
        byz.resize(sn_ids.size());        
        seyz.resize(sn_ids.size());
        pyz.resize(sn_ids.size());        
        curId.resize(sn_ids.size());
        bpsnp.resize(sn_ids.size());
        snpchrom.resize(sn_ids.size());
        rs.resize(sn_ids.size());
        allele1.resize(sn_ids.size());
        allele2.resize(sn_ids.size());
        _X.resize(X.rows(), sn_ids.size());
        
        //#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            for(int t=0;t<smrwk->bxz.size();t++) {
                bxz[t].push_back(smrwk->bxz[t][sn_ids[j]]);
                sexz[t].push_back(smrwk->sexz[t][sn_ids[j]]);
                zxz[t].push_back(smrwk->zxz[t][sn_ids[j]]);
                freq[t].push_back(smrwk->freq[t][sn_ids[j]]);
            }
            byz[j]=(smrwk->byz[sn_ids[j]]);
            seyz[j]=(smrwk->seyz[sn_ids[j]]);
            pyz[j]=(smrwk->pyz[sn_ids[j]]);
            curId[j]=smrwk->curId[sn_ids[j]];
            bpsnp[j]=smrwk->bpsnp[sn_ids[j]];
            snpchrom[j]=smrwk->snpchrom[sn_ids[j]];
            rs[j]=smrwk->rs[sn_ids[j]];
            allele1[j]=smrwk->allele1[sn_ids[j]];
            allele2[j]=smrwk->allele2[sn_ids[j]];
            _X.col(j)=X.col(sn_ids[j]);
        }
        for(int t=0;t<smrwk->bxz.size();t++){
           smrwk->bxz[t].swap(bxz[t]);
           smrwk->sexz[t].swap(sexz[t]);
           smrwk->zxz[t].swap(zxz[t]);
           smrwk->freq[t].swap(freq[t]);
        }
        smrwk->byz.swap(byz);
        smrwk->seyz.swap(seyz);
        smrwk->pyz.swap(pyz);
        smrwk->curId.swap(curId);
        smrwk->bpsnp.swap(bpsnp);
        smrwk->snpchrom.swap(snpchrom);
        smrwk->rs.swap(rs);
        smrwk->allele1.swap(allele1);
        smrwk->allele2.swap(allele2);
        X=_X;
    }
    // end
    void update_freq(char* eqtlFileName, char* frqfile)
    {
        eqtlInfo eqtlinfo;
         if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for updating allele frequency information by the flag --eqtl-summary.");
        
        read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
        ifstream frq(frqfile);
        if (!frq) throw ("Error: can not open the file [" + string(frqfile) + "] to read.");
        cout << "Reading allele frequency information from [" + string(frqfile) + "]." << endl;
        FILE* fsnpfptr=fopen((string(frqfile)+".unmatched.snp.list").c_str(), "w");
        if (!(fsnpfptr)) {
            printf("Error: Failed to open unmatched SNP log file.\n");
        }

        vector<string> vs_buf;
        map<string, int>::iterator iter;
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        int fsnpNum(0);
        int falleleNum(0);
        int hitNum(0);
        while(!frq.eof())
        {
            frq.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, ", \t\n");
                if(col_num!=4) {
                    printf("ERROR: column number is not correct in row %d!\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na" || vs_buf[0]=="0"){
                    printf("ERROR: SNP name is \"NA\" in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[3]=="NA" || vs_buf[3]=="na" ){
                    printf("ERROR: frequency of the effect allele is \"NA\" in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(atof(vs_buf[3].c_str())>=1 || atof(vs_buf[3].c_str())<=0)
                {
                    printf("ERROR: frequency should be between (0 , 1) in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                iter=eqtlinfo._snp_name_map.find(vs_buf[0]);
                if(iter!=eqtlinfo._snp_name_map.end())
                {
                    int snpidx=iter->second;
                    if(eqtlinfo._esi_allele1[snpidx]==vs_buf[1] && eqtlinfo._esi_allele2[snpidx]==vs_buf[2]) {
                        eqtlinfo._esi_freq[snpidx]=atof(vs_buf[3].c_str());
                        hitNum++;
                    }
                    else if(eqtlinfo._esi_allele1[snpidx]==vs_buf[2] && eqtlinfo._esi_allele2[snpidx]==vs_buf[1])
                    {
                        eqtlinfo._esi_freq[snpidx]=1-atof(vs_buf[3].c_str());
                        hitNum++;
                    }
                    else {
                        falleleNum++;
                        string logstr=vs_buf[0]+'\t'+"failedAlleleCheck"+'\n';
                        fputs(logstr.c_str(),fsnpfptr);
                    }
                    
                } else {
                    fsnpNum++;
                    string logstr=vs_buf[0]+'\t'+"failedSNPMatch"+'\n';
                    fputs(logstr.c_str(),fsnpfptr);
                }
                lineNum++;
            }
        }
        fclose(fsnpfptr);
        frq.close();
        cout << lineNum << " SNPs to be included from [" + string(frqfile) + "]." << endl;
        printf("%d SNPs hitted, %d SNPs failed in SNP match and %d SNPs failed in allele check.\n",hitNum,fsnpNum,falleleNum);
        
        string esifile = string(eqtlFileName)+".esi";
        ofstream smr(esifile.c_str());
        if (!smr) throw ("Error: can not open the ESI file " + esifile + " to save!");
        for (int i = 0;i <eqtlinfo._snpNum; i++) {
            smr<<eqtlinfo._esi_chr[i]<<'\t'<<eqtlinfo._esi_rs[i]<<'\t'<<eqtlinfo._esi_gd[i]<<'\t'<<eqtlinfo._esi_bp[i]<<'\t'<<eqtlinfo._esi_allele1[i]<<'\t'<<eqtlinfo._esi_allele2[i]<<'\t'<<eqtlinfo._esi_freq[i]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo._snpNum<<" SNPs have been saved in the file [" + esifile + "]."<<endl;
        
       
    }
    
    void write_besd(string outFileName, eqtlInfo* eqtlinfo)
    {
        filter_probe_null(eqtlinfo); // at the same time, reset the vector _include
        //filter_snp_null(eqtlinfo);
        //update_besd();
        cout<<"\nsaving eQTL data..."<<endl;
        string esdfile = string(outFileName)+".esi";
        ofstream smr(esdfile.c_str());
        if (!smr) throw ("Error: can not open the ESI file " + esdfile + " to save!");
        for (int i = 0;i <eqtlinfo->_snpNum; i++) {
            smr<<eqtlinfo->_esi_chr[i]<<'\t'<<eqtlinfo->_esi_rs[i]<<'\t'<<eqtlinfo->_esi_gd[i]<<'\t'<<eqtlinfo->_esi_bp[i]<<'\t'<<eqtlinfo->_esi_allele1[i]<<'\t'<<eqtlinfo->_esi_allele2[i]<<'\t'<<eqtlinfo->_esi_freq[i]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo->_snpNum<<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
        
        esdfile = string(outFileName)+".epi";
        smr.open(esdfile.c_str());
        if (!smr) throw ("Error: can not open the EPI file " + esdfile + " to save!");
        for (int i = 0;i <eqtlinfo->_include.size(); i++) {
            smr<<eqtlinfo->_epi_chr[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_gd[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_bp[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_gene[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_orien[eqtlinfo->_include[i]]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo->_include.size()<<" probes have been saved in the file [" + esdfile + "]."<<endl;
        
        esdfile = string(outFileName)+".besd";
        FILE * smrbesd;
        smrbesd = fopen (esdfile.c_str(), "wb");
        if(eqtlinfo->_valNum==0)
        {
            uint64_t bsize=(eqtlinfo->_include.size()*eqtlinfo->_snpNum<<1)+1;
            float* buffer=(float*)malloc (sizeof(float)*bsize);
            memset(buffer,0,sizeof(float)*bsize);
            float* ptr=buffer;
            *ptr++=0.0;
            uint64_t pro_num=eqtlinfo->_include.size();
            uint64_t snp_num=eqtlinfo->_snpNum;
            for(int i=0;i<pro_num;i++)
            {
                memcpy(ptr+(i<<1)*snp_num,&eqtlinfo->_bxz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
                memcpy(ptr+((i<<1)+1)*snp_num,&eqtlinfo->_sexz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
            }
            fwrite(buffer,sizeof(float), bsize, smrbesd);
            free(buffer);
        }
        else
        {
            uint64_t colSize=sizeof(uint64_t)*((eqtlinfo->_include.size()<<1)+1);
            uint64_t rowSize=sizeof(uint32_t)*eqtlinfo->_valNum;
            uint64_t valSize=sizeof(float)*eqtlinfo->_valNum;
            uint64_t valNum=eqtlinfo->_valNum;
            uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
            
            char* buffer=(char*)malloc (sizeof(char)*bufsize);
            memset(buffer,0,sizeof(char)*bufsize);
            uint32_t ftype=SPARSE_FILE_TYPE_3F;
            memcpy(buffer,&ftype,sizeof(uint32_t));
            char* wptr=buffer+sizeof(float);
            memcpy(wptr,&valNum,sizeof(uint64_t));
            wptr+=sizeof(uint64_t);
            uint64_t* uptr=(uint64_t*)wptr; *uptr++=0;
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                *uptr++=eqtlinfo->_cols[(eqtlinfo->_include[i]<<1)+1];
                *uptr++=eqtlinfo->_cols[eqtlinfo->_include[i]+1<<1];
            }
            wptr+=colSize;
            memcpy(wptr,&eqtlinfo->_rowid[0],rowSize);
            wptr+=rowSize;
            memcpy(wptr,&eqtlinfo->_val[0],valSize);
            fwrite (buffer,sizeof(char), bufsize, smrbesd);
            free(buffer);

        }

        
        fclose (smrbesd);
        
        cout<<"Effect sizes (beta) and SE for "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
        

    }
    
    void meta(char* outFileName,char* eqtlFileName, char* eqtlFileName2)
    {
        setNbThreads(thread_num);
        
        eqtlInfo etrait;
        eqtlInfo esdata;
        eqtlInfo metadata;
        
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(eqtlFileName2==NULL) throw("Error: please input another eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        vector<int> idx;
        vector<string> cmmnSymbs;
        match_only(esdata._esi_rs, etrait._esi_rs, idx);
        etrait._esi_include=idx;
        for(int i=0;i<idx.size();i++) cmmnSymbs.push_back(etrait._esi_rs[idx[i]]);
        idx.clear();
        match_only(cmmnSymbs,esdata._esi_rs,idx);
        esdata._esi_include=idx;
        cout<<"There are "<<idx.size()<<" SNPs in common."<<endl;
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        idx.clear();
        match_only(esdata._epi_prbID, etrait._epi_prbID, idx);
        etrait._include=idx;
        cmmnSymbs.clear();
        for(int i=0;i<idx.size();i++) cmmnSymbs.push_back(etrait._epi_prbID[idx[i]]);
        idx.clear();
        match_only(cmmnSymbs,esdata._epi_prbID,idx);
        esdata._include=idx;
        cout<<"There are "<<idx.size()<<" Probes in common."<<endl;
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s in the analysis.\n",eqtlFileName2);
            exit(EXIT_FAILURE);
        }
        metadata._cols.push_back(0);
        cout<<"Performing Meta analysis..."<<endl;
        map<string, int> unmatched_rs_map;
        int unmatched_map_size=0;
        for(int i=0;i<etrait._probNum;i++)
        {
            string ref_pid=etrait._epi_prbID[i];
            string alt_pid=esdata._epi_prbID[i];
            if(ref_pid!=alt_pid)
            {
                cout<<"Some bugs here, please help to report!"<<endl;
                exit(1);
            }
            
            vector<float> ref_byz;
            vector<float> ref_seyz;
            vector<string> ref_a1;
            vector<string> ref_a2;
            vector<uint32_t> ref_rowid;
            vector<float> meta_beta;
            vector<float> meta_se;
            vector<uint32_t> meta_rowid;
            
            uint64_t beta_start=etrait._cols[i<<1];
            uint64_t se_start=etrait._cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(uint64_t j=0;j<numsnps;j++)
            {
                uint32_t ge_rowid=etrait._rowid[beta_start+j];
                ref_rowid.push_back(ge_rowid);
                ref_a1.push_back(etrait._esi_allele1[ge_rowid]);
                ref_a2.push_back(etrait._esi_allele2[ge_rowid]);
                ref_byz.push_back(etrait._val[beta_start+j]);
                ref_seyz.push_back(etrait._val[se_start+j]);
            }
            
            vector<float> byz;
            vector<float> seyz;
            vector<string> a1;
            vector<string> a2;
            vector<uint32_t> rowid;
          
            
            uint64_t beta_stt=esdata._cols[i<<1];
            uint64_t se_stt=esdata._cols[1+(i<<1)];
            uint64_t nums=se_stt-beta_stt;
            for(uint64_t j=0;j<nums;j++)
            {
                uint32_t ge_rid=esdata._rowid[beta_stt+j];
                rowid.push_back(ge_rid);
                a1.push_back(esdata._esi_allele1[ge_rid]);
                a2.push_back(esdata._esi_allele2[ge_rid]);
                byz.push_back(esdata._val[beta_stt+j]);
                seyz.push_back(esdata._val[se_stt+j]);
               
            }
            
            vector<int> idx;
            match(ref_rowid,rowid,idx);
            for(int j=0;j<idx.size();j++)
            {
                if(idx[j]==-9) continue;
                uint32_t tmp_rowid=ref_rowid[j];
                string tmp_rs=etrait._esi_rs[tmp_rowid];
                if(tmp_rowid!=rowid[idx[j]])
                {
                    cout<<"Some bugs here, please help to report!"<<endl;
                    exit(1);
                }
                
                string tmp_ref_a1=ref_a1[j];
                string tmp_ref_a2=ref_a2[j];
                string tmp_alt_a1=a1[idx[j]];
                string tmp_alt_a2=a2[idx[j]];
                float se1=ref_seyz[j];
                float se2=seyz[idx[j]];
                float beta1=ref_byz[j];
                float beta2=byz[idx[j]];
                
                
                if (tmp_ref_a1 == tmp_alt_a1 && tmp_ref_a2 == tmp_alt_a2)
                {
                    float tmpSE=se1*se2/sqrt(se1*se1+se2*se2);
                    meta_rowid.push_back(tmp_rowid);
                    meta_se.push_back(tmpSE);
                    meta_beta.push_back((beta1/(se1*se1) + beta2/(se2*se2))*tmpSE*tmpSE);
                }
                else if(tmp_ref_a1 == tmp_alt_a2 && tmp_ref_a2 == tmp_alt_a1)
                {
                    beta2=-beta2;
                    float tmpSE=se1*se2/sqrt(se1*se1+se2*se2);
                    meta_rowid.push_back(tmp_rowid);
                    meta_se.push_back(tmpSE);
                    meta_beta.push_back((beta1/(se1*se1) + beta2/(se2*se2))*tmpSE*tmpSE);
                }
                else {
                    unmatched_rs_map.insert(pair<string, int>(tmp_rs, unmatched_map_size));
                   
                }
            }
        
            for(int j=0;j<meta_beta.size();j++) metadata._val.push_back(meta_beta[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._val.push_back(meta_se[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._rowid.push_back(meta_rowid[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._rowid.push_back(meta_rowid[j]);
            metadata._cols.push_back(meta_beta.size()+metadata._cols[metadata._cols.size()-1]);
            metadata._cols.push_back(meta_beta.size()+metadata._cols[metadata._cols.size()-1]);
           
        }
        metadata._esi_allele1=etrait._esi_allele1;
        metadata._esi_allele2=etrait._esi_allele2;
        metadata._esi_bp=etrait._esi_bp;
        metadata._esi_chr=etrait._esi_chr;
        metadata._esi_gd=etrait._esi_gd;
        metadata._esi_include=etrait._esi_include;
        metadata._esi_rs=etrait._esi_rs;
        metadata._epi_bp=etrait._epi_bp;
        metadata._epi_chr=etrait._epi_chr;
        metadata._epi_gd=etrait._epi_gd;
        metadata._epi_gene=etrait._epi_gene;
        metadata._epi_orien=etrait._epi_orien;
        metadata._epi_prbID=etrait._epi_prbID;
        metadata._include=etrait._include;
        metadata._probNum=etrait._probNum;
        metadata._snpNum=etrait._snpNum;
        metadata._valNum=metadata._val.size();
        
        string unmatchedsnpfname = string(outFileName)+".unmatched.snp.list";
        FILE* unmatchedsnpfile=fopen(unmatchedsnpfname.c_str(), "w");
        if (!(unmatchedsnpfile)) {
            printf("Error: Failed to open unmatchedsnpfile file.\n");
            exit(1);
        }
        for(std::map<string,int>::iterator it=unmatched_rs_map.begin(); it!=unmatched_rs_map.end(); ++it){
            string tmpstr=it->first+'\n';
            fputs(tmpstr.c_str(),unmatchedsnpfile);
        }
        fclose(unmatchedsnpfile);
        write_besd(outFileName, &metadata);
        
    }
    
    void esi_man(eqtlInfo* eqtlinfo,char* snplstName,int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int cis_itvl,const char* prbname)
    {
        string logstr;
        int flags4snp=0;
        if(snplstName != NULL) flags4snp++;
        if(snprs != NULL) flags4snp++;
        if(fromsnprs!=NULL) flags4snp++;
        if(fromsnpkb>=0) flags4snp++;
        if(flags4snp>1)
        {
            logstr="WARNING: Flags for SNPs in this section are mutual exclusive. The priority order (from high to low) is: --extract-snp, --snp-wind, --snp, --from(to)--snp, --from(to)-snp-kb.\n";
            fputs(logstr.c_str(), stdout);
        }
        if(snpchr!=0)
        {
            extract_eqtl_by_chr(eqtlinfo, snpchr);
        }
        else if(chr!=0)
        {
            extract_eqtl_by_chr(eqtlinfo, chr);
        }
        
        if(prbname!=NULL && cis_flag)
        {
            extract_eqtl_snp(eqtlinfo, prbname, cis_itvl, "probe"); // extract cis eQTLs
        }
        else if (snplstName != NULL) extract_eqtl_snp(eqtlinfo, snplstName);
        else if (snpwindFlag)
        {
            if(snprs==NULL)
            {
                logstr="ERROR: please specify the SNP name by --snp when using --snp-wind.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_snp(eqtlinfo, snprs, snpWind, "SNP");
        }
        else if(snprs!=NULL)
        {
            extract_eqtl_single_snp(eqtlinfo, snprs);
        }
        else if(fromsnprs!=NULL)
        {
            if(tosnprs==NULL)
            {
                logstr="ERROR: please specify the SNP name by --to-snp.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_snp(eqtlinfo, fromsnprs, tosnprs);
        }
        else if(fromsnpkb>=0)
        {
            
            if(fromsnpkb>=0 && chr==0 && snpchr==0) {
                logstr="ERROR: please specify the chromosome by --snp-chr or --chr.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            
            if(tosnpkb<0)
            {
                logstr="ERROR: SNP BP can't be negative.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            if(snpchr!=0) extract_eqtl_snp(eqtlinfo, snpchr, fromsnpkb, tosnpkb);
            else if (chr!=0) extract_eqtl_snp(eqtlinfo, chr, fromsnpkb, tosnpkb);
        }
        
    }
    void epi_man(eqtlInfo* eqtlinfo,char* problstName,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename)
    {
        string logstr;
        int flags4prb=0;
        if(problstName != NULL) flags4prb++;
        if(prbname != NULL) flags4prb++;
        if(fromprbname!=NULL) flags4prb++;
        if(fromprbkb>=0) flags4prb++;
        if(genename != NULL) flags4prb++;
        if(flags4prb>1)
        {
            logstr="WARNING: Flags for probes in this section are mutual exclusive. The priority order (from high to low) is: --extract-probe, --gene-list, --probe-wind, --probe, --from(to)--probe, --from(to)-probe-kb, --gene.\n";
            fputs(logstr.c_str(), stdout);
        }
        
        if(prbchr!=0)
        {
            extract_epi_by_chr(eqtlinfo, prbchr);
        }
        else if(chr!=0)
        {
            extract_epi_by_chr(eqtlinfo, chr);
        }
        
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_prob(eqtlinfo, problstName);
            if(genelistName != NULL) extract_prob_by_gene(eqtlinfo, genelistName);
        }
        else if(prbwindFlag)
        {
            if(prbname==NULL)
            {
                logstr="ERROR: Please identify the probe name by --probe when using --probe-wind.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_prob(eqtlinfo, prbname, prbWind);
        }
        else if(prbname!=NULL)
        {
            extract_eqtl_single_probe(eqtlinfo, prbname);
        }
        else if(fromprbname!=NULL)
        {
            if(toprbname==NULL)
            {
                logstr="ERROR: Please identify the probe name by --to-probe.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_prob(eqtlinfo, fromprbname, toprbname);
        }
        else if(fromprbkb>=0)
        {
            if(fromprbkb>=0 && chr==0 && prbchr==0) {
                logstr="ERROR: Please identify the chromosome by --probe-chr or --chr.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            if(toprbkb<0)
            {
                logstr="ERROR: probe BP can't be negative.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_prob(eqtlinfo, prbchr, fromprbkb, toprbkb);
        }
        else if(genename!=NULL)
        {
            
            extract_prob_by_single_gene(eqtlinfo, genename);
        }
        
    }
    
}
