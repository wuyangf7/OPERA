//  OPERA.cpp
//
//  Created by Yang Wu on 04/06/18.
//  Copyright (c) 2021 Yang Wu. All rights reserved.
//

#include "SMR_data.h"
#include "OPERA.h"
#include "stat.hpp"
#include "SMR_data_p1.h"

namespace SMRDATA
{
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
         printf(" %ld SNPs to be included from PLINK BIM files.\n", bdata->_snp_num);
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
         cout << "Reading the PLINK BIM files ..." << endl;
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
         printf(" %ld individuals have been included from the PLINK FAM files.\n", bdata->_indi_num);
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

    void read_numfile(string numFilename, vector<long> &numsum)
    {
        FILE* nfile=fopen(string(numFilename).c_str(),"r");
        if(!nfile) {
            printf("File %s open failed.\n",string(numFilename).c_str());
            exit(EXIT_FAILURE);
        }
        printf("\nReading the number of tests performed for different exposure combinations in stage 2 analysis from %s ......\n", string(numFilename).c_str());
        char Tbuf[MAX_LINE_SIZE];    
        vector<string> numlist;
        fgets(Tbuf, MAX_LINE_SIZE, nfile);
        split_string(Tbuf, numlist, " \t\n");
        int numsize = numlist.size();
        int expoNum = numsize - 1;
        if(numlist[0]!="expoNum") 
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(numFilename).c_str());
            exit(EXIT_FAILURE);
        }
        int line_num = 0; 
        while(fgets(Tbuf, MAX_LINE_SIZE, nfile))
        {   
            split_string(Tbuf, numlist, ", \t\n");
            if(numlist.size() != numsize) {
                printf("ERROR: The %ld row from the input file %s is incomplete! Please check.\n", line_num,  string(numFilename).c_str());
                exit(EXIT_FAILURE);
            }
            if(numlist[0]!="testNum") 
            {
                printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(numFilename).c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=1;i<numlist.size();i++) {
                numsum[i-1]+=stol(numlist[i]);
            }
            line_num++;
        }
        fclose(nfile);
        printf("There are %d line(s) for %ld exposure combinations included in %s.\n",line_num, expoNum, string(numFilename).c_str());
    }

    void read_varfile(string varFilename, vector<string> &sigmasplit)
    {
        // Read prior variance file
        sigmasplit.clear();
        FILE* rfile=fopen(varFilename.c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",varFilename.c_str());
            exit(EXIT_FAILURE);
        }
        printf("\nReading the estimated prior variance from %s ...\n", varFilename.c_str());
        char Tbuf[MAX_LINE_SIZE];
        int line_idx=0;
        vector<string> strlist;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(line_idx == 0 && strlist[0]!="Variance") 
            {
                printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check.\n", varFilename.c_str());
                exit(EXIT_FAILURE);
            }        
            if(line_idx == 0) {
                for(int i=1;i<strlist.size();i++) {
                    sigmasplit.push_back(strlist[i]);
                }                
            }
            line_idx++;
        }
        if(line_idx != 1)
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check.\n", varFilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(rfile);
        printf("There are %d prior variance included in %s.\n",sigmasplit.size(),varFilename.c_str());
    }

    void read_rhofile(string rhoFilename, MatrixXd &rho)
    {
        // Read the estimated residual correlation matrix
        FILE* rfile=fopen(rhoFilename.c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",rhoFilename.c_str());
            exit(EXIT_FAILURE);
        }
        int expoNum = rho.rows();
        printf("\nReading the estimated residual correlations from %s ...\n", rhoFilename.c_str());
        char Tbuf[MAX_LINE_SIZE];
        int line_idx=0;
        vector<string> strlist;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(line_idx == 0 && atof(strlist[0].c_str())!=1 || strlist.size()!=expoNum)
            {
                printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check either the dimension of input matrix is consistent with the number of exposures or the diagnol matrix of input matrix is 1.\n", rhoFilename.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<strlist.size();i++) {
                rho(line_idx,i) = atof(strlist[i].c_str());
            }                
            line_idx++;
        }
        if(line_idx != expoNum)
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 1 analysis of OPERA! Please check.\n", rhoFilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(rfile);
        printf("The estimated correlation matrix of %d exposures are included in %s.\n",line_idx,rhoFilename.c_str());
    }

    void extract_snp_by_chr(bInfo* bdata, int chr)
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

    void extract_prob_opera(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        vector<string> raw_problist;
        for(int i=0;i<eqtlinfo->_include.size();i++) raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
        vector<string> common_probes;
        getUnique(problist);
        set_intersect(problist, raw_problist, common_probes);
        eqtlinfo->_include.clear();
        StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+problstName+"]."<<endl;
    }

    void update_ldata(lociData* ldata)
    {
        long ldata_snpNum = ldata->_include.size();
        vector<int> chr_buf, bp_buf;
        vector<string> rs_buf;
        for (int i = 0; i < ldata_snpNum; i++)
        {
            chr_buf.push_back(ldata->_chr[ldata->_include[i]]);            
            bp_buf.push_back(ldata->_bp[ldata->_include[i]]);
            rs_buf.push_back(ldata->_snp_name[ldata->_include[i]]);            
        }
        ldata->_chr.clear();        
        ldata->_bp.clear();
        ldata->_snp_name.clear();        
        ldata->_chr.swap(chr_buf);
        ldata->_bp.swap(bp_buf);
        ldata->_snp_name.swap(rs_buf);        
        ldata->_include.clear();
        for (int i = 0; i < ldata_snpNum; i++)
        {
            ldata->_include.push_back(i);
        }
    }

    void exclude_eqtl_snp_opera(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist,diffsnps;
        vector<string> mapstr;
        vector<int> tmp, incldidx;
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
        
        set_difference(mapstr, snplist, diffsnps);
        match_only(diffsnps, mapstr, incldidx);
        eqtlinfo->_esi_include.clear();   
        for(int i=0;i<incldidx.size();i++) {
            eqtlinfo->_esi_include.push_back(tmp[incldidx[i]]);
        }
        //StrFunc::set_complement(snplist, mapstr, tmp, eqtlinfo->_esi_include); //sorted
        stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());       
        cout << pre_num-eqtlinfo->_esi_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << eqtlinfo->_esi_include.size() << " SNPs remaining." << endl;
    }
    // assume esdata->_esi_include.size() == esdata->_esi_rs.size();
    long fill_smr_wk_new(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag)
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
                    if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[j]]+9>1e-6)
                    {
                        if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                            {
                                smrwk->bxz.push_back(esdata->_bxz[i][j]);
                                smrwk->sexz.push_back(esdata->_sexz[i][j]);
                                smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                                smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                                smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                                smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                                smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);
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
                            smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                            smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);
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
                if(snpchr==esdata->_epi_chr[i] && abs(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[ge_rowid]]+9>1e-6)
                {
                    if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0)
                    {
                        if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                        {
                            smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                            smrwk->sexz.push_back(esdata->_val[se_start+j]);
                            smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                            smrwk->byz.push_back(gdata->byz[gdata->_include[ge_rowid]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[ge_rowid]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[ge_rowid]]);
                            smrwk->splSize.push_back(gdata->splSize[gdata->_include[ge_rowid]]);
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
                        smrwk->byz.push_back(gdata->byz[gdata->_include[ge_rowid]]);
                        smrwk->seyz.push_back(gdata->seyz[gdata->_include[ge_rowid]]);
                        smrwk->pyz.push_back(gdata->pvalue[gdata->_include[ge_rowid]]);
                        smrwk->splSize.push_back(gdata->splSize[gdata->_include[ge_rowid]]);
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
    // extract SNPs based on esdata->_esi_include.size();
    long fill_smr_wk_include(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
        long maxid =-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<esdata->_esi_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (fabs(esdata->_sexz[i][esdata->_esi_include[j]] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[esdata->_esi_include[j]];
                    int snpchr=esdata->_esi_chr[esdata->_esi_include[j]];
                    if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[j]]+9>1e-6)
                    {
                        if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                            {
                                smrwk->bxz.push_back(esdata->_bxz[i][esdata->_esi_include[j]]);
                                smrwk->sexz.push_back(esdata->_sexz[i][esdata->_esi_include[j]]);
                                smrwk->zxz.push_back(esdata->_bxz[i][esdata->_esi_include[j]]/esdata->_sexz[i][esdata->_esi_include[j]]);
                                smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                                smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                                smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                                smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);
                                smrwk->curId.push_back(esdata->_esi_include[j]);
                                smrwk->rs.push_back(esdata->_esi_rs[esdata->_esi_include[j]]);
                                smrwk->snpchrom.push_back(esdata->_esi_chr[esdata->_esi_include[j]]);
                                smrwk->allele1.push_back(esdata->_esi_allele1[esdata->_esi_include[j]]);
                                smrwk->allele2.push_back(esdata->_esi_allele2[esdata->_esi_include[j]]);
                                if(refSNP!=NULL && esdata->_esi_rs[esdata->_esi_include[j]]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata->_esi_bp[esdata->_esi_include[j]]);
                                if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq.push_back(esdata->_esi_freq[esdata->_esi_include[j]]);
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata->_bxz[i][esdata->_esi_include[j]]/esdata->_sexz[i][esdata->_esi_include[j]]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata->_esi_rs[esdata->_esi_include[j]])+"\t"+ atos(esdata->_esi_chr[esdata->_esi_include[j]])+"\t"+ atos(esdata->_esi_bp[esdata->_esi_include[j]])+"\t"+ atos(esdata->_esi_allele1[esdata->_esi_include[j]])+"\t"+ atos(esdata->_esi_allele2[esdata->_esi_include[j]])+"\t"+ atos(esdata->_esi_freq[esdata->_esi_include[j]])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][esdata->_esi_include[j]])+"\t"+ atos(esdata->_sexz[i][esdata->_esi_include[j]])+"\t"+ dtos(p)+"\n";
                        
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            smrwk->bxz.push_back(esdata->_bxz[i][esdata->_esi_include[j]]);
                            smrwk->sexz.push_back(esdata->_sexz[i][esdata->_esi_include[j]]);
                            smrwk->zxz.push_back(esdata->_bxz[i][esdata->_esi_include[j]]/esdata->_sexz[i][esdata->_esi_include[j]]);
                            smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                            smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);
                            smrwk->curId.push_back(esdata->_esi_include[j]);
                            smrwk->rs.push_back(esdata->_esi_rs[esdata->_esi_include[j]]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[esdata->_esi_include[j]]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[esdata->_esi_include[j]]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[esdata->_esi_include[j]]);
                            if(refSNP!=NULL && esdata->_esi_rs[esdata->_esi_include[j]]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[esdata->_esi_include[j]]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                            else smrwk->freq.push_back(esdata->_esi_freq[esdata->_esi_include[j]]);

                        }
                    }
                }
            }
            
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            vector<int> ge_rowid_common, ge_rowid_tmp, val_idx;
            for(int k=0; k<numsnps; k++) ge_rowid_tmp.push_back(esdata->_rowid[beta_start+k]);
            set_intersect(ge_rowid_tmp,esdata->_esi_include,ge_rowid_common);
            numsnps = ge_rowid_common.size();
            match_only(ge_rowid_common,ge_rowid_tmp,val_idx);
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=ge_rowid_common[j];
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(snpchr==esdata->_epi_chr[i] && abs(esdata->_epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[ge_rowid]+9>1e-6)
                {
                    if(esdata->_epi_start.size()>0 && esdata->_epi_end.size()>0)
                    {
                        if(esdata->_epi_end[i]==-9 || (snpbp>esdata->_epi_end[i] && snpbp<esdata->_epi_start[i]))
                        {
                            smrwk->bxz.push_back(esdata->_val[beta_start+val_idx[j]]);
                            smrwk->sexz.push_back(esdata->_val[se_start+val_idx[j]]);
                            smrwk->zxz.push_back(esdata->_val[beta_start+val_idx[j]]/esdata->_val[se_start+val_idx[j]]);
                            smrwk->byz.push_back(gdata->byz[ge_rowid]);
                            smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                            smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                            smrwk->splSize.push_back(gdata->splSize[ge_rowid]);
                            smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                            smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                            smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                            smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                            smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                            if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                            if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[ge_rowid] / 2);
                            else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                        } else {
                            printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                            double z=(esdata->_bxz[i][ge_rowid]/esdata->_sexz[i][ge_rowid]);
                            double p=pchisq(z*z, 1);
                            string tmp=atos(esdata->_esi_rs[ge_rowid])+"\t"+ atos(esdata->_esi_chr[ge_rowid])+"\t"+ atos(esdata->_esi_bp[ge_rowid])+"\t"+ atos(esdata->_esi_allele1[ge_rowid])+"\t"+ atos(esdata->_esi_allele2[ge_rowid])+"\t"+ atos(esdata->_esi_freq[ge_rowid])+"\t"+ atos(esdata->_epi_prbID[i])+"\t"+ atos(esdata->_epi_chr[i])+"\t"+ atos(esdata->_epi_bp[i])+"\t" + atos(esdata->_epi_gene[i])+"\t"+ atos(esdata->_epi_orien[i])+"\t"+ atos(esdata->_bxz[i][ge_rowid])+"\t"+ atos(esdata->_sexz[i][ge_rowid])+"\t"+ dtos(p)+"\n";
                            
                            printf("%s\n",tmp.c_str());
                        }
                        
                    } else {
                        smrwk->bxz.push_back(esdata->_val[beta_start+val_idx[j]]);
                        smrwk->sexz.push_back(esdata->_val[se_start+val_idx[j]]);
                        smrwk->zxz.push_back(esdata->_val[beta_start+val_idx[j]]/esdata->_val[se_start+val_idx[j]]);
                        smrwk->byz.push_back(gdata->byz[ge_rowid]);
                        smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                        smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                        smrwk->splSize.push_back(gdata->splSize[ge_rowid]);
                        smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                        smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                        if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[ge_rowid] / 2);
                        else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                    }
                  
                }
            }
        }
        
        return maxid;
    }
    // allow _include; not used for multi_heidi_func
    long fill_smr_wk_mlt_include(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
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
                if (fabs(esdata[t_min]._bxz[i][esdata[t_min]._esi_include[j]] + 9) > 1e-6)
                {
                    int snpbp=esdata[t_min]._esi_bp[esdata[t_min]._esi_include[j]];
                    int snpchr=esdata[t_min]._esi_chr[esdata[t_min]._esi_include[j]];
                    if(snpchr==esdata[t_min]._epi_chr[i] && ABS(esdata[t_min]._epi_bp[i]-snpbp)<=cis_itvl && gdata->seyz[gdata->_include[j]]+9>1e-6)
                    {
                        if(esdata[t_min]._epi_start.size()>0 && esdata[t_min]._epi_end.size()>0) //technical eQTLs should be removed
                        {
                            if(esdata[t_min]._epi_end[i]==-9 || (snpbp>esdata[t_min]._epi_end[i] && snpbp<esdata[t_min]._epi_start[i]))
                            {
                                for( int t=0; t<outcoNum; t++)
                                {
                                    smrwk->bxz[t].push_back(esdata[t]._bxz[i][esdata[t_min]._esi_include[j]]);
                                    smrwk->sexz[t].push_back(esdata[t]._sexz[i][esdata[t_min]._esi_include[j]]);
                                    smrwk->zxz[t].push_back(esdata[t]._bxz[i][esdata[t_min]._esi_include[j]]/esdata[t]._sexz[i][esdata[t_min]._esi_include[j]]);
                                    if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                    else smrwk->freq[t].push_back(esdata[t]._esi_freq[esdata[t_min]._esi_include[j]]);
                                }
                                smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                                smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                                smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                                smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);                                    
                                smrwk->curId.push_back(esdata[t_min]._esi_include[j]);
                                smrwk->rs.push_back(esdata[t_min]._esi_rs[esdata[t_min]._esi_include[j]]);
                                smrwk->snpchrom.push_back(esdata[t_min]._esi_chr[esdata[t_min]._esi_include[j]]);
                                smrwk->allele1.push_back(esdata[t_min]._esi_allele1[esdata[t_min]._esi_include[j]]);
                                smrwk->allele2.push_back(esdata[t_min]._esi_allele2[esdata[t_min]._esi_include[j]]);
                                if(refSNP!=NULL && esdata[t_min]._esi_rs[esdata[t_min]._esi_include[j]]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                                smrwk->bpsnp.push_back(esdata[t_min]._esi_bp[esdata[t_min]._esi_include[j]]);
                                
                            } else {
                                printf("Shown below is the technical eQTL and will be excluded from the analysis.\n");
                                double z=(esdata[t_min]._bxz[i][esdata[t_min]._esi_include[j]]/esdata[t_min]._sexz[i][esdata[t_min]._esi_include[j]]);
                                double p=pchisq(z*z, 1);
                                string tmp=atos(esdata[t_min]._esi_rs[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._esi_chr[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._esi_bp[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._esi_allele1[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._esi_allele2[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._esi_freq[esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._epi_prbID[i])+"\t"+ atos(esdata[t_min]._epi_chr[i])+"\t"+ atos(esdata[t_min]._epi_bp[i])+"\t" + atos(esdata[t_min]._epi_gene[i])+"\t"+ atos(esdata[t_min]._epi_orien[i])+"\t"+ atos(esdata[t_min]._bxz[i][esdata[t_min]._esi_include[j]])+"\t"+ atos(esdata[t_min]._sexz[i][esdata[t_min]._esi_include[j]])+"\t"+ dtos(p)+"\n";
                                printf("%s\n",tmp.c_str());
                            }
                            
                        } else {
                            for( int t=0; t<outcoNum; t++)
                            {
                                smrwk->bxz[t].push_back(esdata[t]._bxz[i][esdata[t_min]._esi_include[j]]);
                                smrwk->sexz[t].push_back(esdata[t]._sexz[i][esdata[t_min]._esi_include[j]]);
                                smrwk->zxz[t].push_back(esdata[t]._bxz[i][esdata[t_min]._esi_include[j]]/esdata[t]._sexz[i][esdata[t_min]._esi_include[j]]);
                                if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[bdata->_include[j]] / 2); //for bdata, _include should be used with j, for others, fine.
                                else smrwk->freq[t].push_back(esdata[t]._esi_freq[esdata[t_min]._esi_include[j]]);
                            }
                            smrwk->byz.push_back(gdata->byz[gdata->_include[j]]);
                            smrwk->seyz.push_back(gdata->seyz[gdata->_include[j]]);
                            smrwk->pyz.push_back(gdata->pvalue[gdata->_include[j]]);
                            smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);                            
                            smrwk->curId.push_back(esdata[t_min]._esi_include[j]);
                            smrwk->rs.push_back(esdata[t_min]._esi_rs[esdata[t_min]._esi_include[j]]);
                            smrwk->snpchrom.push_back(esdata[t_min]._esi_chr[esdata[t_min]._esi_include[j]]);
                            smrwk->allele1.push_back(esdata[t_min]._esi_allele1[esdata[t_min]._esi_include[j]]);
                            smrwk->allele2.push_back(esdata[t_min]._esi_allele2[esdata[t_min]._esi_include[j]]);
                            if(refSNP!=NULL && esdata[t_min]._esi_rs[esdata[t_min]._esi_include[j]]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                            smrwk->bpsnp.push_back(esdata[t_min]._esi_bp[esdata[t_min]._esi_include[j]]);
                            
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
                    for(int w=0;w<common_idx.size();w++)
                    ge_rowid_common[w]=ge_rowid_first[common_idx[w]];
                    ge_rowid_first.resize(ge_rowid_common.size());
                    for(int w=0;w<ge_rowid_common.size();w++)
                    ge_rowid_first[w]=ge_rowid_common[w];
                }
                common_idx.clear();
                // extract include SNPs only
                vector<uint32_t> esi_include_tmp;
                for(int m=0; m<esdata[0]._esi_include.size();m++) esi_include_tmp.push_back(esdata[0]._esi_include[m]);
                match_only(ge_rowid_first,esi_include_tmp,common_idx);
                ge_rowid_common.resize(common_idx.size());
                #pragma omp parallel for
                for(int w=0;w<common_idx.size();w++)
                ge_rowid_common[w]=ge_rowid_first[common_idx[w]];
                // loop to extract values
                if(ge_rowid_common.size()!=0) {
                    uint64_t numsnps=ge_rowid_common.size();
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=ge_rowid_common[j];
                        smrwk->byz.push_back(gdata->byz[ge_rowid]);
                        smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                        smrwk->pyz.push_back(gdata->pvalue[ge_rowid]); 
                        smrwk->splSize.push_back(gdata->splSize[ge_rowid]);
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
                                        double z=(esdata[t]._bxz[i][ge_rowid]/esdata[t]._sexz[i][ge_rowid]);
                                        double p=pchisq(z*z, 1);
                                        // string tmp=atos(esdata[t]._esi_rs[j])+"\t"+ atos(esdata[t]._esi_chr[j])+"\t"+ atos(esdata[t]._esi_bp[j])+"\t"+ atos(esdata[t]._esi_allele1[j])+"\t"+ atos(esdata[t]._esi_allele2[j])+"\t"+ atos(esdata[t]._esi_freq[j])+"\t"+ atos(esdata[t]._epi_prbID[i])+"\t"+ atos(esdata[t]._epi_chr[i])+"\t"+ atos(esdata[t]._epi_bp[i])+"\t" + atos(esdata[t]._epi_gene[i])+"\t"+ atos(esdata[t]._epi_orien[i])+"\t"+ atos(esdata[t]._bxz[i][j])+"\t"+ atos(esdata[t]._sexz[i][j])+"\t"+ dtos(p)+"\n";                            
                                        // printf("%s\n",tmp.c_str());
                                    }
                                    
                                } else {
                                        smrwk->bxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]);
                                        smrwk->sexz[t].push_back(esdata[t]._val[se_start+val_idx[j]]);
                                        smrwk->zxz[t].push_back(esdata[t]._val[beta_start+val_idx[j]]/esdata[t]._val[se_start+val_idx[j]]);
                                        if(!heidioffFlag) smrwk->freq[t].push_back(bdata->_mu[ge_rowid] / 2);
                                        else smrwk->freq[t].push_back(esdata[t]._esi_freq[ge_rowid]);                                
                                }                 
                        }
                    }
                }                                          
            }       
        return maxid;
    } 

    void extract_smrwk_opera_new(SMRWK* smrwk,vector<int> &sn_ids,SMRWK* smrwk2)
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
            smrwk2->splSize.push_back(smrwk->splSize[sn_ids[i]]);
        }
    }

    void smr_heidi_func_opera(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg)
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

    void smr_heidi_func_para(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor,map<string, string> &prb_snp, bool targetLstFlg)
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
                
            } else maxid =fill_smr_wk(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, heidioffFlag);
        
            
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
            //pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
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

    int max_zsmr_id_opera(MTSMRWKEXP *smrwk , double p_smr)
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
    void combn_marg_pi1(long expoNum, vector<vector<int>> &combins, vector<vector<int>> &combmarg, vector<vector<int>> &idxmarg, vector<double> &pi1, vector<double> prior)
    {
        vector<vector<int>> idxall;
        for(int i=0;i<expoNum;i++) {
            vector<int> index(2);
            std::iota(index.begin(),index.end(),0);
            idxall.push_back(index);
        }
        permute_vector(idxall, combins);
        long combNum=combins.size();
        
        // list marginal PIP combins
        vector<vector<int>> comborg(combins.size());
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
        // compute pi0 correspoding configuration index
        for(int t=0;t<expoNum;t++) {
            vector<int> ppcidx;
            for(int i=1;i<combmarg.size();i++) {
                if(combmarg[i].size()==(t+1)) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        ppcidx.push_back(idxmarg[i][j]);
                    }
                }
            }
            getUnique(ppcidx);
            for(int k=0;k<ppcidx.size();k++) {
                pi1[t]+=prior[ppcidx[k]];
            }
        }
    }

    void stage2_output_file_format(string &smrfile, FILE* &smr, string &resfile, FILE* &res, string &propfile, FILE* &prop, vector<string> &ppafile, vector<FILE*> &ppa, char* gwasFileName, char* GWAScojosnplstName, bool printcombppaflag, bool printsmrflag, vector<vector<int>> &combmarg)
    {        
        long expoNum = ppafile.size(); 
        string outstr="";
        // header for .smr
        if(gwasFileName!=NULL && printsmrflag) {
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
        }
        // header for all combinations .ppa        
        if(printcombppaflag) {
            resfile = string(outFileName)+"_combo.res";
            res = fopen(resfile.c_str(), "w");
            if (!(res)) {
                printf("ERROR: open error %s\n", resfile.c_str());
                exit(1);
            }
            if(GWAScojosnplstName!=NULL) { outstr="Chr\tGWAS_SNP\tGWAS_bp\t";
            } else { outstr="Chr\t"; }
            int j = 0;
            for(int i=0; i<expoNum;i++)
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
            if(fputs_checked(outstr.c_str(),res))
            {
                printf("ERROR: in writing file %s .\n", resfile.c_str());
                exit(EXIT_FAILURE);
            }
        }        
        // header for .ppa                               
        for(int i=1;i<=expoNum;i++) {
            ppafile[i-1] =  string(outFileName)+"_"+atos(i)+"_expos_assoc.ppa";
            ppa[i-1] = fopen(ppafile[i-1].c_str(), "w");
            if (!(ppa[i-1])) {
                printf("ERROR: open error %s\n", ppafile[i-1].c_str());
                exit(1);
            }
            // output header
            if(GWAScojosnplstName!=NULL) { outstr="Chr\tGWAS_SNP\tGWAS_bp\t";
            } else { outstr="Chr\t"; }
            for(int k=1; k<=i; k++) {
                outstr+="Expo"+atos(k)+"_ID"+'\t'+"Expo"+atos(k)+"_bp"+'\t';
            }
            outstr+="PPA(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\t"; outstr+="p_HEIDI(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\n";
            if(fputs_checked(outstr.c_str(),ppa[i-1]))
            {
                printf("ERROR: error in writing file %s .\n", ppafile[i-1].c_str());
                exit(EXIT_FAILURE);
            }
        }
        // .prop file for output
        if(GWAScojosnplstName!=NULL) {
            propfile = string(outFileName)+".prop";
            prop = fopen(propfile.c_str(), "w");
            if (!(prop)) {
                printf("ERROR: open error %s\n", propfile.c_str());
                exit(1);
            }
            string propstr="Overall\t";
            for(int i=1; i<combmarg.size();i++)
            {
                propstr+="Exposures(";
                for(int j=0;j<combmarg[i].size();j++)
                {
                    propstr+=atos(combmarg[i][j]);
                    if(j<combmarg[i].size()-1) propstr+=",";
                }
                if(i==(combmarg.size()-1)) { propstr+=")\n";
                } else { propstr+=")\t"; } 
            }
            if(fputs_checked(propstr.c_str(),prop))
            {
                printf("ERROR: error in writing file %s .\n", propfile.c_str());
                exit(EXIT_FAILURE);
            }
        }        
    }
    // OPERA main function
    // OPERA 06/08/2022 
    void multiexposurepi_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string sigmastr, double sigma_def, double alpha_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag, bool jointsmrflag, int cis_itvl,int piWind, int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
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
        
        printf("\n");
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());        
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besds.size();i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening xQTL summary data file %s.Please check if the file exist or not. \n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
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

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }
        
        // 3. get prior variance if there are inputs
        long expoNum = besds.size();
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
        vector<eqtlInfo> etrait(expoNum), etrait_sig(expoNum);
        vector<eqtlInfo> esdata(expoNum);
        bInfo bdata;   gwasData gdata;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        //extract the SNP list for exposures
        for(int i=0;i<expoNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp_opera(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        //extract the probe list
        for(int i=0;i<expoNum;i++) {
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
        for(int i=0;i<expoNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the OPERA analysis.\n");
               exit(EXIT_FAILURE);
           }
           // extract probe with a xQTL
           cis_xqtl_probe_include_only(&etrait[i], p_smr, cis_itvl, besds[i]);
        }

        // read the GWAS loci data if exists
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
        }

        // 5. select the indepedent (no overlap) genomic loci for stage 1 analysis
        int exposure_probe_wind = op_wind*1000;
        int cis_itvl_wind = cis_itvl*1000;
        int indwin = piWind*1000;  // the epi_bp are required to be sorted
        vector<vector<int>> includepi(expoNum);
        long minprbnum=etrait[0]._include.size();
        int minexpnum=0;
        // find the molecular trait with minimum num. of probes
        if(GWAScojosnplstName!=NULL) {
            for(int ii=0;ii<ldata._include.size();ii++)
            {
                int probechr=ldata._chr[ii];
                int probebp=ldata._bp[ii];
                int lowerbounder=(probebp-indwin/2)>0?(probebp-indwin/2):0;
                int upperbounder=probebp+indwin/2;
               for(int i=0;i<expoNum;i++)
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

        } else {
            for(int i=0;i<expoNum;i++) {            
                if(etrait[i]._include.size() < minprbnum) {                
                    minprbnum = etrait[i]._include.size();
                    minexpnum = i;
                }
                stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
            }
            // find the num. of independent loci as the num. of no overlap probes for molecular trait with minimum probes        
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
               for(int i=0;i<expoNum;i++)
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
        }        
        // update the _include probes list for etrait
        for(int i=0;i<expoNum;i++) {
            getUnique(includepi[i]);
            printf("There are %ld probes included from exposure %ld.\n",includepi[i].size(),i+1);
            etrait[i]._include.clear();
            for(int l=0;l<includepi[i].size();l++) {
                    etrait[i]._include.push_back(includepi[i][l]);
            }
            stable_sort(etrait[i]._include.begin(),etrait[i]._include.end());
        }
        // select _esi_include SNPs list that are within 2Mb window of the target probes
        vector<vector<int>> slctsnpidx(expoNum);
        for(int i=0;i<expoNum;i++) {
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
        for(int i=0;i<expoNum;i++) {
            getUnique(slctsnpidx[i]);
            printf("There are %ld SNPs included from exposure %ld.\n",slctsnpidx[i].size(),i+1);
            etrait[i]._esi_include.clear();
            for(int m=0;m<slctsnpidx[i].size();m++) {
                    etrait[i]._esi_include.push_back(slctsnpidx[i][m]);
            }
            stable_sort(etrait[i]._esi_include.begin(),etrait[i]._esi_include.end());
        }
        // read the final besd file with updated esi_include and _include
        // for(int i=0;i<expoNum;i++) {
        //    read_besdfile(&etrait[i], string(besds[i])+".besd");
        //    if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
        //    {
        //        printf("ERROR: no data included in the analysis.\n");
        //        exit(EXIT_FAILURE);
        //    }
        // }
        
        // update the xQTL data with _esi_include and _include
        #pragma omp parallel for
        for(int i=0;i<expoNum;i++) {
            e2econvert(&etrait[i], &etrait_sig[i]);
        }
        etrait.clear();

        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata, snplstName);
                update_gwas(&gdata);
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
            allele_check_multi_opt(&bdata, etrait_sig, &gdata);
            //  allele_check_multi(&bdata, etrait_sig, &gdata);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx_opera(&bdata, etrait_sig, &gdata);
            }
        } else
        {
            allele_check_multi(etrait_sig, &gdata);
        }
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata); ngwas = median(gdata.splSize);
        }
        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<expoNum;i++) {
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
        // open the output file for sigma_b estimate
        string outstr_var = "";
        string varfile = string(outFileName)+".var";
        FILE* variter = fopen(varfile.c_str(), "w");
        outstr_var="Variance\t"; 
        
        // 8. compute the pairwise SMR effect for all exposure probes
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        for(int i=0;i<expoNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func_opera(smrrltstmp, NULL, &bdata,&gdata,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0) {                
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                // compute exposure variances with probes at P_SMR < 0.05 && abs(b_SMR) < 0.3
                vector<double> psmrVec, bsmrVec;
                for(int j=0;j<smrrltstmp.size();j++) psmrVec.push_back(smrrltstmp[j].p_SMR);
                double FDRthresh = ControlFDR(psmrVec, 0.01, true);
                // printf("The SMR P value threshold at FDR 0.01 is %s for exposure %ld.\n",atos(FDRthresh).c_str(),i+1);
                for(int j=0;j<smrrltstmp.size();j++) {
                    if(smrrltstmp[j].p_SMR <= FDRthresh) bsmrVec.push_back(smrrltstmp[j].b_SMR);
                }
                if(bsmrVec.size() > 1) {
                    double vartmp = var(bsmrVec);
                    if(vartmp>sigma_def && vartmp<0.9) sigma_b[i] = vartmp;
                }                
            }
        }
        // output the estimated variance
        for(int i=0;i<expoNum;i++) {
            outstr_var+=dtos(sigma_b[i])+"\t";
        }
        outstr_var+="\n";
        if(fputs_checked(outstr_var.c_str(),variter))
        {
            printf("ERROR: in writing file %s .\n", varfile.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(variter);
        // end
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposures with significant instruments (P < 5e-8) are less than the number of specified priors! Please check.\n");
            exit(EXIT_FAILURE);
        }

        // 9. sample a combination of exposure probes from each indepedent loci
        vector<vector<int>> includetmp(expoNum), includetmprho(expoNum);
        vector<int> NAidx(probNum[minexpnum]), NAidxrho(probNum[minexpnum]);
        long iinums = probNum[minexpnum];
        if(GWAScojosnplstName!=NULL) iinums = ldata._include.size();        
        for(int ii=0;ii<iinums;ii++)
        {
            int traitchr, traitbp;
            if(GWAScojosnplstName!=NULL) {
                traitchr = ldata._chr[ii]; traitbp=ldata._bp[ii];
            } else {
                traitchr=smrrlts[minexpnum][ii].ProbeChr; traitbp=smrrlts[minexpnum][ii].Probe_bp;
            }
            int lowerbounder=(traitbp-indwin/2)>0?(traitbp-indwin/2):0;
            int upperbounder=traitbp+indwin/2;
            NAidx[ii] = 0; NAidxrho[ii] = 0; // sample a combination of null exposure probes for estimate rho
            for(int i=0;i<expoNum;i++)
            {
               vector<int> slctprbidx, slctprbidxrho;
               for(int j=0;j<probNum[i];j++)
               {
                    int bptmp=smrrlts[i][j].Probe_bp;
                    if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                        slctprbidx.push_back(j); 
                    }
                    if(smrrlts[i][j].ProbeChr==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder && smrrlts[i][j].p_SMR>=pmecs) {
                        slctprbidxrho.push_back(j);    
                    }
               }
               if(slctprbidx.size()>0) {
                int randomIndex = rand()%slctprbidx.size();
                includetmp[i].push_back(slctprbidx[randomIndex]);
               } else { 
                    includetmp[i].push_back(-1);
                    NAidx[ii] = 1;
               }
               if(slctprbidxrho.size()>0) {
                    int randomtmp = rand()%slctprbidxrho.size();
                    includetmprho[i].push_back(slctprbidxrho[randomtmp]);
                } else { 
                    includetmprho[i].push_back(-1);
                    NAidxrho[ii] = 1;
                }
            }
        }
        // remove these loci with missing exposures
        vector<vector<int>> includesmr(expoNum), includesmrrho(expoNum);
        for(int i=0;i<expoNum;i++) {
            long lnums = includetmp[minexpnum].size();
            if(GWAScojosnplstName!=NULL) lnums = ldata._include.size();
            for(int l=0;l<lnums;l++) {
                if(NAidx[l]!=1) {
                    includesmr[i].push_back(includetmp[i][l]);
                }
            }
            if(GWAScojosnplstName==NULL) lnums = includetmprho[minexpnum].size();
            for(int l=0;l<lnums;l++) {
                if(NAidxrho[l]!=1) {
                    includesmrrho[i].push_back(includetmprho[i][l]);
                }
            }
        }
        vector<vector<SMRRLT>> smrrlts_rho;
        for(int kk=0;kk<includesmrrho[minexpnum].size();kk++) {
            vector<SMRRLT> smrrlts_rhotmp;
            for(int t=0; t<expoNum; t++) {
                long idx = includesmrrho[t][kk];
                smrrlts_rhotmp.push_back(smrrlts[t][idx]);
            }
            smrrlts_rho.push_back(smrrlts_rhotmp);
        }

        printf("\nPerforming joint SMR analysis ...\n");
        // get the joint SMR or SMR effect sizes
        vector<vector<SMRRLT>> smrrlts_joint_all;
        for(int ii=0;ii<includesmr[minexpnum].size();ii++) {
            vector<eqtlInfo> esdatabf(expoNum);
            vector<string> outconamec(expoNum);
            vector<vector<string>> prb_cojolist;
            int findcount = 0;
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
                    esdatabf[t] = esdatatmp;
                    findcount += 1; 
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
            if(findcount == expoNum) {
                vector<SMRRLT> smrrlts_joint(expoNum);
                if(!operasmrflag) {
                    // find the topsnplist and cojosnplist
                    vector<vector<string>> cojolist(expoNum); vector<string> toplist;
                    for(int p=0; p<expoNum; p++) { 
                        int idx = includesmr[p][ii];
                        toplist.push_back(smrrlts[p][idx].SNP);
                    }
                    for(int p=0; p<expoNum; p++) {
                        for(int q=0; q<expoNum; q++) {
                            int idx = includesmr[q][ii];                            
                            // if(q!=p) cojolist[p].push_back(toplist[q]);
                            // only perform COJO on colocalization probes
                            if(q!=p && smrrlts[q][idx].p_HET > 0.01) cojolist[p].push_back(toplist[q]);
                        }
                    }
                    if(targetcojosnplstName!=NULL) {
                        for(int p=0; p<expoNum; p++) {
                            for(int q=0; q<expoNum; q++) {
                                int idx = includesmr[q][ii];
                                if(q!=p) {
                                    for(int m=0; m<prb_cojolist[q].size(); m++) {
                                       // if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q])
                                       if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q] && smrrlts[q][idx].p_HET > 0.01)
                                           cojolist[p].push_back(prb_cojolist[q][m]);
                                    }
                                }
                            }
                        }
                    }                    
                    for(int p=0; p<expoNum; p++) {
                        vector<SMRRLT> smrrlts_joint_tmp;
                        if(!heidioffFlag) multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdatabf[p], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                        if(heidioffFlag) smr_heidi_func_opera(smrrlts_joint_tmp, NULL, &bdata,&gdata,&esdatabf[p],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
                        smrrlts_joint[p] = smrrlts_joint_tmp[0];
                    }
                    for(int p=0; p<expoNum; p++) {
                        int idx = includesmr[p][ii];
                        if(smrrlts_joint[p].p_GWAS <= smrrlts[p][idx].p_GWAS) {
                            smrrlts_joint[p] = smrrlts[p][idx];  
                        }
                    }
                } else { // --opera-smr effect
                    for(int p=0; p<expoNum; p++) {
                        int idx = includesmr[p][ii]; smrrlts_joint[p] = smrrlts[p][idx];
                    }
                }          
                if(smrrlts_joint.size() == expoNum) {
                    smrrlts_joint_all.push_back(smrrlts_joint);
                }            
            }
        }    

        //   compute and output sample correlations
        MatrixXd rho(expoNum, expoNum);
        for(int t=0; t<expoNum; t++) {
            vector<double> zxytmp1;
            for(int ii=0;ii<smrrlts_rho.size();ii++) zxytmp1.push_back(smrrlts_rho[ii][t].b_SMR/smrrlts_rho[ii][t].se_SMR);
            for(int k=0; k<expoNum; k++) {
                vector<double> zxytmp2;
                for(int ii=0;ii<smrrlts_rho.size();ii++) zxytmp2.push_back(smrrlts_rho[ii][k].b_SMR/smrrlts_rho[ii][k].se_SMR);
                if(sampleoverlap) {
                    if(k!=t) { rho(t,k) = cor(zxytmp1, zxytmp2);
                    } else { rho(t,k) = 1; }
                } else {
                    if(k==t) { rho(t,k) = 1;
                    } else { rho(t,k) = 0; }
                }                
            }                
        }
        if(sampleoverlap) {
            string outstr_rho = "";
            string rhofile = string(outFileName)+".rho";
            FILE* rhoiter = fopen(rhofile.c_str(), "w");
            for(int t=0;t<expoNum;t++) {
                outstr_rho = "";
                for(int i=0;i<expoNum;i++) {
                    outstr_rho+=dtos(rho(t,i))+"\t";
                }
                outstr_rho+="\n";        
                if(fputs_checked(outstr_rho.c_str(),rhoiter))
                {
                    printf("ERROR: in writing file %s .\n", rhofile.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            fclose(rhoiter);            
        }
        // end of compute the sample correlations        

        printf("\nThere are %ld independent loci included in the estimation of the global pi.\n",smrrlts_joint_all.size());
        // 10. start the MCMC sampling with the indepedent loci SMR data
        // MCMC iteration start, define variables;//
        int nloops = 10000, nburn_in = 0.2 * nloops, nsamples = nloops - nburn_in;
        MatrixXd Pr(nloops,combNum);
        VectorXd ngamma(combNum);
        VectorXd alpha(combNum);
        VectorXd Pripi(combNum);
        VectorXd Primean(combNum), Prisd(combNum);
        Stat::Dirichlet Prob;
        printf("\nMCMC sampling starts ...\n");
        // initialize the alpha value as 0.1
        double sumalpha = 0;
        for(int i=0;i<combNum;i++) {
            alpha[i] = alpha_def;
            sumalpha += alpha_def;
        }
        // output the posterior samples in log file
        printf("\n%s", logpistr.c_str());
        for(int l=0;l<nloops;l++) {
            if(l==0) { //initialize the Pripi as even prior probability
                for(int i=0;i<combNum;i++) {
                    Pripi[i] = (float)1/combNum;
                }
            }
            // if(l >= 1 && l <= nburn_in) {
            // alpha = Pr.topRows(l).colwise().mean() * sumalpha;
            // }
            // if(l >= nburn_in) {
            // alpha = Pr.topRows(nburn_in).colwise().mean() * sumalpha;
            // }
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
                VectorXd bxy_vect(expoNum), se_vect(expoNum);
                MatrixXd sigma_e_mat(expoNum, expoNum), sigma_v_mat(expoNum, expoNum);
                MatrixXd sigma_b_mat = MatrixXd::Zero(expoNum, expoNum);
                // find the summary-level joint-SMR data
                for(int t=0; t<expoNum; t++) {
                    bxy_vect[t] = smrrlts_joint_all[ii][t].b_SMR;
                    se_vect[t] = smrrlts_joint_all[ii][t].se_SMR;
                }
                sigma_e_mat = (se_vect * se_vect.transpose()).array() * rho.array();
                const double PI = 3.141592653589793238463;
                for(int i=0;i<combNum;i++) {
                    for(int t=0;t<expoNum;t++) {
                        if(combins[i][t] == 1) {
                            sigma_b_mat(t,t) = sigma_b[t];
                        } else { sigma_b_mat(t,t) = 0; }
                    }
                    sigma_v_mat = sigma_e_mat + sigma_b_mat;
                    double exp_part = bxy_vect.transpose() * sigma_v_mat.inverse() * bxy_vect;
                    HH[i] = pow(2*PI,-0.5) * pow(sigma_v_mat.determinant(),-0.5) * exp(-0.5 * exp_part);
                }
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
            // Pripi=Prob.sample(combNum,ngamma * (0.1/alpha_def));
            // Pripi = ngamma/smrrlts_joint_all.size();
            logpistr=atos(l);
            for(int c=0;c<combNum;c++) {
                Pr(l,c)=Pripi[c];
                logpistr+='\t'+atos(Pr(l,c));
            }
            logpistr=logpistr+'\n';
            // print the sampled global pi at each 100 iteration               
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

    void multiexposure_jointsmr_old(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* rhoFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, double sigma_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP,double thresh_PP_out, double thresh_smr,double thresh_gwas,double thresh_heidi, char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, bool printcombppaflag, bool printsmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
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
        
        printf("\n");
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besds.size();i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening xQTL summary data file %s. Please check if the file exist or not. \n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
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

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }
        
        // 3. expoNum = expoNum will be used; get prior variance and PIP header
        long expoNum = besds.size();
        long combNum = pow(2,expoNum);
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        }

        // get the priors (variance and pi)
        vector<string> priorsplit, sigmasplit;
        // double sigma_def = 0.02;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }
        MatrixXd rho(expoNum, expoNum);
        if(rhoFileName==NULL && sampleoverlap) throw("Error: Please input the estimated residual correlation file from the stage 1 analysis by --rho-file.");
        if(rhoFileName!=NULL) {
            read_rhofile(rhoFileName, rho);
        } else {
            for(int k=0; k<expoNum; k++) {
                for(int m=0; m<expoNum; m++) {
                    if(m==k) { rho(m,k) = 1; 
                    } else { rho(m,k) = 0; }
                }
            }            
        }        
        if(piFileName!=NULL) {
            read_pifile(piFileName, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }
        if(sigmaFileName!=NULL) {
            read_varfile(sigmaFileName, sigmasplit);
        } else {
            split_string(sigmastr, sigmasplit);
        }        
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of configurations.");
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

        // illustrate the combinations and compute pi for FPR
        vector<vector<int>> combins, combmarg, idxmarg(combNum);
        vector<double> pi1(expoNum,0), pi0(expoNum,0);
        combn_marg_pi1(expoNum, combins, combmarg, idxmarg, pi1, prior);
        for(int t=0;t<expoNum;t++) { pi0[t] = 1 - pi1[t]; }

        // 4. define global variables and extract the snp and probe data         
        vector<eqtlInfo> etrait(expoNum); 
        vector<eqtlInfo> esdata(expoNum);
        bInfo bdata;  gwasData gdata;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        // extract the SNP list for exposures
        for(int i=0;i<expoNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp_opera(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        // extract the probe list for exposures
        for(int i=0;i<expoNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");            
            if(problstName != NULL) extract_prob_opera(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob_opera(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);            
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);            
        }
        //read the besd
        for(int i=0;i<expoNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the OPERA analysis.\n");
               exit(EXIT_FAILURE);
           }
        }
        if(GWAScojosnplstName==NULL) printf("No input file for gwas loci information. The opera analysis will be performed on the genome-wide scale exposure combinations.\n");
        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata, snplstName);
                update_gwas(&gdata);
            } 
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // read the GWAS cojo signals if file exist
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
            if(prbchr!=0) extract_ldata_by_chr(&ldata,prbchr);
            update_ldata(&ldata);
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
            allele_check_multi_opt(&bdata, etrait, &gdata);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx_opera(&bdata, etrait, &gdata);
            }
        } else
        {
            allele_check_multi(etrait,&gdata);
        }
        // update the SNPs after allele checking
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata); ngwas = median(gdata.splSize);
        }
        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<expoNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        // save the bdata_include in bdata_include_org
        vector<int> bdata_include_org;
        for(int b=0;b<bdata._include.size();b++) {
            bdata_include_org.push_back(bdata._include[b]);
        }
        // 6. open .smr, .res, .ppa, and .prop for writing output
        long itemcount=0,itemcountsmr=0,itercountmlt=0,itercounttest=0;        
        vector<long> itermcount(expoNum,0); vector<double> ppasum(expoNum,0);
        string smrfile, resfile, propfile; vector<string> ppafile(expoNum);
        FILE* smr; FILE* res; FILE* prop; vector<FILE*> ppa(expoNum);
        stage2_output_file_format(smrfile, smr, resfile, res, propfile, prop, ppafile, ppa, gwasFileName, GWAScojosnplstName, printcombppaflag, printsmrflag, combmarg);        
        string outstr="";

        // 7. compute the pairwise SMR effect for all exposure probes
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        float thresh_heidi_filter = 1e-2; if(thresh_heidi_filter > thresh_heidi) thresh_heidi_filter = thresh_heidi;
        map<string, double> hdirlts;
        int ldataLineNum = 0;
        for(int i=0;i<expoNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func_opera(smrrltstmp, NULL, &bdata,&gdata,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                for(int j=0;j<smrrltstmp.size();j++)
                {                    
                    hdirlts.insert(pair<string, double>(smrrltstmp[j].ProbeID,smrrltstmp[j].p_HET));
                    if(GWAScojosnplstName==NULL) {
                        if(smrrltstmp[j].p_SMR<=thresh_smr && smrrltstmp[j].p_HET>=thresh_heidi_filter) {
                            ldata._chr.push_back(smrrltstmp[j].ProbeChr); ldata._bp.push_back(smrrltstmp[j].Probe_bp); // ldata._bp.push_back(smrrltstmp[j].SNP_bp);
                            ldata._snp_name.push_back(smrrltstmp[j].ProbeID); ldata._include.push_back(ldataLineNum);
                            ldataLineNum++;
                        }
                        outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                        if(fputs_checked(outstr.c_str(),smr))
                        {
                            printf("ERROR: in writing file %s .\n", smrfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                }            
            }
        }
        if(GWAScojosnplstName==NULL) {
            for(int k=0;k<probNum.size();k++)
            {   
                itemcount = itemcount + probNum[k]; 
            }
            printf("SMR analysis results of %ld exposure probes have been saved in the file %s .\n",itemcount,smrfile.c_str());
            fclose(smr);
        }

        printf("\nPerforming multi-exposure OPERA analysis (including multi-exposure HEIDI tests) ... \n");
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposure probes with significant instrument are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        }

        // 8. loop with each SMR < 0.05 loci and test all possible combinations at each probe loci;
        double cr=0;  long PPANum = combmarg.size() - 1;
        map<string, long> combname_set; vector<string> combname_str;
        vector<string> locisnps; vector<vector<string>> ppalocisnps(PPANum);
        vector<vector<string>> Namevec, HEIDIvec; vector<vector<float>> PPAvec;
        map<string, int> smrrltsprbs; // output unique .smr results
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
            vector<long> probNumbf(expoNum,0);
            vector<long> expoNumbf; //for missing exposures
                
            int locichr=ldata._chr[ldata._include[ii]];
            string locisnp=ldata._snp_name[ldata._include[ii]];
            int locibp=ldata._bp[ldata._include[ii]];
            int lowerbounder=(locibp-exposure_probe_wind)>0?(locibp-exposure_probe_wind):0;
            int upperbounder=locibp+exposure_probe_wind;
            // find all the probes across exposures within the window
            for(int i=0;i<expoNum;i++)
            {
                int countNum = 0;
                for(int j=0;j<probNum[i];j++)
                {
                    int bptmp=smrrlts[i][j].Probe_bp; // int bptmp=smrrlts[i][j].SNP_bp;
                    if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) 
                    {
                        if(smrrlts[i][j].p_SMR<=thresh_smr && smrrlts[i][j].p_HET>=thresh_heidi_filter && smrrlts[i][j].p_GWAS<=thresh_gwas) {
                            smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                        }                        
                        // output GWAS loci pairwise SMR results                        
                        if(GWAScojosnplstName != NULL) {
                            map<string, int>::iterator prb_smr_pos;
                            prb_smr_pos = smrrltsprbs.find(smrrlts[i][j].ProbeID);
                            if(prb_smr_pos == smrrltsprbs.end()) {
                                itemcountsmr = itemcountsmr + 1;
                                outstr=smrrlts[i][j].ProbeID+'\t'+atos(smrrlts[i][j].ProbeChr)+'\t'+smrrlts[i][j].Gene+'\t'+atos(smrrlts[i][j].Probe_bp)+'\t'+smrrlts[i][j].SNP+'\t'+atos(smrrlts[i][j].SNP_Chr)+'\t'+atos(smrrlts[i][j].SNP_bp)+'\t'+smrrlts[i][j].A1+'\t'+smrrlts[i][j].A2+'\t'+atos(smrrlts[i][j].Freq)+'\t'+atos(smrrlts[i][j].b_GWAS)+'\t'+atos(smrrlts[i][j].se_GWAS)+'\t'+dtos(smrrlts[i][j].p_GWAS)+'\t'+atos(smrrlts[i][j].b_eQTL)+'\t'+atos(smrrlts[i][j].se_eQTL)+'\t'+dtos(smrrlts[i][j].p_eQTL)+'\t'+atos(smrrlts[i][j].b_SMR)+'\t'+atos(smrrlts[i][j].se_SMR)+'\t'+dtos(smrrlts[i][j].p_SMR)+'\t'+(smrrlts[i][j].p_HET >= 0 ? dtos(smrrlts[i][j].p_HET) : "NA") + '\t' + (smrrlts[i][j].nsnp > 0 ? atos(smrrlts[i][j].nsnp+1) : "NA") + '\n';
                                if(fputs_checked(outstr.c_str(),smr))
                                {
                                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                                    exit(EXIT_FAILURE);
                                }
                                smrrltsprbs.insert(pair<string, int> (smrrlts[i][j].ProbeID, itemcountsmr));                            
                            }                            
                        }                    
                    }                   
                }                        
                probNumbf[i] = countNum;                                                
            }

            // skip the GWAS loci without an exposure probe with significant instrument
            int expoNumLoci = 0, missNum = 0;
            vector<int> missexpoNum(expoNum);
            for(int i=0;i<expoNum;i++) {
                if(probNumbf[i]>0) { expoNumLoci+=1; expoNumbf.push_back(i);
                } else { missNum += 1;}
                missexpoNum[i] = missNum;
            }
            if(expoNumLoci == 0) { continue; } // skip the loci with no significant exposures

            // update esdata._esi_include/gdata._include/bdata._include with only SNPs in the cis-window; 
            // use the exposure0 bp as gold standard; gdata, bdata and esdata included SNPs are the same;
            lowerbounder=(locibp-cis_itvl*1000)>0?(locibp-cis_itvl*1000):0;
            upperbounder=locibp+cis_itvl*1000;
            // check if the number of included SNPs are equal
            bool snpNuminequal = 0;
            for(int i=0;i<expoNum;i++) {
                if(esdata[i]._snpNum != gdata.snpNum) snpNuminequal = 1;
            }
            if(bdata_include_org.size() != gdata.snpNum) snpNuminequal = 1;
            if(snpNuminequal) printf("\nWarning: the numbers of included SNPs between data sets after allele checking are inconsistent\n");
            // select the included SNPs in the region
            for(int i=0;i<expoNum;i++) esdata[i]._esi_include.clear();
            gdata._include.clear(); bdata._include.clear();
            for(int j=0;j<esdata[0]._snpNum;j++) {
                int bptmp=esdata[0]._esi_bp[j];
                if(esdata[0]._esi_chr[j]==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                   for(int i=0;i<expoNum;i++) {
                        esdata[i]._esi_include.push_back(j); 
                   } 
                   gdata._include.push_back(j);
                   bdata._include.push_back(bdata_include_org[j]);
                }
            }
            if(bdata._include.size() == 0)  { continue; }               

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
            vector<long> idxcomb_smrrltsbf(expoNumLoci,-9), idxcomb_smrrltsbf_last(expoNumLoci,-9);
            vector<eqtlInfo> esdatabf(expoNumLoci);
            // vector<vector<string>> Namevec, HEIDIvec; vector<vector<float>> PPAvec;
            for(int cc=0; cc<combines.size(); cc++)
            {
                vector<vector<string>> prb_cojolist;
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);                    
                vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum,0),PIP(combNum,0);
                vector<string> outconamec(expoNum), outcogenec(expoNum); vector<long> outcobpc(expoNum);
                MatrixXd lh(2,expoNum);          
                VectorXd bxy_vect(expoNumLoci), se_vect(expoNumLoci);
                MatrixXd sigma_e_mat(expoNumLoci, expoNumLoci), sigma_v_mat(expoNumLoci, expoNumLoci), rho_mat(expoNumLoci, expoNumLoci);
                MatrixXd sigma_b_mat = MatrixXd::Zero(expoNumLoci, expoNumLoci);          
                // get the target combination name
                long postmp = 0; string combname; vector<string> Nametmp;
                if(GWAScojosnplstName!=NULL) {
                    Nametmp.push_back(atos(locichr)); Nametmp.push_back(locisnp); Nametmp.push_back(atos(locibp));
                } else { Nametmp.push_back(atos(locichr)); }
                for(int t=0; t<expoNum; t++)
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
                    Nametmp.push_back(outconamec[t]); Nametmp.push_back(atos(outcobpc[t]));
                    if(t<(expoNum-1)) combname.append("\t");
                }
                // skip the tested combination
                map<string, long>::iterator comb_pos;
                comb_pos = combname_set.find(combname);
                if(comb_pos != combname_set.end()) { 
                    continue; 
                } else {
                    combname_set.insert(pair<string, long> (combname, itercounttest));
                    combname_str.push_back(combname);    
                }
                // find the probe in smrrltsbf and esdata                    
                postmp = 0; idxcomb_smrrltsbf.clear(); idxcomb_smrrltsbf.resize(expoNumLoci,-9);
                if(operasmrflag) jointsmrflag = false;
                int findcount = 0;
                for(int t=0; t<expoNum; t++)
                {   
                    if(probNumbf[t] > 0) {
                        int t_new = t - missexpoNum[t];
                        long idxtmp = combines[cc][t] + postmp;
                        idxcomb_smrrltsbf[t_new] = idxtmp;
                        postmp = postmp + probNumbf[t];
                        if(!heidioffFlag || jointsmrflag) {
                            // if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new]) {
                                esdata[t]._include.clear();
                                map<string, int>::iterator itt;
                                eqtlInfo esdatatmp;
                                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                                if(itt != esdata[t]._probe_name_map.end()) {
                                    esdata[t]._include.push_back(itt->second);
                                    findcount = findcount + 1;
                                    // e2econvert(&esdata[t], &esdatatmp);                                
                                    // esdatabf[t_new] = esdatatmp;
                                }
                            // }                            
                            if(targetcojosnplstName!=NULL) {
                                // find the target probe COJO signals
                                map<string, vector<string>>::iterator prb_pos;
                                prb_pos = prb_cojosnps.find(outconamec[t]);
                                vector<string> navector; navector.push_back("");
                                if(prb_pos!=prb_cojosnps.end()) {
                                    prb_cojolist.push_back(prb_pos->second);
                                    //cojolist.push_back(prb_pos->second);
                                } else { prb_cojolist.push_back(navector); }
                            }
                        }                        
                    }                                                
                }
                // idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;
                
                //if(esdatabf.size()!=0 && esdatabf.size()!= expoNumLoci) { continue; }
                if(findcount!= expoNumLoci) { continue; }

                // 8.1 perform joint-SMR analysis or extract the SMR effect                
                vector<SMRRLT> smrrlts_joint(expoNumLoci);
                if(!operasmrflag) { // compute the joint SMR effect                                            
                    // find the topsnplist and cojosnplist
                    vector<vector<string>> cojolist(expoNumLoci); vector<string> toplist;
                    for(int p=0; p<expoNumLoci; p++) { 
                        toplist.push_back(smrrltsbf[idxcomb_smrrltsbf[p]].SNP);
                    }
                    for(int p=0; p<expoNumLoci; p++) {
                        for(int q=0; q<expoNumLoci; q++) {
                            // if(q!=p) cojolist[p].push_back(toplist[q]);
                            // only perform COJO on colocalization probes
                            if(q!=p && smrrltsbf[idxcomb_smrrltsbf[q]].p_HET > thresh_heidi) cojolist[p].push_back(toplist[q]);
                        }
                    }
                    if(targetcojosnplstName!=NULL) {
                        for(int p=0; p<expoNumLoci; p++) {
                            for(int q=0; q<expoNumLoci; q++) {
                                if(q!=p) {
                                    for(int m=0; m<prb_cojolist[q].size(); m++) {
                                       // if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q])
                                       if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q] && smrrltsbf[idxcomb_smrrltsbf[q]].p_HET > thresh_heidi)
                                           cojolist[p].push_back(prb_cojolist[q][m]);
                                    }
                                }
                            }
                        }
                    }
                    // if(targetcojosnplstName!=NULL) {
                    //  multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata, esdatabf, ngwas, prb_cojolist, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    // } else { multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata, esdatabf, ngwas, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor); }
                    // for(int t=0; t<expoNum; t++) {
                    //     if(probNumbf[t] > 0) {
                    //         int p = t - missexpoNum[t];
                    //         vector<SMRRLT> smrrlts_joint_tmp;
                    //         multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdatabf[t], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    //         smrrlts_joint.push_back(smrrlts_joint_tmp[0]);
                    //     }
                    // }
                    #pragma omp parallel for
                    for(int t=0; t<expoNum; t++) {
                        if(probNumbf[t] > 0) {
                            int p = t - missexpoNum[t];
                            vector<SMRRLT> smrrlts_joint_tmp;
                            multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdata[t], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor);
                            smrrlts_joint[p] = smrrlts_joint_tmp[0];
                        }
                    }
                    // for(int p=0; p<expoNumLoci; p++) {
                    //     vector<SMRRLT> smrrlts_joint_tmp;
                    //     multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdatabf[p], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    //     smrrlts_joint.push_back(smrrlts_joint_tmp[0]);
                    // } 
                    // replace the joint effect when top SNP has been changed or p < p_smr due to extract common SNP
                    for(int p=0; p<smrrlts_joint.size(); p++) {                            
                        //if(smrrlts_joint[p].p_eQTL >= p_smr) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        //if(smrrlts_joint[p].p_eQTL >= p_smr && smrrlts_joint.size() == expoNumLoci) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        if(smrrlts_joint[p].SNP != smrrltsbf[idxcomb_smrrltsbf[p]].SNP || smrrlts_joint[p].b_eQTL != smrrltsbf[idxcomb_smrrltsbf[p]].b_eQTL || smrrlts_joint[p].p_GWAS <= smrrltsbf[idxcomb_smrrltsbf[p]].p_GWAS) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                    }
                } else { // compute the marginal SMR effect
                    // for(int s=0; s<expoNumLoci; s++) {
                    //     smrrlts_joint.push_back(smrrltsbf[idxcomb_smrrltsbf[s]]);
                    // }
                    #pragma omp parallel for
                    for(int t=0; t<expoNum; t++) {
                        if(probNumbf[t] > 0) {
                            int p = t - missexpoNum[t];
                            smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        }
                    }
                }
                // skip no joint SMR effect due to no common SNPs
                if(smrrlts_joint.size() != expoNumLoci) { continue; }                
                //get the bxy, sigma_b and sigma_e from joint-SMR
                int k_joint = 0;
                for(int t=0; t<expoNum; t++)
                {
                    if(probNumbf[t] > 0) {
                        int t_new = t - missexpoNum[t];
                        bxy_vect[t_new] = smrrlts_joint[k_joint].b_SMR;
                        se_vect[t_new] = smrrlts_joint[k_joint].se_SMR;
                        k_joint = k_joint + 1;
                    } //else {
                        //bxy_vect[t] = 0; se_vect[t] = 0;
                    //}
                }
                for(int t=0; t<expoNum; t++)
                {
                    int t_new = t - missexpoNum[t];
                    for(int m=0; m<expoNum; m++) {
                        int m_new = m - missexpoNum[m];
                        if(probNumbf[t] > 0 & probNumbf[m] > 0) {
                            rho_mat(t_new,m_new) = rho(t,m);
                        }
                    }
                }                
                sigma_e_mat = (se_vect * se_vect.transpose()).array() * rho_mat.array();
                // 8.2 compute PIP
                const double PI = 3.141592653589793238463;
                for(int i=0;i<combNum;i++) {
                    bool misscomb = false;
                    for(int t=0;t<expoNum;t++) {
                        int t_new = t - missexpoNum[t];
                        if(probNumbf[t] > 0) {
                            if(combins[i][t] == 1) {
                                sigma_b_mat(t_new,t_new) = sigma_b[t];
                            } else { sigma_b_mat(t_new,t_new) = 0; }
                        }
                        if(probNumbf[t] == 0 && combins[i][t] == 1) misscomb = true;
                    }
                    if(!misscomb) {
                        sigma_v_mat = sigma_e_mat + sigma_b_mat;
                        double exp_part = bxy_vect.transpose() * sigma_v_mat.inverse() * bxy_vect;
                        HH[i] = pow(2*PI,-0.5) * pow(sigma_v_mat.determinant(),-0.5) * exp(-0.5 * exp_part);
                    } else { HH[i] = 0; }
                }
                float POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i] = HH[i]*prior[i];
                    POall+=PO[i];
                }
                if(POall > 0) {
                    for(int i=0;i<combNum;i++) {
                        PP[i] = PO[i]/POall;
                    }
                } else { continue; }
                for(int i=0;i<combmarg.size();i++) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        PIP[i] += PP[idxmarg[i][j]];
                    }
                }
                // output the probe information in _combo.res;
                if(GWAScojosnplstName!=NULL) {
                    outstr=atos(locichr)+'\t'+locisnp+'\t'+atos(locibp)+'\t';
                } else { outstr=atos(locichr)+'\t'; }                    
                for(int t=0; t<expoNum; t++) {
                    outstr+=outconamec[t]+'\t'+atos(outcobpc[t])+'\t';
                }                
                for(int i=0;i<PIP.size();i++) { 
                    ostringstream streamObj;
                    streamObj << std::setprecision(8);
                    streamObj << PIP[i];                      
                    outstr = outstr + streamObj.str() +'\t';
                    // outstr = outstr + atos(PIP[i])+'\t';                        
                }
                bool sigflag = false, sigoutflag = false; vector<int> sigcomb;
                
                for(int i=1;i<combmarg.size();i++) {
                    if(PIP[i]>=thresh_PP) { sigflag = true; sigcomb.push_back(i); }
                    if(PIP[i]>=thresh_PP_out) { sigoutflag = true;}
                }
                // get Namevec, PPAvec, HEIDIvec for section 9;
                vector<float> PPAtmp;
                for(int i=1;i<PIP.size();i++) PPAtmp.push_back(PIP[i]);
                PPAvec.push_back(PPAtmp);                    
                Namevec.push_back(Nametmp);
                // 8.3 multi-exposure HEIDI test
                vector<vector<SMRRLT>> smrrltsheidi(combmarg.size()); vector<string> HEIDItmp;                    
                if(! heidioffFlag && sigflag) {
                    // extract the esdatabf for multi-HEIDI analyses
                    postmp = 0; idxcomb_smrrltsbf.clear(); idxcomb_smrrltsbf.resize(expoNumLoci,-9);
                    for(int t=0; t<expoNum; t++)
                    {   
                        if(probNumbf[t] > 0) {
                            int t_new = t - missexpoNum[t];
                            long idxtmp = combines[cc][t] + postmp;
                            idxcomb_smrrltsbf[t_new] = idxtmp;
                            postmp = postmp + probNumbf[t];
                            if(!heidioffFlag || jointsmrflag) {
                                if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new]) {
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
                            }                    
                        } 
                    }
                    idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;
                    // perform multi-HEIDI analyses for exposure combinations with PPA > thresh
                    for(int h=0;h<sigcomb.size();h++) {
                        vector<eqtlInfo> esdataheidi;
                        int tmpidx = sigcomb[h];
                        string prbname;
                        for(int c=0;c<combmarg[tmpidx].size();c++) {
                            int combidx = combmarg[tmpidx][c]-1;
                            for(int k=0;k<expoNumLoci;k++) {
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
                            multi_heidi_func(smrrltsheidi[tmpidx], NULL, &bdata, &gdata, esdataheidi, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                            hdirlts.insert(pair<string, double>(prbname,smrrltsheidi[tmpidx][0].p_HET));
                        }
                    }                    
                    // output heidi pvalue
                    for(int i=1;i<combmarg.size();i++) {
                        if(smrrltsheidi[i].size()>0) {
                            if(i<(combmarg.size()-1)) {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\t' ;
                            } else {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\n';}
                            HEIDItmp.push_back(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA");
                        } else {
                            if(i<(combmarg.size()-1)) {outstr=outstr + "NA" + '\t';
                            } else {outstr=outstr + "NA" + '\n';}
                            HEIDItmp.push_back("NA");
                        }
                    }                
                } else {
                    for(int i=1;i<combmarg.size();i++) {
                        if(i<(combmarg.size()-1)) { outstr=outstr + "NA" + '\t';
                        } else { outstr=outstr + "NA" + '\n'; }
                        HEIDItmp.push_back("NA");
                    }
                }
                itercounttest += 1;
                HEIDIvec.push_back(HEIDItmp);
                if(sigoutflag && printcombppaflag) {
                    itercountmlt += 1;
                    if(fputs_checked(outstr.c_str(),res))
                    {
                        printf("ERROR: in writing file %s .\n", resfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
            }                            
        }
        if(PPAvec.size() != Namevec.size() || PPAvec.size() != HEIDIvec.size()) 
        { 
            printf("\nWarning: the size of name, PPA and HEIDI vectors are inconsistent\n"); 
            exit(EXIT_FAILURE);
        }
        // 9 take average PPA output significant results
        for(int i=1;i<=expoNum;i++) {
            for(int p=0;p<PPANum;p++) {
                if(combmarg[p+1].size() == i) {                
                    // scan the results and find max and unique ppa value
                    map<string, long> outidx; // prbname + lineidx
                    map<string, double> pparlts; // prbname + sum(ppa)
                    map<string, long> pparlts_count; // prbname + ppa num count
                    for(long l=0;l<PPAvec.size();l++) {
                        string prbname;
                        if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[l][1]); // add the loci snp in
                        } else { prbname.append(Namevec[l][0]); }
                        for(int j=0;j<combmarg[p+1].size();j++) {
                            int expoIdx = combmarg[p+1][j];
                            if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[l][2*(expoIdx+1) - 1]);
                            } else { prbname.append(Namevec[l][2*expoIdx - 1]);}
                        }                        
                        map<string, double>::iterator itmp;
                        itmp = pparlts.find(prbname);
                        if(itmp == pparlts.end()) { // not found in pparlts
                            pparlts.insert(pair<string, double>(prbname, PPAvec[l][p]));
                            outidx.insert(pair<string, long>(prbname, l));
                            pparlts_count.insert(pair<string, long>(prbname, 1));
                        } else { // found in pparlts                            
                            double pparltssum;
                            pparltssum = PPAvec[l][p] + itmp->second;
                            itmp->second = pparltssum;
                            map<string, long>::iterator itmp0;
                            itmp0 = outidx.find(prbname);                                
                            if(itmp0 != outidx.end() && HEIDIvec[l][p]!="NA") { 
                                itmp0->second = l;
                            }
                            map<string, long>::iterator itmp1;
                            itmp1 = pparlts_count.find(prbname);                                
                            if(itmp1 != pparlts_count.end()) {
                                itmp1->second = itmp1->second + 1;
                            }
                        }
                    }
                    // output string
                    map<string, long>::iterator it;
                    for (it = outidx.begin(); it != outidx.end(); it++) {
                        long ii = it->second;
                        // chr or gwas loci information
                        if(GWAScojosnplstName!=NULL) { outstr = Namevec[ii][0] + '\t'+ Namevec[ii][1] + '\t' + Namevec[ii][2] + '\t';
                        } else { outstr = Namevec[ii][0] + '\t'; }
                        string prbname;
                        if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[ii][1]); // add the loci snp in
                        } else { prbname.append(Namevec[ii][0]); }
                        for(int j=0;j<combmarg[p+1].size();j++) {
                            int expoIdx = combmarg[p+1][j];
                            if(GWAScojosnplstName!=NULL) {
                                prbname.append(Namevec[ii][2*(expoIdx+1) - 1]);
                            } else {
                                prbname.append(Namevec[ii][2*expoIdx - 1]);                        
                            }
                        }
                        // take the mean ppa
                        map<string, double>::iterator itmp;
                        itmp = pparlts.find(prbname);
                        map<string, long>::iterator itmp1;
                        itmp1 = pparlts_count.find(prbname);
                        double ppamean = 0;
                        if(itmp != pparlts.end() && itmp1 != pparlts_count.end()) {
                            ppamean = itmp->second/itmp1->second;
                        }
                        if(ppamean >= thresh_PP && HEIDIvec[ii][p]!="NA" && stod(HEIDIvec[ii][p]) >= thresh_heidi) {
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(GWAScojosnplstName!=NULL) { outstr+= Namevec[ii][2*(expoIdx+1) - 1] + '\t'+ atos(Namevec[ii][2*(expoIdx+1)]) + '\t';
                                ppalocisnps[p].push_back(Namevec[ii][1]); locisnps.push_back(Namevec[ii][1]);
                                } else { outstr+= Namevec[ii][2*expoIdx - 1] + '\t'+ atos(Namevec[ii][2*expoIdx]) + '\t'; }
                            }
                            outstr+=atos(ppamean) +'\t'+ HEIDIvec[ii][p]+'\n';
                            if(fputs_checked(outstr.c_str(), ppa[i-1]))
                            {
                            printf("ERROR: in writing file %s .\n", ppafile[i-1].c_str());
                            exit(EXIT_FAILURE);
                            }
                            // ppasum += (1 - stof(PPAvec[ii][p]));
                            ppasum[i-1] += (1 - ppamean);
                            itermcount[i-1]+=1;
                        }                                                       
                    }
                }
            }
        }
        // 10. count the actual number of test for different number of exposures
        vector<int> testNum;
        for(int i=1;i<=expoNum;i++) {
            vector<string> combname_expo;
            for(int p=1;p<=combmarg.size();p++) {
                if(combmarg[p].size() == i) {
                    for(int c=0;c<combname_str.size();c++) {
                        vector<string> combname_tmp;
                        split_string(combname_str[c],combname_tmp);
                        string name_buf;
                        for(int k=0;k<combmarg[p].size();k++) {
                            int expoIdx = combmarg[p][k] - 1;
                            name_buf.append(combname_tmp[expoIdx]);
                        }
                        combname_expo.push_back(name_buf);        
                    }                    
                }
            }
            getUnique(combname_expo);
            testNum.push_back(combname_expo.size());
        }
        // 11. estimated FDR, FPR
        float fdr, fpr;
        for(int f=1;f<=expoNum;f++) {
            double nullNum = pi0[f-1] * testNum[f-1];
            printf("\nPPA results for %ld combinatorial associations between %ld exposure(s) and 1 outcome have been extracted and saved in the file %s.\n",itermcount[f-1],f,ppafile[f-1].c_str());
            if(itermcount[f-1] > 0) {
                fdr = ppasum[f-1]/itermcount[f-1];
                printf("The estimated FDR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fdr).c_str(),f);
            }
            if(nullNum > 0) {
                fpr = ppasum[f-1]/nullNum;
                printf("The estimated FPR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fpr).c_str(),f);
            }    
        }
        // 12. output proportion of GWAS loci explained by combinations
        if(GWAScojosnplstName!=NULL) {
            getUnique(locisnps); 
            vector<double> ppaprop(PPANum,0); vector<int> ppacount(PPANum,0);
            double gwasprop = 0; vector<int> matchsnpallidx, matchsnpidx;
            StrFunc::match_only(locisnps, ldata._snp_name, matchsnpallidx);
            if(ldata._include.size() > 0) gwasprop = (double)matchsnpallidx.size()/(ldata._include.size());            
            printf("\nThere are %s%% GWAS loci were detected to be associated with at least one xQTL data.\n",atos(gwasprop*100).c_str());
            for(int p=0;p<PPANum;p++) {
                getUnique(ppalocisnps[p]); matchsnpidx.clear();
                StrFunc::match_only(ppalocisnps[p], ldata._snp_name, matchsnpidx);
                if(matchsnpidx.size() > 0) {
                    ppaprop[p] = (double)matchsnpidx.size()/(ldata._include.size());
                    ppacount[p] = matchsnpidx.size();
                }
            }            
            // .prop output
            string propstr=""; propstr += atos(gwasprop); propstr += "\t";
            for(int p=0;p<PPANum;p++) {
                if(p < (PPANum-1)) { propstr += atos(ppaprop[p]) + "\t";
                } else { propstr += atos(ppaprop[p]) + "\n"; }
            }
            propstr += atos(matchsnpallidx.size()); propstr += "\t";
            for(int p=0;p<PPANum;p++) {
                if(p < (PPANum-1)) { propstr += atos(ppacount[p]) + "\t";
                } else { propstr += atos(ppacount[p]) + "\n"; }
            }
            if(fputs_checked(propstr.c_str(), prop))
            {
                printf("ERROR: in writing file %s .\n", propfile.c_str());
                exit(EXIT_FAILURE);
            }
        }        
        if(GWAScojosnplstName!=NULL) {
            fclose(smr);
            printf("\nPairwise SMR and HEIDI analyses for %ld exposure probes have been saved in the file %s.\n",itemcountsmr,smrfile.c_str());    
        }
        if(printcombppaflag) {
            fclose(res);
            printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s.\n",itercounttest,expoNum,itercountmlt,resfile.c_str());
        } else {
            printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\n",itercounttest,expoNum);
        }                
        //end
    }

    void multiexposure_jointsmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* rhoFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, double sigma_def, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,double thresh_PP,double thresh_PP_out, double thresh_smr,double thresh_gwas,double thresh_heidi, char* refSNP, bool heidioffFlag, bool jointsmrflag, bool operasmrflag, bool printcombppaflag, bool printsmrflag, int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* GWAScojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
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
        
        printf("\n");
        read_msglist(eqtlsmaslstName, besds,"xQTL summary file names");
        if(besds.size()<1) {
            printf("Less than 1 BESD files list in %s.\n",eqtlsmaslstName);
            exit(EXIT_FAILURE);
        }
        printf("%ld xQTL summary file names are included.\n",besds.size());
        FILE** fptrs = (FILE**)malloc(sizeof(FILE*) * besds.size());
        for(int i=0;i<besds.size();i++) {
            string besdFileName=besds[i]+".besd";
            fptrs[i]=fopen(besdFileName.c_str(),"rb");
            if(fptrs[i]==NULL) {
                printf("ERROR: in opening xQTL summary data file %s. Please check if the file exist or not. \n", besdFileName.c_str());
                exit(EXIT_FAILURE);
            }
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

        // check multiple bfiles input
        if(mbFileName!=NULL) {
            read_msglist(mbFileName, multi_bfiles,"PLINK bed file names");
            if(multi_bfiles.size()<1) {
                printf("Less than 1 PLINK bed file list in %s.\n",mbFileName);
                exit(EXIT_FAILURE);
            }
            printf("%ld PLINK genotype files are included.\n",multi_bfiles.size());
        }
        
        // 3. expoNum = expoNum will be used; get prior variance and PIP header
        long expoNum = besds.size();
        long combNum = pow(2,expoNum);
        printf("There are %ld exposure(s) and 1 outcome included in the OPERA analysis.\n",expoNum);
        if(expoNum < 2) {
            printf("\nWARNING: The program can not perform the OPERA analsyis with joint SMR effect because there is only one exposure included.\nThe SMR effect will be used for OPERA analysis.\n");
            operasmrflag = true;
        }

        // get the priors (variance and pi)
        vector<string> priorsplit, sigmasplit;
        // double sigma_def = 0.02;
        if(sigmastr.size() == 0) {
            for(int i=0;i<expoNum;i++) {
                sigmastr+=atos(sigma_def);
                if(i < (expoNum - 1))  sigmastr+=","; 
            }            
        }
        MatrixXd rho(expoNum, expoNum);
        if(rhoFileName==NULL && sampleoverlap) throw("Error: Please input the estimated residual correlation file from the stage 1 analysis by --rho-file.");
        if(rhoFileName!=NULL) {
            read_rhofile(rhoFileName, rho);
        } else {
            for(int k=0; k<expoNum; k++) {
                for(int m=0; m<expoNum; m++) {
                    if(m==k) { rho(m,k) = 1; 
                    } else { rho(m,k) = 0; }
                }
            }            
        }        
        if(piFileName!=NULL) {
            read_pifile(piFileName, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }
        if(sigmaFileName!=NULL) {
            read_varfile(sigmaFileName, sigmasplit);
        } else {
            split_string(sigmastr, sigmasplit);
        }        
        if(sigmasplit.size()!=expoNum)
            throw("Error: The number of input prior variances is not consistent with the number of input exposures.");
        if(priorsplit.size()!=combNum)
            throw("Error: The number of input prior probabilities is not consistent with the number of configurations.");
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

        // illustrate the combinations and compute pi for FPR
        vector<vector<int>> combins, combmarg, idxmarg(combNum);
        vector<double> pi1(expoNum,0), pi0(expoNum,0);
        combn_marg_pi1(expoNum, combins, combmarg, idxmarg, pi1, prior);
        for(int t=0;t<expoNum;t++) { pi0[t] = 1 - pi1[t]; }

        // 4. define global variables and extract the snp and probe data         
        vector<eqtlInfo> etrait(expoNum); 
        vector<eqtlInfo> esdata(expoNum);
        bInfo bdata;  gwasData gdata;
        map<string, string> prb_snp;
        bool heidiFlag=false, targetLstflg=false;
        
        printf("\nReading the xQTL summary data file ...\n");
        if((!heidioffFlag && bFileName == NULL && mbFileName == NULL) || (jointsmrflag && bFileName == NULL && mbFileName == NULL)) throw("Error: please input Plink file for SMR analysis by either the flag --bfile or --mbfile.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe works when the probes are used as either exposures dataset or outcomes.\n"<<endl;
      
        // extract the SNP list for exposures
        for(int i=0;i<expoNum;i++) {
            read_esifile(&etrait[i], string(besds[i])+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&etrait[i], snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp_opera(&etrait[i], snplst2exclde);
            if(snpchr!=0) extract_eqtl_by_chr(&etrait[i], snpchr);
        }
        // extract the probe list for exposures
        for(int i=0;i<expoNum;i++) {
            read_epifile(&etrait[i], string(besds[i])+".epi");            
            if(problstName != NULL) extract_prob_opera(&etrait[i], problstName);
            if(problst2exclde != NULL) exclude_prob(&etrait[i], problst2exclde);
            if(eproblstName != NULL ) extract_prob_opera(&etrait[i], eproblstName);
            else if(eprobe != NULL) extract_eqtl_single_probe(&etrait[i], eprobe);
            if(eproblst2exclde != NULL) exclude_prob(&etrait[i], eproblst2exclde);
            else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&etrait[i], eprobe2rm);            
            if(prbchr!=0) extract_epi_by_chr(&etrait[i],prbchr);            
        }
        //read the besd
        for(int i=0;i<expoNum;i++) {
           read_besdfile(&etrait[i], string(besds[i])+".besd");
           if(etrait[i]._rowid.empty() && etrait[i]._bxz.empty())
           {
               printf("ERROR: no data included in the OPERA analysis.\n");
               exit(EXIT_FAILURE);
           }
        }
        if(GWAScojosnplstName==NULL) printf("No input file for gwas loci information. The opera analysis will be performed on the genome-wide scale exposure combinations.\n");
        // read gwas
        if(gwasFileName!=NULL) {
            read_gwas_data(&gdata, gwasFileName);
            if (snplstName!= NULL) {
                extract_gwas_snp(&gdata, snplstName);
                update_gwas(&gdata);
            } 
        }
        // read the cojo independent SNPs for each probe
        vector<string> cojoprbs; map<string, vector<string>> prb_cojosnps;
        if(targetcojosnplstName!=NULL) {
            read_prb_cojo_snplist(targetcojosnplstName, cojoprbs, prb_cojosnps);
        }
        // read the GWAS cojo signals if file exist
        lociData ldata;
        if(GWAScojosnplstName!=NULL) {
            read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);
            if(prbchr!=0) extract_ldata_by_chr(&ldata,prbchr);
            update_ldata(&ldata);
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
            allele_check_multi_opt(&bdata, etrait, &gdata);
            if(bFileName!=NULL) read_bedfile(&bdata, string(bFileName)+".bed");
            if(mbFileName!=NULL) read_multi_bedfiles(&bdata, multi_bfiles, snp_name_per_chr);
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx_opera(&bdata, etrait, &gdata);
            }
        } else
        {
            allele_check_multi(etrait,&gdata);
        }
        // update the SNPs after allele checking
        double ngwas = 0.0;
        if(gwasFileName!=NULL)  {
            update_gwas(&gdata); ngwas = median(gdata.splSize);
        }
        // update the SNPs after allele checking
        #pragma omp parallel for
        for(int i=0;i<expoNum;i++) {
            e2econvert(&etrait[i], &esdata[i]);
        }
        // save the bdata_include in bdata_include_org
        vector<int> bdata_include_org;
        for(int b=0;b<bdata._include.size();b++) {
            bdata_include_org.push_back(bdata._include[b]);
        }
        // 6. open .smr, .res, .ppa, and .prop for writing output
        long itemcount=0,itemcountsmr=0,itercountmlt=0,itercounttest=0;        
        vector<long> itermcount(expoNum,0); vector<double> ppasum(expoNum,0);
        string smrfile, resfile, propfile; vector<string> ppafile(expoNum);
        FILE* smr; FILE* res; FILE* prop; vector<FILE*> ppa(expoNum);
        stage2_output_file_format(smrfile, smr, resfile, res, propfile, prop, ppafile, ppa, gwasFileName, GWAScojosnplstName, printcombppaflag, printsmrflag, combmarg);        
        string outstr="";

        // 7. compute the pairwise SMR effect for all exposure probes
        printf("\nPerforming SMR analysis ...\n");
        vector<vector<SMRRLT>> smrrlts;
        vector<long> probNum;
        vector<string> outconame;
        vector<int> outcochr;
        vector<string> outcogene;
        vector<int> outcobp;
        int exposure_probe_wind=op_wind*1000;
        float thresh_heidi_filter = 1e-2; if(thresh_heidi_filter > thresh_heidi) thresh_heidi_filter = thresh_heidi;
        map<string, double> hdirlts;
        int ldataLineNum = 0;
        for(int i=0;i<expoNum;i++)
        {
            vector<SMRRLT> smrrltstmp;
            if(esdata[i]._val.size()>0) {
                smr_heidi_func_opera(smrrltstmp, NULL, &bdata,&gdata,&esdata[i],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor,prb_snp,targetLstflg);
            }
            if(smrrltstmp.size()>0)
            {
                probNum.push_back(smrrltstmp.size());
                smrrlts.push_back(smrrltstmp);
                for(int j=0;j<smrrltstmp.size();j++)
                {                    
                    hdirlts.insert(pair<string, double>(smrrltstmp[j].ProbeID,smrrltstmp[j].p_HET));
                    if(GWAScojosnplstName==NULL) {
                        if(!heidioffFlag) {
                            if(smrrltstmp[j].p_SMR<=thresh_smr && smrrltstmp[j].p_HET>=thresh_heidi_filter) {
                                ldata._chr.push_back(smrrltstmp[j].ProbeChr); ldata._bp.push_back(smrrltstmp[j].SNP_bp); // ldata._bp.push_back(smrrltstmp[j].Probe_bp);
                                ldata._snp_name.push_back(smrrltstmp[j].ProbeID); ldata._include.push_back(ldataLineNum);
                                ldataLineNum++;
                            }
                        } else {
                            if(smrrltstmp[j].p_SMR<=thresh_smr) {
                                ldata._chr.push_back(smrrltstmp[j].ProbeChr); ldata._bp.push_back(smrrltstmp[j].SNP_bp); // ldata._bp.push_back(smrrltstmp[j].Probe_bp);
                                ldata._snp_name.push_back(smrrltstmp[j].ProbeID); ldata._include.push_back(ldataLineNum);
                                ldataLineNum++;
                            }
                        }                        
                        if(printsmrflag) {
                            outstr=smrrltstmp[j].ProbeID+'\t'+atos(smrrltstmp[j].ProbeChr)+'\t'+smrrltstmp[j].Gene+'\t'+atos(smrrltstmp[j].Probe_bp)+'\t'+smrrltstmp[j].SNP+'\t'+atos(smrrltstmp[j].SNP_Chr)+'\t'+atos(smrrltstmp[j].SNP_bp)+'\t'+smrrltstmp[j].A1+'\t'+smrrltstmp[j].A2+'\t'+atos(smrrltstmp[j].Freq)+'\t'+atos(smrrltstmp[j].b_GWAS)+'\t'+atos(smrrltstmp[j].se_GWAS)+'\t'+dtos(smrrltstmp[j].p_GWAS)+'\t'+atos(smrrltstmp[j].b_eQTL)+'\t'+atos(smrrltstmp[j].se_eQTL)+'\t'+dtos(smrrltstmp[j].p_eQTL)+'\t'+atos(smrrltstmp[j].b_SMR)+'\t'+atos(smrrltstmp[j].se_SMR)+'\t'+dtos(smrrltstmp[j].p_SMR)+'\t'+(smrrltstmp[j].p_HET >= 0 ? dtos(smrrltstmp[j].p_HET) : "NA") + '\t' + (smrrltstmp[j].nsnp > 0 ? atos(smrrltstmp[j].nsnp+1) : "NA") + '\n';
                            if(fputs_checked(outstr.c_str(),smr))
                            {
                                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                                exit(EXIT_FAILURE);
                            }
                        }                        
                    }
                }            
            } else { probNum.push_back(0); }
        }
        if(GWAScojosnplstName==NULL && printsmrflag) {
            for(int k=0;k<probNum.size();k++)
            {   
                itemcount = itemcount + probNum[k]; 
            }
            printf("SMR analysis results of %ld exposure probes have been saved in the file %s .\n",itemcount,smrfile.c_str());
            fclose(smr);
        }

        printf("\nPerforming multi-exposure OPERA analysis (including multi-exposure HEIDI tests) ... \n");
        if(probNum.size()!=expoNum) {
            throw("ERROR: The number of exposure probes with significant instrument are less than the number of specified priors.\n");
            exit(EXIT_FAILURE);
        }

        // 8. loop with each SMR < 0.05 loci and test all possible combinations at each probe loci;
        double cr=0;  long PPANum = combmarg.size() - 1;
        map<string, long> combname_set; vector<string> combname_str;
        vector<string> locisnps; vector<vector<string>> ppalocisnps(PPANum);
        // vector<vector<string>> Namevec, HEIDIvec; vector<vector<float>> PPAvec;
        // Bayesian average across loci
        vector<vector<string>> Namevecall; vector<string> HEIDIvecall; 
        vector<float> PPASumall; vector<long> PPACountall; vector<int> PPANumall;
        map<string, int> smrrltsprbs; // output unique .smr results
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
            vector<long> probNumbf(expoNum,0);
            vector<long> expoNumbf; //for missing exposures
                
            int locichr=ldata._chr[ldata._include[ii]];
            string locisnp=ldata._snp_name[ldata._include[ii]];
            int locibp=ldata._bp[ldata._include[ii]];
            int lowerbounder=(locibp-exposure_probe_wind)>0?(locibp-exposure_probe_wind):0;
            int upperbounder=locibp+exposure_probe_wind;
            // find all the probes across exposures within the window
            for(int i=0;i<expoNum;i++)
            {
                int countNum = 0;
                for(int j=0;j<probNum[i];j++)
                {
                    int bptmp=smrrlts[i][j].SNP_bp; // int bptmp=smrrlts[i][j].Probe_bp;
                    if(smrrlts[i][j].ProbeChr==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) 
                    {
                        if(!heidioffFlag) {
                            if(smrrlts[i][j].p_SMR<=thresh_smr && smrrlts[i][j].p_HET>=thresh_heidi_filter && smrrlts[i][j].p_GWAS<=thresh_gwas) {
                                smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                            } 
                        } else {
                            if(smrrlts[i][j].p_SMR<=thresh_smr && smrrlts[i][j].p_GWAS<=thresh_gwas) {
                                smrrltsbf.push_back(smrrlts[i][j]); countNum = countNum + 1;
                            }   
                        }                                             
                        // output GWAS loci pairwise SMR results                        
                        if(GWAScojosnplstName != NULL) {
                            map<string, int>::iterator prb_smr_pos;
                            prb_smr_pos = smrrltsprbs.find(smrrlts[i][j].ProbeID);
                            if(prb_smr_pos == smrrltsprbs.end()) {
                                if(printsmrflag) {
                                    itemcountsmr = itemcountsmr + 1;
                                    outstr=smrrlts[i][j].ProbeID+'\t'+atos(smrrlts[i][j].ProbeChr)+'\t'+smrrlts[i][j].Gene+'\t'+atos(smrrlts[i][j].Probe_bp)+'\t'+smrrlts[i][j].SNP+'\t'+atos(smrrlts[i][j].SNP_Chr)+'\t'+atos(smrrlts[i][j].SNP_bp)+'\t'+smrrlts[i][j].A1+'\t'+smrrlts[i][j].A2+'\t'+atos(smrrlts[i][j].Freq)+'\t'+atos(smrrlts[i][j].b_GWAS)+'\t'+atos(smrrlts[i][j].se_GWAS)+'\t'+dtos(smrrlts[i][j].p_GWAS)+'\t'+atos(smrrlts[i][j].b_eQTL)+'\t'+atos(smrrlts[i][j].se_eQTL)+'\t'+dtos(smrrlts[i][j].p_eQTL)+'\t'+atos(smrrlts[i][j].b_SMR)+'\t'+atos(smrrlts[i][j].se_SMR)+'\t'+dtos(smrrlts[i][j].p_SMR)+'\t'+(smrrlts[i][j].p_HET >= 0 ? dtos(smrrlts[i][j].p_HET) : "NA") + '\t' + (smrrlts[i][j].nsnp > 0 ? atos(smrrlts[i][j].nsnp+1) : "NA") + '\n';
                                    if(fputs_checked(outstr.c_str(),smr))
                                    {
                                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                                        exit(EXIT_FAILURE);
                                    }
                                    smrrltsprbs.insert(pair<string, int> (smrrlts[i][j].ProbeID, itemcountsmr));  
                                }                                                          
                            }                            
                        }                    
                    }                   
                }                        
                probNumbf[i] = countNum;                                                
            }

            // skip the GWAS loci without an exposure probe with significant instrument
            int expoNumLoci = 0, missNum = 0;
            vector<int> missexpoNum(expoNum);
            for(int i=0;i<expoNum;i++) {
                if(probNumbf[i]>0) { expoNumLoci+=1; expoNumbf.push_back(i);
                } else { missNum += 1;}
                missexpoNum[i] = missNum;
            }
            if(expoNumLoci == 0) { continue; } // skip the loci with no significant exposures

            // update esdata._esi_include/gdata._include/bdata._include with only SNPs in the cis-window; 
            // use the exposure0 bp as gold standard; gdata, bdata and esdata included SNPs are the same;
            lowerbounder=(locibp-cis_itvl*1000)>0?(locibp-cis_itvl*1000):0;
            upperbounder=locibp+cis_itvl*1000;
            // check if the number of included SNPs are equal
            bool snpNuminequal = 0;
            for(int i=0;i<expoNum;i++) {
                if(esdata[i]._snpNum != gdata.snpNum) snpNuminequal = 1;
            }
            if(bdata_include_org.size() != gdata.snpNum) snpNuminequal = 1;
            if(snpNuminequal) printf("\nWarning: the numbers of included SNPs between data sets after allele checking are inconsistent\n");
            // select the included SNPs in the region
            for(int i=0;i<expoNum;i++) esdata[i]._esi_include.clear();
            gdata._include.clear(); bdata._include.clear();
            for(int j=0;j<esdata[0]._snpNum;j++) {
                int bptmp=esdata[0]._esi_bp[j];
                if(esdata[0]._esi_chr[j]==locichr && bptmp>=lowerbounder && bptmp<=upperbounder) {
                   for(int i=0;i<expoNum;i++) {
                        esdata[i]._esi_include.push_back(j); 
                   } 
                   gdata._include.push_back(j);
                   bdata._include.push_back(bdata_include_org[j]);
                }
            }
            if(bdata._include.size() == 0)  { continue; }               

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
            vector<long> idxcomb_smrrltsbf(expoNumLoci,-9), idxcomb_smrrltsbf_last(expoNumLoci,-9);
            vector<eqtlInfo> esdatabf(expoNumLoci);
            vector<vector<string>> Namevec, HEIDIvec; vector<vector<float>> PPAvec;
            for(int cc=0; cc<combines.size(); cc++)
            {
                vector<vector<string>> prb_cojolist;
                vector<float> bxy(expoNum), sigma_e(expoNum), c(expoNum);                    
                vector<float> Pr(combNum),HH(combNum),PO(combNum),PP(combNum,0),PIP(combNum,0);
                vector<string> outconamec(expoNum), outcogenec(expoNum); vector<long> outcobpc(expoNum);
                MatrixXd lh(2,expoNum);          
                VectorXd bxy_vect(expoNumLoci), se_vect(expoNumLoci);
                MatrixXd sigma_e_mat(expoNumLoci, expoNumLoci), sigma_v_mat(expoNumLoci, expoNumLoci), rho_mat(expoNumLoci, expoNumLoci);
                MatrixXd sigma_b_mat = MatrixXd::Zero(expoNumLoci, expoNumLoci);          
                // get the target combination name
                long postmp = 0; string combname; vector<string> Nametmp;
                if(GWAScojosnplstName!=NULL) {
                    Nametmp.push_back(atos(locichr)); Nametmp.push_back(locisnp); Nametmp.push_back(atos(locibp));
                } else { Nametmp.push_back(atos(locichr)); }
                for(int t=0; t<expoNum; t++)
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
                    Nametmp.push_back(outconamec[t]); Nametmp.push_back(atos(outcobpc[t]));
                    if(t<(expoNum-1)) combname.append("\t");
                }
                // skip the tested combination
                map<string, long>::iterator comb_pos;
                comb_pos = combname_set.find(combname);
                if(comb_pos != combname_set.end()) { 
                    continue; 
                } else {
                    combname_set.insert(pair<string, long> (combname, itercounttest));
                    combname_str.push_back(combname);    
                }
                // find the probe in smrrltsbf and esdata                    
                postmp = 0; idxcomb_smrrltsbf.clear(); idxcomb_smrrltsbf.resize(expoNumLoci,-9);
                if(operasmrflag) jointsmrflag = false;
                int findcount = 0;
                for(int t=0; t<expoNum; t++)
                {   
                    if(probNumbf[t] > 0) {
                        int t_new = t - missexpoNum[t];
                        long idxtmp = combines[cc][t] + postmp;
                        idxcomb_smrrltsbf[t_new] = idxtmp;
                        postmp = postmp + probNumbf[t];
                        if(!heidioffFlag || jointsmrflag) {
                            // if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new]) {
                                esdata[t]._include.clear();
                                map<string, int>::iterator itt;
                                eqtlInfo esdatatmp;
                                itt = esdata[t]._probe_name_map.find(outconamec[t]);
                                if(itt != esdata[t]._probe_name_map.end()) {
                                    esdata[t]._include.push_back(itt->second);
                                    findcount = findcount + 1;
                                    // e2econvert(&esdata[t], &esdatatmp);                                
                                    // esdatabf[t_new] = esdatatmp;
                                }
                            // }                            
                            if(targetcojosnplstName!=NULL) {
                                // find the target probe COJO signals
                                map<string, vector<string>>::iterator prb_pos;
                                prb_pos = prb_cojosnps.find(outconamec[t]);
                                vector<string> navector; navector.push_back("");
                                if(prb_pos!=prb_cojosnps.end()) {
                                    prb_cojolist.push_back(prb_pos->second);
                                    //cojolist.push_back(prb_pos->second);
                                } else { prb_cojolist.push_back(navector); }
                            }
                        }                        
                    }                                                
                }
                // idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;
                
                //if(esdatabf.size()!=0 && esdatabf.size()!= expoNumLoci) { continue; }
                if(findcount!= expoNumLoci) { continue; }

                // 8.1 perform joint-SMR analysis or extract the SMR effect                
                vector<SMRRLT> smrrlts_joint(expoNumLoci);
                if(!operasmrflag) { // compute the joint SMR effect                                            
                    // find the topsnplist and cojosnplist
                    vector<vector<string>> cojolist(expoNumLoci); vector<string> toplist;
                    for(int p=0; p<expoNumLoci; p++) { 
                        toplist.push_back(smrrltsbf[idxcomb_smrrltsbf[p]].SNP);
                    }
                    for(int p=0; p<expoNumLoci; p++) {
                        for(int q=0; q<expoNumLoci; q++) {
                            // if(q!=p) cojolist[p].push_back(toplist[q]);
                            // only perform COJO on colocalization probes
                            if(q!=p && smrrltsbf[idxcomb_smrrltsbf[q]].p_HET > thresh_heidi) cojolist[p].push_back(toplist[q]);
                        }
                    }
                    if(targetcojosnplstName!=NULL) {
                        for(int p=0; p<expoNumLoci; p++) {
                            for(int q=0; q<expoNumLoci; q++) {
                                if(q!=p) {
                                    for(int m=0; m<prb_cojolist[q].size(); m++) {
                                       // if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q])
                                       if(prb_cojolist[q][m]!="" && prb_cojolist[q][m]!=toplist[q] && smrrltsbf[idxcomb_smrrltsbf[q]].p_HET > thresh_heidi)
                                           cojolist[p].push_back(prb_cojolist[q][m]);
                                    }
                                }
                            }
                        }
                    }
                    // if(targetcojosnplstName!=NULL) {
                    //  multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata, esdatabf, ngwas, prb_cojolist, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    // } else { multi_joint_smr_func(smrrlts_joint, NULL, &bdata, &gdata, esdatabf, ngwas, cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor); }
                    // for(int t=0; t<expoNum; t++) {
                    //     if(probNumbf[t] > 0) {
                    //         int p = t - missexpoNum[t];
                    //         vector<SMRRLT> smrrlts_joint_tmp;
                    //         multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdatabf[t], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    //         smrrlts_joint.push_back(smrrlts_joint_tmp[0]);
                    //     }
                    // }
                    #pragma omp parallel for
                    for(int t=0; t<expoNum; t++) {
                        if(probNumbf[t] > 0) {
                            int p = t - missexpoNum[t];
                            vector<SMRRLT> smrrlts_joint_tmp;
                            // multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdata[t], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor);
                            if(!heidioffFlag) multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdata[t], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor);
                            if(heidioffFlag) smr_heidi_func_opera(smrrlts_joint_tmp, NULL, &bdata,&gdata,&esdata[t],cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap=false,pmecs,minCor,prb_snp,targetLstflg);
                            smrrlts_joint[p] = smrrlts_joint_tmp[0];
                        }
                    }
                    // for(int p=0; p<expoNumLoci; p++) {
                    //     vector<SMRRLT> smrrlts_joint_tmp;
                    //     multi_joint_smr_func_v2(smrrlts_joint_tmp, NULL, &bdata, &gdata, esdatabf[p], ngwas, cojolist[p], cis_itvl, heidioffFlag,refSNP,p_hetero,ld_top,m_hetero,p_smr,threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                    //     smrrlts_joint.push_back(smrrlts_joint_tmp[0]);
                    // } 
                    // replace the joint effect when top SNP has been changed or p < p_smr due to extract common SNP
                    for(int p=0; p<smrrlts_joint.size(); p++) {                            
                        //if(smrrlts_joint[p].p_eQTL >= p_smr) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        //if(smrrlts_joint[p].p_eQTL >= p_smr && smrrlts_joint.size() == expoNumLoci) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        if(smrrlts_joint[p].SNP != smrrltsbf[idxcomb_smrrltsbf[p]].SNP || smrrlts_joint[p].b_eQTL != smrrltsbf[idxcomb_smrrltsbf[p]].b_eQTL || smrrlts_joint[p].p_GWAS <= smrrltsbf[idxcomb_smrrltsbf[p]].p_GWAS) smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                    }
                } else { // compute the marginal SMR effect
                    // for(int s=0; s<expoNumLoci; s++) {
                    //     smrrlts_joint.push_back(smrrltsbf[idxcomb_smrrltsbf[s]]);
                    // }
                    #pragma omp parallel for
                    for(int t=0; t<expoNum; t++) {
                        if(probNumbf[t] > 0) {
                            int p = t - missexpoNum[t];
                            smrrlts_joint[p] = smrrltsbf[idxcomb_smrrltsbf[p]];
                        }
                    }
                }
                // skip no joint SMR effect due to no common SNPs
                if(smrrlts_joint.size() != expoNumLoci) { continue; }                
                //get the bxy, sigma_b and sigma_e from joint-SMR
                int k_joint = 0;
                for(int t=0; t<expoNum; t++)
                {
                    if(probNumbf[t] > 0) {
                        int t_new = t - missexpoNum[t];
                        bxy_vect[t_new] = smrrlts_joint[k_joint].b_SMR;
                        se_vect[t_new] = smrrlts_joint[k_joint].se_SMR;
                        k_joint = k_joint + 1;
                    } //else {
                        //bxy_vect[t] = 0; se_vect[t] = 0;
                    //}
                }
                for(int t=0; t<expoNum; t++)
                {
                    int t_new = t - missexpoNum[t];
                    for(int m=0; m<expoNum; m++) {
                        int m_new = m - missexpoNum[m];
                        if(probNumbf[t] > 0 & probNumbf[m] > 0) {
                            rho_mat(t_new,m_new) = rho(t,m);
                        }
                    }
                }                
                sigma_e_mat = (se_vect * se_vect.transpose()).array() * rho_mat.array();
                // 8.2 compute PIP
                const double PI = 3.141592653589793238463;
                for(int i=0;i<combNum;i++) {
                    bool misscomb = false;
                    for(int t=0;t<expoNum;t++) {
                        int t_new = t - missexpoNum[t];
                        if(probNumbf[t] > 0) {
                            if(combins[i][t] == 1) {
                                sigma_b_mat(t_new,t_new) = sigma_b[t];
                            } else { sigma_b_mat(t_new,t_new) = 0; }
                        }
                        if(probNumbf[t] == 0 && combins[i][t] == 1) misscomb = true;
                    }
                    if(!misscomb) {
                        sigma_v_mat = sigma_e_mat + sigma_b_mat;
                        double exp_part = bxy_vect.transpose() * sigma_v_mat.inverse() * bxy_vect;
                        HH[i] = pow(2*PI,-0.5) * pow(sigma_v_mat.determinant(),-0.5) * exp(-0.5 * exp_part);
                    } else { HH[i] = 0; }
                }
                float POall = 0;
                for(int i=0;i<combNum;i++) {
                    PO[i] = HH[i]*prior[i];
                    POall+=PO[i];
                }
                if(POall > 0) {
                    for(int i=0;i<combNum;i++) {
                        PP[i] = PO[i]/POall;
                    }
                } else { continue; }
                for(int i=0;i<combmarg.size();i++) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        PIP[i] += PP[idxmarg[i][j]];
                    }
                }
                // output the probe information in _combo.res;
                if(GWAScojosnplstName!=NULL) {
                    outstr=atos(locichr)+'\t'+locisnp+'\t'+atos(locibp)+'\t';
                } else { outstr=atos(locichr)+'\t'; }                    
                for(int t=0; t<expoNum; t++) {
                    outstr+=outconamec[t]+'\t'+atos(outcobpc[t])+'\t';
                }                
                for(int i=0;i<PIP.size();i++) { 
                    ostringstream streamObj;
                    streamObj << std::setprecision(8);
                    streamObj << PIP[i];                      
                    outstr = outstr + streamObj.str() +'\t';
                    // outstr = outstr + atos(PIP[i])+'\t';                        
                }
                bool sigflag = false, sigoutflag = false; vector<int> sigcomb;
                
                for(int i=1;i<combmarg.size();i++) {
                    if(PIP[i]>=thresh_PP) { sigflag = true; sigcomb.push_back(i); }
                    if(PIP[i]>=thresh_PP_out) { sigoutflag = true;}
                }
                // get Namevec, PPAvec, HEIDIvec for section 8.4;
                vector<float> PPAtmp;
                for(int i=1;i<PIP.size();i++) PPAtmp.push_back(PIP[i]);
                PPAvec.push_back(PPAtmp);                    
                Namevec.push_back(Nametmp);
                // 8.3 multi-exposure HEIDI test
                vector<vector<SMRRLT>> smrrltsheidi(combmarg.size()); vector<string> HEIDItmp;                    
                if(! heidioffFlag && sigflag) {
                    // extract the esdatabf for multi-HEIDI analyses
                    postmp = 0; idxcomb_smrrltsbf.clear(); idxcomb_smrrltsbf.resize(expoNumLoci,-9);
                    for(int t=0; t<expoNum; t++)
                    {   
                        if(probNumbf[t] > 0) {
                            int t_new = t - missexpoNum[t];
                            long idxtmp = combines[cc][t] + postmp;
                            idxcomb_smrrltsbf[t_new] = idxtmp;
                            postmp = postmp + probNumbf[t];
                            if(!heidioffFlag || jointsmrflag) {
                                if(idxcomb_smrrltsbf[t_new] != idxcomb_smrrltsbf_last[t_new]) {
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
                            }                    
                        } 
                    }
                    idxcomb_smrrltsbf_last = idxcomb_smrrltsbf;
                    // perform multi-HEIDI analyses for exposure combinations with PPA > thresh
                    for(int h=0;h<sigcomb.size();h++) {
                        vector<eqtlInfo> esdataheidi;
                        int tmpidx = sigcomb[h];
                        string prbname;
                        for(int c=0;c<combmarg[tmpidx].size();c++) {
                            int combidx = combmarg[tmpidx][c]-1;
                            for(int k=0;k<expoNumLoci;k++) {
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
                            multi_heidi_func(smrrltsheidi[tmpidx], NULL, &bdata, &gdata, esdataheidi, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero, p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor);
                            hdirlts.insert(pair<string, double>(prbname,smrrltsheidi[tmpidx][0].p_HET));
                        }
                    }                    
                    // output heidi pvalue
                    for(int i=1;i<combmarg.size();i++) {
                        if(smrrltsheidi[i].size()>0) {
                            if(i<(combmarg.size()-1)) {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\t' ;
                            } else {outstr=outstr+(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA") + '\n';}
                            HEIDItmp.push_back(smrrltsheidi[i][0].p_HET >= 0 ? dtos(smrrltsheidi[i][0].p_HET) : "NA");
                        } else {
                            if(i<(combmarg.size()-1)) {outstr=outstr + "NA" + '\t';
                            } else {outstr=outstr + "NA" + '\n';}
                            HEIDItmp.push_back("NA");
                        }
                    }                
                } else {
                    for(int i=1;i<combmarg.size();i++) {
                        if(i<(combmarg.size()-1)) { outstr=outstr + "NA" + '\t';
                        } else { outstr=outstr + "NA" + '\n'; }
                        HEIDItmp.push_back("NA");
                    }
                }
                itercounttest += 1;
                HEIDIvec.push_back(HEIDItmp);
                if(sigoutflag && printcombppaflag) {
                    itercountmlt += 1;
                    if(fputs_checked(outstr.c_str(),res))
                    {
                        printf("ERROR: in writing file %s .\n", resfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
            }
            if(PPAvec.size() != Namevec.size() || PPAvec.size() != HEIDIvec.size()) 
            { 
                printf("\nWarning: the size of name, PPA and HEIDI vectors are inconsistent\n"); 
                exit(EXIT_FAILURE);
            }
            // 8.4 take average PPA across combination at a loci
            for(int i=1;i<=expoNum;i++) {
                for(int p=0;p<PPANum;p++) {
                    if(combmarg[p+1].size() == i) {                
                        // scan the results and find max and unique ppa value
                        map<string, long> outidx; // prbname + lineidx
                        map<string, double> pparlts; // prbname + sum(ppa)
                        map<string, long> pparlts_count; // prbname + ppa num count
                        for(long l=0;l<PPAvec.size();l++) {
                            string prbname;
                            if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[l][1]); // add the loci snp in
                            } else { prbname.append(Namevec[l][0]); }
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[l][2*(expoIdx+1) - 1]);
                                } else { prbname.append(Namevec[l][2*expoIdx - 1]);}
                            }                        
                            map<string, double>::iterator itmp;
                            itmp = pparlts.find(prbname);
                            if(itmp == pparlts.end()) { // not found in pparlts
                                pparlts.insert(pair<string, double>(prbname, PPAvec[l][p]));
                                outidx.insert(pair<string, long>(prbname, l));
                                pparlts_count.insert(pair<string, long>(prbname, 1));
                            } else { // found in pparlts                            
                                double pparltssum;
                                pparltssum = PPAvec[l][p] + itmp->second;
                                itmp->second = pparltssum;
                                map<string, long>::iterator itmp0;
                                itmp0 = outidx.find(prbname);
                                if(!heidioffFlag) {
                                    if(itmp0 != outidx.end() && HEIDIvec[l][p]!="NA") { 
                                        itmp0->second = l;
                                    } 
                                } else {
                                    if(itmp0 != outidx.end()) {
                                        itmp0->second = l;
                                    }
                                }                                
                                map<string, long>::iterator itmp1;
                                itmp1 = pparlts_count.find(prbname);                                
                                if(itmp1 != pparlts_count.end()) {
                                    itmp1->second = itmp1->second + 1;
                                }
                            }
                        }
                        // output string
                        map<string, long>::iterator it;
                        for (it = outidx.begin(); it != outidx.end(); it++) {
                            long ii = it->second;
                            string prbname;
                            if(GWAScojosnplstName!=NULL) { prbname.append(Namevec[ii][1]); // add the loci snp in
                            } else { prbname.append(Namevec[ii][0]); }
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(GWAScojosnplstName!=NULL) {
                                    prbname.append(Namevec[ii][2*(expoIdx+1) - 1]);
                                } else {
                                    prbname.append(Namevec[ii][2*expoIdx - 1]);                        
                                }
                            }
                            // take the mean ppa
                            map<string, double>::iterator itmp;
                            itmp = pparlts.find(prbname);
                            map<string, long>::iterator itmp1;
                            itmp1 = pparlts_count.find(prbname);
                            // double ppamean = 0;
                            // if(itmp != pparlts.end() && itmp1 != pparlts_count.end()) {
                            //     ppamean = itmp->second/itmp1->second;
                            // }
                            // chr or gwas loci information
                            vector<string> Namevectmp;
                            if(GWAScojosnplstName!=NULL) { Namevectmp.push_back(Namevec[ii][0]); Namevectmp.push_back(Namevec[ii][1]); Namevectmp.push_back(Namevec[ii][2]);
                            } else { Namevectmp.push_back(Namevec[ii][0]); }
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(GWAScojosnplstName!=NULL) { Namevectmp.push_back(Namevec[ii][2*(expoIdx+1) - 1]); Namevectmp.push_back(Namevec[ii][2*(expoIdx+1)]);
                                } else { Namevectmp.push_back(Namevec[ii][2*expoIdx - 1]); Namevectmp.push_back(Namevec[ii][2*expoIdx]); }
                            }
                            if(itmp != pparlts.end()) {
                                PPASumall.push_back(itmp->second);
                            }
                            if(itmp1 != pparlts_count.end()) {
                                PPACountall.push_back(itmp1->second);
                            }                                                        
                            Namevecall.push_back(Namevectmp);                            
                            HEIDIvecall.push_back(HEIDIvec[ii][p]);
                            PPANumall.push_back(p);                       
                        }
                    }
                }
            }

        }
        if(PPASumall.size() != PPACountall.size() || Namevecall.size() != PPANumall.size() || PPASumall.size() != HEIDIvecall.size()) 
        { 
            printf("\nWarning: the size of name, PPA and HEIDI vectors are inconsistent in across loci Bayesian average \n"); 
            exit(EXIT_FAILURE);
        }
        // 9. Bayesian average for duplicated combinations across loci
        map<string, long> outidx; // prbname + lineidx
        map<string, double> pparlts; // prbname + sum(ppa)
        map<string, long> pparlts_count; // prbname + ppa num count
        for(int l=0; l<Namevecall.size(); l++) {
            string prbname;
            for(int j=0;j<Namevecall[l].size();j++) {
                prbname.append(Namevecall[l][j]);
            }
            map<string, double>::iterator itmp;
            itmp = pparlts.find(prbname);
            if(itmp == pparlts.end()) { // not found in pparlts
                pparlts.insert(pair<string, double>(prbname, PPASumall[l]));
                outidx.insert(pair<string, long>(prbname, l));
                pparlts_count.insert(pair<string, long>(prbname, PPACountall[l]));
            } else { // found in pparlts                            
                double pparltssum;
                pparltssum = PPASumall[l] + itmp->second;
                itmp->second = pparltssum;
                map<string, long>::iterator itmp0;
                itmp0 = outidx.find(prbname);                                
                if(!heidioffFlag) {
                    if(itmp0 != outidx.end() && HEIDIvecall[l]!="NA") { 
                        itmp0->second = l;
                    } 
                } else {
                    if(itmp0 != outidx.end()) {
                        itmp0->second = l;
                    } 
                }                
                map<string, long>::iterator itmp1;
                itmp1 = pparlts_count.find(prbname);                                
                if(itmp1 != pparlts_count.end()) {
                    itmp1->second =  PPACountall[l] + itmp1->second;
                }
            }
        }
        // output string; loop with unique lines
        map<string, long>::iterator it;
        for (it = outidx.begin(); it != outidx.end(); it++) {
            long ii = it->second;                        
            string prbname;
            for(int j=0;j<Namevecall[ii].size();j++) {
                prbname.append(Namevecall[ii][j]);
            }
            int expofile = 0; 
            if(GWAScojosnplstName!=NULL) { expofile = (Namevecall[ii].size() - 3)/2;
            } else { expofile = (Namevecall[ii].size() - 1)/2; }
            // take the mean ppa
            map<string, double>::iterator itmp;
            itmp = pparlts.find(prbname);
            map<string, long>::iterator itmp1;
            itmp1 = pparlts_count.find(prbname);
            double ppamean = 0;
            if(itmp != pparlts.end() && itmp1 != pparlts_count.end()) {
                ppamean = itmp->second/itmp1->second;
            }
            outstr="";
            if(!heidioffFlag) {
                if(ppamean >= thresh_PP && HEIDIvecall[ii]!="NA" && stod(HEIDIvecall[ii]) >= thresh_heidi) {
                    for(int j=0;j<Namevecall[ii].size();j++) {
                        outstr += Namevecall[ii][j] + '\t';
                    }
                    if(GWAScojosnplstName!=NULL) {
                        ppalocisnps[PPANumall[ii]].push_back(Namevecall[ii][1]); locisnps.push_back(Namevecall[ii][1]);
                    }                
                    outstr+=atos(ppamean) +'\t'+ HEIDIvecall[ii]+'\n';
                    if(fputs_checked(outstr.c_str(), ppa[expofile-1]))
                    {
                    printf("ERROR: in writing file %s .\n", ppafile[expofile-1].c_str());
                    exit(EXIT_FAILURE);
                    }
                    // ppasum += (1 - stof(PPAvec[ii][p]));
                    ppasum[expofile-1] += (1 - ppamean);
                    itermcount[expofile-1]+=1;
                }
            } else {
                if(ppamean >= thresh_PP) {
                    for(int j=0;j<Namevecall[ii].size();j++) {
                        outstr += Namevecall[ii][j] + '\t';
                    }
                    if(GWAScojosnplstName!=NULL) {
                        ppalocisnps[PPANumall[ii]].push_back(Namevecall[ii][1]); locisnps.push_back(Namevecall[ii][1]);
                    }                
                    outstr+=atos(ppamean) +'\t'+ HEIDIvecall[ii]+'\n';
                    if(fputs_checked(outstr.c_str(), ppa[expofile-1]))
                    {
                    printf("ERROR: in writing file %s .\n", ppafile[expofile-1].c_str());
                    exit(EXIT_FAILURE);
                    }
                    // ppasum += (1 - stof(PPAvec[ii][p]));
                    ppasum[expofile-1] += (1 - ppamean);
                    itermcount[expofile-1]+=1;
                }
            }
            
        }
        // 10. count the actual number of test for different number of exposures
        vector<int> testNum;
        for(int i=1;i<=expoNum;i++) {
            vector<string> combname_expo;
            for(int p=1;p<=combmarg.size();p++) {
                if(combmarg[p].size() == i) {
                    for(int c=0;c<combname_str.size();c++) {
                        vector<string> combname_tmp;
                        split_string(combname_str[c],combname_tmp);
                        string name_buf;
                        for(int k=0;k<combmarg[p].size();k++) {
                            int expoIdx = combmarg[p][k] - 1;
                            name_buf.append(combname_tmp[expoIdx]);
                        }
                        combname_expo.push_back(name_buf);        
                    }                    
                }
            }
            getUnique(combname_expo);
            testNum.push_back(combname_expo.size());
        }
        // 11. estimated FDR, FPR
        float fdr, fpr;
        for(int f=1;f<=expoNum;f++) {
            double nullNum = pi0[f-1] * testNum[f-1];
            printf("\nPPA results for %ld combinatorial associations between %ld exposure(s) and 1 outcome have been extracted and saved in the file %s.\n",itermcount[f-1],f,ppafile[f-1].c_str());
            if(itermcount[f-1] > 0) {
                fdr = ppasum[f-1]/itermcount[f-1];
                printf("The estimated FDR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fdr).c_str(),f);
            }
            if(nullNum > 0) {
                fpr = ppasum[f-1]/nullNum;
                printf("The estimated FPR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fpr).c_str(),f);
            }    
        }
        // 12. output proportion of GWAS loci explained by combinations
        if(GWAScojosnplstName!=NULL) {
            getUnique(locisnps); 
            vector<double> ppaprop(PPANum,0); vector<int> ppacount(PPANum,0);
            double gwasprop = 0; vector<int> matchsnpallidx, matchsnpidx;
            StrFunc::match_only(locisnps, ldata._snp_name, matchsnpallidx);
            if(ldata._include.size() > 0) gwasprop = (double)matchsnpallidx.size()/(ldata._include.size());            
            printf("\nThere are %s%% GWAS loci were detected to be associated with at least one xQTL data.\n",atos(gwasprop*100).c_str());
            for(int p=0;p<PPANum;p++) {
                getUnique(ppalocisnps[p]); matchsnpidx.clear();
                StrFunc::match_only(ppalocisnps[p], ldata._snp_name, matchsnpidx);
                if(matchsnpidx.size() > 0) {
                    ppaprop[p] = (double)matchsnpidx.size()/(ldata._include.size());
                    ppacount[p] = matchsnpidx.size();
                }
            }            
            // .prop output
            string propstr=""; propstr += atos(gwasprop); propstr += "\t";
            for(int p=0;p<PPANum;p++) {
                if(p < (PPANum-1)) { propstr += atos(ppaprop[p]) + "\t";
                } else { propstr += atos(ppaprop[p]) + "\n"; }
            }
            propstr += atos(matchsnpallidx.size()); propstr += "\t";
            for(int p=0;p<PPANum;p++) {
                if(p < (PPANum-1)) { propstr += atos(ppacount[p]) + "\t";
                } else { propstr += atos(ppacount[p]) + "\n"; }
            }
            if(fputs_checked(propstr.c_str(), prop))
            {
                printf("ERROR: in writing file %s .\n", propfile.c_str());
                exit(EXIT_FAILURE);
            }
        }        
        if(GWAScojosnplstName!=NULL && printsmrflag) {
            fclose(smr);
            printf("\nPairwise SMR and HEIDI analyses for %ld exposure probes have been saved in the file %s.\n",itemcountsmr,smrfile.c_str());    
        }
        if(printcombppaflag) {
            fclose(res);
            printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\nPosterior probability and HEIDI results of %ld combinations have been saved in the file %s.\n",itercounttest,expoNum,itercountmlt,resfile.c_str());
        } else {
            printf("\nOPERA analyses for %ld combinations between %ld exposures and 1 outcome completed.\n",itercounttest,expoNum);
        }                
        //end
    }

    void multioutcomesmr(char* outFileName, char* bFileName, char* mbFileName, char* piFileName, char* sigmaFileName, char* eqtlFileName, char* eqtlsmaslstName, char* gwasFileName, double maf,string priorstr,string sigmastr, char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, int opt_hetero,  char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,int snpchr,int prbchr,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double threshpsmrest, bool new_het_mth, bool opt, double ld_min,bool cis2all, bool sampleoverlap, double pmecs, int minCor, char* targetcojosnplstName, char* snpproblstName,double afthresh,double percenthresh)
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
        if(snplst2exclde != NULL) exclude_eqtl_snp_opera(&esdata, snplst2exclde);
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
            if(snplst2exclde != NULL) exclude_eqtl_snp_opera(&etrait[i], snplst2exclde);
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
                    update_geIndx_opera(&bdata, etrait, &gdata1, &esdata);
                } else {update_geIndx_opera(&bdata, etrait, &esdata);}
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
                smr_heidi_func_opera(smrrlts1, NULL, &bdata, &gdata, &esdata, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);                
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
            smr_heidi_func_opera(smrrltsgwas, NULL, &bdata,&gdata1,&esdata, cis_itvl, heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest,new_het_mth,opt,ld_min,opt_hetero,sampleoverlap,pmecs,minCor,prb_snp,targetLstflg);
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

    void ppafile_summary_old1(char* ppaFilename, char* numFilename,char* piFilename,string priorstr,char* outFileName, double thresh_PP, double thresh_heidi)
    {
        // Read posterior probabilities for associations file
        if(ppaFilename==NULL) {
            printf("There is no .ppa file input from the stage 2 analysis, please specify it with --ppa-file.\n");
            exit(EXIT_FAILURE);
        }
        if(numFilename==NULL) {
            printf("There is no .num file input from the stage 2 analysis, please specify it with --num-file.\n");
            exit(EXIT_FAILURE);
        }
        if(piFilename==NULL && priorstr=="") {
            printf("There is no prior probability input from the stage 1 analysis, please specify it with --prior-pi-file or --prior-pi.\n");
            exit(EXIT_FAILURE);
        }
        if(outFileName==NULL) {
            printf("There is no output file name specified, please check!");
            exit(EXIT_FAILURE);
        }
        // ppa file
        FILE* rfile=fopen(string(ppaFilename).c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        printf("\nReading the posterior probabilities for associations (PPA) from %s ......\n", string(ppaFilename).c_str());
        char Tbuf[MAX_LINE_SIZE];    
        vector<string> strlist;
        int expoNum, PPANum;
        vector<int> cidx_PPA, cidx_HEIDI;
        bool gwaslociflag = false;
        // header line
        fgets(Tbuf, MAX_LINE_SIZE, rfile);
        split_string(Tbuf, strlist, " \t\n");
        int strsize = strlist.size();
        if(strlist[0]!="Chr") 
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<strlist.size();i++) {
            bool PPAflag = has_prefix(strlist[i], "PPA");
            bool HEIDIflag = has_prefix(strlist[i], "p_HEIDI");
            bool tmpflag = has_prefix(strlist[i], "GWAS_SNP");
            if(PPAflag) cidx_PPA.push_back(i);
            if(HEIDIflag) cidx_HEIDI.push_back(i);
            if(tmpflag) gwaslociflag = true;
        }
        if(cidx_PPA.size()==0 || cidx_HEIDI.size()==0)
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        PPANum = cidx_PPA.size() - 1;
        if(gwaslociflag) { expoNum = (cidx_PPA[0] - 3)/2;
        } else { expoNum = (cidx_PPA[0] - 1)/2; }
        vector<vector<string>> PPAvec, HEIDIvec, Namevec;
        
        // read prior pi file
        vector<string> priorsplit; vector<double> prior;
        if(piFilename!=NULL) {
            read_pifile(piFilename, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }
        // read_pifile(piFilename, priorsplit);
        if(priorsplit.size()!=(PPANum+1))
            throw("Error: The number of input prior probabilities is not consistent with the number of test configurations.");
        for(int t=0; t<(PPANum+1); t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }

        // read num file
        vector<long> numsum(expoNum, 0);
        read_numfile(numFilename, numsum);

        // scan the ppa file
        long line_idx = 0;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            vector<string> PPAtmp, HEIDItmp, Nametmp;
            split_string(Tbuf, strlist, ", \t\n");
            if(strlist.size() != strsize) {
                printf("ERROR: The %ld row from the input file %s is incomplete! Please check.\n", line_idx,  string(ppaFilename).c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<strlist.size();i++) {
                if(i < cidx_PPA[0]) {
                    Nametmp.push_back(strlist[i]);
                } else if(i >= cidx_PPA[1] && i < cidx_HEIDI[0]) {
                    PPAtmp.push_back(strlist[i]);
                } else if(i >= cidx_HEIDI[0]) {
                    HEIDItmp.push_back(strlist[i]);
                }
            }
            Namevec.push_back(Nametmp);
            PPAvec.push_back(PPAtmp);
            HEIDIvec.push_back(HEIDItmp);
            line_idx++;        
        }
        fclose(rfile);
        printf("There are %d association combinations between %ld exposures and 1 outcome included in %s.\n",line_idx, expoNum, string(ppaFilename).c_str());
        // get the corresponding marginal index
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
        // compute pi0 correspoding configuration index
        vector<double> pi1(expoNum,0), pi0(expoNum,0);
        for(int t=0;t<expoNum;t++) {
            vector<int> ppcidx;
            for(int i=1;i<combmarg.size();i++) {
                if(combmarg[i].size()==(t+1)) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        ppcidx.push_back(idxmarg[i][j]);
                    }
                }
            }
            getUnique(ppcidx);
            for(int k=0;k<ppcidx.size();k++){
                pi1[t]+=prior[ppcidx[k]];
            }
        }
        for(int t=0;t<expoNum;t++) {
            pi0[t] = 1 - pi1[t];
        }
        // extract the significant results and remove duplicated
        // loop with the expoNum      
        for(int i=1;i<=expoNum;i++) {
            double nullNum = pi0[i-1] * numsum[i-1];
            // open the output file
            string ppafile0 = string(outFileName)+"_"+atos(i)+"_exposures_ppa.summary";
            FILE* ppa0;
            ppa0 = fopen(ppafile0.c_str(), "w");
            if (!(ppa0)) {
                printf("ERROR: open error %s\n", ppafile0.c_str());
                exit(1);
            }
            // output header
            string outstr;
            if(gwaslociflag) { outstr="Chr\tGWAS_SNP\tGWAS_bp\t";
            } else { outstr="Chr\t"; }
            for(int k=1; k<=i; k++) {
                outstr+="Expo"+atos(k)+"_ID"+'\t'+"Expo"+atos(k)+"_bp"+'\t';
            }
            outstr+="PPA(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\t"; outstr+="p_HEIDI(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\n";
            if(fputs_checked(outstr.c_str(),ppa0))
            {
                printf("ERROR: error in writing file %s .\n", ppafile0.c_str());
                exit(EXIT_FAILURE);
            }
            long itermcount = 0;
            double ppasum = 0, fdr = 0, fpr = 0;
            for(int p=0;p<PPANum;p++) {
                if(combmarg[p+1].size() == i) {                
                    // scan the results and find max and unique ppa value
                    map<string, long> outidx; // prbname + lineidx
                    map<string, double> pparlts; // prbname + sum(ppa)
                    map<string, long> pparlts_count; // prbname + ppa num count
                    for(long l=0;l<line_idx;l++) {
                        //if(HEIDIvec[l][p]!="NA" && stof(PPAvec[l][p]) >= thresh_PP && stod(HEIDIvec[l][p]) >= thresh_heidi) {
                            string prbname;
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(gwaslociflag) { prbname.append(Namevec[l][2*(expoIdx+1) - 1]);
                                } else { prbname.append(Namevec[l][2*expoIdx - 1]);}
                            }                        
                            map<string, double>::iterator itmp;
                            itmp = pparlts.find(prbname);
                            if(itmp == pparlts.end()) { // not found in pparlts
                                pparlts.insert(pair<string, double>(prbname, stof(PPAvec[l][p])));
                                outidx.insert(pair<string, long>(prbname, l));
                                pparlts_count.insert(pair<string, long>(prbname, 1));
                            } else { // found in pparlts                            
                                double pparltssum;
                                pparltssum = stof(PPAvec[l][p]) + itmp->second;
                                itmp->second = pparltssum;
                                map<string, long>::iterator itmp0;
                                itmp0 = outidx.find(prbname);
                                if(itmp0 != outidx.end() && HEIDIvec[l][p]!="NA") {                                
                                //if(itmp0 != outidx.end()) {
                                    itmp0->second = l;
                                }
                                map<string, long>::iterator itmp1;
                                itmp1 = pparlts_count.find(prbname);                                
                                if(itmp1 != pparlts_count.end()) {
                                    itmp1->second = itmp1->second + 1;
                                }
                            }
                        //}
                    }
                    // output string
                    map<string, long>::iterator it;
                    for (it = outidx.begin(); it != outidx.end(); it++) {
                        long ii = it->second;
                        // file header
                        if(gwaslociflag) { outstr = Namevec[ii][0] + '\t'+ Namevec[ii][1] + '\t' + Namevec[ii][2] + '\t';
                        } else { outstr = Namevec[ii][0] + '\t'; }
                        string prbname;
                        for(int j=0;j<combmarg[p+1].size();j++) {
                            int expoIdx = combmarg[p+1][j];
                            if(gwaslociflag) { prbname.append(Namevec[ii][2*(expoIdx+1) - 1]);
                            } else { prbname.append(Namevec[ii][2*expoIdx - 1]);}
                        }
                        // take the mean ppa
                        map<string, double>::iterator itmp;
                        itmp = pparlts.find(prbname);
                        map<string, long>::iterator itmp1;
                        itmp1 = pparlts_count.find(prbname);
                        double ppamean = 0;
                        if(itmp != pparlts.end() && itmp1 != pparlts_count.end()) {
                            ppamean = itmp->second/itmp1->second;
                        }
                        if(ppamean >= thresh_PP && HEIDIvec[ii][p]!="NA" && stod(HEIDIvec[ii][p]) >= thresh_heidi) {
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(gwaslociflag) { outstr+= Namevec[ii][2*(expoIdx+1) - 1] + '\t'+ atos(Namevec[ii][2*(expoIdx+1)]) + '\t';
                                } else { outstr+= Namevec[ii][2*expoIdx - 1] + '\t'+ atos(Namevec[ii][2*expoIdx]) + '\t'; }
                            }
                            // outstr+=PPAvec[ii][p] +'\t'+ HEIDIvec[ii][p]+'\n';
                            outstr+=atos(ppamean) +'\t'+ HEIDIvec[ii][p]+'\n';
                            if(fputs_checked(outstr.c_str(), ppa0))
                            {
                            printf("ERROR: in writing file %s .\n", ppafile0.c_str());
                            exit(EXIT_FAILURE);
                            }
                            // ppasum += (1 - stof(PPAvec[ii][p]));
                            ppasum += (1 - ppamean);
                            itermcount+=1;
                        }                                                       
                    }
                }
            }
            fclose(ppa0);
            printf("\nPPA results for %ld combinatorial associations between %ld exposure(s) and 1 outcome have been extracted and saved in the file %s.\n",itermcount,i,ppafile0.c_str());
            if(itermcount > 0) {
                fdr = ppasum/itermcount;
                printf("The estimated FDR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fdr).c_str(),i);
            }
            if(nullNum > 0) {
                fpr = ppasum/nullNum;
                printf("The estimated FPR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fpr).c_str(),i);
            }
        }
    }

    void ppafile_summary(char* ppaFilename,char* GWAScojosnplstName,char* numFilename,char* piFilename,string priorstr,char* outFileName, double thresh_PP, double thresh_heidi)
    {
        // Read posterior probabilities for associations file
        if(ppaFilename==NULL) {
            printf("There is no .ppa file input from the stage 2 analysis, please specify it with --ppa-file.\n");
            exit(EXIT_FAILURE);
        }
        if(numFilename==NULL) {
            printf("There is no .num file input from the stage 2 analysis, please specify it with --num-file.\n");
            exit(EXIT_FAILURE);
        }
        if(piFilename==NULL && priorstr=="") {
            printf("There is no prior probability input from the stage 1 analysis, please specify it with --prior-pi-file or --prior-pi.\n");
            exit(EXIT_FAILURE);
        }
        if(outFileName==NULL) {
            printf("There is no output file name specified, please check!");
            exit(EXIT_FAILURE);
        }
        // ppa file
        FILE* rfile=fopen(string(ppaFilename).c_str(),"r");
        if(!rfile) {
            printf("File %s open failed.\n",string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        printf("\nReading the posterior probabilities for associations (PPA) from %s ......\n", string(ppaFilename).c_str());
        char Tbuf[MAX_LINE_SIZE];    
        vector<string> strlist, strheader;
        long expoNum, PPANum;
        vector<int> cidx_PPA, cidx_HEIDI;
        bool gwaslociflag = false;
        // header line
        fgets(Tbuf, MAX_LINE_SIZE, rfile);
        split_string(Tbuf, strlist, " \t\n");
        long strsize = strlist.size();
        if(strlist[0]!="Chr") 
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<strlist.size();i++) {
            bool PPAflag = has_prefix(strlist[i], "PPA");
            bool HEIDIflag = has_prefix(strlist[i], "p_HEIDI");
            bool tmpflag = has_prefix(strlist[i], "GWAS_SNP");
            if(PPAflag) cidx_PPA.push_back(i);
            if(HEIDIflag) cidx_HEIDI.push_back(i);
            if(tmpflag) gwaslociflag = true;
        }
        strheader = strlist;
        if(cidx_PPA.size()==0 || cidx_HEIDI.size()==0)
        {
            printf("ERROR: The input file %s doesn't follow the output file format from the stage 2 analysis of OPERA! Please check.\n", string(ppaFilename).c_str());
            exit(EXIT_FAILURE);
        }
        if(gwaslociflag && GWAScojosnplstName==NULL) {
            printf("Note: The OPERA analysis was performed based on GWAS loci only. If the proportion of GWAS loci explained by at least one xQTL data is interested, please input the GWAS loci file by --extract-GWAS-loci.\n", string(ppaFilename).c_str());
        }
        if(!gwaslociflag && GWAScojosnplstName!=NULL) {
            printf("Warning: The OPERA analysis was performed genome-wide scale rather than on GWAS loci only. If the proportion of GWAS loci explained by at least one xQTL data is interested, please run OPERA analysis on the GWAS loci region by adding loci file --extract-GWAS-loci.\n", string(ppaFilename).c_str());
        }
        PPANum = cidx_PPA.size() - 1;
        if(gwaslociflag) { expoNum = (cidx_PPA[0] - 3)/2;
        } else { expoNum = (cidx_PPA[0] - 1)/2; }
                        
        // read prior pi file
        vector<string> priorsplit; vector<double> prior;
        if(piFilename!=NULL) {
            read_pifile(piFilename, priorsplit);
        } else {
            split_string(priorstr,priorsplit);
        }
        // read_pifile(piFilename, priorsplit);
        if(priorsplit.size()!=(PPANum+1))
            throw("Error: The number of input prior probabilities is not consistent with the number of test configurations.");
        for(int t=0; t<(PPANum+1); t++)
        {
            prior.push_back(atof(priorsplit[t].c_str()));
            if(prior[t]<0 || prior[t]>1) throw("Error: --prior-pi. Prior probability values should be betweeen 0 and 1.");
        }
        // read num file
        vector<long> numsum(expoNum, 0);
        read_numfile(numFilename, numsum);
        // read gwas loci file
        lociData ldata;
        if(GWAScojosnplstName!=NULL) read_GWAS_cojo_snplist(&ldata, GWAScojosnplstName);

        vector<vector<string>> PPAvec, HEIDIvec, Namevec;
        // scan the ppa file
        long line_idx = 0;
        while(fgets(Tbuf, MAX_LINE_SIZE, rfile))
        {
            vector<string> PPAtmp, HEIDItmp, Nametmp;
            split_string(Tbuf, strlist, ", \t\n");
            if(strlist.size() != strsize) {
                printf("ERROR: The %ld row from the input file %s is incomplete! Please check.\n", line_idx,  string(ppaFilename).c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<strlist.size();i++) {
                if(i < cidx_PPA[0]) {
                    Nametmp.push_back(strlist[i]);
                } else if(i >= cidx_PPA[1] && i < cidx_HEIDI[0]) {
                    PPAtmp.push_back(strlist[i]);
                } else if(i >= cidx_HEIDI[0]) {
                    HEIDItmp.push_back(strlist[i]);
                }
            }
            Namevec.push_back(Nametmp);
            PPAvec.push_back(PPAtmp);
            HEIDIvec.push_back(HEIDItmp);
            line_idx++;        
        }
        fclose(rfile);
        printf("There are %d association combinations between %ld exposures and 1 outcome included in %s.\n",line_idx, expoNum, string(ppaFilename).c_str());
        // get the corresponding marginal index
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
        // compute pi0 correspoding configuration index
        vector<double> pi1(expoNum,0), pi0(expoNum,0);
        for(int t=0;t<expoNum;t++) {
            vector<int> ppcidx;
            for(int i=1;i<combmarg.size();i++) {
                if(combmarg[i].size()==(t+1)) {
                    for(int j=0;j<idxmarg[i].size();j++) {
                        ppcidx.push_back(idxmarg[i][j]);
                    }
                }
            }
            getUnique(ppcidx);
            for(int k=0;k<ppcidx.size();k++){
                pi1[t]+=prior[ppcidx[k]];
            }
        }
        for(int t=0;t<expoNum;t++) {
            pi0[t] = 1 - pi1[t];
        }
        // extract the significant results and remove duplicated
        // loop with the expoNum 
        vector<string> locisnps; vector<vector<string>> ppalocisnps(PPANum);
        for(int i=1;i<=expoNum;i++) {
            double nullNum = pi0[i-1] * numsum[i-1];
            // open the output file
            string ppafile0 = string(outFileName)+"_"+atos(i)+"_exposures_ppa.summary";
            FILE* ppa0;
            ppa0 = fopen(ppafile0.c_str(), "w");
            if (!(ppa0)) {
                printf("ERROR: open error %s\n", ppafile0.c_str());
                exit(1);
            }
            // output header
            string outstr;
            if(gwaslociflag) { outstr="Chr\tGWAS_SNP\tGWAS_bp\t";
            } else { outstr="Chr\t"; }
            for(int k=1; k<=i; k++) {
                outstr+="Expo"+atos(k)+"_ID"+'\t'+"Expo"+atos(k)+"_bp"+'\t';
            }
            outstr+="PPA(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\t"; outstr+="p_HEIDI(";
            for(int k=1; k<=i; k++) {
                outstr+=atos(k);
                if(k<i) outstr+=",";
            }
            outstr+=")\n";
            if(fputs_checked(outstr.c_str(),ppa0))
            {
                printf("ERROR: error in writing file %s .\n", ppafile0.c_str());
                exit(EXIT_FAILURE);
            }
            long itermcount = 0; 
            double ppasum = 0, fdr = 0, fpr = 0;
            for(int p=0;p<PPANum;p++) {
                if(combmarg[p+1].size() == i) {                
                    // scan the results and find max and unique ppa value
                    map<string, long> outidx; // prbname + lineidx
                    map<string, double> pparlts; // prbname + sum(ppa)
                    map<string, long> pparlts_count; // prbname + ppa num count
                    for(long l=0;l<line_idx;l++) {
                        //if(HEIDIvec[l][p]!="NA" && stof(PPAvec[l][p]) >= thresh_PP && stod(HEIDIvec[l][p]) >= thresh_heidi) {
                            string prbname;
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(gwaslociflag) { prbname.append(Namevec[l][2*(expoIdx+1) - 1]);
                                } else { prbname.append(Namevec[l][2*expoIdx - 1]);}
                            }                        
                            map<string, double>::iterator itmp;
                            itmp = pparlts.find(prbname);
                            if(itmp == pparlts.end()) { // not found in pparlts
                                pparlts.insert(pair<string, double>(prbname, stof(PPAvec[l][p])));
                                outidx.insert(pair<string, long>(prbname, l));
                                pparlts_count.insert(pair<string, long>(prbname, 1));
                            } else { // found in pparlts                            
                                double pparltssum;
                                pparltssum = stof(PPAvec[l][p]) + itmp->second;
                                itmp->second = pparltssum;
                                map<string, long>::iterator itmp0;
                                itmp0 = outidx.find(prbname);                                
                                if(itmp0 != outidx.end() && HEIDIvec[l][p]!="NA") {
                                    itmp0->second = l;
                                }
                                map<string, long>::iterator itmp1;
                                itmp1 = pparlts_count.find(prbname);                                
                                if(itmp1 != pparlts_count.end()) {
                                    itmp1->second = itmp1->second + 1;
                                }
                            }
                        //}
                    }
                    // output string
                    map<string, long>::iterator it;
                    for (it = outidx.begin(); it != outidx.end(); it++) {
                        long ii = it->second;
                        // chr or gwas loci information
                        if(gwaslociflag) { outstr = Namevec[ii][0] + '\t'+ Namevec[ii][1] + '\t' + Namevec[ii][2] + '\t'; 
                        } else { outstr = Namevec[ii][0] + '\t'; }
                        string prbname;
                        for(int j=0;j<combmarg[p+1].size();j++) {
                            int expoIdx = combmarg[p+1][j];
                            if(gwaslociflag) {
                                prbname.append(Namevec[ii][2*(expoIdx+1) - 1]);
                            } else {
                                prbname.append(Namevec[ii][2*expoIdx - 1]);                        
                            }
                        }
                        // take the mean ppa
                        map<string, double>::iterator itmp;
                        itmp = pparlts.find(prbname);
                        map<string, long>::iterator itmp1;
                        itmp1 = pparlts_count.find(prbname);
                        double ppamean = 0;
                        if(itmp != pparlts.end() && itmp1 != pparlts_count.end()) {
                            ppamean = itmp->second/itmp1->second;
                        }
                        if(ppamean >= thresh_PP && HEIDIvec[ii][p]!="NA" && stod(HEIDIvec[ii][p]) >= thresh_heidi) {
                            for(int j=0;j<combmarg[p+1].size();j++) {
                                int expoIdx = combmarg[p+1][j];
                                if(gwaslociflag) { outstr+= Namevec[ii][2*(expoIdx+1) - 1] + '\t'+ atos(Namevec[ii][2*(expoIdx+1)]) + '\t';
                                ppalocisnps[p].push_back(Namevec[ii][1]); locisnps.push_back(Namevec[ii][1]);
                                } else { outstr+= Namevec[ii][2*expoIdx - 1] + '\t'+ atos(Namevec[ii][2*expoIdx]) + '\t'; }
                            }
                            // outstr+=PPAvec[ii][p] +'\t'+ HEIDIvec[ii][p]+'\n';
                            outstr+=atos(ppamean) +'\t'+ HEIDIvec[ii][p]+'\n';
                            if(fputs_checked(outstr.c_str(), ppa0))
                            {
                            printf("ERROR: in writing file %s .\n", ppafile0.c_str());
                            exit(EXIT_FAILURE);
                            }
                            // ppasum += (1 - stof(PPAvec[ii][p]));
                            ppasum += (1 - ppamean);
                            itermcount+=1;
                        }                                                       
                    }
                }
            }
            fclose(ppa0);
            
            printf("\nPPA results for %ld combinatorial associations between %ld exposure(s) and 1 outcome have been extracted and saved in the file %s.\n",itermcount,i,ppafile0.c_str());
            if(itermcount > 0) {
                fdr = ppasum/itermcount;
                printf("The estimated FDR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fdr).c_str(),i);
            }
            if(nullNum > 0) {
                fpr = ppasum/nullNum;
                printf("The estimated FPR is %s for combinatorial associations between %ld exposure(s) and 1 outcome.\n",atos(fpr).c_str(),i);
            }
        }
        // output proportion of GWAS loci explained by combinations
        if(GWAScojosnplstName!=NULL) {
            getUnique(locisnps); ldata._include.clear(); double gwasprop = 0; vector<double> ppaprop(PPANum,0);
            StrFunc::match_only(locisnps, ldata._snp_name, ldata._include);
            gwasprop = ldata._include.size()*100/(ldata._snp_name.size());
            printf("\nThere are %s%% GWAS loci were detected to be associated with at least one xQTL data.\n",atos(gwasprop).c_str());
            for(int p=0;p<PPANum;p++) {
                getUnique(ppalocisnps[p]); ldata._include.clear();
                StrFunc::match_only(ppalocisnps[p], ldata._snp_name, ldata._include);
                if(ldata._include.size() > 0) ppaprop[p] = (double)ldata._include.size()/(ldata._snp_name.size());
            }
            // open .prop file for output
            string propfile = string(outFileName)+".prop";
            FILE* prop;
            prop = fopen(propfile.c_str(), "w");
            if (!(prop)) {
                printf("ERROR: open error %s\n", propfile.c_str());
                exit(1);
            }
            // .prop file header
            string propstr="";
            for(int i=1; i<=PPANum; i++) {
                if(i!=PPANum) {propstr+=strheader[cidx_PPA[i]]+'\t';
                } else {propstr+=strheader[cidx_PPA[i]]+'\n';}
            }
            if(fputs_checked(propstr.c_str(),prop))
            {
                printf("ERROR: error in writing file %s .\n", propfile.c_str());
                exit(EXIT_FAILURE);
            }
            // .prop output
            propstr="";
            for(int p=0;p<PPANum;p++) {
                if(p < (PPANum - 1)) { propstr += atos(ppaprop[p]) + "\t";
                } else { propstr += atos(ppaprop[p]) + "\n"; }
            }
            if(fputs_checked(propstr.c_str(), prop))
            {
                printf("ERROR: in writing file %s .\n", propfile.c_str());
                exit(EXIT_FAILURE);
            }
        }
        
    }

    void multi_joint_smr_func_old(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, double ngwas, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
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
            // get all the summary data        
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0 || smrwk.zxz.size()!=expoNum) {
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
            extract_smrwk_opera(&smrwk,max_ids,&smrwk_joint);
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
            // get all the summary data        
            long maxid = fill_smr_wk_mlt(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0 || smrwk.zxz.size()!=expoNum) {
                SMRRLT currlt;
                smrrlts.push_back(currlt);
                continue;
            }
            // compute the mean ypy
            //vector<int> splSizeidx;
            //match_only(smrwk.rs, gdata->snpName, splSizeidx);
            vector<double> ypy(smrwk.byz.size());
            double ypymedian = 0.0;
            for(int s=0; s<smrwk.byz.size(); s++) {
                uint32_t snpsplSize = smrwk.splSize[s];
                double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*snpsplSize;
                ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (snpsplSize-1) + pow(smrwk.byz[s],2));
                // double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*ngwas;
                // ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (ngwas-1) + pow(smrwk.byz[s],2));
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
                CommFunc::getUnique(all_ids[t]);
            }
            CommFunc::getUnique(cojo_ids);            
            // only include COJO SNPs for summary data
            MTSMRWKEXP smrwk_joint;
            extract_smrwk_opera(&smrwk,cojo_ids,&smrwk_joint);
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
            // find the joint SNP sample size
            vector<double> njointsplSize;
            for(int m=0; m<smrwk_joint.splSize.size();m++)
            {
                njointsplSize.push_back(smrwk_joint.splSize[m]);
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
                run_joint_effect_func(max_pos[t],cond_pos[t],_byz,_seyz,_byz_adj[t],_seyz_adj[t],_X,ngwas,ypymedian,njointsplSize);
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
                currlt.b_GWAS = byzt;
                currlt.se_GWAS = seyzt;
                currlt.p_GWAS = pchisq(byzt*byzt/(seyzt*seyzt), 1);
                currlt.b_eQTL = bxzt;
                currlt.se_eQTL = sexzt;
                currlt.p_eQTL = pchisq(bxzt*bxzt/(sexzt*sexzt), 1);
                currlt.b_SMR = bxy_val;
                currlt.se_SMR = sexy_val;
                currlt.p_SMR = pxy_val;
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
            //vector<int> splSizeidx;
            //match_only(smrwk.rs, gdata->snpName, splSizeidx);
            vector<double> ypy(smrwk.byz.size());
            double ypymedian = 0.0;
            for(int s=0; s<smrwk.byz.size(); s++) {
                uint32_t snpsplSize = smrwk.splSize[s];
                double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*snpsplSize;
                ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (snpsplSize-1) + pow(smrwk.byz[s],2));
                //double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*ngwas;
                //ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (ngwas-1) + pow(smrwk.byz[s],2));
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
            extract_smrwk_opera(&smrwk,cojo_ids,&smrwk_joint);
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
            // find the joint SNP sample size
            vector<double> njointsplSize;
            for(int m=0; m<smrwk_joint.splSize.size();m++)
            {
                njointsplSize.push_back(smrwk_joint.splSize[m]);
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
                run_joint_effect_func(max_pos[t],cond_pos[t],_byz,_seyz,_byz_adj[t],_seyz_adj[t],_X,ngwas,ypymedian,njointsplSize);
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
                currlt.b_GWAS = byzt;
                currlt.se_GWAS = seyzt;
                currlt.p_GWAS = pchisq(byzt*byzt/(seyzt*seyzt), 1);
                currlt.b_eQTL = bxzt;
                currlt.se_eQTL = sexzt;
                currlt.p_eQTL = pchisq(bxzt*bxzt/(sexzt*sexzt), 1);
                currlt.b_SMR = bxy_val;
                currlt.se_SMR = sexy_val;
                currlt.p_SMR = pxy_val;
                smrrlts.push_back(currlt);
            }
        }
    }
    void multi_joint_smr_func_v2(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, eqtlInfo &esdata, double ngwas, vector<string> &cojolist, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {
        MatrixXd _X;
        SMRWK smrwk;
        int expoNum = 1;
        uint64_t probNum = esdata._include.size();
        // compute the average GWAS sample size
        // double ngwas = 0.0;
        // ngwas = median(gdata->splSize);
        
        cis_itvl=cis_itvl*1000;
        for(int ii=0;ii<probNum;ii++)
        {
            int i=esdata._include[0];
            int probechr=esdata._epi_chr[i];            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx = i;
            smrwk.cur_chr = probechr;        
            // get all the summary data
            long maxid = fill_smr_wk_include(bdata, gdata, &esdata, &smrwk, refSNP, cis_itvl, false);
            if(smrwk.byz.size()==0) {
                SMRRLT currlt;
                smrrlts.push_back(currlt);
                continue;
            }
            // compute the mean ypy
            //vector<int> splSizeidx;
            //match_only(smrwk.rs, gdata->snpName, splSizeidx);
            vector<double> ypy(smrwk.byz.size());
            double ypymedian = 0.0;
            for(int s=0; s<smrwk.byz.size(); s++) {
                uint32_t snpsplSize = smrwk.splSize[s];
                double DJ = 2*smrwk.freq[s]*(1-smrwk.freq[s])*snpsplSize;
                ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (snpsplSize-1) + pow(smrwk.byz[s],2));
                //double DJ = 2*smrwk.freq[0][s]*(1-smrwk.freq[0][s])*ngwas;
                //ypy[s] = DJ * (pow(smrwk.seyz[s],2) * (ngwas-1) + pow(smrwk.byz[s],2));
            }
            ypymedian = median(ypy);
            // find the top and COJO SNPs index in smrwk
            int max_ids;
            vector<int> cojo_ids, all_ids;
            maxid = max_abs_id(smrwk.zxz);
            max_ids = maxid;
            all_ids.push_back(maxid);
            cojo_ids.push_back(maxid);
            for(int j=0; j<cojolist.size(); j++) {
                int cojoidx = find(smrwk.rs.begin(),smrwk.rs.end(),cojolist[j])-smrwk.rs.begin();
                if(cojoidx!=smrwk.rs.size()) {
                    all_ids.push_back(cojoidx);
                    cojo_ids.push_back(cojoidx);
                }
            }
            CommFunc::getUnique(all_ids);
            CommFunc::getUnique(cojo_ids);            
            // only include COJO SNPs for summary data
            SMRWK smrwk_joint;
            extract_smrwk_opera_new(&smrwk,cojo_ids,&smrwk_joint);
            bool minus_2p = true;
            make_XMat(bdata,smrwk_joint.curId, _X, minus_2p);
            // find the top and conditional SNP positions in cojo_ids
            int max_pos;
            vector<int> cond_pos;
            max_pos = find(cojo_ids.begin(),cojo_ids.end(),max_ids) - cojo_ids.begin();
            for(int k=0; k<cojo_ids.size(); k++) cond_pos.push_back(k);
            long tmp_ids = find(cojo_ids.begin(),cojo_ids.end(),maxid)-cojo_ids.begin();
            long cond_ids = find(cond_pos.begin(),cond_pos.end(),tmp_ids)-cond_pos.begin();
            cond_pos.erase(cond_pos.begin()+cond_ids);
            // check the top and COJO SNP list
            string maxSNPs = smrwk_joint.rs[max_pos];
            vector<string> condSNPs;
            for(int m=0; m<cond_pos.size(); m++) {
                condSNPs.push_back(smrwk_joint.rs[cond_pos[m]]);
            }
            // find the joint SNP sample size
            vector<double> njointsplSize;
            for(int m=0; m<smrwk_joint.splSize.size();m++)
            {
                njointsplSize.push_back(smrwk_joint.splSize[m]);
            }
            // get the summary data
            VectorXd _bxz(cojo_ids.size()),_sexz(cojo_ids.size()),_zsxz(cojo_ids.size());
            VectorXd _byz(cojo_ids.size()),_seyz(cojo_ids.size()),_byz_adj(cojo_ids.size()),_seyz_adj(cojo_ids.size());
            for(int j=0;j<cojo_ids.size();j++)
            {
                _bxz[j]=smrwk_joint.bxz[j];
                _sexz[j]=smrwk_joint.sexz[j];
                _zsxz[j]=smrwk_joint.zxz[j];
                _byz[j]=smrwk_joint.byz[j];
                _seyz[j]=smrwk_joint.seyz[j];
            }
            // run joint GWAS effect            
             run_joint_effect_eigenVar_func(max_pos,cond_pos,_byz,_seyz,_byz_adj,_seyz_adj,_X,ngwas,ypymedian,njointsplSize);
            //run_joint_effect_shrink_NoLDprune_func(max_pos,cond_pos,_byz,_seyz,_byz_adj,_seyz_adj,_X,ngwas,ypymedian,njointsplSize);
            // compute the joint bxy
            SMRRLT currlt;
            // compute and output joint SMR effect
            double byzt=_byz_adj[max_pos], seyzt=_seyz_adj[max_pos];
            double bxzt=_bxz[max_pos], sexzt=_sexz[max_pos];
            double bxy_val = byzt / bxzt;
            double sexy_val = -9;
            sexy_val = sqrt((seyzt * seyzt * bxzt * bxzt + sexzt * sexzt * byzt * byzt) / (bxzt * bxzt * bxzt * bxzt));
            double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
            double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
            currlt.ProbeID = esdata._epi_prbID[0];
            currlt.ProbeChr = esdata._epi_chr[0];
            currlt.Probe_bp = esdata._epi_bp[0];
            currlt.SNP = maxSNPs;
            currlt.b_GWAS = byzt;
            currlt.se_GWAS = seyzt;
            currlt.p_GWAS = pchisq(byzt*byzt/(seyzt*seyzt), 1);
            currlt.b_eQTL = bxzt;
            currlt.se_eQTL = sexzt;
            currlt.p_eQTL = pchisq(bxzt*bxzt/(sexzt*sexzt), 1);
            currlt.b_SMR = bxy_val;
            currlt.se_SMR = sexy_val;
            currlt.p_SMR = pxy_val;
            smrrlts.push_back(currlt);
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
        float ldr2_top_thresh = 0.81, ldr2_thresh = 0.9;
        vector<int> sn_ids, snp_ids;
        VectorXd ld_v;
        ld_calc_o2m(ld_v,0,Xtmp);
        for(int i=0;i<ld_v.size();i++)
        {
           double ldr2tmp=ld_v(i)*ld_v(i);
           if(ldr2tmp < ldr2_top_thresh) sn_ids.push_back(i);
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

    void run_joint_effect_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp)
    {
        int nld = X.rows(), jointsize_tmp = condpos.size() + 1;
        MatrixXd Xtmp(nld,jointsize_tmp);
        Xtmp.col(0) = X.col(maxpos);
        for(int s=1; s<jointsize_tmp; s++) {
            Xtmp.col(s) = X.col(condpos[s-1]);
        }
        // compute the GWAS zscore and p value
        VectorXd zyz, pyz(jointsize_tmp);
        zyz = byz.array()/seyz.array();
        for(int q=0;q<jointsize_tmp;q++)
        {
            pyz[q] = pchisq(zyz[q]*zyz[q],1);
        }
        // LD pruning with top > 0.9
        float ldr2_top_thresh = 0.81, ldr2_thresh = 0.81, gwas_lrt = 0.05, gwas_crt = 5e-8;//gwas_crt = 1e-100;
        vector<int> sn_ids, snp_ids;        
        VectorXd ld_v;
        ld_calc_o2m(ld_v,0,Xtmp);
        for(int i=0;i<ld_v.size();i++)
        {
           double ldr2tmp=ld_v(i)*ld_v(i);
           if(ldr2tmp < ldr2_top_thresh && pyz[i] <= gwas_lrt) sn_ids.push_back(i);
        }
        // initialize the adjusted value as marginal effect
        byz_adj = byz; seyz_adj = seyz;
        // perform joint analysis when condition matched
        if(sn_ids.size() > 0 && pyz[maxpos] <= gwas_lrt && pyz[maxpos] >= gwas_crt) {
        //if(sn_ids.size() > 0) {
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
            vector<int> jointpos; vector<double> ncleansnp;
            jointpos.push_back(maxpos);
            ncleansnp.push_back(njointsnp[maxpos]);
            if(snp_ids.size() > 0) {
                for(int s=0; s<snp_ids.size(); s++) {
                    jointpos.push_back(condpos[snp_ids[s]-1]);
                    ncleansnp.push_back(njointsnp[condpos[snp_ids[s]-1]]);
                }
                int jointsize = jointpos.size();
                VectorXd byz_joint(jointsize), seyz_joint(jointsize), byz_tmp(jointsize), seyz_tmp(jointsize);
                // run joint analysis on bzy
                MatrixXd X_joint(X.rows(),jointsize), D(jointsize,jointsize), V(jointsize,jointsize), Vinv(jointsize,jointsize), nV(jointsize,jointsize);
                double sigma_c=0.0;
                for(int s=0; s<jointsize; s++) {
                    byz_tmp[s] = byz[jointpos[s]];
                    seyz_tmp[s] = seyz[jointpos[s]];
                    X_joint.col(s) = X.col(jointpos[s]);
                }
                for(int p=0; p<jointsize; p++) {
                    for(int q=0; q<jointsize; q++) {
                        nV(p,q) = Min(ncleansnp[p],ncleansnp[q]);
                    }
                }
                V = (X_joint.transpose() * X_joint).array() * nV.array() * (1/(double)nld);
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
                sigma_c = (ypy - byz_joint.transpose() * D * byz_tmp)/(ncleansnp[0] - jointsize);
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

    void run_joint_effect_eigenVar_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp)
    {
        int nld = X.rows(), jointsize_tmp = condpos.size() + 1;
        MatrixXd Xtmp(nld,jointsize_tmp);
        Xtmp.col(0) = X.col(maxpos);
        for(int s=1; s<jointsize_tmp; s++) {
            Xtmp.col(s) = X.col(condpos[s-1]);
        }
        // compute the GWAS zscore and p value
        VectorXd zyz, pyz(jointsize_tmp);
        zyz = byz.array()/seyz.array();
        for(int q=0;q<jointsize_tmp;q++)
        {
            pyz[q] = pchisq(zyz[q]*zyz[q],1);
        }

        // LD pruning with top > 0.9
        // float ldr2_top_thresh = 0.81, ldr2_thresh = 0.81, gwas_lrt = 0.05, gwas_crt = 5e-8;//gwas_crt = 1e-100;
        float ldr2_top_thresh = 0.9, ldr2_thresh = 0.9, gwas_lrt = 0.05, gwas_crt = 5e-8;//gwas_crt = 1e-100;
        // float ldr2_top_thresh = 0.9, ldr2_thresh = 0.9, gwas_lrt = 0.99, gwas_crt = 1e-100;
        vector<int> sn_ids, snp_ids;        
        VectorXd ld_v;
        ld_calc_o2m(ld_v,0,Xtmp);
        for(int i=0;i<ld_v.size();i++)
        {
           double ldr2tmp=ld_v(i)*ld_v(i);
           if(ldr2tmp < ldr2_top_thresh && pyz[i] <= gwas_lrt) sn_ids.push_back(i);
        }
        // initialize the adjusted value as marginal effect
        byz_adj = byz; seyz_adj = seyz;
        // perform joint analysis when condition matched
        if(sn_ids.size() > 0 && pyz[maxpos] <= gwas_lrt && pyz[maxpos] >= gwas_crt) {
        //if(sn_ids.size() > 0) {
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
            vector<int> jointpos; vector<double> ncleansnp;
            jointpos.push_back(maxpos);
            ncleansnp.push_back(njointsnp[maxpos]);
            if(snp_ids.size() > 0) {
                for(int s=0; s<snp_ids.size(); s++) {
                    jointpos.push_back(condpos[snp_ids[s]-1]);
                    ncleansnp.push_back(njointsnp[condpos[snp_ids[s]-1]]);
                }
                int jointsize = jointpos.size();
                VectorXd byz_joint(jointsize), seyz_joint(jointsize), byz_tmp(jointsize), seyz_tmp(jointsize);
                // run joint analysis on bzy
                MatrixXd X_joint(X.rows(),jointsize), D(jointsize,jointsize), V(jointsize,jointsize), Rinv(jointsize,jointsize), nV(jointsize,jointsize);
                double sigma_c=0.0;
                for(int s=0; s<jointsize; s++) {
                    byz_tmp[s] = byz[jointpos[s]];
                    seyz_tmp[s] = seyz[jointpos[s]];
                    X_joint.col(s) = X.col(jointpos[s]);
                }
                for(int p=0; p<jointsize; p++) {
                    for(int q=0; q<jointsize; q++) {
                        nV(p,q) = Min(ncleansnp[p],ncleansnp[q]);
                    }
                }
                MatrixXd R;
                cor_calc(R, X_joint);
                V = (X_joint.transpose() * X_joint).array() * nV.array() * (1/(double)nld);
                D = V.diagonal().asDiagonal();
                // do eigen decomposition
                SelfAdjointEigenSolver<MatrixXd> eigensolver(R);
                MatrixXd evec = eigensolver.eigenvectors();
                VectorXd eval = eigensolver.eigenvalues();
                vector<int> eigen_index, posidx;
                // float eigen_thresh = 1e-8;
                float eigen_var_thresh = 0.995, eigenvardiff_min = 1;
                for(int i=0; i<eval.size(); i++) {
                    if(eval[i] > 0) posidx.push_back(i);
                }
                int poslength = posidx.size();
                vector<float> eigensum(poslength), eigenvardiff(poslength);
                float eigenplus = 0;
                for(int i=0; i<poslength; i++) {
                    eigenplus += eval[i];
                    eigensum[i] = eigenplus;
                }
                for(int i=0; i<poslength; i++) {
                    eigenvardiff[i] = Abs(eigensum[i]/eigensum[poslength - 1] - eigen_var_thresh);
                    if(eigenvardiff[i] <= eigenvardiff_min) {
                        eigenvardiff_min = eigenvardiff[i];
                        eigen_index.push_back(i);
                    }
                }
                MatrixXd evec_trunc(evec.rows(),eigen_index.size());
                VectorXd eval_trunc(eigen_index.size());
                for(int s=0; s<eigen_index.size(); s++) {
                    evec_trunc.col(s) = evec.col(eigen_index[s]);
                    eval_trunc(s) = eval(eigen_index[s]);
                }
                Rinv = evec_trunc * eval_trunc.asDiagonal().inverse() * evec_trunc.transpose();
                // Vinv = V.inverse();
                // compute joint effect and se
                byz_joint = Rinv * byz_tmp;
                sigma_c = (ypy - byz_joint.transpose() * D * byz_tmp)/(ncleansnp[0] - jointsize);
                if(sigma_c > 0) {
                    seyz_joint = (Rinv * D.inverse()).diagonal().cwiseSqrt() * sqrt(sigma_c);
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

    void run_joint_effect_shrink_NoLDprune_func(int maxpos, vector<int> &condpos, VectorXd &byz, VectorXd &seyz, VectorXd &byz_adj, VectorXd &seyz_adj, MatrixXd &X, double ngwas, double ypy, vector<double> njointsnp)
    {
        int nld = X.rows(), jointsize = condpos.size() + 1;
        // compute the GWAS zscore and p value
        VectorXd zyz, pyz(jointsize);
        zyz = byz.array()/seyz.array();
        for(int q=0;q<jointsize;q++)
        {
            pyz[q] = pchisq(zyz[q]*zyz[q],1);
        }
        // initialize the adjusted value as marginal effect
        byz_adj = byz; seyz_adj = seyz;
        // run joint analysis on bzy
        vector<int> jointpos; vector<double> ncleansnp;
        jointpos.push_back(maxpos);
        ncleansnp.push_back(njointsnp[maxpos]);
        for(int s=1; s<jointsize; s++) {
            jointpos.push_back(condpos[s-1]);
            ncleansnp.push_back(njointsnp[condpos[s-1]]);
        }
        VectorXd byz_joint(jointsize), seyz_joint(jointsize), byz_tmp(jointsize), seyz_tmp(jointsize);
        MatrixXd X_joint(X.rows(),jointsize), D(jointsize,jointsize), V(jointsize,jointsize), Rinv(jointsize,jointsize), nV(jointsize,jointsize);
        double sigma_c=0.0;
        for(int s=0; s<jointsize; s++) {
            byz_tmp[s] = byz[jointpos[s]];
            seyz_tmp[s] = seyz[jointpos[s]];
            X_joint.col(s) = X.col(jointpos[s]);
        }
        for(int p=0; p<jointsize; p++) {
            for(int q=0; q<jointsize; q++) {
                nV(p,q) = Min(ncleansnp[p],ncleansnp[q]);
            }
        }
        MatrixXd R; float lambda = 0.1, df = 0;
        cor_calc(R, X_joint);
        // get the degree of freedom
        SelfAdjointEigenSolver<MatrixXd> eigensolver(R);
        VectorXd eval = eigensolver.eigenvalues();
        vector<int> posidx;
        for(int i=0; i<eval.size(); i++) {if(eval[i] > 0) posidx.push_back(i);}
        int poslength = posidx.size();
        for(int i=0; i<poslength; i++) { df += eval[i]/(eval[i] + lambda); };
        V = (X_joint.transpose() * X_joint).array() * nV.array() * (1/(double)nld);
        D = V.diagonal().asDiagonal();
        // add the shrinkage lambda;
        for(int i = 0; i < R.rows(); ++i) R(i, i) += lambda;
        Rinv = R.inverse();
        byz_joint = Rinv * byz_tmp; 
        sigma_c = (ypy - byz_joint.transpose() * D * byz_tmp)/(ncleansnp[0] - df);
        if(sigma_c > 0) {
            seyz_joint = (Rinv * R * Rinv * D.inverse()).diagonal().cwiseSqrt() * sqrt(sigma_c);
            // seyz_joint = (Rinv * D.inverse()).diagonal().cwiseSqrt() * sqrt(sigma_c);
        } else {
            seyz_joint = seyz_tmp;
        }
        // push back the joint effect
        for(int b=0; b<jointsize; b++) {
            byz_adj[jointpos[b]] = byz_joint[b];
            seyz_adj[jointpos[b]] = seyz_joint[b];
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
            extract_smrwk_opera(&smrwk,max_ids,&smrwk_condi);
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
    void multi_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, vector<gwasData> &gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero, double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
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
             //if(opt) maxid=max_zsmr_id_opera(&smrwk, p_smr);
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
        int m_hetero_tmp = m_hetero > (expoNum + 1) ? m_hetero : (expoNum + 1);
        // int m_hetero_tmp = m_hetero + expoNum;
        m_hetero = m_hetero_tmp;
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
                if(opt) maxid=max_zsmr_id_opera(&smrwk, p_smr);
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
    // multiple exposures HEIDI with sample overlap
    void multi_heidi_func_so(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, vector<eqtlInfo> &esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest, bool new_heidi_mth, bool opt, double ld_min,int opt_hetero, bool sampleoverlap, double pmecs, int minCor)
    {    
        uint64_t probNum = esdata[0]._include.size();
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        MTSMRWKEXP smrwk;
        long expoNum = esdata.size();
        int m_hetero_tmp = m_hetero > expoNum ? m_hetero : expoNum;
        // int m_hetero_tmp = m_hetero + expoNum;
        m_hetero = m_hetero_tmp;
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
            if(smrwk.bxz[0].size() == 0) {
              currlt.p_HET=-9;
              currlt.nsnp=-9;
              smrrlts.push_back(currlt);
              continue;
            }
        
            Map<VectorXd> ei_bxz(&smrwk.bxz[0][0],smrwk.bxz[0].size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0][0],smrwk.sexz[0].size());
            MatrixXd theta(expoNum, expoNum);
            if(sampleoverlap)
            {
                // printf("Estimating the correlation ...\n");
                double z2mecs=qchisq(pmecs,1);
                double zmecs=sqrt(z2mecs);
                for(int t=0;t<expoNum;t++)
                {
                    for(int j=0;j<expoNum;j++)   
                    {
                        vector<double> zxz1, zxz2;
                        for(int k=0;k<smrwk.bxz[0].size();k++)
                        {                     
                            double z1=smrwk.bxz[t][k]/smrwk.sexz[t][k];
                            double z2=smrwk.bxz[j][k]/smrwk.sexz[j][k];
                            if(abs(z1)<zmecs && abs(z2)<zmecs )
                            {
                                zxz1.push_back(z1);
                                zxz2.push_back(z2);
                            }
                        }
                        if(zxz1.size()< minCor) {
                            // printf("WARNING: less than %d common SNPs obtained from the cis-region of probe %s at a p-value threshold %5.2e.\n ", minCor,probename.c_str(), pmecs);
                            // printf("probe %s is skipped.\n ", probename.c_str());
                            // continue;
                            currlt.p_HET=-9;
                            currlt.nsnp=-9;
                            smrrlts.push_back(currlt);
                            continue;
                        }
                        else
                        {
                            theta(t,j) = cor(zxz1,zxz2);
                            // printf("The estimated correlation between exposure and outcome%d is %f.\n",t+1,theta[t]);
                            // if(j==t) {theta(t,j) = cor(zxz1,zxz2);
                            // } else { theta(t,j) = 0;}
                        }
                    }
                }
            }
            zsxz=ei_bxz.array()/ei_sexz.array();
            if(refSNP==NULL) {
                if(opt) maxid=max_zsmr_id_opera(&smrwk, p_smr);
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
                        pdev = multi_heidi_test_new_v2_so(bdata, &smrwk, ld_top,  thresh_heidi,  m_hetero, nsnp, ld_min ,opt_hetero,sampleoverlap, theta);
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
        update_snidx_opera(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        MTSMRWK smrwk_heidi;
        extract_smrwk_opera(smrwk,sn_ids,&smrwk_heidi);
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
        update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
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
        update_snidx_opera(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
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
       // vector<int> tmp_ids;
       // match_only(maxSNPs,cmmnSNPs,tmp_ids);
        
       // for(int i=0;i<smrwk->zxz[0].size();i++)
       // {
       //     if(smrwk->zxz[0][i]*smrwk->zxz[0][i]-threshold>1e-6) sn_ids.push_back(i);
       // }
     
       //  printf("%ld SNPs left after filtering.\n",sn_ids.size());
       // update_snidx_opera(smrwk,sn_ids,MAX_NUM_LD,"LD pruning");
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk_opera(smrwk,sn_ids,&smrwk_heidi);
       // long maxid_heidi=max_abs_id(smrwk_heidi.zxz[0]);
        make_XMat(bdata,smrwk_heidi.curId, _X);
       // //printf("Removing SNPs with LD r-squared between top-SNP %s > %f or < %f...\n",smrwk_heidi.rs[maxid_heidi].c_str(),ldr2_top,ld_min);
       // ld_calc_o2m(ld_v,maxid_heidi,_X);

       // if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
       //     sn_ids.clear();
       //     for(int i=0;i<smrwk_heidi.zxz[0].size();i++)
       //     {
       //         if(i!= maxid_heidi)
       //         {
       //             double ldr2tmp=ld_v(i)*ld_v(i);
       //             if((ldr2tmp<ldr2_top) && (ldr2tmp > ld_min)) sn_ids.push_back(i);
       //         }
       //         else{
       //             sn_ids.push_back(i);
       //         }
       //     }
       // }
       // //printf("%ld SNPs are removed and %ld SNPs are retained.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
       // if(sn_ids.size() < m_hetero) {
       //     //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", sn_ids.size(), m_hetero);
       //     return -9;
       // }
       for(int i=0;i<sn_ids.size();i++) sn_ids[i]=i;
       update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
      // maxid_heidi=max_abs_id(smrwk_heidi.zxz[0]);
      //  printf("Removing one of each pair of remaining SNPs with LD r-squared > %f...\n",ldr2_top);
      // int m = (int)smrwk_heidi.bxz[0].size();
      // vector<int> rm_ID1;
       MatrixXd C;
       cor_calc(C, _X);
      // double ld_top=sqrt(ldr2_top);
      // if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
      // //printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are retained.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
      // if(m-rm_ID1.size() < m_hetero) {

      //     //printf("INFO: HEIDI test is skipped because the number of SNPs (%ld) is smaller than %d.\n", m-rm_ID1.size(), m_hetero);
      //     return -9;
      // }
      // //Create new index
      // sn_ids.clear();
      // int qi=0;
      // for (int i=0 ; i<m ; i++) {
      //     if (rm_ID1.size() == 0) sn_ids.push_back(i);
      //     else {
      //         if (qi<rm_ID1.size() && rm_ID1[qi] == i) qi++; //Skip removed snp
      //         else sn_ids.push_back(i);
      //     }
      // }
      // update_snidx_opera(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
       
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
        extract_smrwk_opera(smrwk,sn_ids,&smrwk_heidi);
        make_XMat(bdata,smrwk_heidi.curId, _X);
        for(int i=0;i<sn_ids.size();i++) sn_ids[i]=i;
        //update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
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

    // multi-exposure HEIDI with new pruning with sample overlap
    double multi_heidi_test_new_v2_so(bInfo* bdata,MTSMRWKEXP* smrwk,double ldr2_top, double threshold, int m_hetero,long &nsnp, double ld_min, int opt_hetero, bool sampleoverlap, MatrixXd theta)
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
        extract_smrwk_opera(smrwk,sn_ids,&smrwk_heidi);
        make_XMat(bdata,smrwk_heidi.curId, _X);
        for(int i=0;i<sn_ids.size();i++) sn_ids[i]=i;
        //update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
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
        pdev = bxy_mltheter_v2_so(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp, theta);
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
        // update_snidx_opera(smrwk_tmp,t,sn_ids,MAX_NUM_LD,"LD pruning");
        vector<long> maxidvect(expoNum);
        sort(sn_ids.begin(),sn_ids.end());
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk_opera(smrwk_tmp,sn_ids,&smrwk_heidi);
        for(int i=0;i<smrwk_heidi.byz.size();i++) sn_ids[i] = i;
        for(int t=0; t<expoNum; t++) {
            VectorXd ld_v; MatrixXd _Xtmp;
            long maxid = max_abs_id(smrwk_heidi.zxz[t]);
            maxidvect[t] = maxid;
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

        // 3. return NA if LD between top SNPs < 0.8;  
        double ldr2_between_top = 0.8;
        bool topLDpass = 1;
        for(int i=0; i<expoNum; i++) {
            VectorXd ld_v; MatrixXd _Xtmp;
            long maxid = max_abs_id(smrwk_heidi.zxz[i]);
            make_XMat(bdata, smrwk_heidi.curId, _Xtmp);
            ld_calc_o2m(ld_v,maxid,_Xtmp);
            for(int j=0; j<expoNum; j++) {
                if(j !=  i) {
                    if(abs(ldr2_top-1)>1e-6 || ld_min>0) {
                        long maxidx_j = maxidvect[j];                            
                        double ldr2tmp=ld_v(maxidx_j)*ld_v(maxidx_j);
                        if(ldr2tmp < ldr2_between_top) topLDpass = 0;
                    }
                }                
            }
        }
        // return NA if fail
        if(!topLDpass) {
           slctsnps.push_back("NA");
           return;
        }

        // 4. pairwise LD pruning
        make_XMat(bdata,smrwk_heidi.curId, _X);
        update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
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
        update_snidx_opera(smrwk_tmp,t,sn_ids,MAX_NUM_LD,"LD pruning");
        sort(sn_ids.begin(),sn_ids.end());
        
        MTSMRWKEXP smrwk_heidi;
        extract_smrwk_opera(smrwk_tmp,sn_ids,&smrwk_heidi);
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

        update_smrwk_x_opera(&smrwk_heidi,sn_ids,_X);
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
        
       // update_snidx_opera(&smrwk_heidi,sn_ids,opt_hetero,"HEIDI test");
        
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

    void est_cov_dxy_so(MatrixXd &covdxy,VectorXd &vardev,vector<int> maxid, VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz,vector<VectorXd> &_sexz, MatrixXd &_LD_heidi, MatrixXd theta)
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
                    MatrixXd cov1 = theta(k,j)*_LD_heidi.array()*(_sexz[k]*_sexz[j].transpose()).array()*boutco.array()/((_bxz[k]*_bxz[j].transpose()).array()*(_bxz[k]*_bxz[j].transpose()).array()).array();
                    MatrixXd cov2 = _LD_heidi.array()*seoutco.array()/(_bxz[k]*_bxz[j].transpose()).array();
                    // if(j==k) {covbxy = cov1 + cov2;
                    // } else {covbxy = cov2;}
                    covbxy = cov1 + cov2;

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
    // multi_heidi with sample overlap
    float bxy_mltheter_v2_so(VectorXd &_byz, vector<VectorXd> &_bxz, VectorXd &_seyz, vector<VectorXd> &_sexz, vector<VectorXd> &_zsxz, MatrixXd &_LD_heidi, long* snp_num, MatrixXd theta)
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
            
            est_cov_dxy_so(covdxy,vardev,maxid,_byz,_bxz,_seyz,_sexz,_LD_heidi,theta);
            
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
        //sort(gsnp.begin(), gsnp.end()); sort(bsnp.begin(), bsnp.end());
        //set_intersection(gsnp.begin(), gsnp.end(), bsnp.begin(), bsnp.end(), back_inserter(slctsnp));
        set_intersect(bsnp, gsnp, slctsnp);
        for(int t=0; t<etraitNum; t++) {
            esnp.clear(); esnp.resize(etrait[t]._esi_include.size());
            #pragma omp parallel for    
            for(int i=0; i<etrait[t]._esi_include.size(); i++)
                esnp[i]=etrait[t]._esi_rs[etrait[t]._esi_include[i]];
            cmmsnp.clear();
            //sort(slctsnp.begin(), slctsnp.end()); sort(esnp.begin(), esnp.end());
            //set_intersection(slctsnp.begin(), slctsnp.end(), esnp.begin(), esnp.end(), back_inserter(cmmsnp));
            set_intersect(slctsnp, esnp, cmmsnp);
            slctsnp.clear(); slctsnp.resize(cmmsnp.size());
            #pragma omp parallel for
            for(int i=0; i<cmmsnp.size(); i++)
                slctsnp[i]=cmmsnp[i];        
        }
        // reorder the SNP rsIDs based on bp positions using bfile
        // #pragma omp parallel for
        // for(int i=0; i<bdata->_include.size(); i++)
        //     bsnp[i]=bdata->_snp_name[bdata->_include[i]];
        // vector<int> bidtmp;
        // match_only(slctsnp, bsnp, bidtmp);        
        // if(bidtmp.empty()) throw("Error: no common SNPs found between data! Please check your LD reference data.");
        // sort(bidtmp.begin(), bidtmp.end());
        // slctsnp.clear(); slctsnp.resize(bidtmp.size());
        // for(int i=0; i<bidtmp.size(); i++) {
        //     slctsnp[i] = bsnp[bidtmp[i]];
        // }        
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
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata)
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
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, gwasData* gdata, eqtlInfo* esdata)
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
    void update_geIndx_opera(bInfo* bdata, vector<eqtlInfo> &etrait, eqtlInfo* esdata)
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
    long fill_smr_wk_mlt_old3(bInfo* bdata,gwasData* gdata,vector<eqtlInfo> &esdata,MTSMRWKEXP* smrwk,const char* refSNP,int cis_itvl,bool heidioffFlag)
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
    // updated; need to update the function to extract the workspace for cis variants
    // assume _esi_include.size() = snp.size();
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
                                smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);                                    
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
                            smrwk->splSize.push_back(gdata->splSize[gdata->_include[j]]);                            
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
                        smrwk->splSize.push_back(gdata->splSize[gdata->_include[ge_rowid]]);
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
                        }
                    }

                }                                          
            }       
        return maxid;
    } 

    void update_snidx_opera(MTSMRWK* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
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
    void update_snidx_opera(MTSMRWKEXP* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
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
    void update_snidx_opera(MTSMRWKEXP* smrwk,int t, vector<int> &sn_ids,int max_snp_slct, string forwhat)
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
    void extract_smrwk_opera(MTSMRWK* smrwk,vector<int> &sn_ids,MTSMRWK* smrwk2)
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
    void extract_smrwk_opera(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MTSMRWKEXP* smrwk2)
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
            smrwk2->splSize.push_back(smrwk->splSize[sn_ids[i]]);
        }
    }
    void update_smrwk_x_opera(MTSMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X)
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
    void update_smrwk_x_opera(MTSMRWKEXP* smrwk,vector<int> &sn_ids,MatrixXd &X)
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
}
