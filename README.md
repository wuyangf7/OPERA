# OPERA
This software tool implements the OPERA (Omics PlEiotRopic Association) method to test for combinatorial pleiotropic associations of molecular phenotypes (e.g., expression level of a gene and DNA methylation level at CpG sites) with a complex trait of interest using summary-level data from GWAS and molecular QTL studies. OPERA is a Bayesian generalization of the SMR &amp; HEIDI approach to a multi-omics model, where the molecular phenotypes are considered as exposures and only the complex trait is considered as the outcome. This tool can therefore be used to prioritize molecular phenotypes that mediate the genetic effects for complex trait and further provide mechanistic interpretation of the GWAS signal.

## Credits
* Yang Wu developed the software tool with the support from Ting Qi, Jian Zeng, and Jian Yang.
* Yang Wu, Jian Zeng and Jian Yang developed the OPERA methods. 
* The software is programmed and maintained by Yang Wu.

## Questions and Help Requests
If you have any bug reports or questions, please send an email to Yang Wu (y.wu2@uq.edu.au) and Jian Zeng (j.zeng@uq.edu.au) at Institute for Molecular Bioscience, The University of Queensland, and Jian Yang (jian.yang@westlake.edu.cn) at School of Life Sciences, Westlake University.

## Citations
Wu Y., Qi T., Wray N.R., Visscher P.M., Zeng J. & Yang J. (2021) Joint analysis of GWAS and multi-omics QTL summary statistics reveals a large fraction of GWAS signals shared with molecular phenotypes. Under review.

## Installation
To install OPERA, you can download the [opera_Linux.zip](https://github.com/wuyangf7/OPERA/blob/main/opera_Linux.zip) package, which contains a standalone (i.e., statically linked) 64-bit Linux executable file opera_Linux. We strongly recommend using this static executable because it is well-optimized and no further installation is required.

For compiling yourself, you should clone this repository and use the MAKEFILE,  
```
git clone https://github.com/wuyangf7/OPERA
cd OPERA
make
```
There are dependencies on your local MKL, BOOST and EIGEN Libraries.

# Tutorial
The OPERA analysis consists of three steps. OPERA first estimates the global frequencies of each possible association pattern, then computes the posterior probability supporting each configuration by weighting the data likelihood with the estimated global frequencies. We compute the posterior probability of associations (PPA) for any combinatorial sites by summing up the posterior probability of configurations where the site combination present. For results passed the association test with a PPA threshold (e.g., 0.9 by default), OPERA performs the heterogeneity test to reject the associations that are due to linkage. 

We collected and prepared multiple public available molecular QTL data for users to perform the OPERA analysis with their specific complex trait of interest, which is available for download [here](https://cnsgenomics.com/software/smr/#DataResource). For illustration purpose, we also provide the demonstration [data](https://github.com/wuyangf7/OPERA/tree/main/demo) to run opera analysis with command line below. 

## Run OPERA for stage 1 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --mbfile mybdatalist --estimate-pi --out myopera

* --besd-flist reads a file to get the full paths of the multiple xQTL BESD files. The input format follows that for the SMR analysis (https://cnsgenomics.com/software/smr/#DataManagement). 
* --gwas-summary reads summary-level data from GWAS. The input format follows that for GCTA-COJO analysis (http://cnsgenomics.com/software/gcta/#COJO).
* --mbfile reads a text file with each row representing a PLINK binary file (e.g., one for each chromosome) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files. The chromosome-wide .fam files are required to contain the same individuals. Note that stage 1 analysis requires the genome-wide LD reference to estimate the global parameters. If the genome-wide genotype data in PLINK format (containing the genotype data for each chromosome) is ready, please switch the flag --bfile to read the standard PLINK binary file as LD reference. 
* --out saves the estimation of prior probabilities from the OPERA stage 1 analysis in .pi file and the estimation of prior variances in .var file (text format, see example from demo data below).

```
Posteriors      Pi1(0:0)        Pi2(0:1)        Pi3(1:0)        Pi4(1:1)
Mean    0.605441        0.168855        0.0743364       0.151368
SD      0.20765 0.186273        0.134037        0.16549
```
```
Variance        2.000000e-02    2.000000e-02
```
The .pi output includes the posterior mean and SD from MCMC for each possible configuration. The .var output includes the prior variance for each xQTL data. The posterior samples from MCMC are also printed in the log file, for example, 
```
Iteration       Pi1(0:0)        Pi2(0:1)        Pi3(1:0)        Pi4(1:1)
0       0.0557406       0.319719        0.0834392       0.541101
100     0.163492        0.422161        0.413673        0.000673509
200     0.661116        0.338883        2.4066e-07      6.7008e-15
300     0.587679        0.412321        5.49276e-11     3.1432e-09
400     0.814359        0.101207        0.0586814       0.0257529
500     0.430485        8.37885e-05     0.0929275       0.476504
600     0.282069        0.657061        0.017055        0.0438149
...
```
Columns are iteration numbers and posterior samples for each configuration from the MCMC.  

### Other parameters for stage 1 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --sample-overlap --estimate-pi --extract-snp mysnplist --prior-var 0.02,0.02 --pi-wind 100 --out myopera --thread-num 3

* --bfile reads individual-level SNP genotype data (in PLINK binary format) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files. 
* --extract-snp specifies a snplist (e.g., [Hapmap3 SNP list](https://github.com/wuyangf7/OPERA/blob/main/demo/hapmap.snplist)) to be extracted across LD reference, xQTL data and GWAS summary data, and used for stage 1 analysis. 
* –-prior-var specifies the prior variance of the non-zero mediated effects for each molecular trait on the complex trait. It can be computed by the variance of the estimated SMR effects of each molecular trait on complex trait at a significance level (e.g., FDR < 0.05) adjusting for the estimation errors.
* --sample-overlap specifies the flag to let OPERA consider the between-study correlations due to overlapping samples. OPERA will automatically output the estimated correlations in .rho file (text format, see below example). 
```
1.000000e+00	0.200000e+00
0.200000e+00	1.000000e+00
```
* --opera-smr turns on the flag of using the estimated smr effect rather than the estimated joint smr effect to run the stage 1 analysis.  
* --pi-wind defines a window centered on the molecular phenotype with smallest number of sites to select no overlap independent loci, e.g., 200 Kb (default). 
* --thread-num specifies the number of OpenMP threads for parallel computing.

## Run OPERA for stage 2 analysis and heterogeneity analysis
> opera --besd-flist mylist --extract-gwas-loci myloci --chr 7 --gwas-summary mygwas.ma --bfile mydata --prior-pi-file myopera.pi --prior-var-file myopera.var --out myopera_chr7

Note: Only the cis-SNPs of each exposure site are used, so the stage 2 analysis can be performed for each chromosome seperately. The genome-wide analysis results can be combined through shell command below, e.g., for association between single exposure and outcome,
> awk 'NR==1 || FNR>1' myopera_chr*_1_expos_assoc.ppa > myopera_1_expos_assoc.ppa  

* --extract-gwas-loci extracts a subset of genomic region with GWAS COJO or LD clumpped loci (i.e., 2Mb region by default) to perform opera analysis. The input file format is
```
Chr     SNP     bp
7       rs1859788       99971834
7       rs7810606       143108158
```
* --chr specifies the target chromosome for SNPs and exposure sites to run chromosome-wide opera stage 2 analysis.
* --bfile reads individual-level SNP genotype data (in PLINK binary format) from a reference sample for LD estimation. 
* --prior-pi-file reads the prior probabilities estimated from the stage 1 analysis (i.e., the posterior Mean from stage 1 analysis).  
* --prior-var-file reads the prior variances estimated from the stage 1 analysis.  
* --out saves the marginal or joint PPA and multi-exposure HEIDI test P-values for each association hypothesis passed significance threshold (e.g., PPA > 0.9 and P-HEIDI > 0.01 default) in .ppa file (text format, see below example). Results are seperated as different number of exposure combinations associated with the outcome, e.g., single exposure association, 

```
Chr	GWAS_SNP	GWAS_bp	Expo1_ID	Expo1_bp	PPA(1)	p_HEIDI(1)
7	rs1859788	99971834	cg08582801	99588335	0.929121	8.152293e-02
7	rs1859788	99971834	cg13210467	99775443	0.999054	1.625980e-01
7	rs1859788	99971834	cg19116668	99932089	1	1.125248e-01
7	rs1859788	99971834	cg26429636	99573747	0.915179	5.926392e-02
...
```
and associations between 2 exposure(s) and 1 outcome.
```
Chr	GWAS_SNP	GWAS_bp	Expo1_ID	Expo1_bp	Expo2_ID	Expo2_bp	PPA(1,2)	p_HEIDI(1,2)
7	rs1859788	99971834	ENSG00000085514	99981436	cg13210467	99775443	0.998746	1.787164e-02
7	rs1859788	99971834	ENSG00000146828	100444536	cg01869186	100423987	0.906761	1.641383e-01
7	rs1859788	99971834	ENSG00000146828	100444536	cg12616177	100434510	0.902426	1.044594e-01
...
```
Columns are chromosome, rs ID for GWAS SNP, physical position for GWAS SNP, probe ID for the 1st exposure, probe position for the 1st exposure, probe ID for the 2nd exposure, probe position for the 2nd exposure, PPA for the 1st and 2nd exposures joint association, and p-value from HEIDI for the 1st and 2nd exposures joint association.  

Note: If gwas loci file is not specified, the program will scan all the possible combinations between exposure sites for the whole genome. 

The program will automatically ouput the proportion and number of GWAS loci associated with different combination of xQTL data in .prop file (see example below).    
```
Overall	Exposures(1)	Exposures(2)	Exposures(1,2)
0.5	0	0.5	0.5
1	0	1	1
```
The output also includes the estimated false discovery rate (FDR) and false positive rate (FPR) for any combinatorial associations, which are printed in the log file, for example,
```
PPA results for 4 combinatorial associations between 1 exposure(s) and 1 outcome have been extracted and saved in the file /Users/uqywu16/Desktop/mtSMR/software/demo/myopera_1_expos_assoc.ppa.
The estimated FDR is 0.0391616 for combinatorial associations between 1 exposure(s) and 1 outcome.
The estimated FPR is 0.0168247 for combinatorial associations between 1 exposure(s) and 1 outcome.

PPA results for 3 combinatorial associations between 2 exposure(s) and 1 outcome have been extracted and saved in the file /Users/uqywu16/Desktop/mtSMR/software/demo/myopera_2_expos_assoc.ppa.
The estimated FDR is 0.0640222 for combinatorial associations between 2 exposure(s) and 1 outcome.
The estimated FPR is 0.00503074 for combinatorial associations between 2 exposure(s) and 1 outcome.

There are 50% GWAS loci were detected to be associated with at least one type of xQTL data.
```
Note: we suggest a PPA threshold of 0.9 to roughly control the FDR below 0.05. However, if more strigent FDR or FPR is required and interested, OPERA can acheive this by increasing the PPA threshold (e.g., 0.995). 

OPERA also automatically outputs the results from the SMR analysis of molecular phenotypes and complex trait in the .smr file
 
```
probeID	ProbeChr	Gene	Probe_bp	topSNP	topSNP_chr	topSNP_bp	A1	A2	Freq	b_GWAS	se_GWAS	p_GWAS	b_eQTL	se_eQTL	p_eQTL	b_SMR	se_SMR	p_SMR	p_HEIDI	nsnp_HEIDI
ENSG00000238109	7	AC004893.10	98596857	rs2283015	7	98600839	G	C	0.0847203	-0.00193708	0.00365864	5.964895e-01	0.711699	0.0382071	1.926562e-77	-0.00272178	0.00514278	5.966378e-01	7.575220e-01	20
ENSG00000146833	7	TRIM4	99495902	rs2571997	7	99514417	A	C	0.410541	0.00614763	0.00213414	3.969059e-03	-0.690266	0.00715596	0.000000e+00	-0.00890618	0.00309315	3.985252e-03	8.417932e-05	20
ENSG00000166526	7	ZNF3	99670913	rs67110214	7	99633739	C	G	0.28583	0.00401872	0.0024395	9.948536e-02	-0.0968669	0.0088226	4.801025e-28	-0.041487	0.0254659	1.032880e-01	4.811370e-05	20
cg03856969	7	AK001533	98600799	rs4729505	7	98601025	C	T	0.0844313	-0.0024541	0.00373076	5.106651e-01	-0.431667	0.0598434	5.462031e-13	0.00568517	0.00867854	5.124136e-01	6.233012e-01	20
cg10154880	7	AK001533	98603502	rs17720576	7	98616657	A	G	0.0865118	-0.0018608	0.00371972	6.168967e-01	0.489588	0.05912	1.218873e-16	-0.00380076	0.00761151	6.175377e-01	8.987915e-01	20
...
```
Columns are probe ID, probe chromosome, gene name, probe position, SNP name, SNP chromosome, SNP position, the effect (coded) allele, the other allele, frequency of the effect allele (estimated from the reference samples), effect size from GWAS, SE from GWAS, p-value from GWAS, effect size from eQTL study, SE from eQTL study, p-value from eQTL study, effect size from SMR, SE from SMR, p-value from SMR, p-value from HEIDI test, and number of SNPs used in the HEIDI test.

The heterogeneity test (i.e., multi-exposure HEIDI) will be automatically performed for any combinatorial associations passed a PPA threshold (0.9 as default). If the heterogeneity test is not interested, it can be turned off by specifying --heidi-off.


### Other parameters for stage 2 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --chr 7 --extract-exposure-probe myexposure --sample-overlap --rho-file myopera.rho --outcome-wind 1000 --thresh-ppa 0.5 --thresh-smr 0.05 --thresh-heidi 0.01 --print-combo-ppa-res --extract-target-cojo-snps mycojo --extract-gwas-loci myloci --prior-pi 0.8,0.09,0.09,0.02 --prior-var 0.02,0.02 --out myopera --thread-num 3
* --prior-pi specifies the prior probalities for stage 2 analysis.
* --prior-var specifies the prior variance for stage 2 analysis (seperated by comma).
* --rho-file specifies the estimated between-study correlations due to sample overlap from the stage 1 analysis. See the input format for the .rho file above. 
* --extract-exposure-probe	extracts a subset of exposure sites for opera analysis.
* --outcome-wind specifies the window around each GWAS loci for stage 2 analysis, or the window around each site across exposures for stage 2 analysis when GWAS loci were not specified, e.g., 1Mb in either direction (default). 
* --extract-target-cojo-snps specifies full COJO SNP list for each site of molecular phenotype as the conditional SNP set to compute the joint SMR effect for other exposure sites. The input file format is
```
ENSG00000196367 rs182325057,rs12532598,rs17638906,rs62472014,rs219843,rs183732601,rs150746244
ENSG00000238109 rs2283015,rs142345619,rs117749026,rs117295696,rs139511767
ENSG00000242687 rs34631688,rs187375676,rs149211972,rs219813,rs6976207,rs7789895,rs10264067,rs150746244
```
* --thresh-ppa specifies significance threshold of PPA to perform heterogeneity test and output the significant results. 
* --thresh-smr specifies significance threshold of SMR to perform the OPERA analysis, e.g., 0.05 (default).
* --thresh-heidi specifies significance threshold of single-exposure HEIDI test to perform the OPERA analysis, e.g., 0.01 (default).
* --opera-smr turns on the flag of runing OPERA analysis using the estimated SMR effect rather than estimated joint-SMR effect. 
* --thread-num specifies the number of OpenMP threads for parallel computing. 
* --print-combo-ppa-res saves the ppa for each possible association hypothesis across all combinations in .res file if specified. 
```
Chr	GWAS_SNP	GWAS_bp	Expo1_ID	Expo1_bp	Expo2_ID	Expo2_bp	PPA(0)	PPA(1)	PPA(2)	PPA(1,2)	p_HEIDI(1)	p_HEIDI(2)	p_HEIDI(1,2)
7	rs1859788	99971834	ENSG00000085514	99981436	cg27486786	99074036	0.4934397	0.37980783	0.42463091	0.2978785	NA	NA	NA
7	rs1859788	99971834	ENSG00000085514	99981436	cg05832509	99515080	0.5101769	0.36135745	0.41064838	0.28218272	NA	NA	NA
7	rs1859788	99971834	ENSG00000085514	99981436	cg26429636	99573747	0.14851999	0.56356138	0.83166325	0.54374462	NA	NA	NA
7	rs1859788	99971834	ENSG00000085514	99981436	cg08582801	99588335	0.13705845	0.57444245	0.84428662	0.55578756	NA	NA	NA
7	rs1859788	99971834	ENSG00000085514	99981436	cg13210467	99775443	2.3377886e-06	0.99889582	0.99984795	0.9987461	5.470704e-02	1.625980e-01	1.787164e-02
7	rs1859788	99971834	ENSG00000085514	99981436	cg19116668	99932089	1.4435816e-11	0.54813665	1	0.54813665	NA	1.125248e-01	NA
...
```
Columns are chromosome, rs ID for GWAS SNP, physical position for GWAS SNP, probe ID for the 1st exposure, probe position for the 1st exposure, probe ID for the 2nd exposure, probe position for the 2nd exposure, PPA for no associations, PPA for the 1st exposure marginal association, PPA for the 2nd exposure marginal association, PPA for the 1st and 2nd exposures joint association, p-value from HEIDI for the 1st exposure association, p-value from HEIDI for the 2nd exposure association, and p-value from HEIDI for the 1st and 2nd exposures joint association. Missing value is represented by "NA". 

## Command line to generate a data file for xQTL plot
> opera --besd-flist mylist1 --gwas-summary mygwas.ma --bfile mydata --beqtl-summary xQTL1 --out myplot --plot --probe ENSG00000085514 --probe-wind 500 --gene-list glist-hg19
* --beqtl-summary reads the BESD file for one type of xQTL (e.g.,eQTL).
* --besd-flist reads a file to get the full paths of the other xQTL BESD files.
* --probe specifies a probe site (e.g., often a gene) as locus center to generate the plot.
* --probe-wind specifies a window to extract the exposure site in the window
* --gene-list specifies a gene range list, which is available for download [here](https://github.com/wuyangf7/OPERA/blob/main/demo/glist-hg19).

To generate the loucs omics plot with extracted file, please see the omics SMR plot page (https://plot.cnsgenomics.com/omicsplot/).

We also provide an [R scirpt](https://github.com/wuyangf7/OPERA/blob/main/demo/plot/plot_OmicsSMR_xQTL.r) to plot the omics SMR plot as presented in Wu et al.. Please see the demo plot below.  

## R commands to draw the plots
```
source("./plot/plot_OmicsSMR_xQTL.r") 
SMRData = ReadomicSMRData("./plot/myplot.ENSG00000085514.txt")
omicSMRLocusPlot(data=SMRData,esmr_thresh=1e-4,msmr_thresh=1e-4,eprobeNEARBY="ENSG00000085514",mprobeNEARBY=c("cg13210467"),trait_name="Test",funcAnnoFile="./plot/funcAnno.RData")
```
* esmr_thresh and msmr_thresh are the threshold for first xQTL and other xQTLs, respectively
* eprobeNEARBY specifies the eQTL association pattern for specific gene
* mprobeNEARBY specifies the xQTL association pattern for other specific sites
* funcAnnoFile reads the epigenome annotation data from the Roadmap Epigenomics Consortium, which can be downloaded [here](https://www.dropbox.com/s/0rf4hmdjzostuv2/funcAnno.RData?dl=0).

OPERA shares the same data management function and flags with the SMR software, for a full list of option reference, please see [here](https://cnsgenomics.com/software/smr/#OptionsReference). 
