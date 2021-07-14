# OPERA
This software tool implements the OPERA (Omics PlEiotRopic Association) method to test for combinatorial pleiotropic associations of molecular phenotypes (e.g., expression level of a gene and DNA methylation level at CpG sites) with a complex trait of interest using summary-level data from GWAS and molecular QTL studies. OPERA is a Bayesian generalization of the SMR &amp; HEIDI approach to a multi-omics model, where the molecular phenotypes are considered as exposures and only the complex trait is considered as the outcome. This tool can therefore be used to prioritize molecular phenotypes that mediate the genetic effects for complex trait and further provide mechanistic interpretation of the GWAS signal.

## Credits
* Yang Wu developed the software tool with the support from Ting Qi, Jian Zeng, and Jian Yang.
* Yang Wu, Jian Zeng and Jian Yang developed the OPERA methods. 
* The software is programmed and maintained by Yang Wu.

## Questions and Help Requests
If you have any bug reports or questions, please send an email to Yang Wu (y.wu2@uq.edu.au) and Jian Zeng (j.zeng@uq.edu.au) at Institute for Molecular Bioscience, The University of Queensland, and Jian Yang (jian.yang@westlake.edu.cn) at School of Life Sciences, Westlake University.

## Citations
Wu Y., Qi T., Wray N.R., Visscher P.M., Zeng J. & Yang J. (2021) Joint analysis of multi-omics data reveals molecular mechanisms at GWAS loci. bioRxiv.

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
Our OPERA analysis consists of three steps. OPERA first estimates the global frequencies of each possible association patterns, then computes the posterior probability supporting each configuration by weighting the data likelihood with the estimated global frequencies. We compute the posterior probability of associations (PPA) for any combinatorial sites by summing up the posterior probability of configurations where the site combination present. For results passed the association test with a PPA threshold, OPERA performs the heterogeneity test to reject the associations that are due to linkage. 

We collected and prepared multiple public available molecular QTL data for users to perform the OPERA analysis with their specific complex trait of interest, which is available for download [here](https://cnsgenomics.com/software/smr/#DataResource). For illustration purpose, we also provide the demonstration [data](https://github.com/wuyangf7/OPERA/tree/main/demo) to run opera analysis with command line below. 

## Run OPERA for stage 1 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --estimate-pi --out myopera --thread-num 3

* --besd-flist reads a file to get the full paths of the multiple xQTL BESD files. The input format follows that for the SMR analysis (https://cnsgenomics.com/software/smr/#DataManagement). 
* --gwas-summary reads summary-level data from GWAS. The input format follows that for GCTA-COJO analysis (http://cnsgenomics.com/software/gcta/#COJO).
* --bfile reads individual-level SNP genotype data (in PLINK binary format) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files.
* --out saves the estimation of prior proportions from the OPERA stage 1 analysis in .pi file (text format, see below example).

```
Iteration	Pi1(0:0)	Pi2(0:1)	Pi3(1:0)	Pi4(1:1)
0	0.206741	0.618764	0.163905	0.0105903
1	0.57161	0.390242	0.0271588	0.0109884
2	0.535602	0.191737	0.272659	1.04125e-06
3	0.821147	0.0184099	0.160443	1.01198e-08
4	0.312696	0.583291	0.0777026	0.0263109
5	0.83348	0.118327	0.000677142	0.0475152
6	0.611609	0.365378	0.00324214	0.019771
...
```
Columns are iteration numbers from MCMC and posterior samples for each configuration from the MCMC.  


## Other parameters for stage 1 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --estimate-pi --prior-sigma 0.02,0.02 --pi-wind 100 --out myopera --thread-num 3 

* –-prior-sigma specifies the estimated variance of the non-zero mediated effects for each molecular trait on the complex trait.  It can be computed by the variance of the estimated SMR effects at the nominal significance level (i.e., 0.05) adjusting for the estimation errors, e.g., 0.02 (default).  
* --pi-wind defines a window centered on the molecular phenotype with smallest number of sites to select no overlap independent loci, e.g., 100 (default). 

## Run OPERA for stage 2 analysis and heterogeneity analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --prior-pi 0.8,0.09,0.09,0.02 --out myopera --thread-num 3

* --prior-pi the estimated global proportions of each configuration from the stage 1 analysis. 
* --out saves the PPA and multi-exposure HEIDI test P-values for each possible association hypotheses in .ppa file (text format, see below example), and saves the results from the pairwise SMR analysis in .smr file

```
Chr	Expo1_ID	Expo1_Gene	Expo1_bp	Expo2_ID	Expo2_Gene	Expo2_bp	PPA(0)	PPA(1)	PPA(2)	PPA(1:2)	HEIDI(1)	HEIDI(2)	HEIDI(1:2)
7	ENSG00000238109	AC004893.10	98596857	cg19636519	GJC3	99541626	0.180167	0.00832795	0.818992	0.0074865	NA	1.248951e-05	NA
7	ENSG00000238109	AC004893.10	98596857	cg08582801	AZGP1P1	99588335	0.135318	0.00809727	0.864083	0.00749766	NA	1.184727e-01	NA
7	ENSG00000238109	AC004893.10	98596857	cg07693238	AZGP1P1	99595437	0.197766	0.00824852	0.801311	0.00732488	NA	2.194922e-02	NA
7	ENSG00000085514	PILRA	99981436	cg19116668	PILRB	99932089	4.09321e-15	0.999997	1	0.999997	1.036234e-01	7.709418e-02	1.604312e-01
...
```
Columns are chromosome,probe ID for the 1st exposure, gene name for the 1st exposure, probe position for the 1st exposure, probe ID for the 2nd exposure, gene name for the 2nd exposure, probe position for the 2nd exposure, PPA for no associations, PPA for the 1st probe association, PPA for the 2nd probe association, PPA for the 1st and 2nd probes joint association, p-value from HEIDI for the 1st probe association, p-value from HEIDI for the 2nd probe association, and p-value from HEIDI for the 1st and 2nd probes joint association.  
Missing Value is represented by "NA".

```
ProbeID    Probe_Chr   Gene    Probe_bp    SNP SNP_Chr SNP_bp  A1  A2  Freq    b_GWAS  se_GWAS p_GWAS  b_eQTL  se_eQTL p_eQTL  b_SMR   se_SMR  p_SMR   p_HEIDI nsnp_HEIDI
prb01    1   Gene1   1001    rs01    1   1011    C   T   0.95    -0.024  0.0063  1.4e-04 0.36    0.048   6.4e-14 -0.0668 0.0197  6.8e-04 NA  NA
prb02    1   Gene2   2001    rs02    1   2011    G   C   0.0747  0.0034  0.0062  5.8e-01 0.62    0.0396  2e-55   0.0055  0.01    5.8e-01 4.17e-01    28
...
```
Columns are probe ID, probe chromosome, gene name, probe position, SNP name, SNP chromosome, SNP position, the effect (coded) allele, the other allele, frequency of the effect allele (estimated from the reference samples), effect size from GWAS, SE from GWAS, p-value from GWAS, effect size from eQTL study, SE from eQTL study, p-value from eQTL study, effect size from SMR, SE from SMR, p-value from SMR, p-value from HEIDI (HEterogeneity In Depedent Instruments) test, and number of SNPs used in the HEIDI test.

* The heterogeneity test (i.e., multi-exposure HEIDI) will be automatically performed for any combinatorial associations passed a PPA threshold (0.8 as default). If the heterogeneity test is not interested, it can be turned off by specifying --heidi-off.

## Other parameters for stage 2 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --extract-exposure-probe myexposure --outcome-wind 1000 --thresh-PP 0.5 --thresh-SMR 0.05 --extract-target-cojo-snps mycojo --extract-GWAS-loci myloci --prior-pi 0.8,0.09,0.09,0.02 --prior-sigma 0.02,0.02 --out myopera --thread-num 3 

* --extract-exposure-probe	extracts a subset of exposure sites for analysis
* --outcome-wind specifies the window around each GWAS loci for stage 2 analysis, e.g., 500 (default). 
* --extract-GWAS-loci extracts a subset of GWAS loci for analysis
* --extract-target-cojo-snps specifies full COJO SNP list for each site of molecular phenotype as the target to compute the joint SMR effect.
* --thresh-PP specifies significance threshold of PPA to perform heterogeneity test and output the significant results. 
* --thresh-SMR specifies significance threshold of SMR to perform the OPERA analysis.

OPERA shares the same data management function and flags with the SMR software, for a full list of option reference, please see [here](https://cnsgenomics.com/software/smr/#OptionsReference). 













