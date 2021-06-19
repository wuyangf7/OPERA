# OPERA
This software tool implements the OPERA (Omics PlEiotRopic Association) method to test for combinatorial pleiotropic associations of molecular phenotypes (e.g., expression level of a gene and DNA methylation level at CpG sites) with a complex trait of interest using summary-level data from GWAS and molecular QTL studies. OPERA is a Bayesian generalization of the SMR &amp; HEIDI approach8 to a multi-omics model, where the molecular phenotypes are considered as exposures and only the complex trait is considered as the outcome. This tool can therefore be used to prioritize molecular phenotypes that mediate the genetic effects for complex trait and provide mechanistic interpretation of the GWAS signal.

# Credits
Yang Wu developed the software tool with the support from Ting Qi, Jian Zeng, and Jian Yang.
Yang Wu, Jian Zeng and Jian Yang developed the OPERA methods. 
The software is programmed and maintained by Yang Wu.

# Questions and Help Requests
If you have any bug reports or questions, please send an email to Yang Wu (y.wu2@uq.edu.au), Jian Zeng (j.zeng@uq.edu.au), and Jian Yang (jian.yang@westlake.edu.cn) at Institute for Molecular Bioscience, The University of Queensland.

# Citations
Wu Y., Qi T., Wray N.R., Visscher P.M., Zeng J. & Yang J. (2021) Joint analysis of multi-omics data reveals molecular mechanisms at GWAS loci. bioRxiv.

# Installation
To install OPERA, you can download the opera_Linux.zip package, which contains a standalone (i.e., statically linked) 64-bit Linux executable file opera_Linux. We strongly recommend using this static executable because it is well-optimized and no further installation is required.

For compiling yourself, you should clone this repository and use the MAKEFILE,  
```
git clone https://github.com/Yves-CHEN/DENTIST
cd DENTIST
make
```
There are dependencies on your local MKL, BOOST and EIGEN Libraries.

# Tutorial
Our OPERA analysis consists of three steps. OPERA first estimates the global frequencies of each possible association patterns, then computes the posterior probability for supporting each configuration by weighting the data likelihood with the estimated prior proportions. We compute the posterior probability of associations (PPA) for any combinatorial sites by summing up the posterior probability of configurations where the site combination present. For results passed the association test with PPA threshold, OPERA performs the heterogeneity analysis test to reject associations that are due to linkage. 

# run OPERA for stage 1 analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --estimate-pi –prior-sigma 0.02,0.02 --out myopera --thread-num 3

* --besd-flist reads a file to get the full paths of the multiple xQTL BESD files. The input format follows that for the SMR analysis (https://cnsgenomics.com/software/smr/#DataManagement). 
* --gwas-summary reads summary-level data from GWAS. The input format follows that for GCTA-COJO analysis (http://cnsgenomics.com/software/gcta/#COJO).
* --bfile reads individual-level SNP genotype data (in PLINK binary format) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files.
* –prior-sigma the estimated variance of the non-zero mediated effects for each molecular trait on the complex trait.  It can be computed by the variance of the estimated SMR effects at the nominal significance level (i.e., 0.05) adjusting for the estimation errors. 
* --out saves the estimation of prior proportions from the OPERA stage 1 analysis in .pi file (text format, see below example).

# run OPERA for stage 2 analysis and heterogeneity analysis
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --prior-pi 0.8,0.09,0.09,0.02 –prior-sigma 0.02,0.02 --out myopera --thread-num 3

* --prior-pi the estimated global proportions of each configuration from the stage 1 analysis. 
* --out saves the posterior probability of association for each possible association hypotheses in .ppa file (text format, see below example).

* The heterogeneity test (i.e., multi-exposure HEIDI) will be automatically performed for any combinatorial associations passed a PPA threshold (0.8 as default). Users can change the PPA threshold by --thresh-PP to define significant associations to perform further heterogeneity test. The heterogeneity test can be turned off by specifying --heidi-off if not interested.

# Other parameters
> opera --besd-flist mylist --gwas-summary mygwas.ma --bfile mydata --extract-exposure-probe myexposure --outcome-wind 1000 --joint-smr --extract-target-cojo-snps mycojo --extract-GWAS-loci myloci --prior-pi 0.8,0.09,0.09,0.02 –prior-sigma 0.02,0.02 --out myopera --thread-num 3 

* --extract-exposure-probe	extracts a subset of exposure probes for analysis
* --outcome-wind specify the window around each GWAS loci
* --extract-GWAS-loci extracts a subset of GWAS loci for analysis












