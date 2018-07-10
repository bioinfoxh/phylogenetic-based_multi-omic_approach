# phylogenetic-based_multi-omic_approach

## Prerequisite
The codes for phylogenetic-based multi-omic approach are implemented in R and using the R packages mvSLOUCH and permute. Before running the codes, please install the R environment and the mvSLOUCH package first. We recommend to install them following the instructions of CRAN:
https://cran.r-project.org/

## Implementation 
### 1. Create a folder “Data” in the home directory of the codes 
### 2. Put the input datasets in the folder “Data”
The input datasets include three files: i) the information of gene expression and DNA methylation of normal samples in a csv file; ii) the information of gene expression and DNA methylation of patient samples in a csv file; iii) the nexus file of phlogenetic tree of the gene family. 

In the "Data folder", we have provided the example data for the phylogenetic constrained multi-omic analysis includes the gene expression and gene body DNA methylation of TGFbeta gene family in asthmatics:
"asthmaticsNormal.csv" includes the means and variances of gene expression and gene body methylation in normal samples;
"asthmaticsPatient.csv" includes the means and variances of gene expression and gene body methylation in patient samples;
"tree.nj.nexus" provides the phylogenetic tree of the gene family. 

### 3. Run the analysis

The script "toRun.R" performs the whole phylogenetic constrained multi-omic analysis using the up-mentioned example datasets as input and the corresponding result files will be automatically generated in the "Data" folder. 

To run your own analysis, just replace the input datasets with your own files and change the prefix "asthmatics" of the files to the disease name that you use. Please also remember to change the corresponding disease name in the "toRun.R": disease<-"asthmatics. After that, just run toRun.R in R environment in the home directory of the codes.

