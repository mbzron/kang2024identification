## Reproduction of Kang <em>et al</em>. (2024) "[Identification macrophage signatures in prostate cancer by single-cell sequencing and machine learning.](https://link.springer.com/article/10.1007/s00262-024-03633-5)"

This is a reproduction of the analysis described in the study

<b>Kang, Zhen, Yu-Xuan Zhao, Ren Shun Qian Qiu, Dong-Ning Chen, Qing-Shui Zheng, Xue-Yi Xue, Ning Xu, and Yong Wei. "Identification macrophage signatures in prostate cancer by single-cell sequencing and machine learning." <em>Cancer Immunology, Immunotherapy</em> 73, no. 3 (2024): 41. </b>

The original study published its analysis script as a .pdf file. This repo was created for learning purposes and (hopefully) to make the analysis more easily reproducible.

## Data
Input data was downloaded from [www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193337](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193337) (GSE193337_RAW.tar).

To organize the data into the required directory structure, cd into the same directory as GSE193337_RAW and run

```. createDataDirs.sh```

## Install dependencies

Install R package dependencies with

```Rscript install_requirements.r```

or, in an R session, run

```source("install_requirements.r")```.

Install linux dependencies with

```. install_linux_requirements.sh```

## Debugging Installation of FGSEA during clusterProfiler install

Installation of the R-package 'clusterProfiler' may throw errors relating to failed compilation of c++ programs during installation of the dependent package, FGSEA. 

I.e.,

```BiocManager::install("clusterProfiler")```

eventually spits out a large error message culminating in

```ERROR: compilation failed for package ‘fgsea’```.

This was apparently caused by use of the wrong c++ compiler (gnu++11) at some point during the FGSEA installation process.

To fix this, I directly modified /etc/R/Makeconf to force use of a different version of the compiler, gnu++14, by changing

```CXX11STD = -std=gnu++11``` 
to 
```CXX11STD = -std=gnu++14```.

(After the installation completed, I reverted /etc/R/Makeconf to its previous state, with 'CXX11STD = -std=gnu++11'.)

## Run Analysis

To run the full workflow, at the command line, run

```Rscript kang2024_master```

OR, start an R session (e.g., by entering "R" at the command line), and run the line

```> source("kang2024_master")```.

### Parameters

The script "kang2024_master" contains two parameters: ```use_existing_files``` and ```data_subset_proportion```. These parameters are designed to help run the analysis on personal computer or laptop. (Without using these parameters, I encountered problems running the full analysis on a Dell 5310 laptop with 32 GB RAM and x8 intel CORE i7 processors.)

Set ```use_existing_files <- TRUE``` to make use of existing data-dump files (variously written as *.RData, *.rds, *.txt) created during prior executions of the script to avoiding re-running the corresponding analyses.

Set ```data_subset_proportion``` in the range 0 < data_subset_proportion < 1 to run the Tumor/normal classification package 'copykat' on a subset of the data to reduce the computational workload.