## Reproduction of Kang <em>et al</em>. (2024) "[Identification macrophage signatures in prostate cancer by single-cell sequencing and machine learning.](https://link.springer.com/article/10.1007/s00262-024-03633-5)"

This is a reproduction of the analysis described in the study

<b>Kang, Zhen, Yu-Xuan Zhao, Ren Shun Qian Qiu, Dong-Ning Chen, Qing-Shui Zheng, Xue-Yi Xue, Ning Xu, and Yong Wei. "Identification macrophage signatures in prostate cancer by single-cell sequencing and machine learning." <em>Cancer Immunology, Immunotherapy</em> 73, no. 3 (2024): 41. </b>

The original study published its analysis script as a .pdf file. This repo was created for learning purposes and (hopefully) to make the analysis more easily reproducible.

## Data
Input data was downloaded from [www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193337](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193337) (GSE193337_RAW.tar).

To organize the data into the directory structure expected by the analysisR-script, cd into the same directory as GSE193337_RAW and run

```. createDataDirs.sh```

## Install dependencies

Install R package dependencies with

```source("install_requirements.r")```

Install linux dependencies with

```. install_linux_requirements.sh```

## Debugging Installation of FGSEA during clusterProfiler install

Installation of the R-package 'clusterProfiler' threw errors that seemed to relate to failed compilation of c++ programs during installation of the dependent R-package, FGSEA. 

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

```> source("kang2024_master")```