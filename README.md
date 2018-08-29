# ISnorm
ISnorm is a method implemented in R for normalizing single-cell RNA sequencing (scRNA-seq) data. We will demonstrate how to normalize scRNA-seq data using ISnorm in this tutorial.

## Installation
ISnorm requires R 3.4.3 or higher, which is available in http://www.r-project.org, and R package dbscan, which can be installed by running following command in R terminal:
```{r }
install.packages("dbscan")
```
The source code of ISnorm can be found in the file `source/ISnorm_function.R`. Put it into your workding directory.<br>
We also provide one example dataset from Klein et al. 2015, containing scRNA-seq data of 933 mouse embryonic stem cells. You should also put it into your work directory to run the scripts in this tutorial.

## Normalization
Let us import ISnorm with required R packages and read the example dataset:
```{r }
library(dbscan)
library(parallel)
source("ISnorm_function.R")
mat<-as.matrix(read.csv(file="GSM1599494_ES_d0_main.csv",sep=",",header=F,row.names=1))
```
The example dataset is a UMI count matrix. But generally the inputs can be in many forms, including un-normalized matrix such as UMI count, reads count and transcripts count, or normalized matrix such as rpm , tpm and fpkm (see our article for more details).<br>
The first step of ISnorm is to calculate pairwise distance between genes:
```{r }
gene_dis<-calculate.dis(mat=mat,detection_rate=0.9,ncore=4)
```
The function `calculate.dis`
