# ISnorm
ISnorm is a method implemented in R for normalizing single-cell RNA sequencing (scRNA-seq) data. We will demonstrate how to normalize scRNA-seq data using ISnorm in this tutorial.

## Installation
ISnorm requires R 3.4.3 or higher, which is available in http://www.r-project.org, and R package dbscan, which can be installed by running following command in R terminal:
```{r }
install.packages("dbscan")
```
The source code of ISnorm can be found in the file `source/ISnorm_function.R`. Put it into your work directory.<br>
We also provide one example dataset from Klein et al. 2015, containing scRNA-seq data of 933 mouse embryonic stem cells. You should also put it into your work directory to run the scripts below.

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
The function `calculate.dis` requires three parameters. The parameter `mat` specifies the input data, which is a numeric matrix containing expression values, with each row representing one gene and each column representing one cell; The parameter `detection_rate` specifies threshold to filter genes; `detection_rate=0.9` means genes without at least 90% cells having nonzero expression will not be included in further analyis. The parameter `ncore` specifies the number of cores used to calculate the distance. A symmetrical matrix containing the distance between genes will be returned.
This is the most time-consuming step in ISnorm. Utilizing multiple cores can reduce the running time. Also you can check the number of genes included in downstream analysis (eg. the number of genes with detection rate larger than 0.9):
```{r }
sum(apply(mat,1,function(x) sum(x>0)/ncol(mat))>0.9)
```
If you have a dataset with low sparsity (eg. more than 5000 genes included in dowstream analysis), you can set `detection_rate` to a higher value such as 0.95 to filter more genes, which can help reduce the running time. In all the cases we have tested, this won't change the results. If you have few genes included in downstream analysis (eg. fewer than 100 genes), you can set `detection_rate` to a lower value such as 0.85 to include more genes. But this is not recommanded as it simply means the quality of your data is poor.


apply(mat,1,function(x) sum(x>0)/ncol(mat))>0.9
