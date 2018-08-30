# ISnorm tutorial
ISnorm is a method implemented in R for normalizing single-cell RNA sequencing (scRNA-seq) data by a set of constantly expressed genes across all cells (internal spike-in genes, IS genes). We will demonstrate how to normalize scRNA-seq data using ISnorm in this tutorial.

## Installation
ISnorm requires R 3.4.3 or higher version, which is available in http://www.r-project.org, and R package `dbscan`, which can be installed by running following command in R terminal:
```{r }
install.packages("dbscan")
```
The source code of ISnorm can be found in the file `source/ISnorm_function.R`. Put it into your work directory.<br>
We will use one dataset from [Klein *et al.*](https://linkinghub.elsevier.com/retrieve/pii/S0092867415005000) as example, which is available at GEO database under accession number [GSE65525](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525). The file `GSM1599494_ES_d0_main.csv` contains scRNA-seq data of 933 mouse embryonic stem cells. You should also put it into your work directory to run the scripts below.

## Normalization
Let us import ISnorm and read the example dataset:
```{r }
source("ISnorm_function.R")
mat<-as.matrix(read.csv(file="GSM1599494_ES_d0_main.csv",sep=",",header=F,row.names=1))
```
The example dataset is a UMI count matrix. But generally the inputs can be in many forms, including un-normalized matrix such as UMI count, reads count and transcripts count, or normalized matrix such as rpm , tpm and fpkm (see our article for more details).<br>
The first step of ISnorm is to calculate pairwise distance between genes:
```{r }
gene_dis<-calculate.dis(mat=mat,detection_rate=0.9,ncore=4)
```
The function `calculate.dis` returns a symmetrical matrix containing the distance between genes. It requires 3 parameters. The parameter `mat` specifies expression matrix, which is a numeric matrix containing expression values, with each row representing one gene and each column representing one cell. The parameter `detection_rate` specifies threshold to filter genes; `detection_rate=0.9` means genes without at least 90% cells having nonzero expression will not be included in further analyis. The parameter `ncore` specifies the number of cores used to calculate the distance.<br><br>
This is the most time-consuming step in ISnorm. Utilizing multiple cores can reduce the running time. If you have a dataset with low sparsity (eg. more than 5000 genes included after filtering), you can set `detection_rate=0.95` to filter more genes, which can help reduce the running time. But we recommand `detection_rate=0.9` as it works well for all the datasets we've tested.<br><br>
Next we use DBscan algorithm to predict IS genes:
```{r }
spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5,solution=100)
```
The function `dbscan.pick` returns a list with each element containing a set of candidate IS geneset. It require 3 parameters. The parameter `dis` specifies the output from `calculate.dis`. The parameter `ngene` specifies a series of expected number of IS genes and the parameter `solution` specifies the increasing rate of scanning radius. See our article for detailed description of these two parameters. You do not need to change them as they works well for almost all datasets.<br><br>
We normalize the matrix with each candidate set:
```{r }
candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=4)
```
The function `candidate.norm` requires 3 parameters. The parameter `mat` specifies the expression matrix. The parameter `spike_candidate` specifies the output from `dbscan.pick`. The parameter `ncore` specifies the number of cores used.<br><br>

The function `candidate.norm` returns a list containing the normalization results for each candidate set. The output `candidate_res$sf` is a numeric matrix containing the size factors, with each row representing one cell and each column representing size factors estimated by one candidate set. The output `candidate_res$inst` is a numeric matrix containing the instability scores, with each row representing one cell and each column representing instability scores estimated by one candidate set. Instability score can be used to measure the reliability of each candidate set (see our article for more details). The output `candidate_res$spike` is the same as `spike_candidate`.<br><br>
You can check the results using following commands:
```{r }
sapply(candidate_res$spike,length)  ##check the number of genes in each candidate set
boxplot(candidate_res$inst)  ##draw a boxplot to see the instability score of cells for each candidate set
apply(candidate_res$inst,2,mean)  ##check the average instability score for each candidate set
```
An appropriate IS geneset can be chosen manually from the information above. We've also developed a method to choose the best IS geneset automatically:
```{r }
ISnorm_res<-opt.candidate(mat=mat,candidate_res=candidate_res,threshold=0.1,switch_check=2)
```
The function `opt.candidate` requires 4 parameters. The parameter `mat` specifies the expression matrix. The parameter `candidate_res` specifies the output from `candidate.norm`. The parameter `threshold` specifies the threshold for instability socre (see our article for more details). The results of ISnorm is sensitive to `threshold`. We recommand `threshold=0.1` as it will provide correct results in most of the dataset we've tested. The parameter `switch_check` specifies the number of candidate set to check for switch pattern. By `switch_check=2`, ISnorm will check the next 2 candidate sets after the chosen set to see whether genes from the chosen set are also included in them. If genes from the chosen set are not included in any of these 2 candidate sets, ISnorm will print a message suggesting the results may be not reliable (see our article for more details). But this parameter won't change the results.<br><br>

The function `opt.candidate` returns a list with 4 elements. The output `ISnorm_res$normalized` contains the normalized matrix. The output `ISnorm_res$size_factor` contains the size factor for each cell. These two outputs can be used as inputs for other scRNA-seq analysis. The output `ISnorm_res$ISgenes` contains the name of IS genes used for normalization. The output `ISnorm_res$inst_cell` contains the instability score for each cell.
