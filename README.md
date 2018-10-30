# ISnorm tutorial
Li Lin<br>

ISnorm is a method implemented in R for normalizing single-cell RNA sequencing (scRNA-seq) data by a set of constantly expressed genes across all cells (internal spike-in genes, IS genes). We will demonstrate how to normalize scRNA-seq data using ISnorm in this tutorial.

## Installation
ISnorm requires R 3.4.3 or higher version, which is available in http://www.r-project.org, and R package `dbscan`, which can be installed by running following command in R terminal:
```{r }
install.packages("dbscan")
```
The source code of ISnorm can be found in the file `source/ISnorm_function.R`. Put it into your work directory.<br>
We will use one dataset from [Klein *et al.*](https://linkinghub.elsevier.com/retrieve/pii/S0092867415005000) as example, which is available at GEO database under accession number [GSE65525](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525). Extract `GSM1599494_ES_d0_main.csv` from the tar file to your work directory. This file contains scRNA-seq data of 933 mouse embryonic stem cells and is needed to run the scripts below.

## Normalization
Import the source code of ISnorm and read the example dataset:
```{r }
source("ISnorm_function.R")
mat<-as.matrix(read.csv(file="GSM1599494_ES_d0_main.csv",sep=",",header=F,row.names=1))
```
The example dataset is a UMI count matrix. But generally the inputs can be in many forms, including un-normalized matrix such as UMI count, reads count and transcripts count, or normalized matrix such as rpm , tpm and fpkm (see our article for more details).<br>
The first step of ISnorm is to calculate pairwise distance between genes:
```{r }
gene_dis<-calculate.dis(mat=mat,detection_rate=0.9,ncore=4)
```
The function `calculate.dis` returns a symmetrical matrix containing the distance between genes. It requires 3 parameters. The parameter `mat` specifies expression matrix, which is a numeric matrix containing expression values, with each row representing one gene and each column representing one cell. The parameter `detection_rate` specifies threshold to filter genes; `detection_rate=0.9` means genes without at least 90% cells having nonzero expression will not be included in further analyis. The parameter `ncore` specifies the number of cores used in parallel.<br><br>
This is the most time-consuming step in ISnorm and generally will take several minutes. Utilizing multiple cores can reduce the running time. The running time mainly depends on the number of genes retained after filtering. You can manually check it using the following command:
```{r }
sum(apply(mat,1,function(x) sum(x>0)/length(x))>0.9)
```
If the running time is too long, you can set `detection_rate=0.95` or other higher values to filter more genes. But we recommand `detection_rate=0.9` as it works well for all the datasets we've tested.<br><br>
Next we use DBscan algorithm to predict IS genes:
```{r }
spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5,solution=100)
```
The function `dbscan.pick` returns a list with each element containing a set of candidate IS geneset. It require 3 parameters. The parameter `dis` specifies the output from `calculate.dis`. The parameter `ngene` specifies a series of expected number of IS genes and the parameter `solution` specifies the increasing rate of scanning radius. See our article for detailed description of these two parameters. Generally there is no need to change them as they works well for almost all datasets.<br><br>
Then we normalize the matrix with each candidate set:
```{r }
candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=4)
```
The function `candidate.norm` requires 3 parameters. The parameter `mat` specifies the expression matrix. The parameter `spike_candidate` specifies the output from `dbscan.pick`. The parameter `ncore` specifies the number of cores used.<br><br>

The function `candidate.norm` returns a list containing the normalization results for each candidate set. The output `candidate_res$sf` is a numeric matrix containing the size factors, with each row representing one cell and each column representing size factors estimated by one candidate set. Normalized matrix can be obtained by dividing all the counts of one cell by its size factor. The output `candidate_res$inst` is a numeric matrix containing the instability scores, with each row representing one cell and each column representing instability scores estimated by one candidate set. Instability score can be used to measure the reliability of each candidate set (see our article for more details). The output `candidate_res$spike` is the same as `spike_candidate`.<br><br>
You can briefly check the reliability of each candidate set using following commands:
```{r }
sapply(candidate_res$spike,length)  ##check the number of genes in each candidate set
boxplot(candidate_res$inst)  ##draw a boxplot to see the instability score of cells for each candidate set
apply(candidate_res$inst,2,mean)  ##check the average instability score for each candidate set
```
An appropriate IS geneset can be chosen manually based on the results from `candidate_res` (see our article for more details). For now we've developed a method to choose the best IS geneset based on F-test:
```{r }
ISnorm_res<-opt.candidate(mat=mat,candidate_res=candidate_res,baseline_threshold=0.1,p_value=0.05,switch_check=2)
```
The function `opt.candidate` requires 5 parameters. The parameter `mat` specifies the expression matrix. The parameter `candidate_res` specifies the output from `candidate.norm`. The parameter `baseline_threshold` specifies the threshold of instability score used to choose the baseline geneset (see our article for more details). The parameter `p_value` specifies the threshold of p value for F-test.<br><br> 
The parameter `switch_check` specifies the parameter used for reliability check. By `switch_check=2`, ISnorm will examine the next two candidate genesets after the chosen geneset. If they share no common genes, ISnorm will print a warning message (in this case, the results of ISnorm is not strongly supported and may be sensitive to the value of `baseline_threshold`, see our article for more details). But this parameter won't change the results.<br><br>

The function `opt.candidate` returns a list with 5 elements. The output `ISnorm_res$normalized` contains the normalized matrix. The output `ISnorm_res$size_factor` contains the size factor for each cell. These two outputs can be used as inputs for other scRNA-seq analysis. The output `ISnorm_res$ISgenes` contains the name of IS genes used for normalization. The output `ISnorm_res$inst_cell` contains the instability score for each cell. The output `ISnorm_res$picked` contains the index of the optimized candidate geneset.
