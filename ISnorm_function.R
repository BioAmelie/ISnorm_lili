library(dbscan)
library(parallel)

calculate.dor<-function(xx,yy){
  retain<-xx!=0&yy!=0
  ratio<-xx[retain]/yy[retain]
  output<-sqrt(sum(log2(ratio/median(ratio))^2)/(sum(retain)-1))
  return(output)
}  ##calculate DoR between two vector

estimate.ref<-function(mat,ref_gene){
  mat<-mat[ref_gene,]
  cell_sum<-colSums(mat)
  mat<-mat[,cell_sum!=0]
  cell_sum<-cell_sum[cell_sum!=0]
  ref_expr<-apply(sweep(mat,2,cell_sum/median(cell_sum),FUN="/"),1,function(x) mean(x[x!=0]))
  ref_expr<-sort(ref_expr)
  return(ref_expr)
}  ##calculate reference expression given internal spike-in genes


calculate.dis<-function(mat,detection_rate=0.9,ncore=1,exclude=""){
  mat2<-mat[apply(mat,1,function(x) sum(x>0)/length(x))>detection_rate,]
  mat2<-mat2[!rownames(mat2)%in%exclude,]
  gene_dis<-matrix(0,nrow(mat2),nrow(mat2))
  cl <- makePSOCKcluster(ncore)

  for(i in 2:nrow(mat2)-1){
    gene_dis[i,i:nrow(mat2)]<-parRapply(cl=cl,mat2[i:nrow(mat2),],calculate.dor,yy=mat2[i,])
    gene_dis[i:nrow(mat2),i]<-gene_dis[i,i:nrow(mat2)]
  }
  stopCluster(cl)
  rownames(gene_dis)<-rownames(mat2)
  colnames(gene_dis)<-rownames(mat2)
  return(gene_dis)
}  ##calculate pairwise DoR distance between genes

cluster.count<-function(x){
  output<-numeric()
  for(i in 1:max(x)){
    output[i]<-sum(x==i)
  }
  return(output)
}  ##count cluster size

dbscan.pick<-function(dis,ngene=(1:floor(nrow(dis)/25))*5,solution=100){
  genelist<-rownames(dis)
  ngene<-ngene[ngene<=length(genelist)]
  ngene<-sort(ngene)
  dis<-as.dist(dis)
  xx<-as.vector(dis)
  xx<-xx[xx!=0]
  step<-min(xx)
  ngene_index<-1
  output<-list()
  while(ngene_index<=length(ngene)){
    cluster<-0
    while(max(cluster.count(cluster))<ngene[ngene_index]){
      step<-step*10
      cluster<-dbscan(dis,eps=step)$cluster
    }
    step<-step/10
    cluster<-0
    i<-1
    while(i<=solution&ngene_index<=length(ngene)){
      i<-i+1
      cluster<-dbscan(dis,eps=step+step*9/solution*i)$cluster
      if(max(cluster.count(cluster))<ngene[ngene_index]) next
      ngene_index2<-max(which(max(cluster.count(cluster))>=ngene))
      while(ngene_index<=ngene_index2){
        output[[ngene_index]]<-genelist[cluster==which.max(cluster.count(cluster))]
        ngene_index<-ngene_index+1
      }
    }
  }
  output<-output[c(T,!sapply(2:length(output),function(x) identical(output[[x-1]],output[[x]])))]
  return(output)
}  ##pick internal spike-in genes


estimate.sf<-function(cell,ref_expr){
  cell<-cell[names(ref_expr)]
  if(sum(cell!=0)<=1)
    return(c(size_factor=NA,var=NA,ngene=sum(cell!=0)))
  sf<-cell[cell!=0]/ref_expr[cell!=0]
  ngene<-length(sf)
  sf_es<-median(sf)
  var<-sum(log10(sf/sf_es)^2)/(length(sf)-1)
  output<-c(sf_es,var,length(sf))
  names(output)<-c("size_factor","var","ngene")
  return(output)
}  ##calculate size factor for one cell

candidate.norm<-function(mat,spike_candidate,ncore=1){
  cl <- makeCluster(ncore)
  candidate_ref<-parLapply(cl,spike_candidate,estimate.ref,mat=mat)
  temp<-lapply(candidate_ref,function(x) apply(mat,2,estimate.sf,ref=x))
  stopCluster(cl)
  candidate_res<-list(sf=sapply(temp,function(x) x[1,]),inst=sqrt(sapply(temp,function(x) x[2,])),ngene=sapply(temp,function(x) x[3,]))
  candidate_res$ref<-candidate_ref
  candidate_res$spike<-spike_candidate
  
  return(candidate_res)
}  ##normalization by spike candidate


opt.candidate<-function(mat,candidate_res,baseline_threshold=0.1,p_value=0.05){
  instability<-apply(candidate_res$inst,2,function(x) mean(x[!is.na(x)]))
  if(sum(instability<baseline_threshold)==0){
    baseset<-1
  }
  else{
    baseset<-max((1:length(candidate_res$spike))[instability<0.1])
  }
  inst<-candidate_res$inst
  ngene<-candidate_res$ngene
  pmat<-sapply(baseset:ncol(inst),function(x){
    pf((inst[,baseset]/inst[,x])^2,df1=ngene[,baseset]-1,df2=ngene[,x]-1)
    })
  picked<-max(which(apply(pmat,2,function(x) (sum(x[!is.na(x)]<p_value)/length(x[!is.na(x)]))<p_value)))+baseset-1
  cat("Candidate set",picked,"is chosen.\n")
  
  expr<-sweep(mat,2,candidate_res$sf[,picked],FUN="/")
  inst_cell<-candidate_res$inst[,picked]
  ISgenes<-candidate_res$spike[[picked]]
  return(list(normalized=expr,size_factor=candidate_res$sf[,picked],ISgenes=ISgenes,inst_cell=inst_cell,picked=picked))
}


pois.pvalue<-function(xx,calculate.dor,seed=1,nrep=1000){
  set.seed(seed)
  inst<-xx[length(xx)]
  lambda<-xx[1:(length(xx)-1)]
  smat<-sapply(lambda,rpois,n=nrep)
  null_dis<-apply(smat,1,calculate.dor,yy=lambda)
  return(sum(null_dis>inst)/nrep)
}

cal.pvalue<-function(candidate_res,spike_index=1,seed=1,nrep=1000,ncore=4){
  lambda<-sapply(candidate_res$sf[,spike_index],function(x) x*candidate_res$ref[[spike_index]])
  temp<-rbind(lambda,candidate_res$inst[,spike_index])
  cl <- makeCluster(ncore)
  pvalue<-parApply(cl=cl,temp,2,pois.pvalue,calculate.dor=calculate.dor,seed=seed,nrep=nrep)
  stopCluster(cl)
  return(pvalue)
  
}
