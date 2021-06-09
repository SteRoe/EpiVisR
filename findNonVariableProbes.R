betaFileName <- "f:/roeder/methylation_transposed2.csv" #config$betaFileName
beta <- fread(betaFileName,stringsAsFactors=FALSE,header=TRUE,sep="\t")

beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))

beta_wo_outliers<-removeOutliers3IQR(as.matrix(beta))
beta_wo_outliers<-as.data.frame(beta_wo_outliers[[1]])
beta<-beta_wo_outliers

beta<-removeMultiModelCpGFromBeta(beta,MultiModProbes)
assign("beta",beta,envir=globalenv())
beta.t<-t(beta)
assign("beta.t",beta.t,envir=globalenv())

library(parallel)
library(doParallel)
library(foreach)

#find probes with delta methylation less than 0.05 (5%)
findNonVariableProbes <- function(){
  DList = list()

  foreach (i=1:nrow(beta)) %do% {
    delta = max(beta[i,]) - min(beta[i,])
    if (delta < 0.05) {
      DList = append(DList,rownames(beta)[i])
    }
  }
  return (DList)
}
a=findNonVariableProbes()
saveRDS(a,file="f:/roeder/nonVariableProbes.RDS")
fwrite(a,file="f:/roeder/nonVariableProbes.csv", sep="\t",dec=".")
