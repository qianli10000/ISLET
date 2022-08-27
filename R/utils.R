###function to obtain sum of squares
ss<-function(x){
  a<-sum(x^2)
  return(a)
}

###function to obtain sum of squares
colss<-function(x){
    a<-colSums(x^2)
    return(a)
}

###function to make the design matrix [A] for random effect
#updated on 05/31/2022 to reflect the change in ID order
#user should sort their data by subject ID
#makea <- function(onectprop, ind_id = sub_id, datX = X, aNU = NU){
makea <- function(onectprop, ind_id, datX , aNU){
   lp <- split(onectprop, ind_id)
  a1 <- matrix(0, nrow=nrow(datX), ncol=aNU)
#  ct_sub=table(sub_id)[as.character(unique(sub_id))]
#  lp=lp[names(ct_sub)]
  chk <- unique(ind_id) #chk should have the length of NU
  lp<-lp[as.character(chk)]
  count <- rep(0, length(chk))
  for(i in seq_len(aNU)){
    tmp <- sum(ind_id == chk[i])
    count[i] <- tmp
  }

  for(i in seq_len(aNU)){
    s <- 1+sum(count[0:(i-1)])
    e <- sum(count[seq_len(i)])
    a1[s:e, i] <- lp[[i]]
  }
  return(a1)
}


LRT<-function(llk_f, llk_0, df){
    test.stat<-llk_f-llk_0
    p<-pchisq(as.numeric(test.stat), df, lower.tail=FALSE)
    return(p)
}


#clexp<-function(){
#    #clusterExport(cl, list('colss'))
#    clusterEvalQ(cl,  {
#        library(Matrix)})
#}

inputSet <- setClass("inputSet", slots=c(exp_case="data.frame",
                                                 exp_ctrl="data.frame",
                                                 X="Matrix",
                                                 A="Matrix",
                                                 K="numeric",
                                                 NS="integer",
                                                 NU="numeric",
                                                 case_num="numeric",
                                                 ctrl_num="numeric",
                                                 CT="character",
                                                 SubjectID="numeric",
                                                 type="character"
)
)

setMethod("show", "inputSet",
          function(object){
              cat("First couple of elements from cases and controls:", "\n")
              print(object@exp_case[seq_len(3), seq_len(6)])
              print(object@exp_ctrl[seq_len(3), seq_len(6)])
              cat("Design matrices hidded.", "\n")
              cat("Total cell type number:", "\n")
              print(object@K)
              cat("Cell type categories:", "\n")
              print(object@CT)
              cat("Total sample number and subject number:", "\n")
              print(c(object@NS, object@NU))
              cat("Total case number and ctrl number:","\n")
              print(c(object@case_num, object@ctrl_num))
              cat("First several subject ID for the samples:", "\n")
              print(object@SubjectID[seq_len(10)])
              cat("Data preparation type (intercept/slope):", "\n")
              print(object@type)
          }
)


outputSol <- setClass("outputSol", slots=c(case.ind.ref="list",
                                             ctrl.ind.ref="list",
                                             mLLK="numeric"
                                             )
                      )


setMethod("show", "outputSol",
          function(object){
              cat("For cell type 1, the first several individual-specific reference from the case group:", "\n")
              print(object@case.ind.ref[[1]][seq_len(3), seq_len(5)])
              cat(paste0(nrow(object@case.ind.ref[[1]]), " by ", ncol(object@case.ind.ref[[1]]),
                         " (Gene by Individual) reference panel is repeated for all ", length(object@case.ind.ref), " cell types." ), "\n")
              cat("For cell type 1, the first several individual-specific reference from the control group:", "\n")
              print(object@ctrl.ind.ref[[1]][seq_len(3), seq_len(5)])
              cat(paste0(nrow(object@ctrl.ind.ref[[1]]), " by ", ncol(object@ctrl.ind.ref[[1]]),
                         " (Gene by Individual) reference panel is repeated for all ", length(object@ctrl.ind.ref), " cell types." ), "\n")
              cat("First several Log-likelihoods stored for downstream test:", "\n")
              print(object@mLLK[seq_len(3)])
          }
)
