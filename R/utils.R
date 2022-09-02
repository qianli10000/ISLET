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
