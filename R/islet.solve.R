

###function to run ISLET, using parallel computing
#ipc is the index of parallel computing for
islet.solve<-function(input){

    if(.Platform$OS.type == "unix") {
    ## do some parallel computation under Unix
      G = nrow(input$exp_case)
      res = bplapply(seq_len(G), islet.est.bp, datuse = input)
  }
  else {
    ## This will be windows
    ## Use serial param or do not use any parallel functions, just use ‘lapply’
    ## result should be of the same “type” from both the if and else statements.
    Yall=as.matrix(cbind(input$exp_case, input$exp_ctrl))
    #make Yall a list
    Yall.list <- split(t(Yall), rep(seq_len(ncol(t(Yall))), each = nrow(t(Yall))))
    #res = lapply(X = Yall.list, FUN = islet.est.win, datuse = input)
    nworkers=min(detectCores()-1,15)
    cl <- makeCluster(nworkers)
#    clusterExport(cl,list('ss'))
#    clusterEvalQ(cl, {
#      require(Matrix)})
    res =parLapply(cl, X=Yall.list, islet.est.win, datuse = input)
    stopCluster(cl)
  }

  return(res)
}

