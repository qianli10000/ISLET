###function to read-in and check data from case and control: observed expression, proportion, and sample-to-subject relationship

dataPrepSlope<-function(dat_se){
    message("Begin: working on data preparation as the input for ISLET algorithm.")
#    if (missing(case_dat_se) || missing(ctrl_dat_se))
#        stop("SummarizedExperiment objects from both groups are needed.")
#    if (length(intersect(colData(case_dat_se)[,1], colData(ctrl_dat_se)[,1]))>0)
#        stop("Subject IDs across case group and control group must be unique.")
#    if (ncol(colData(case_dat_se))!=ncol(colData(ctrl_dat_se)))
#        stop("Case group and control group must have the same number of cell types.")

    #subject id between cases and ctrls should also be unique, check and implement this later
    #check for negative values, implement later
    if (!is(dat_se, "SummarizedExperiment"))
        stop("The input dataset must be a SummarizedExperiment object.")
    if (length(unique(colData(dat_se)$group))!=2)
        stop("There must be two groups (case/ctrl) in the input SummarizedExperiment object.")
    if (unique(colData(dat_se)$group)[1]!="case" || unique(colData(dat_se)$group)[2]!="ctrl")
        stop("The names for the two groups in comparison should be
             labeled as `case` and `ctrl` in the input SummarizedExperiment object.")

    idx <- which(colData(dat_se)$group == "case")

    case_dat_se <- SummarizedExperiment(assays=list(counts=assays(dat_se)$counts[, idx]),
                                       colData=colData(dat_se)[idx, -1])
    ctrl_dat_se <- SummarizedExperiment(assays=list(counts=assays(dat_se)$counts[, -idx]),
                                       colData=colData(dat_se)[-idx, -1])




    #K = number of cell types
    K <- ncol(colData(case_dat_se))-2

    #N1 = number of samples for group 1
    N1 <- ncol(assays(case_dat_se)$counts)
    #N1 = number of samples for group 2
    N2 <- ncol(assays(ctrl_dat_se)$counts)
    #NS = total number of Samples for group 1&2
    NS <- N1 + N2
    #NU = total number of Unique subjects for group 1&2
    caseUN <- length(unique(colData(case_dat_se)[, 1]))
    ctrlUN <- length(unique(colData(ctrl_dat_se)[, 1]))
    NU <- caseUN + ctrlUN
    case_age<-colData(case_dat_se)[, 2] ## The first two columns are: subject ID and sample age.
    ctrl_age<-colData(ctrl_dat_se)[, 2]


    X_sub1 <- as.matrix(rbind(colData(case_dat_se)[, -(seq_len(2))], colData(ctrl_dat_se)[, -(seq_len(2))]))
    X_sub2 <- rbind(matrix(1, nrow = N1, ncol = K), matrix(0, nrow = N2, ncol = K))*X_sub1
    X_age <- c(case_age, ctrl_age)
    X_sub3 <- X_sub1*X_age
    X_sub4 <- rbind(matrix(1, nrow = N1, ncol = K), matrix(0, nrow = N2, ncol = K))*X_sub1*X_age
    ## This is difference in slope between two groups
    X_sub4 <- X_sub4
    X_0 <- cbind(X_sub1, X_sub2, X_sub3, X_sub4)
    X_list <- lapply(1, function(x){return(X_0)})
    X <- bdiag(X_list)

    #obtain a vector of unique subject IDs, for all, to use later
    sub_id <- c(colData(case_dat_se)[, 1], colData(ctrl_dat_se)[, 1])

    propm <- as.matrix(rbind(colData(case_dat_se)[, -(seq_len(2))], colData(ctrl_dat_se)[, -(seq_len(2))]))
    # propd = apply(propm, MARGIN = 2, makea, sub_id = sub_id, X = X, NU = NU, simplify = F)
    propd <- apply(X = propm, MARGIN = 2, FUN = makea,
                   ind_id = sub_id, datX = X, aNU = NU, simplify = FALSE)

    A_0 <- do.call(cbind, propd)
    #A_list=lapply(1,function(x){return(A_0)})
    A<-bdiag(A_0)

    CT<-colnames(propm)

    datuse <- inputSet(exp_case = assays(case_dat_se)$counts,
                  exp_ctrl = assays(ctrl_dat_se)$counts,
                  X = X,
                  A = A,
                  K = K,
                  NS = NS,
                  NU = NU,
                  case_num = caseUN,
                  ctrl_num = ctrlUN,
                  CT = CT,
                  SubjectID=sub_id,
                  type='slope'
                  )
    message("Complete: data preparation for ISLET.")
    return(datuse)
}

