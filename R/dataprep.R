
###function to read-in and check data from case and control: observed expression, proportion, and sample-to-subject relationship
dataprep<-function(case_obs, ctrl_obs, case_prop, ctrl_prop, case_subject, ctrl_subject){
  message("Begin: working on data preparation as the input for ISLET algorithm.")
  if (missing(case_obs))
    stop("Observed value matrix for case (group 1) must be provided. E.x. (gene by sample) gene expression data matrix for group 1")
  if (missing(ctrl_obs))
    stop("Observed value matrix for control (group 2) must be provided. E.x. (gene by sample) gene expression data matrix for group 2")
  if (missing(case_prop))
    stop("Cell-type proportion matrix for case (group 1) must be provided. E.x. (sample by cell-type) proportion data matrix for group 1")
  if (missing(ctrl_prop))
    stop("Cell-type proportion matrix for control (group 2) must be provided. E.x. (sample by cell-type) proportion data matrix for group 2")
  if (missing(case_subject))
    stop("Sample to subject relationship, for case (group 1), must be provided.")
  if (missing(ctrl_subject))
    stop("Sample to subject relationship, for case (group 1), must be provided.")

  #subject id between cases and ctrls should also be unique, check and implement this later
  #check for negative values, implement later

  #K = number of cell types
  K = ncol(case_prop)

  #N1 = number of samples for group 1
  N1 = ncol(case_obs)
  #N1 = number of samples for group 2
  N2 = ncol(ctrl_obs)
  #NS = total number of Samples for group 1&2
  NS = N1 + N2
  #NU = total number of Unique subjects for group 1&2
  NU = length(unique(case_subject$subject_id)) + length(unique(ctrl_subject$subject_id))



  X_sub1 = rbind(case_prop, ctrl_prop)
  X_sub2 = rbind(matrix(1, nrow = N1, ncol = K), matrix(0, nrow = N2, ncol = K))*X_sub1
#  X_sub2 = X_sub2[,1:para]
  X_0 = cbind(X_sub1, X_sub2)
  X_list = lapply(1,function(x){return(X_0)})
  X = bdiag(X_list)

  #obtain a vector of unique subject IDs, for all, to use later
  sub_id = c(case_subject$subject_id, ctrl_subject$subject_id)

  propm = rbind(case_prop, ctrl_prop)
 # propd = apply(propm, MARGIN = 2, makea, sub_id = sub_id, X = X, NU = NU, simplify = F)
  propd = apply(X = propm, MARGIN = 2, FUN = makea, ind_id = sub_id, datX = X, aNU = NU, simplify = FALSE)

  A_0 = do.call(cbind, propd)
  #A_list=lapply(1,function(x){return(A_0)})
  A=bdiag(A_0)

  datuse = list(exp_case = case_obs,
                exp_ctrl = ctrl_obs,
                X = X,
                A = A,
                K = K,
                NS = NS,
                NU = NU,
                case_num = length(unique(case_subject$subject_id)),
                ctrl_num = length(unique(ctrl_subject$subject_id))
                )
  message("Complete: data preparation for ISLET.")
  return(datuse)
}

