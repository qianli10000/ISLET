
###function to implement EM algorithm by block of genes in ISLET algorithm
#function here for unix and windows, using lapply, no parallel computing
#Y is a GxN matrix for gene expression
islet.solve.block<-function(Y, datuse){
    #exp_case = as.matrix(datuse@exp_case)
    #exp_ctrl = as.matrix(datuse@exp_ctrl)
    X <- datuse@X
    A <- datuse@A
    K <- datuse@K
    NU <- datuse@NU
    NS <- datuse@NS
#    para<-datuse@para

    #initialization of parameters parameter estimation storage
    B_est<-NULL
    Sig0_est <- NULL
    SigU_est <- NULL
    E_U_est <- NULL
    llk <- NULL
    ##
    Y<-t(Y)
    G<-ncol(Y)
    #  Y=log2(Y+1)

    ####1. Initialization of parameters
    #1.1 cell type profiles AND csDE B parameters
    #B_0 = solve(X,Y)
    B_0 <- Matrix::tcrossprod( Matrix::tcrossprod(solve( Matrix::crossprod(X)), X), t(Y))

    #1.2 error terms
    sig <- mean((Y-X%*%B_0)^2)
    #sig <- 20

    #1.3 missing values
    U_0 <- rep(0, NU*K)

    B_t <- B_0
    #sig_t = rep(sig, 7)
    U_t <- U_0
    sig0_t <- rep(sig, G)
    sigK_t <- rep(sig, K)

    iem <- 1
    diff1 <- 100
    diff2 <- 100
    pp <- 1


    #Sig_U = diag(rep(sigK_t, each = NU))
    Sig_p<-lapply(seq_len(G), function(x, A, sig0_t, sigK_t, NU, Y, X, B_t){
        invSig_U<-Matrix::bdiag(diag(rep(1/sigK_t, each=NU)))
        Sig<-solve( Matrix::crossprod(A)/sig0_t[x]+invSig_U)
        U<- Matrix::tcrossprod( Matrix::tcrossprod(Sig, A),
                                BiocGenerics::t(Y[, x] - Matrix::tcrossprod(X,
                               BiocGenerics::t(B_t[, x])))
                               )/sig0_t[x]
        return(list(Sig=Sig, U=U))
    }, A, sig0_t, sigK_t, NU, Y, X, B_t)
    E_Up<-do.call(cbind, lapply(Sig_p, function(x)x$U))


    while(diff2>0.01 & iem<30){
#        cat("iteration=", iem, "\n")
        iem <- iem + 1
        ####2. E-step
        #observed data COV(Y) = V

        #V = A%*%Sig_U%*%t(A) + diag(rep(sig0_t, 5*600))

        #2.1 E[U|Y]: missing data [U|Y] given observed data
        #invV = solve(V)
        # E_U = mu_p = t(Sig_U) %*% t(A) %*% invV %*% (Y - X %*% B_t)
        # Sig_p = Sig_U - crossprod(Sig_U,t(A)) %*% invV %*% A %*% Sig_U

        # Estimate from last iteration

        E_U <- E_Up
        mu_p <- E_Up
        E_U_frame <- as.data.frame(as.matrix(E_U))
        #2.2 E[t(S)S|Y]
        E_StS <- lapply(seq_len(G), function(x, A, Sig_p, mu_p, X, B_t, Y){
           sum( Matrix::diag(Matrix::tcrossprod( Matrix::tcrossprod(A, Sig_p[[x]]$Sig), A))) +
                sum(( Matrix::tcrossprod(A, BiocGenerics::t(mu_p[, x])) +
                          Matrix::tcrossprod(X, BiocGenerics::t(B_t[, x])) - Y[, x])^2)},
           A, Sig_p, mu_p, X, B_t, Y)
        E_StS <-unlist(E_StS)

        #2.3 E[U_k^T U_k|Y]
        mutra_split <- lapply(Sig_p, function(x){
            sig_p<-split(diag(x$Sig), ceiling(seq_len(NU*K)/NU))
            tra<-unlist(lapply(sig_p, sum))
            return(tra)
        })
        mu_split <- split(E_U_frame, ceiling(seq_along(E_U_frame[, 1])/NU))

        E_UkTUk <- do.call('cbind', mutra_split) + do.call('rbind', lapply(mu_split, colss))

        ####3. M-step
        #3.1 B
        B_tp <-  Matrix::tcrossprod( Matrix::tcrossprod(solve( Matrix::crossprod(X)), X),
                                     BiocGenerics::t(Y- Matrix::tcrossprod(A, BiocGenerics::t(E_U))) )

        #make correction in case B[1:K]<0 or B_tp[(K+1):2K]<0
        #important to bound the estimation to positive values

        # B_tp[1:K,]=ifelse(B_tp[1:K,]<0,0,B_tp[1:K,])
        # B_tp[-(1:K),]=ifelse(B_tp[1:K,]+B_tp[-(1:K),]<0,-B_tp[1:K,],B_tp[-(1:K),])


        #3.2 sigma_0^2
        sig0_tp <- E_StS/(NS)

        #3.3 sigma_k^2
        sigK_tp <- E_UkTUk/(NU)

        ####4. Stopping criteria
        diff1 <- sum(abs(B_tp - B_t)) + abs(sig0_tp - sig0_t) + sum(abs(sigK_tp - sigK_t))

        n1 <- sum(abs(B_tp - B_t))/length(B_tp)
        n2 <- sum(abs(B_tp))/length(B_tp)
        pp <- n1/n2
#        cat("B_sum_val=", n2, "\n")
#        cat("B_change_val=", n1, "\n")
#        cat("B_change_prop=", pp*100,"% \n")

        ####5. Update params
        B_t<-B_tp
        sig0_t <- sig0_tp
        sigK_t <- sigK_tp

        Sig_p<-lapply(seq_len(G), function(x, A, sig0_t, sigK_t, NU, Y, X, B_t){
            Sig_U<- Matrix::bdiag(diag(rep(sigK_t[, x], each=NU)))
            invSig_U<-Matrix::bdiag(diag(rep(1/sigK_t[, x], each=NU)))
            Sig<-solve( Matrix::crossprod(A)/sig0_t[x]+invSig_U)
            U<- Matrix::tcrossprod(Matrix::tcrossprod(Sig, A),
                                   t(as.matrix(Y[, x] -
                                    Matrix::tcrossprod(X, BiocGenerics::t(B_t[, x]))))
                                   )/sig0_t[x]
            return(list(Sig_U=Sig_U, Sig=Sig, U=U))
        }, A, sig0_t, sigK_t, NU, Y, X, B_t)
        Sig_p_all<-NULL
        E_Up<-do.call(cbind, lapply(Sig_p, function(x)x$U))

        diff2 <- sum(abs(E_Up - E_U))/(length(E_U)*mean(colMeans(Y))^2)
#        cat("Random effect diff2=", diff2, "\n")

    }
    # Estimate of fixed effect
    B_est<-cbind(B_est, B_t)
    # Estimate of random effect
    E_U_est<-cbind(E_U_est, E_Up)
    # Estimate of variance Sigma_U,  Sigma_0
    Sig0_est <- cbind(Sig0_est, sig0_t)
    SigU_est <- cbind(SigU_est, sigK_t)

    #calculate LLK
    llk<-lapply(seq_len(G), function(x){
        Sig<- Matrix::tcrossprod( Matrix::tcrossprod(A, Sig_p[[x]]$Sig_U),  A)+
            Matrix::bdiag(diag(sig0_t[x], nrow = nrow(A)))
        l<- Matrix::determinant(Sig)$modulus+
            Matrix::tcrossprod(Matrix::crossprod(Y[, x]-
                      Matrix::tcrossprod(X, BiocGenerics::t(B_t[, x])), solve(Sig)),
                      BiocGenerics::t(Y[, x]-
                           Matrix::tcrossprod(X, BiocGenerics::t(B_t[, x]))) )
        return(-as.numeric(l))
    })

    llk<-unlist(llk)


    #compile return list
    case.m <- B_est[seq_len(K), ]+B_est[K+seq_len(K), ]
    ctrl.m <- B_est[seq_len(K), ]

    #(2) the individual value for case and control, for all cell types. 2 matrices of NU by K.
    rel <- split(as.data.frame(as.matrix(E_U_est)), ceiling(seq_along(E_U_est[, 1])/NU))


    case.indv <- lapply(seq_len(K), function(k){rel[[k]][seq_len(datuse@case_num), ] +
            matrix(rep(case.m[k, ], each = datuse@case_num), nrow=datuse@case_num)})
    ctrl.indv <- lapply(seq_len(K), function(k){rel[[k]][-seq_len(datuse@case_num), ] +
            matrix(rep(ctrl.m[k, ], each = datuse@ctrl_num), nrow=datuse@ctrl_num)})
    names(case.indv) <- names(rel)
    names(ctrl.indv) <- names(rel)

    #(3) Variance for K cell types. 1 vector of length K.
    #'SigU_est' is already to be rendered.

    #(4) Variance for grand residuals. 1 scalar.
    #'Sig0_est' is already to be rendered.

    #(5) the model likelihood. 1 scalar.
    #'llk' is already to be rendered.

    #compile return list
    rval <- list(
        case.m = case.m,
        ctrl.m = ctrl.m,
        case.indv = case.indv,
        ctrl.indv = ctrl.indv,
        var.k = SigU_est,
        var.0 = Sig0_est,
        LLK = llk)
    return(rval)

    message("Complete: parameter estimation from ISLET is complete.")
}

