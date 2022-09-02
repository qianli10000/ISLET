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




outputSol <- setClass("outputSol", slots=c(case.ind.ref="list",
                                           ctrl.ind.ref="list",
                                           mLLK="numeric"
                                           )
                      )

