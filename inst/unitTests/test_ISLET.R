##Unit Test using RUnit

#test dataPrep
test_dataPrep <- function() {
    data("GE600")
    s1 = dataPrep(dat_se=GE600_se)

    checkEquals(s1@type, "intercept")
    checkEquals(s1@K, 6)
    checkEquals(s1@case_num, 100)
    checkEquals(s1@ctrl_num, 100)
    checkEquals(s1@NS, 600)
    checkEquals(s1@NU, 200)

    checkTrue(length(unique(s1@SubjectID))==200)

 #   checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}


#test dataPrepSlope
test_dataPrepSlope <- function() {
    data("GE600age")
    s1 = dataPrepSlope(dat_se=GE600age_se)

    checkEquals(s1@type, "slope")
    checkEquals(s1@K, 6)
    checkEquals(s1@case_num, 100)
    checkEquals(s1@ctrl_num, 100)
    checkEquals(s1@NS, 600)
    checkEquals(s1@NU, 200)

    checkTrue(length(unique(s1@SubjectID))==200)

    #   checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}

