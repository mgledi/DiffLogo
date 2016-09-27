library(testthat)

pwm1 = matrix( c(0.6, 0.2, 0.1, 0.1,
                 0.3, 0.3, 0.2, 0.2,
                 0.1, 0.1, 0.2, 0.6,
                 0.15, 0.35, 0.25, 0.25,
                 0.85, 0.05, 0.05, 0.05,
                 0.2, 0.1, 0.3, 0.4), 4, 6)
rownames(pwm1) = c("A", "C", "G", "T")
pwm1_shifted = cbind(pwm1[,2:6], pwm1[,1])
pwm1_revcomp = pwm1[4:1, 6:1]
rownames(pwm1_revcomp) = c("A", "C", "G", "T")
pwm1_revcomp_shifted = cbind(cbind(pwm1_revcomp[,3:6], pwm1_revcomp[,1:1]),
                                                       pwm1_revcomp[,1:1])

ACTG_pwm = matrix( c(1 ,0 ,0, 0,
                     0 ,1 ,0, 0,
                     0, 0 ,0 ,1,
                     0, 0, 1, 0), 4, 4)
TGAC_pwm =  cbind(matrix(c(0.25,0.25,0.25,0.25), 4, 1),
                  ACTG_pwm[4:1, 4:1],
                  matrix(c(0.25,0.25,0.25,0.25), 4, 1) )

test_that("SameLengthMatrixesAlign", {
    alignment = localPwmAlignment(pwm1, pwm1)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, pwm1_shifted)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 1)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1_shifted, pwm1)
    expect_equal(alignment$vector[[1]]$shift, 1)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, pwm1_revcomp)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')
    expect_equal(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1_revcomp, pwm1)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')

    alignment = localPwmAlignment(pwm1, pwm1_revcomp_shifted)
    expect_equal(alignment$vector[[1]]$shift, 1)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
})

short_pwm_inside = pwm1[,2:4]
short_pwm_shifted = cbind(pwm1[,2:4], pwm1[,1])

test_that("DifferentLengthMatrixesAlign", {
    alignment = localPwmAlignment(pwm1, short_pwm_inside)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 3)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(short_pwm_inside, pwm1)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, short_pwm_shifted)
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 1)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(short_pwm_shifted, pwm1)
    expect_equal(alignment$vector[[1]]$shift, 1)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, revCompPwm(short_pwm_shifted))
    expect_equal(alignment$vector[[1]]$shift, 0)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 1)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(ACTG_pwm, TGAC_pwm)
    expect_equal(alignment$vector[[1]]$shift, 1)
    expect_equal(alignment$vector[[1]]$direction, 'forward')
    expect_equal(alignment$vector[[2]]$shift, 0)
    expect_equal(alignment$vector[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)
})

test_that("createDiffLogoObjRespectsAlignPwmObj", {
    diffLogoObj = createDiffLogoObject(pwm1, pwm1_revcomp)
    expect_gt(diffLogoObj$distance, 0)

    diffLogoObj = createDiffLogoObject(pwm1, pwm1_revcomp, align_pwms=T)
    expect_equal(diffLogoObj$distance, 0, tolerance=1e-5)
});

test_that("createDiffLogoObjAcceptsExtendedPwms", {
    pwm_list = list(pwm1, pwm1_revcomp_shifted)
    alignment = multipleLocalPwmsAlignment(pwm_list)
    extended_pwms = extendPwmsFromAlignmentVector(pwm_list,
                                                  alignment$alignment$vector)
    diffLogoObj = createDiffLogoObject(extended_pwms[[1]],
                                       extended_pwms[[2]], align_pwms=T)
    expect_equal(is.null(diffLogoObj), F)
});

test_that("extendPwmsFromAlignmentVector", {
    alignment = localPwmAlignment(pwm1, short_pwm_shifted)
    extended_pwms_alignment = extendPwmsFromAlignmentVector(
        list(pwm1, short_pwm_shifted),
        alignment$vector)
    expect_equal(alignment$divergence,
                 localPwmAlignment(pwm1, short_pwm_shifted)$divergence,
                 tolerance=1e-5)
    expect_equal(ncol(extended_pwms_alignment$pwms[[1]]),
                 ncol(extended_pwms_alignment$pwms[[2]]))
});

test_that("pwmDistanceMatrixForTwoPwms", {
    two_pwms_distance_matrix = pwmsDistanceMatrix(list(pwm1, short_pwm_shifted),
                                                  diagonal_value=999,
                                                  bottom_default_value=-1)
    expect_equal(ncol(two_pwms_distance_matrix), 2)
    expect_equal(nrow(two_pwms_distance_matrix), 2)

    expect_equal(two_pwms_distance_matrix[[1, 2]],
                 localPwmAlignment(pwm1, short_pwm_shifted)$divergence)
    expect_equal(two_pwms_distance_matrix[[2, 1]],
                 -1)
    expect_equal(two_pwms_distance_matrix[[1, 1]], 999)
    expect_equal(two_pwms_distance_matrix[[2, 2]], 999)
});

test_that("pwmDistanceMatrixForTreePwms", {
    three_pwms_distance_matrix = pwmsDistanceMatrix(list(pwm1,
                                                         short_pwm_shifted,
                                                         ACTG_pwm),
                                                    diagonal_value = Inf,
                                                    bottom_default_value = Inf)
    expect_equal(ncol(three_pwms_distance_matrix), 3)
    expect_equal(nrow(three_pwms_distance_matrix), 3)
    expect_equal(three_pwms_distance_matrix[[1, 2]],
                 localPwmAlignment(pwm1, short_pwm_shifted)$divergence)
    expect_equal(three_pwms_distance_matrix[[1, 3]],
                 localPwmAlignment(pwm1, ACTG_pwm)$divergence)
    expect_equal(three_pwms_distance_matrix[[2, 3]],
                 localPwmAlignment(short_pwm_shifted, ACTG_pwm)$divergence)

    expect_equal(three_pwms_distance_matrix[[2, 1]], Inf)
    expect_equal(three_pwms_distance_matrix[[3, 1]], Inf)
    expect_equal(three_pwms_distance_matrix[[3, 2]], Inf)

    expect_equal(three_pwms_distance_matrix[[1, 1]], Inf)
    expect_equal(three_pwms_distance_matrix[[2, 2]], Inf)
    expect_equal(three_pwms_distance_matrix[[3, 3]], Inf)
});

test_that("minimumMatrixRowCol", {
    mat = matrix(c(5,3,2,
                   2,2,0,
                   3,1,2), ncol=3, byrow=T)
    min_row_col = minimumMatrixElementIndexes(mat)
    expect_equal(min_row_col[[1]], 2)
    expect_equal(min_row_col[[2]], 3)

    mat = matrix(c(5,7,0,
                   1,2,4,
                   3,4,2), ncol=3, byrow=T)
    min_row_col = minimumMatrixElementIndexes(mat)
    expect_equal(min_row_col[[1]], 1)
    expect_equal(min_row_col[[2]], 3)
});

test_that("addLastTreeNodeToDistanceMatrix", {
    mat = matrix(c(0,0,
                   0,0), ncol=2, byrow=T)
    tree_nodes = list(list("pwms"=list(pwm1)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward"))
                            ,"divergence"=2))
                     ,list("pwms"=list(pwm1, ACTG_pwm)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward")
                                          ,list("shift"=1,
                                                "direction"="reverse"))
                            ,"divergence"=3))
                     ,list("pwms"=list(pwm1_revcomp_shifted)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward"))
                            ,"divergence"=4)))
    new_distance_matrix = addLastTreeNodeToDistanceMatrix(mat, tree_nodes)
    expect_equal(new_distance_matrix[[3, 1]], Inf) 
    expect_equal(new_distance_matrix[[1, 3]], localPwmAlignment(
                                                pwm1,
                                                pwm1_revcomp_shifted)$divergence)
});

test_that("joinTwoNodesInAlignmentTree", {
    mat = matrix(c(Inf, 2,1,
                   Inf,Inf,2,
                   Inf,Inf,Inf), ncol=3, byrow=T)
    tree_nodes = list(list("pwms"=list(pwm1)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward"))
                            ,"divergence"=2))
                     ,list("pwms"=list(pwm1, ACTG_pwm)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward")
                                          ,list("shift"=1,
                                                "direction"="reverse"))
                            ,"divergence"=3))
                     ,list("pwms"=list(pwm1_revcomp_shifted)
                          ,"pwms_alignment"=list(
                             "vector"=list(list("shift"=0,
                                                "direction"="forward"))
                            ,"divergence"=4)))
    joined = joinTwoNodesInAlignmentTree(mat, tree_nodes, 3)
    expect_equal(length(joined), 2)
});


test_that("multipleLocalPwmsAlignment", {
    alignment = multipleLocalPwmsAlignment(list(pwm1,
                                                pwm1_revcomp_shifted,
                                                ACTG_pwm))

    pwms_list = list(pwm1, pwm1_revcomp_shifted, ACTG_pwm)
    names(pwms_list) = sapply(1:length(pwms_list), paste)
    alignment = multipleLocalPwmsAlignment(pwms_list)
});


test_that("diffLogoTablePlotsWithoutAligning", {
    pwm_list = list(pwm1, pwm1_revcomp_shifted, ACTG_pwm)
    alignment = multipleLocalPwmsAlignment(pwm_list)
    extended_pwms = extendPwmsFromAlignmentVector(pwm_list,
                                                  alignment$alignment$vector)
    diffLogoTable(extended_pwms, list("", "", ""))
    diffLogoTable(list(pwm1, pwm1_revcomp_shifted, ACTG_pwm), names=1:3,
                  configuration=list("multiple_align_pwms"=T))
    diffLogoTable(list(pwm1, pwm1_revcomp_shifted, ACTG_pwm), names=1:3,
                  configuration=list(enableClustering=T, multiple_align_pwms=T))
});
