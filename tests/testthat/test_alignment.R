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

test_that("SameLengthMatrixesAlign", {
    alignment = localPwmAlignment(pwm1, pwm1)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, pwm1_shifted)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 1)
    expect_equal(alignment[[2]]$direction, 'forward')
    #expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1_shifted, pwm1)
    expect_equal(alignment[[1]]$shift, 1)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'forward')
    #expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, pwm1_revcomp)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'reverse')

    alignment = localPwmAlignment(pwm1_revcomp, pwm1)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'reverse')

    alignment = localPwmAlignment(pwm1, pwm1_revcomp_shifted)
    expect_equal(alignment[[1]]$shift, 1)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'forward')
})

short_pwm_inside = pwm1[,2:4]
short_pwm_shifted = cbind(pwm1[,2:4], pwm1[,1])

test_that("DifferentLengthMatrixesAlign", {
    alignment = localPwmAlignment(pwm1, short_pwm_inside)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 3)
    expect_equal(alignment[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(short_pwm_inside, pwm1)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, short_pwm_shifted)
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 1)
    expect_equal(alignment[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(short_pwm_shifted, pwm1)
    expect_equal(alignment[[1]]$shift, 1)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 0)
    expect_equal(alignment[[2]]$direction, 'forward')
    expect_gte(alignment$divergence, 0)

    alignment = localPwmAlignment(pwm1, revCompPwm(short_pwm_shifted))
    expect_equal(alignment[[1]]$shift, 0)
    expect_equal(alignment[[1]]$direction, 'forward')
    expect_equal(alignment[[2]]$shift, 1)
    expect_equal(alignment[[2]]$direction, 'reverse')
    expect_gte(alignment$divergence, 0)
})

test_that("createDiffLogoObjRespectsAlignPwmObj", {
    diffLogoObj = createDiffLogoObject(pwm1, pwm1_revcomp)
    expect_gt(diffLogoObj$distance, 0)

    diffLogoObj = createDiffLogoObject(pwm1, pwm1_revcomp, align_pwms=T)
    expect_eq(diffLogoObj$distance, 0)
});
