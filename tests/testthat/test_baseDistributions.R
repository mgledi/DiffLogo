library(DiffLogo)
context("Base Distributions")

a=c(1,0,0,0)
b=c(0,1,0,0)
c=c(0.5,0.5,0,0)
d=c(0,0,0.5,0.5)
e=c(0,0.5,0.5,0)
f=c(.25,.25,.25,.25)


test_that("sumOfAbsProbabilityDifferences", {
   expect_equal(normalizedDifferenceOfProbabilities(a,a), c(.25,.25,.25,.25))
   expect_equal(normalizedDifferenceOfProbabilities(a,b), c(.5,-.5,0,0))
   expect_equal(normalizedDifferenceOfProbabilities(b,a), c(-.5,.5,0,0))
   expect_equal(normalizedDifferenceOfProbabilities(c,d), c(.25,.25,-.25,-.25))
   expect_equal(normalizedDifferenceOfProbabilities(c,e), c(.5,0,-.5,0))
   expect_equal(normalizedDifferenceOfProbabilities(e,f), c(-.25,.25,.25,-.25))
})

test_that("differenceOfICs", {
   expect_equal(differenceOfICs(a,a), c(.25,.25,.25,.25))
   expect_equal(differenceOfICs(a,b), c(.5,-.5,0,0))
   expect_equal(differenceOfICs(b,a), c(-.5,.5,0,0))
   expect_equal(differenceOfICs(c,d), c(.25,.25,-.25,-.25))
   expect_equal(differenceOfICs(c,e), c(.5,0,-.5,0))
   expect_equal(differenceOfICs(e,f), c(.0,.5,.5,0))
   expect_equal(differenceOfICs(f,e), c(.0,-.5,-.5,0))
})