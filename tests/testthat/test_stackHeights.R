library(DiffLogo)
context("Stack Heights")

a=c(1,0,0,0)
b=c(0,1,0,0)
c=c(0.5,0.5,0,0)
d=c(0,0,0.5,0.5)
e=c(0,0.5,0.5,0)
f=c(.25,.25,.25,.25)
test_that("sumOfAbsProbabilityDifferences", {
   expect_equal(sumOfAbsProbabilityDifferences(a,a)$height, 0)
   expect_equal(sumOfAbsProbabilityDifferences(a,b)$height, 2)
   expect_equal(sumOfAbsProbabilityDifferences(b,a)$height, 2)
   expect_equal(sumOfAbsProbabilityDifferences(c,d)$height, 2)
   expect_equal(sumOfAbsProbabilityDifferences(d,c)$height, 2)
   expect_equal(sumOfAbsProbabilityDifferences(c,e)$height, 1)
   expect_equal(sumOfAbsProbabilityDifferences(e,c)$height, 1)
})

test_that("sumOfAbsICDifferences", {
   expect_equal(sumOfAbsICDifferences(a,a)$height, 0)
   expect_equal(sumOfAbsICDifferences(a,b)$height, 4)
   expect_equal(sumOfAbsICDifferences(b,a)$height, 4)
   expect_equal(sumOfAbsICDifferences(c,d)$height, 2)
   expect_equal(sumOfAbsICDifferences(d,c)$height, 2)
   expect_equal(sumOfAbsICDifferences(c,e)$height, 1)
   expect_equal(sumOfAbsICDifferences(e,c)$height, 1)
})

test_that("lossOfAbsICDifferences", {
   expect_equal(lossOfAbsICDifferences(a,a)$height, 0)
   expect_equal(lossOfAbsICDifferences(a,b)$height, 200)
   expect_equal(lossOfAbsICDifferences(b,a)$height, 200)
   expect_equal(lossOfAbsICDifferences(c,d)$height, 200)
   expect_equal(lossOfAbsICDifferences(d,c)$height, 200)
   expect_equal(lossOfAbsICDifferences(c,e)$height, 100)
   expect_equal(lossOfAbsICDifferences(e,c)$height, 100)
})

test_that("shannonDivergence", {
   expect_equal(shannonDivergence(a,a)$height, 0)
   expect_equal(shannonDivergence(a,b)$height, 1)
   expect_equal(shannonDivergence(b,a)$height, 1)
   expect_equal(shannonDivergence(c,d)$height, 1)
   expect_equal(shannonDivergence(d,c)$height, 1)
   expect_equal(shannonDivergence(c,e)$height, 0.5)
   expect_equal(shannonDivergence(e,c)$height, 0.5)
})

test_that("informationContent", {
       expect_equal(informationContent(a)$height, 2)
       expect_equal(informationContent(c)$height, 1)
       expect_equal(informationContent(f)$height, 0)
})
  