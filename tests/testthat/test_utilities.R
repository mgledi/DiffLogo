library(DiffLogo)
context("Utilities")


expectation = matrix(nrow = DNA$size, ncol = 1)
colnames(expectation) = 1
rownames(expectation) = DNA$chars
expectation[,1]=c(.4,.3,.2,.1)

expectation2 = matrix(nrow = DNA$size, ncol = 1)
colnames(expectation2) = 1
rownames(expectation2) = DNA$chars
expectation2[,1]=(c(4,3,2,1)+2)/(10+4*2)


input = c("A","A","A","A","C","C","C","G","G","T");
test_that("getPwmFromAlignment", {
   expect_equal(getPwmFromAlignment(input,DNA,0), expectation)
   expect_equal(getPwmFromAlignment(input,DNA,2), expectation2)
})
