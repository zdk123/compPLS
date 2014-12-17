context("partial least squares discriminant analysis")


irisdat <- clr(iris[,1:4], 1)
y       <- iris[,5]
mod <- plsDA(irisdat, y)

test_that("plsDA returns correct class", {

    expect_that(mod, is_a('plsda'))
    expect_that(mod, is_a('mvr'))

})
