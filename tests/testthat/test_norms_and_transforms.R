context("norms and transforms")

test_that("total sum norms add up", {
    rvec <- exp(rnorm(10))
    expect_equal(sum(norm_to_total(rvec)), 1)
    expect_equal(sum(norm_pseudo(rvec)), 1)
})


test_that("clr vector output is expected", {
    testvec1    <- 1:10
    testvecrand <- runif(10, 0, 1)
    # test clr function and various equivilent definitions of the clr
    expect_equal(clr(testvec1), log(testvec1/prod(testvec1)^(1/length(testvec1))))
    expect_equal(clr(testvec1), log(testvec1/exp(mean(log(testvec1)))))
    expect_equal(clr(testvec1), log(testvec1) - mean(log(testvec1)))

    expect_equal(clr(testvecrand), log(testvecrand/prod(testvecrand)^(1/length(testvecrand))))
    expect_equal(clr(testvecrand), log(testvecrand/exp(mean(log(testvecrand)))))
    expect_equal(clr(testvecrand), log(testvecrand) - mean(log(testvecrand)))

    ## test presence of zero component
    expect_equal(clr(c(0, testvecrand))[1], 0)
    
    ## test that a composition doesn't need to be closed
    expect_equal(clr(testvec1), clr(norm_to_total(testvec1)))
    expect_equal(clr(testvecrand), clr(norm_to_total(testvecrand)))
})

test_that("clr matrix output is expected", {

    testmat1    <- matrix(1:110, 10)
    testmatrand <- matrix(runif(110, 0, 1), 10)

    expect_equal(clr(testmat1, 1), t(apply(testmat1, 1, clr.default)))
    expect_equal(clr(testmat1, 2),   apply(testmat1, 2, clr.default))

    expect_equal(clr(testmatrand, 1), t(apply(testmatrand, 1, clr.default)))
    expect_equal(clr(testmatrand, 2),   apply(testmatrand, 2, clr.default))

    ## test zero counts
    expect_equal(clr(rbind(0, testmatrand, 2))[1,], rep(0, ncol(testmatrand)))

    ## test that a composition doesn't need to be closed
    expect_equal(clr(testmat1, 2),    clr(apply(testmat1, 2, norm_to_total), 2))
    expect_equal(clr(testmatrand, 2), clr(apply(testmatrand, 2, norm_to_total), 2))

})


test_that("clr data.frame output is expected", {

    testdf <- as.data.frame(testmatrand)
    colnames(testmatrand) <- colnames(testdf)
    expect_equal(clr(testdf, 2), clr(testmatrand, 2))
    expect_equal(clr(testdf, 1), clr(testmatrand, 1))
})



test_that("alr vector output is expected", {
    testvec1    <- 1:10
    testvecrand <- runif(10, 0, 1)
    # test alr function and various equivilent definitions of the alr
    expect_equal(alr(testvec1, divcomp=1), log(testvec1/testvec1[1])[-1])
    expect_equal(alr(testvec1, divcomp=1), log(testvec1)[-1] - log(testvec1[1]))
    expect_equal(alr(testvec1, divcomp=5), log(testvec1/testvec1[5])[-5])
    expect_equal(alr(testvec1, divcomp=5), log(testvec1)[-5] - log(testvec1[5]))

    expect_equal(alr(testvecrand, divcomp=1), log(testvecrand/testvecrand[1])[-1])
    expect_equal(alr(testvecrand, divcomp=1), log(testvecrand)[-1] - log(testvecrand[1]))
    expect_equal(alr(testvecrand, divcomp=5), log(testvecrand/testvecrand[5])[-5])
    expect_equal(alr(testvecrand, divcomp=5), log(testvecrand)[-5] - log(testvecrand[5]))



    ## test presence of zero component
    expect_equal(alr(c(0, testvecrand), divcomp=2)[1], 0)
    
    ## test that a composition doesn't need to be closed
    expect_equal(alr(testvec1), alr(norm_to_total(testvec1)))
    expect_equal(alr(testvecrand), alr(norm_to_total(testvecrand)))
})

test_that("alr matrix output is expected", {

    testmat1    <- matrix(1:110, 10)
    testmatrand <- matrix(runif(110, 0, 1), 10)

    expect_equal(alr(testmat1, 1), t(apply(testmat1, 1, alr.default)))
    expect_equal(alr(testmat1, 2),   apply(testmat1, 2, alr.default))

    expect_equal(alr(testmatrand, 1), t(apply(testmatrand, 1, alr.default)))
    expect_equal(alr(testmatrand, 2),   apply(testmatrand, 2, alr.default))

    ## test zero counts
    expect_equal(alr(cbind(0, testmatrand), 1, divcomp=2)[,1], rep(0, nrow(testmatrand)))
    expect_equal(alr(rbind(0, testmatrand), 2, divcomp=2)[1,], rep(0, ncol(testmatrand)))

    ## test that a composition doesn't need to be closed
    expect_equal(alr(testmat1, 2),    alr(apply(testmat1, 2, norm_to_total), 2))
    expect_equal(alr(testmatrand, 2), alr(apply(testmatrand, 2, norm_to_total), 2))

})


test_that("alr data.frame output is expected", {

    testdf <- as.data.frame(testmatrand)
    colnames(testmatrand) <- colnames(testdf)
    expect_equal(alr(testdf, 2), alr(testmatrand, 2))
    expect_equal(alr(testdf, 1), alr(testmatrand, 1))
})



