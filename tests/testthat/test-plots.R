library(mitology)
library(testthat)

# test_that("mitoTreeHeatmap work", {
#     n <- length(names(MCgs)) * 5
#     rmatrix <- matrix(rnorm(n, 0), ncol = 5)
#     rownames(rmatrix) <- names(MCgs)
#     colnames(rmatrix) <- paste0("Sample_", seq_len(5))
#     pres <- mitoTreeHeatmap(data = rmatrix, database = "MitoCarta")
#     expect_true(is(pres, "ggtree"))
#     expect_true(is(pres, "gg"))
#     expect_true(is(pres, "ggplot"))
# })

test_that("mitoHeatmap work", {
    n <- length(names(MCgs)) * 5
    rmatrix <- matrix(rnorm(n, 0), ncol = 5)
    rownames(rmatrix) <- names(MCgs)
    colnames(rmatrix) <- paste0("Sample_", seq_len(5))
    pres <- mitoHeatmap(data = rmatrix, database = "MitoCarta")
    expect_true(is(pres, "Heatmap"))
})
