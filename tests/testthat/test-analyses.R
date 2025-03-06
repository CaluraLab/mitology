library(mitology)
library(testthat)

test_that("getGeneSets work", {
    dbnames <- c("MitoCarta", "Reactome", "GO-CC", "GO-BP")
    pname <- sample(dbnames, 1)
    
    res <- getGeneSets(database = pname)
    
    expect_true(is(res, "list"))
    expect_type(res[[1]], "character")
})
