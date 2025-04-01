library(mitology)
library(testthat)

test_that(".DBgeneset works", {
    dbnames <- c("MitoCarta", "Reactome", "GO-CC", "GO-BP")
    dbname <- sample(dbnames, 1)
    res <- .DBgeneset(database = dbname)
    expect_true(is(res, "list"))
    expect_type(res[[1]], "character")
})

test_that(".geneIDtrans works", {
    dbnames <- c("MitoCarta", "Reactome", "GO-CC", "GO-BP")
    dbname <- sample(dbnames, 1)
    idnames <- c("SYMBOL", "ENTREZID", "ENSEMBL")
    idname <- sample(idnames, 1)
    res <- .DBgeneset(database = dbname)
    res <- .geneIDtrans(nametype = idname, genes = res, database = dbname)
    expect_true(is(res, "list"))
    expect_type(res[[1]], "character")
})

test_that(".DBtree works", {
    dbnames <- c("MitoCarta", "Reactome", "GO-CC", "GO-BP")
    dbname <- sample(dbnames, 1)
    res <- .DBtree(database = dbname)
    expect_true(is(res, "phylo"))
})

test_that(".plotParams works", {
    dbnames <- c("MitoCarta", "Reactome", "GO-CC", "GO-BP")
    dbname <- sample(dbnames, 1)
    res <- .plotParams(database = dbname)
    expect_true(is(res, "list"))
    expect_type(res[[1]], "character")
    expect_type(res[[2]], "double")
})
