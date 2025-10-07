
.consistencyCheck <- function(
        database = "MitoCarta", nametype = "ENSEMBL", objectType = "list", 
        sections = FALSE, labelNames = "sections") {
    if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
        stop("Database has to be MitoCarta, Reactome, GO-CC or GO-BP")}
    if (!(nametype %in% c("SYMBOL", "ENTREZID", "ENSEMBL"))) {
        stop("The name of genes must be either SYMBOL, ENTREZID or ENSEMBL")}
    if(!(objectType %in% c("list", "dataframe"))){
        stop("Database has to be list or dataframe")}
    if(!is.logical(sections)){ stop("sections must be logical") }
    if(labelNames != "leaves" & labelNames != "sections"){
        stop("Label names must be leaves or sections")}
}

.DBgeneset <- function(database){
    if (database == "MitoCarta") {MCgs
    } else if (database == "Reactome") {RTgs
    } else if (database == "GO-CC") {CCgs
    } else if (database == "GO-BP") {BPgs}}

#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
.geneIDtrans <- function(nametype, genes, database) {
    if (nametype != "ENSEMBL") { 
        if(nametype == "SYMBOL" & database == "MitoCarta"){
            converted_genes <- lapply(genes, function(g){ names(g) })
        } else {
            converted_genes <- lapply(genes, function(g){ mapIds(
                org.Hs.eg.db, keys = g, column = nametype, 
                keytype = "ENSEMBL", multiVals = "first") })}
        names(converted_genes) <- names(genes)
        return(converted_genes)
    } else { genes }
}

.DBtree <- function(database){
    if (database == "MitoCarta") {MCtree
    } else if (database == "Reactome") {RTtree
    } else if (database == "GO-CC") {CCtree
    } else if (database == "GO-BP") {BPtree}}

.plotParams <- function(database){
    if (database == "MitoCarta") {
        DBcol <- c("#823cef", "#2f0967"); hjust <- c(1,0,1,1,0,0,0)
    } else if (database == "Reactome") {
        DBcol <- c("#36ABC9", "#1d5b6b")
        hjust <- c(0,0,0,1,0,0,0,1,1,.5,0,1,1,0,0,1,0,1,0,0)
    } else if (database == "GO-CC") {
        DBcol <- c("#ECC300", "#866f00"); hjust <- c(0,1,0)
    } else if (database == "GO-BP") {
        DBcol <- c("#34CB80", "#1b6a43"); hjust <- c(1,1,rep(0, 10))}
    return(list(DBcol = DBcol, hjust = hjust))}

.matchArguments <- function(dots, defaults) {
    defaults[names(defaults) %in% names(dots)] <- NULL
    c(defaults, dots)
}

# convert DB names
.mitologyName <- function(names, database){
    if(database == "Reactome"){
        names <- gsub(" ", "_", names)
        names <- gsub("(,|\\(|\\))", "", names)
        names <- gsub("\\/", "-", names)
    } else if(database == "GO"){
        names <- gsub(" |,", "_", names)}
    return(names)
}