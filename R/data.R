#' Example expression data.
#'
#' This is an example dataset containing gene expression values (in normalized
#' counts) of 40 ovarian cancer (OVC) patients extracted
#' from the Cancer Genome Atlas (TCGA) database.
#' This dataset should be used only with example purpose.
#' RNA sequencing OVC data were retrieved using
#' \code{\link[curatedTCGAData]{curatedTCGAData}} package. Data were then
#' normalized with the \code{\link[EDASeq]{betweenLaneNormalization}} function.
#' To lighten the dataset, the \code{\link[signifinder]{consensusOVSign}}
#' function was computed, which return 4 different scores, one for each OVC
#' subtype (Chen et al, 2018, Clinical Cancer Research) and the 10 samples
#' with the highest scores were selected for each subgroup.
#' Further, only the mitochondrial genes included in mitology were kept.
#' Finally, the log fold change of the IMR versus the PRO samples were computed.
#' Further details in mitology/inst/scripts/howToGenerateOvse.Rmd.
#'
#' @docType data
#'
#' @usage data(ovse)
#'
"ovse"


#' Mitochondrial-related genes
#'
#' Here are listed all the mitochondrial-related genes collected in mitology.
#' Columns are:
#' \itemize{
#'  \item the ENSEMBL id
#'  \item the SYMBOL id
#'  \item the databases from where the gene was collected (IMPI, MitoCarta3.0, 
#' MSeqDR and/or GO-CC)
#'  \item the localization class:
#'  \itemize{
#' \item “Ancillary” for genes encoding ancillary mitochondrial proteins with 
#' no evidence of mitochondrial localization but reported to affect 
#' mitochondrial function or morphology;
#' \item “Mitochondrial only” for genes encoding proteins with evidence of 
#' mitochondrial localization and not in additional cellular compartments;
#' \item “Multi-localized” for genes encoding proteins with evidence of 
#' mitochondrial localization and in additional cellular compartments.
#' }}
#' 
#' @docType data
#' 
#' @usage data(MitoGenesDB)
#' 
"MitoGenesDB"
