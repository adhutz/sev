
gat<-GeneAnnoTable(PanelWidth=8L)
got<-GOTable(PanelWidth=8L)
rst <- RowDataTable(PanelWidth = 12L)

#' Explore summarized experiments
#' Datasets collected in the data folder of the github repo are available.
#'
#' @param se_name Name of summarized experiment
#'
#' @return Returns an iSEE object
#'
#' @import iSEE
#' @export
#'
#' @examples
#' isee_data(se = sev::data(test)) ## runs iSEE to explore the selected test dataset
isee_data <- function(se) {
  iSEE::iSEE(se)
}

#' Explores your own summarized experiment
#'
#' @param se summarized experiment object
#'
#' @return an iSEE app
#' @export
isee_mini <- function(se) {
  iSEE::iSEE(se, initial=list(got, rst, VolcanoPlot(PanelWidth=6L), gat, FeatureSetTable(PanelWidth=6L)))
}

#' List available datasets
#'
#' @return list of available summarized experiments
#' @export
#'
#' @examples
#' list_data()[1] ## returns the first file of the data folder
#'
#' isee_data(list_data()[1]) ## Starts iSEE for the first dataset
list_data <- function() {
  gsub(".rds","",list.files("data"))
}

#' Install Bioconductor dependencies
#'
#' @return None
#' @export
#'
sev_depend <- function(){
  BiocManager::install("iSEE")
  BiocManager::install("iSEEu")
}

