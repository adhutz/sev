
gat<-GeneAnnoTable(PanelWidth=8L)
got<-GOTable(PanelWidth=8L)
rst <- RowDataTable(PanelWidth = 12L)

#' Explore created summarized experiments uploaded in the data folder of my github repo.
#' Run list_data() to see a list of all available datasets.
#'
#' @param se_name Name of summarized experiment
#'
#' @return Returns an iSEE object
#' @export
#'
#' @examples
#' isee_data("test") ## runs iSEE to explore the selected dataset
isee_data <- function(se_name = "test") {
   iSEE(readRDS(paste0("data/", se_name, ".rds")))
}

#' Title
#'
#' @param se summarized experiment object
#'
#' @return an iSEE app
#' @export
isee_mini <- function(se) {
  iSEE(se, initial=list(got, rst, VolcanoPlot(PanelWidth=6L), gat, FeatureSetTable(PanelWidth=6L)))
}


#' Title
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

sev::isee_data()
