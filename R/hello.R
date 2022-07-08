library(iSEE)
library(iSEEu)

#' Title
#'
#' @param se_name Name of summarized experiment
#'
#' @return an iSEE app
#' @export
#'
#' @examples
#' isee_data("test")
isee_data <- function(se_name = "test") {
   iSEE(readRDS(paste0("data/", se_name, ".rds")))
}

isee_mini <- function(se) {

  gat<-GeneAnnoTable(PanelWidth=8L)
  got<-GOTable(PanelWidth=8L)
  rst <- RowDataTable(RowSelectionSource="GOTable1")

  iSEE(se, initial=list(got, rst, VolcanoPlot(PanelWidth=6L), gat, FeatureSetTable(PanelWidth=6L)))
}

isee_mini(readRDS(paste0("data/", "test", ".rds")))
