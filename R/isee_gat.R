library(iSEE)
setClass("GeneAnnoTable", contains="RowDataTable",
         slots=c(
           IDColumn="character_OR_NULL",
           IDType="character",
           Organism="character",
           AnnoBoxOpen="logical"
         )
)

allowable <- c("ENSEMBL", "SYMBOL", "ENTREZID")
setValidity2("GeneAnnoTable", function(object) {
  msg <- character(0)

  if (!is.null(val <- object[["IDColumn"]])) {
    msg <- .validStringError(msg, object, "IDColumn")
  }

  msg <- .validStringError(msg, object, "Organism")

  msg <- .allowableChoiceError(msg, object, "IDType", allowable)

  msg <- .validLogicalError(msg, object, "AnnoBoxOpen")

  if (length(msg)) {
    return(msg)
  }
  TRUE
})

setMethod("initialize", "GeneAnnoTable", function(.Object, IDColumn="gene_names",
                                                  Organism="org.Hs.eg.db", IDType="SYMBOL", AnnoBoxOpen=FALSE, ...)
{
  callNextMethod(.Object, IDColumn=IDColumn, IDType=IDType,
                 Organism=Organism, AnnoBoxOpen=AnnoBoxOpen, ...)
})

GeneAnnoTable <- function(...) {
  new("GeneAnnoTable", ...)
}

setMethod(".fullName", "GeneAnnoTable", function(x) "Annotated gene table")

setMethod(".panelColor", "GeneAnnoTable", function(x) "#AA1122")

setMethod(".defineOutput", "GeneAnnoTable", function(x, ...) {
  panel_name <- .getEncodedName(x)
  tagList(
    callNextMethod(), # Re-using RowDataTable's definition.
    uiOutput(paste0(panel_name, "_annotation")),
    hr()
  )
})

setMethod(".defineInterface", "GeneAnnoTable", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  c(
    list(
      collapseBox(
        paste0(panel_name, "_AnnoBoxOpen"),
        title="Annotation parameters",
        open=x[["AnnoBoxOpen"]],
        selectInput(paste0(panel_name, "_IDColumn"),
                    label="ID-containing column:",
                    choices=colnames(rowData(se)),
                    selected=x[["IDColumn"]]
        ),
        selectInput(paste0(panel_name, "_IDType"),
                    label="ID type:",
                    choices=allowable,
                    selected=x[["IDType"]]
        ),
        selectInput(paste0(panel_name, "_Organism"),
                    label="Organism",
                    choices=c("org.Hs.eg.db", "org.Mm.eg.db"),
                    selected=x[["Organism"]]
        )
      )
    ),
    callNextMethod()
  )
})

setMethod(".createObservers", "GeneAnnoTable",
          function(x, se, input, session, pObjects, rObjects)
          {
            callNextMethod()

            plot_name <- .getEncodedName(x)

            .createUnprotectedParameterObservers(plot_name,
                                                 fields=c("IDColumn", "Organism", "IDType"),
                                                 input=input, pObjects=pObjects, rObjects=rObjects)
          })

setMethod(".renderOutput", "GeneAnnoTable", function(x, se, ..., output, pObjects, rObjects) {
  callNextMethod() # Re-using RowDataTable's output rendering.

  panel_name <- .getEncodedName(x)
  output[[paste0(panel_name, "_annotation")]] <- renderUI({
    .trackSingleSelection(panel_name, rObjects)
    instance <- pObjects$memory[[panel_name]]

    rowdata_col <- instance[["IDColumn"]]
    selectedGene <- instance[["Selected"]]
    if (!is.null(rowdata_col)) {
      selectedGene <- rowData(se)[selectedGene,rowdata_col]
    }

    keytype <- instance[["IDType"]]
    selgene_entrez <- NA
    if (keytype!="ENTREZID") {
      ORG <- instance[["Organism"]]
      if (require(ORG, character.only=TRUE, quietly=TRUE)) {
        orgdb <- get(ORG)
        selgene_entrez <- try(mapIds(orgdb, selectedGene, "ENTREZID", keytype),
                              silent=TRUE)
      }
    } else {
      selgene_entrez <- selectedGene
    }

    if (is.na(selgene_entrez) || is(selgene_entrez, "try-error")) {
      return(NULL)
    }

    fullinfo <- rentrez::entrez_summary("gene", selgene_entrez)
    link_pubmed <- paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',
                          selgene_entrez,
                          '" target="_blank">Click here to see more at the NCBI database</a>')

    mycontent <- paste0("<b>",fullinfo$name, "</b><br/><br/>",
                        fullinfo$description,"<br/><br/>",
                        ifelse(fullinfo$summary == "","",paste0(fullinfo$summary, "<br/><br/>")),
                        link_pubmed)

    HTML(mycontent)
  })
})

