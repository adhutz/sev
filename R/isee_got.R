
setClass("GOTable", contains="Panel",
         slots=c(
           IDType="character",
           Organism="character",
           Selected="character",
           Search="character",
           SearchColumns="character"
         )
)

allowable.ids <- c("ENSEMBL", "SYMBOL", "ENTREZID")
allowable.org <- c("org.Mm.eg.db", "org.Hs.eg.db")

setValidity2("GOTable", function(object) {
  msg <- character(0)

  msg <- .allowableChoiceError(msg, object, "Organism", allowable.org)

  msg <- .allowableChoiceError(msg, object, "IDType", allowable.ids)

  msg <- .singleStringError(msg, object, c("Selected", "Search"))

  if (length(msg)) {
    return(msg)
  }
  TRUE
})

setMethod("initialize", "GOTable", function(.Object,
                                            Organism="org.Mm.eg.db", IDType="SYMBOL",
                                            Selected="", Search="", SearchColumns=character(0), ...)
{
  callNextMethod(.Object, IDType=IDType, Organism=Organism,
                 Selected=Selected, Search=Search,
                 SearchColumns=SearchColumns, ...)
})

GOTable <- function(...) new("GOTable", ...)

setMethod(".fullName", "GOTable", function(x) "Gene ontology table")

setMethod(".panelColor", "GOTable", function(x) "#BB00FF")

setMethod(".defineOutput", "GOTable", function(x, ...) {
  panel_name <- .getEncodedName(x)
  tagList(DT::dataTableOutput(panel_name))
})

setMethod(".defineDataInterface", "GOTable", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  list(
    selectInput(paste0(panel_name, "_IDType"),
                label="ID type:",
                choices=allowable.ids,
                selected=x[["IDType"]]
    ),
    selectInput(paste0(panel_name, "_Organism"),
                label="Organism",
                choices=allowable.org,
                selected=x[["Organism"]]
    )
  )
})

setMethod(".hideInterface", "GOTable", function(x, field) {
  if (field %in% "SelectionBoxOpen") {
    TRUE
  } else {
    callNextMethod()
  }
})

setMethod(".generateOutput", "GOTable", function(x, se, ..., all_memory, all_contents) {
  envir <- new.env()
  commands <- c("require(GO.db);",
                "tab <- select(GO.db, keys=keys(GO.db), columns='TERM');",
                "rownames(tab) <- tab$GOID;",
                "tab$GOID <- NULL;")
  eval(parse(text=commands), envir=envir)
  list(
    commands=list(commands),
    contents=list(table=envir$tab, available=nrow(se)),
    varname="tab"
  )
})

setMethod(".createObservers", "GOTable",
          function(x, se, input, session, pObjects, rObjects)
          {
            callNextMethod()

            panel_name <- .getEncodedName(x)

            .createUnprotectedParameterObservers(panel_name,
                                                 fields=c("Organism", "IDType"),
                                                 input=input, pObjects=pObjects, rObjects=rObjects)

            # Observer for the DataTable row selection:
            select_field <- paste0(panel_name, "_rows_selected")
            multi_name <- paste0(panel_name, "_", iSEE:::.flagMultiSelect)
            observeEvent(input[[select_field]], {
              chosen <- input[[select_field]]
              if (length(chosen)==0L) {
                chosen <- ""
              } else {
                chosen <- rownames(pObjects$contents[[panel_name]]$table)[chosen]
              }

              previous <- pObjects$memory[[panel_name]][["Selected"]]
              if (chosen==previous) {
                return(NULL)
              }
              pObjects$memory[[panel_name]][["Selected"]] <- chosen
              .requestActiveSelectionUpdate(panel_name, session, pObjects, rObjects, update_output=FALSE)
            }, ignoreInit=TRUE, ignoreNULL=FALSE)

            # Observer for the search field:
            search_field <- paste0(panel_name, "_search")
            observeEvent(input[[search_field]], {
              search <- input[[search_field]]
              if (identical(search, pObjects$memory[[panel_name]][["Search"]])) {
                return(NULL)
              }
              pObjects$memory[[panel_name]][["Search"]] <- search
            })

            # Observer for the column search fields:
            colsearch_field <- paste0(panel_name, "_search_columns")
            observeEvent(input[[colsearch_field]], {
              search <- input[[colsearch_field]]
              if (identical(search, pObjects$memory[[panel_name]][["SearchColumns"]])) {
                return(NULL)
              }
              pObjects$memory[[panel_name]][["SearchColumns"]] <- search
            })
          })

setMethod(".renderOutput", "GOTable", function(x, se, ..., output, pObjects, rObjects) {
  callNextMethod()

  panel_name <- .getEncodedName(x)
  output[[panel_name]] <- DT::renderDataTable({
    t.out <- .retrieveOutput(panel_name, se, pObjects, rObjects)
    full_tab <- t.out$contents$table

    param_choices <- pObjects$memory[[panel_name]]
    chosen <- param_choices[["Selected"]]
    search <- param_choices[["Search"]]
    search_col <- param_choices[["SearchColumns"]]
    search_col <- lapply(search_col, FUN=function(x) { list(search=x) })

    # If the existing row in memory doesn't exist in the current table, we
    # don't initialize it with any selection.
    idx <- which(rownames(full_tab)==chosen)[1]
    if (!is.na(idx)) {
      selection <- list(mode="single", selected=idx)
    } else {
      selection <- "single"
    }

    DT::datatable(
      full_tab, filter="top", rownames=TRUE,
      options=list(
        search=list(search=search, smart=FALSE, regex=TRUE, caseInsensitive=FALSE),
        searchCols=c(list(NULL), search_col), # row names are the first column!
        scrollX=TRUE),
      selection=selection
    )
  })
})

setMethod(".multiSelectionDimension", "GOTable", function(x) "row")

setMethod(".multiSelectionCommands", "GOTable", function(x, index) {
  orgdb <- x[["Organism"]]
  type <- x[["IDType"]]
  c(
    sprintf("require(%s);", orgdb),
    sprintf("selected <- tryCatch(select(%s, keys=%s, keytype='GO',
    column=%s)$SYMBOL, error=function(e) character(0));",
            orgdb, deparse(x[["Selected"]]), deparse(type)),
    "selected <- intersect(selected, rownames(se));"
  )
})

setMethod(".multiSelectionActive", "GOTable", function(x) {
  if (x[["Selected"]]!="") {
    x[["Selected"]]
  } else {
    NULL
  }
})

setMethod(".multiSelectionClear", "GOTable", function(x) {
  x[["Selected"]] <- ""
  x
})

setMethod(".multiSelectionAvailable", "GOTable", function(x, contents) {
  contents$available
})
