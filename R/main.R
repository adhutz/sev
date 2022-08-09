
gat<-GeneAnnoTable(PanelWidth=8L)
got<-GOTable(PanelWidth=8L)
rst <- RowDataTable(PanelWidth = 12L)


#' isee_mini()
#' Explores your own summarized experiment
#'
#' @param se summarized experiment object
#'
#' @return an iSEE app
#' @importFrom iSEE iSEE
#' @export
isee_mini <- function(se) {
  iSEE::iSEE(se, initial=list(got, rst, VolcanoPlot(PanelWidth=6L), gat, FeatureSetTable(PanelWidth=6L)))
}


#' Make_mini()
#' Reduces the columns of rowData to the minimum. Names and IDs of genes and
#' proteins are retained together with t-test results including fold change,
#' p-value, adjusted p-value and if the result is significant (adj. p-value <0.05).
#'
#' @param se summarized Experiment
#'
#' @return se summarized Experiment with reduced rowData
#' @importFrom dplyr ends_with
#' @export
make_mini <- function(se) {
  rowData(se) <- rowData(se) %>%
    as.data.frame() %>%
    select(gene_names,
           protein_ids,
           majority_protein_ids,
           protein_names,
           dplyr::ends_with("_diff|_p.val|p.adj|significant"))

  return(se)
}

#' filter_perseus()
#' Filters proteins with too many missing values. Parameters mimic the options
#' given by Perseus.
#'
#' @param se summarized experiment
#' @param perc_na fraction of missing values that is tolerated
#' @param filter_mode "each_group": perc_na in every group; "one_group": perc_na
#' in at least one group; "total": perc_na across all samples.
#'
#' @return summarized experiment without entries with too many missing values
#' @importFrom dplyr group_by summarize filter ungroup select mutate n
#' @importFrom BiocGenerics unique
#' @export
filter_perseus<-function(se, perc_na = 0.33, filter_mode = "each_group"){

  if(filter_mode %in% c("one_group", "each_group", "total")){

      # 1. Step: create long table with groups

      data_long <- get_df_long(se)

      # 2. Step: Apply filtering

      if(filter_mode=="one_group"){

        filtered_terms<-data_long%>%
          group_by(name, condition)%>%
          summarize(number_of_na=sum(is.na(intensity.x)), group_size=n()) %>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          ungroup()%>%
          select(name)%>%
          BiocGenerics::unique()

      }

      if(filter_mode=="each_group"){
        filtered_terms<-data_long%>%
          group_by(name, condition)%>%
          summarize(number_of_na=sum(is.na(intensity.x)), group_size=n())%>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          mutate(n_data=length(unique(data_long$condition)))%>%
          group_by(name)%>%dplyr::filter(n()==n_data)%>%
          select(name)%>%unique()

      }

      if(filter_mode=="total"){

        filtered_terms<-data_long%>%
          group_by(sample, name)%>%
          summarize(m=mean(intensity.x))%>%ungroup()%>%
          group_by(name)%>%
          summarise(number_of_na=sum(is.na(m)), group_size=n())%>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          select(name)%>%unique()

      }

      # 3. Step: return summarized experiment with only remaining genes
      return(se[filtered_terms %>% unlist(),])
  }else{
    cat("Unknown filter_mode. Select one of: \"one_group\", \"each_group\", \"total\".")
  }
}


#'Impute_perseus()
#'Imputes missing values in a summarized experiment similar to Perseus. Missing values are
#'replaced by drawing from a left shifted gaussian distribution that is calculated based on either
#'individual samples or all measurements.
#' @param se summarized experiment
#' @param width width of gaussian distribution
#' @param downshift shift of gaussian distribution
#' @param per_col if distribution should be calculated per column (sample)
#' or accross all samples
#'
#' @return summarized experiment without missing values
#' @importFrom dplyr group_by summarize filter ungroup select mutate_all mutate
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @export
#'
impute_perseus = function(se, width = 0.3, downshift = 1.8, per_col=T) {

  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  
  # 1. transform to long and set lfq_imputed = TRUE for all missing values
  df <- se %>%
    get_df_long() %>%
    select(name, label, condition, intensity.x) %>%
    ungroup()%>%
    mutate(lfq_imputed = ifelse(is.na(.$intensity.x), TRUE, FALSE))

  # 2. impute with the selected method
  if(per_col){
    cols<-split(df, df$label)

    for (cols_ in names(cols)){
      set.seed(1)
      temp = cols[[cols_]]$intensity.x
      temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
      temp.mean = mean(temp, na.rm = TRUE) -
        downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
      n.missing = sum(is.na(temp))
      cols[[cols_]]$intensity.x[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    }

    df <- do.call(rbind, cols)

  }

  else{
    set.seed(1)
    temp = df$intensity.x
    temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
    temp.mean = mean(temp, na.rm = TRUE) -
      downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
    n.missing = sum(is.na(temp))
    df$lfq_intensity_norm[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  }

  # Calculate table with all values and another assay describing which values were imputed
  imp_val <- df  %>%
    select(label, intensity.x, name) %>%
    mutate(label = gsub("lfq_intensity_", "", label)) %>%
    pivot_wider(names_from = "label", values_from = "intensity.x") %>%
    as.data.frame() %>%
    select(-name)

  imp_yn <- df  %>%
    select(label, lfq_imputed, name) %>%
    pivot_wider(names_from = "label", values_from = "lfq_imputed") %>%
    as.data.frame() %>%
    select(-name)

  colnames(imp_yn) <- paste0(se$ID,"_imputed")

  # Add as assay
  se <- add_assay(se, imp_val %>% as.matrix(), name = "imputed_perseus", withDimnames = FALSE)
  assays(se, withDimnames = FALSE)$imputed <- imp_yn %>% dplyr::mutate_all(~ ifelse(.x, 1,0)) %>% as.matrix()

  # Add rowData
  rowData(se) <- cbind(rowData(se), imp_yn)
  
  return(se)

}

###################################

#' se_volcano()
#' Create a volcano plot for a summarized experiment and the selected contrast.
#'
#' @param se summarized experiment
#' @param contrast_ contrast
#'
#' @return volcano plot
#' @importFrom ggplot2 ggplot geom_point scale_color_manual theme_bw labs lims
#' @importFrom ggrepel geom_text_repel
#'

se_volcano<-function(se, contrast_){
  LFC <- rowData(se)[,paste0(contrast_, "_diff")]
  p <- rowData(se)[,paste0(contrast_, "_p.val")]%>%-log2(.)
  p.adj <- rowData(se)[,paste0(contrast_, "_p.adj")]%>%-log2(.)
  sig <- rowData(se)[,paste0(contrast_, "_significant")]

  change<-ifelse(sig & (LFC>0), "Increase",
                 ifelse(sig & (LFC < 0), "Decrease","None"))

  df <- data.frame(LFC,p,sig,change, "SYMBOL"=rowData(se)$gene_names)

  #Calculate the axis limits
  ex_x<-max(LFC)%>%ceiling(.)
  ex_y<-max(p)%>%ceiling(.)
  plot<-ggplot(data=df, aes(labels=SYMBOL,x=LFC, y=p))+
    geom_point(data=df, pos="identity", aes(color=change), size=2, alpha=0.6)+
    scale_color_manual(values=c("Increase"="#ED553B","None"="grey","Decrease"="#20639B"))+
    geom_text_repel(data=subset(df, sig==TRUE), aes(x=LFC, y=p, label=SYMBOL), max.overlaps = Inf)+
    #geom_vline(xintercept = FC_cut, linetype="dashed")+
    #geom_vline(xintercept = -FC_cut, linetype="dashed")+
    theme_bw()+
    labs(x="Log2 fold-change", y="-log2 (p-value)", color="Change", title=paste("Volcano plot -", contrast_))+
    lims(x=c(-ex_x, ex_x), y=c(0, ex_y))

  return(plot)

}

#' Split_genes
#'
#' @param table dataframe
#' @param colname column that needs to be split (e.g., gene_names, protein_ids)
#' @param keep_all if TRUE, all entries are retained in additional rows. If FALSE,
#' only the first name/id is kept.
#' @importFrom splitstackshape cSplit
#' @return table with split entries
split_genes <- function(table, colname="gene_names", keep_all=FALSE){
  require(splitstackshape)
  if (keep_all){
    return(cSplit(table, colname, direction = "long", sep=";"))
  }
  else{
    table[[colname]]<-gsub(";.*", "", table[[colname]])
    return(table)
  }
}

#' se_read_in()
#' Reads in MaxQuant or other output and creates a summarized experiment.
#'
#' @param file path to proteinGroups or similar file
#' @param gene_column name of gene_name column after janitor
#' @param protein_column name of protein column after janitor
#' @param sep character describing the separator between sample name and replicate number (e.g. "_rep_", "_r")
#' @param filt character vector of column names used for filtering. Genes with a "+" will be dropped.
#' @param keep_all_proteins if FALSE, only first protein per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param keep_all_genes if FALSE, only first gene per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param experimental_design dataframe with information regarding samples. If not specified,
#' sample names and groups are read automatically (start with "lfq_intensity"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @export
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names
#' @examples
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#'
se_read_in <- function(file, gene_column = "gene_names", protein_column = "protein_ids", sep="_rep_",
                       filt = c("reverse", "only_identified_by_site", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){

  #ProteinGroups
  data<-read.delim(file, sep="\t")

  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()

  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_ids,
           orig_gene_names = gene_names) %>%
    split_genes(colname = "protein_ids", keep_all = keep_all_proteins) %>%
    split_genes(colname = "gene_names", keep_all = keep_all_genes)

  #Filter false and low quality hits
  data <- data %>% filter(!if_any(filt, ~ .x == "+"))

  #Make gene_names unique
  data_unique <- make_unique(data, "gene_names", "protein_ids", delim=";")

  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){

    LFQ_labels <- colnames(data)[grep("^lfq_intensity_", colnames(data_unique))]

    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("lfq_intensity_", "", LFQ_labels),
                                    condition=gsub(paste0("lfq_intensity_|",sep, "[0-9].*"), "", LFQ_labels),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE))
  }else{
    LFQ_labels<-experimental_design$label
  }

  data_se <- make_se(data_unique, which(colnames(data_unique) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}

#' add_randna()
#' Add column to rowData that specifies for each row if values are missing at random or
#' not. The method is very rudimentary:
#' 1. Missing values are replaced by zeros, measured values by 1.
#' 2. ANOVA between conditions is performed to determine significant differences
#' of the mean.
#' 3. Rows with p-values < 0.05 are termed missing not at random.
#'
#' @param se
#'
#' @return se with added rowData column
#' @importFrom HybridMTest row.oneway.anova
#' @importFrom fdrtool gcmlcm
#' @export
add_randna <- function(se){

  # replace na with 0, values with 1
  dummy <- as.matrix(ifelse(is.na(assay(se)), 0,1))

  #ANOVA
  oa <- HybridMTest::row.oneway.anova(dummy, colData(se)$condition)

  #Add information regarding mode of missingness (MAR = TRUE)
  rowData(se)$randna<-ifelse(is.na(oa$pval), FALSE, ifelse(oa$pval<0.05, FALSE, TRUE))

  return(se)
}

#' se_GOE()
#' Performs GO-term enrichment via clusterprofiler for all contrasts of a summarized experiment.
#' Results are added to the metadata.
#'
#' @param se with t-test results (e.g. via DEP test_diff())
#' @param col_names name of columns on that enrichment is performed. If left empty, all columns ending with "_diff" are selected.
#' @param simplify Should similar terms be merged (sim > 0.7)
#' @param ont ontology
#' @param keyType Type of identifier
#' @param OrgDb Organism
#' @param pvalueCutoff significance threshold
#'
#' @return se with GO-term enrichment in metadata
#' @export
#'
#' @importFrom clusterProfiler simplify gseGO
se_GOE <- function(se, col_names =c(), simplify = TRUE, ont="all", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pvalueCutoff = 1){

  tab_ <- rowData(se)%>%as.data.frame()

  if(length(col_names) == 0){
    col_names <- colnames(tab_)[grep("_diff",colnames(tab_))]
  }

  m<-list()
  GO_raw <- list()

  for(comp_ in col_names){

    l <- tab_[,comp_]
    names(l)<-tab_$gene_names
    l<-sort(l, decreasing=TRUE)

    #Calculate GO-term enrichment
    temp <- clusterProfiler::gseGO(l, ont="all", keyType = "SYMBOL", OrgDb = OrgDb, pvalueCutoff = 1)

    #Reduce similar terms

    if(simplify == TRUE){
      tryCatch(
        expr = {

          temp <- temp %>% clusterProfiler::simplify(., 0.7, by="NES", select_fun=max)

        },

        error = function(e){
          print("Not enough terms for reduction ")
        }

      )
    }
    #Split genes to list that can be compared to rownames of se
    gene_sets <- temp$core_enrichment %>% strsplit(.,split="/")

    #Assign GO-term names
    names(gene_sets)<-temp$Description

    #Get unique gene names to allow filtering of se by GO-terms
    for(names_ in names(gene_sets)){
      gene_sets[[names_]]<- tab_ %>% filter(gene_names %in% gene_sets[[names_]])%>%rownames()
    }

    gene_sets <- CharacterList(gene_sets)

    mcols(gene_sets) <- temp

    metadata(se)$GO_enrichment[[comp_]] <- gene_sets

  }

  return(se)
}


#' se_to_isee()
#' Transform summarized experiment. The result is a summarized experiment
#' of the type "SingleCellExperiment" to allow dimensional reduction plots.
#' GO-enrichment is registered from the metadata as a feature set collection and
#' p-values and differences are registered to be plotted via volcano plots.
#' #'
#' @param se summarized experiment
#'
#' @return se summarized experiment with registered features
#' @export
#' @importFrom iSEEu registerFeatureSetCollections registerPValuePatterns registerLogFCPatterns
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom stats prcomp
se_to_isee <- function(se, PValuePatterns = "p.val", LogFCPatterns = "_diff"){
  se <- as(se, "SingleCellExperiment")

  if("GO_enrichment" %in% names(metadata(se))){
  # Add GO-term enrichment as FeatureSet to allow plotting
  se <- registerFeatureSetCollections(se, metadata(se)$GO_enrichment)
  }
  
  if(!any(is.na(assay(se)))){
  # Calculate pca data
  pca_data <- prcomp(t(assay(se)), rank=50)

  # Add PCA data to the experiment object
  reducedDims(se)<-list(PCA=pca_data$x)
  }
  
  #Register columns that contain p-values and log differences (for volcano mainly)
  se<-registerPValuePatterns(se, PValuePatterns)
  se<-registerLogFCPatterns(se, LogFCPatterns)

}


#' impute_DEP()
#' Utilizes the DEP::impute() function to replace missing values. In addition, original raw data
#' is retained in an additional assay.
#'
#' @param se with missing data
#' @param ... parameters for DEP::impute()
#'
#' @return se with imputed data in the main assay and raw data in another assay
#' @export
#'
impute_DEP <- function(se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
                               "man", "min", "zero", "mixed", "nbavg"), ...){
  
  if(fun == "mixed"){
    
    temp <- DEP::impute(se, fun = fun, randna = rowData(se)$randna, ...)
    
  }else{
    
    temp <- DEP::impute(se, fun = fun, ...)
    
  }
  
  se <- add_assay(se, assay(temp), "imputed_DEP")
  
  assays(se, withDimnames = FALSE)$imputed <- assay(se) %>% as.data.frame() %>% dplyr::mutate_all(~ ifelse(is.na(.x), 1,0)) %>% as.matrix()
  
  return(se)
  
}

#' add_assay()
#' 
#' Convinience function to add a new assay directly as the main assay. The old assay is reatined. 
#' @param se Summarized experiment
#' @param assay Assay to add as the main assay (first in the assays() list) 
#' @param name Name of the added assay
#'
#' @return se with added assay. Old main assay is retained. 
#' @export
#'
add_assay <- function(se, new_assay, name, withDimnames = TRUE){
  temp <- assay(se)
  temp_n <- names(assays(se))[1]
  
  assay(se, withDimnames = withDimnames) <- new_assay

  names(assays(se))[1] <- name
  assays(se)[[temp_n]] <- temp
  
  return(se)  
}


#' to_pdf()
#' 
#' Convenience function to print plots as pdf document. 
#'
#' @param .data plot object
#' @param filename String specifying filename. Can also specify path. 
#' @param w width in inches
#' @param h height in inches
#'
#' @return plot
#' @export
#'
to_pdf<-function(.data, filename, w=7,h=7){
  
  pdf(  sub("(.*)(?<=/)(.*)", paste0("\\1",gsub(":","-",format(Sys.time(), "%Y%m%d_")), "\\2", ".pdf"), filename, perl=TRUE ), # File name
        width = w, height = h, # Width and height in inches
        bg = "white",          # Background color
        colormodel = "srgb",    # Color model (cmyk is required for most publications)
        onefile=T)          # Paper size
  
  # Creating a plot
  plot(.data)
  # Closing the graphical device
  dev.off() 
}
