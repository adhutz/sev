
gat<-GeneAnnoTable(PanelWidth=8L)
got<-GOTable(PanelWidth=8L)
rst <- RowDataTable(PanelWidth = 12L)


#' Create iSEE app with all custom panels
#'
#' Includes GeneAnnoTable, GOTable and RowDataTable. 
#'
#' @param se summarized experiment object
#'
#' @return an iSEE app
#' @importFrom iSEE iSEE
#' @export
isee_mini <- function(se) {
  iSEE::iSEE(se, initial=list(got, rst, VolcanoPlot(PanelWidth=6L), gat, FeatureSetTable(PanelWidth=6L)))
}


#' Minimize results
#' 
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

#' Filter on occurence 
#' 
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
          summarize(number_of_na=sum(is.na(intensity)), group_size=n()) %>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          ungroup()%>%
          select(name)%>%
          BiocGenerics::unique()

      }

      if(filter_mode=="each_group"){
        filtered_terms<-data_long%>%
          group_by(name, condition)%>%
          summarize(number_of_na=sum(is.na(intensity)), group_size=n())%>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          mutate(n_data=length(unique(data_long$condition)))%>%
          group_by(name)%>%dplyr::filter(n()==n_data)%>%
          select(name)%>%unique()

      }

      if(filter_mode=="total"){

        filtered_terms<-data_long%>%
          group_by(sample, name)%>%
          summarize(m=mean(intensity))%>%ungroup()%>%
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


#'Impute as with perseus
#'
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


  # 1. transform to long and set lfq_imputed = TRUE for all missing values
  df <- se %>%
    get_df_long() %>%
    select(sample, condition, intensity, name) %>%
    ungroup()%>%
    mutate(lfq_imputed = ifelse(is.na(intensity), TRUE, FALSE))

  # 2. impute with the selected method
  if(per_col){
    cols<-split(df, df$sample)

    for (cols_ in names(cols)){
      set.seed(1)
      temp = cols[[cols_]]$intensity
      temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
      temp.mean = mean(temp, na.rm = TRUE) -
        downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
      n.missing = sum(is.na(temp))
      cols[[cols_]]$intensity[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    }

    df <- do.call(rbind, cols)

  }

  else{
    set.seed(1)
    temp = df$intensity
    temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
    temp.mean = mean(temp, na.rm = TRUE) -
      downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
    n.missing = sum(is.na(temp))
    df$lfq_intensity_norm[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  }

  # Calculate table with all values and another assay describing which values were imputed
  imp_val <- df  %>%
    select(sample, intensity, name) %>%
    pivot_wider(names_from = "sample", values_from = "intensity") %>%
    as.data.frame() %>%
    select(-name)

  imp_yn <- df  %>%
    select(sample, lfq_imputed, name) %>%
    pivot_wider(names_from = "sample", values_from = "lfq_imputed") %>%
    as.data.frame() %>%
    select(-name)

  colnames(imp_yn) <- paste0(se$sample,"_imputed")

  # Add as assay
  se <- add_assay(se, imp_val %>% as.matrix(), name = "imputed_perseus", withDimnames = FALSE)
  assays(se, withDimnames = FALSE)$imputed <- imp_yn %>% dplyr::mutate_all(~ ifelse(.x, 1,0)) %>% as.matrix()

  # Add rowData
  rowData(se) <- cbind(rowData(se), imp_yn)
  
  return(se)

}

###################################

#' Colored volcano plot
#' 
#' Create a volcano plot for a summarized experiment and the selected contrast.
#'
#' @param se summarized experiment
#' @param contrast_ contrast
#' @param id_col id column that contains targets
#' @param target_names targets to label
#' @param label_sign if TRUE, significant genes are labeled
#' @param label_targets if TRUE, target genes are labeled
#' @param max.overlaps maximal overlaps for labels. Set to Inf to label everything
#' 
#' @return volcano plot
#' @export
#' 
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @import ggplot2

se_volcano <- function(se, contrast, id_col = "gene_names", target_names = c(""), label_sign = TRUE, label_targets = TRUE, max.overlaps = Inf){
  pval <- paste0(contrast, "_p.val")
  padj <- paste0(contrast, "_p.adj")
  diff <- paste0(contrast, "_diff")
  sign <- paste0(contrast, "_significant")
  
  plot_data <- rowData(se) %>% as.data.frame() %>% 
    select(all_of(c(id_col, "name")), starts_with(contrast)) %>% 
    mutate(target_gene = factor(ifelse(!!sym(id_col) %in% target_names, "Target", "Non-Target"), levels = c( "Non-Target", "Target"))) %>%
    mutate(sign = !!sym(sign)) %>%
    mutate(target_sign = factor(ifelse(!!sym(id_col) %in% target_names, "Target", 
                                       ifelse(!!sym(sign) == TRUE, "Significant", "None")), levels = c("None", "Significant", "Target"))) %>%
    mutate(alpha = ifelse(target_sign == "None", "None", "goi")) %>%
    plyr::arrange(target_gene, sign)
  
  p1 <- ggplot(plot_data, aes_string(x = diff, y = paste0("-log10(", pval, ")"))) +
    
    # Plot points based on another cutoff
    geom_point(aes(fill = target_sign, alpha = alpha, shape = sign), size = 2, stroke = 0.1, show.legend = T) +
    scale_fill_manual(values = c("Significant" = "firebrick", "Target" = "darkorange", "None" = "grey90"), drop = FALSE) +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), drop = FALSE) +
    scale_alpha_manual(values = c("None" = 0.2, "goi" = 0.7), drop = FALSE) +
    guides(alpha = FALSE,
           fill = guide_legend( 
             override.aes=list(shape = 21, size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    
    labs(title = contrast,
         fill = "Targets",
         shape = "Significant",
         x = "Log2(Fold Change)",
         y = "-Log10(P-value)") +
    
    theme_bw()
  
  if(label_sign){
    if(label_targets){
      p1 <- p1 + geom_text_repel(data = subset(plot_data, target_gene != "Target" & sign), aes(label = name), size = 3, color = "firebrick", min.segment.length = 0, max.overlaps = max.overlaps)
    } else{
      p1 <- p1 + geom_text_repel(data = subset(plot_data, sign), aes(label = name), size = 3, color = "firebrick", min.segment.length = 0, max.overlaps = max.overlaps)
    }
  }
  
  if(label_targets){
    p1 <- p1 + geom_label_repel(data = subset(plot_data, target_gene == "Target"), aes(label = name), size = 3, color = "darkorange", min.segment.length = 0, max.overlaps = max.overlaps)
    
  }
  
  return(p1)
}

#' Handling of columns with multiple entries
#'
#' Split columns with multiple entries (need to be separated by a semicolon) into multiple rows with one entry each.
#' @param table dataframe
#' @param colname column that needs to be split (e.g., gene_names, protein_ids)
#' @param keep_all if TRUE, all entries are retained in additional rows. If FALSE,
#' only the first name/id is kept.
#' @param sep separator
#' @importFrom splitstackshape cSplit
#' @return table with split entries
split_genes <- function (table, colname = "gene_names", keep_all = FALSE, sep = ";") 
{
  require(splitstackshape)
  if (keep_all) {
    return(cSplit(table, colname, direction = "long", sep = sep))
  }
  else {
    table[[colname]] <- gsub(";.*", "", table[[colname]])
    return(table)
  }
}

#' Create se from MaxQuant
#' 
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
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all if_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names
#' @examples
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#'
se_read_in <- function(file, gene_column = "gene_names", protein_column = "protein_ids", sep="_rep_",
                       filt = c("reverse", "only_identified_by_site", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){

  #ProteinGroups
  data <- read.delim(file, sep="\t")

  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()

  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_column,
           orig_gene_names = gene_column) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes) %>%
    dplyr::rename("perseus_intensity" = "intensity")

  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ .x == ""))
  
  #Make gene_names unique
  data_unique <- make_unique(data, gene_column, protein_column, delim=";")

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


#' Explore type of missingness
#' 
#' Add column to rowData that specifies for each row if values are missing at random or
#' not (MAR vs MNAR). The method is very rudimentary:
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


#' GO term enrichment
#' 
#' Performs GO-term enrichment via clusterprofiler for all contrasts of a summarized experiment.
#' Results are added to the metadata.
#'
#' @param se with t-test results (e.g. via DEP2 test_diff())
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


#' Make se compatible to iSEE 
#' 
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


#' Imputation via DEP2 while keeping original data
#' 
#' Utilizes the DEP2::impute() function to replace missing values. In addition, original raw data
#' is retained in an additional assay.
#'
#' @param se with missing data
#' @param ... parameters for DEP2::impute()
#'
#' @return se with imputed data in the main assay and raw data in another assay
#' @export
#'
impute_DEP2 <- function(se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
                               "man", "min", "zero", "mixed", "nbavg"), ...){
  
  assays(se, withDimnames = FALSE)$imputed <- assay(se) %>% as.data.frame() %>% dplyr::mutate_all(~ ifelse(is.na(.x), 1,0)) %>% as.matrix()
  
  if(fun == "mixed"){
    
    temp <- DEP2::impute(se, fun = fun, randna = rowData(se)$randna, ...)
    
  }else{
    
    temp <- DEP2::impute(se, fun = fun, ...)
    
  }
  
  se <- add_assay(se, assay(temp), "imputed_DEP2")
  
  
  return(se)
  
}

#' Add assay to se
#' 
#' Convinience function to add a new assay directly as the main assay. The old assay is retained. 
#' @param se Summarized experiment
#' @param assay Assay to add as the main assay (first in the assays() list) 
#' @param name Name of the added assay
#'
#' @return se with added assay. Old main assay is retained. 
#' @export
#'
add_assay <- function(se, new_assay, name, withDimnames = TRUE){
  
    if(name %in% names(assays(se))){
      names(assays(se))[names(assays(se)) == name] <- paste0(name, "_old")
    }
    temp <- assay(se)
    temp_n <- names(assays(se))[1]
    assay(se, withDimnames = withDimnames) <- new_assay
    assays(se, withDimnames = withDimnames)[[1]] <- new_assay
    names(assays(se))[1] <- name
    
    # add old assay
    assays(se)[[temp_n]] <- temp
    
    return(se)
  }


#' Plot to pdf
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



#' Handling of pig derived data
#' 
#' Pig proteins are poorly annotated and have a different "style" when it comes to fasta files. Thus, the output of MaxQuant 
#' is not directly usable. This function takes in a MaxQuant file from an experiment using pig derived proteins and a fasta
#' file and returns MaxQuant like files with correct protein names and gene names. In addition, columns are added that specify
#' the organism in case proteins from multiple organisms are expected. In this case, the fasta input must contain entries for 
#' the other organisms as well in the "pig-style". 
#'
#' @param proteingroups String specifying path to ProteinGroups file
#' @param peptides String specifying path to ProteinGroups file
#' @param fasta String specifying path to ProteinGroups file
#' @param mult_org Logical specifying if fasta and experiment contain multiple organisms.
#' @param obj Logical specifying if results should be returned as one named list. If FALSE, results are saved as individual files. 
#' @importFrom dplyr mutate select
#' @return Either a named list containing all dataframes with correct protein names and gene names or locations where results have been stored.
#' @export
fix_maxq_pig <- function(proteingroups, peptides, fasta, mult_org = FALSE, obj = FALSE){
  
  # Read fasta file and create additional file containing only headers
  system2(command = "grep", args = c("\"^>\"", fasta), stdout = paste0(gsub("\\..*", "", fasta), "_fasta_headers.txt"))
  
  # Read newly created header file
  fasta_headers_ <- read.delim(paste0(gsub("\\..*", "", fasta), "_fasta_headers.txt"), header = FALSE, col.names = "fasta")
  
  # Create annotation from fasta headers
  fasta_headers <- fasta_headers_ %>% 
    mutate(
      uniprot = sub("^>.*?\\|(.*?)\\|.*", "\\1", fasta, perl = TRUE), # extract UNIPROT ID (Between "|")
      
      uniprot_name = ifelse(grepl("\\|.*\\|[A-Z0-9]*_[A-Z0-9]* ", fasta), # Extract uniprot name if present
                            sub("^>.*?\\|.*\\|(.*?) .*", "\\1", fasta, perl = TRUE), 
                            NA),
      protein_name = ifelse(grepl("\\|.*\\|[A-Z0-9]*_[A-Z0-9]*", fasta), # Extract protein name
                            sub("^>.*?\\|.*?\\|.*? (.*?) OS.*", "\\1", fasta, perl = TRUE), # If uniprot_name is present 
                            sub("^>.*?\\|.*?\\|(.*?)\\ OS.*", "\\1", fasta, perl = TRUE)), #Else
      organism = ifelse(grepl("OS=", fasta), 
                        sub(".*OS=(.*? [a-z]*).*", "\\1", fasta, perl = TRUE), #Extract Organism ID if present. If no OS is specified, protein_names won't be returned
                        NA),
      organism_id = ifelse(grepl("OX=", fasta), 
                           sub(".* OX=([0-9]*).*", "\\1", fasta, perl = TRUE), #Extract Organism ID
                           NA),
      gene_name = ifelse(grepl("GN=.*? ", fasta), 
                         sub(".*GN=(.*?) .*", "\\1", fasta, perl = TRUE), 
                         ifelse(grepl("GN=.*?", fasta), 
                                sub(".*GN=.*?", "\\1", fasta, perl = TRUE), 
                                NA)
      )
    )
  
  # For peptides
  peptides_first<-read.delim(file = peptides) %>%
    sev:::split_genes(., "Proteins", FALSE) %>%
    merge(fasta_headers, by.x = "Proteins", by.y = "uniprot", all.x = T) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "organism_id") 
      else c("gene_name", "protein_name", "uniprot_name")
    )
  
  peptides_all<-read.delim(file = peptides) %>%
    sev:::split_genes(., "Proteins", TRUE) %>%
    merge(fasta_headers, by.x="Proteins", by.y = "uniprot", all.x = T) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "organism_id") 
      else c("gene_name", "protein_name", "uniprot_name")
    )
  
  # For protein groups
  protein_groups <- read.delim(file = proteingroups)
  
  protein_groups_first <- sev:::split_genes(protein_groups, "Protein.IDs", FALSE) %>%
    merge(fasta_headers, by.x = "Protein.IDs", by.y = "uniprot", all.x = TRUE) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name,
           Fasta.headers = fasta) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "fasta", "organism_id") 
      else c("gene_name", "protein_name", "fasta", "uniprot_name")
    )
  
  protein_groups_all <- sev:::split_genes(protein_groups, "Protein.IDs", TRUE) %>%
    merge(fasta_headers, by.x = "Protein.IDs", by.y = "uniprot", all.x = TRUE) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name,
           Fasta.headers = fasta) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "fasta", "organism_id") 
      else c("gene_name", "protein_name", "fasta", "uniprot_name")
    )
  
  if(obj){
    return(list("fasta_headers" = fasta_headers, 
                "peptides_first" = peptides_first, 
                "peptides_all" = peptides_all, 
                "proteins_groups_first" = protein_groups_first,
                "protein_groups_all" = protein_groups_all, 
                "fasta_file" = fasta)
    )
  }
  
  else{
    write.table(fasta_headers, format(Sys.time(), "%Y%m%d_fasta_headers.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(peptides_first, format(Sys.time(), "%Y%m%d_peptides_first.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(peptides_all, format(Sys.time(), "%Y%m%d_peptides_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(protein_groups_first, format(Sys.time(), "%Y%m%d_protein_groups_first.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(protein_groups_all, format(Sys.time(), "%Y%m%d_protein_groups_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    return(
      paste0("Results are accessible in the folder:\n\"", getwd(), "\".\n\nFile names:\n", 
             format(Sys.time(), "%Y%m%d_fasta_headers.txt\n"),
             format(Sys.time(), "%Y%m%d_peptides_first.txt\n"),
             format(Sys.time(), "%Y%m%d_peptides_all.txt\n"),
             format(Sys.time(), "%Y%m%d_protein_groups_first.txt\n"),
             format(Sys.time(), "%Y%m%d_protein_groups_all.txt\n"),
             "\n\nInput fasta file was:\n", 
             fasta
      ) %>% 
        cat(.)
    ) 
  }
}


#' Strind-db network for goi
#' 
#' Creates a network with known protein-protein interactions and maps expression values to nodes.
#'
#' @param genes List or vector of gene symbols
#' @param species ID of species (human = 9606)
#' @param expression_data dataframe containing a column "gene_names" and additional columns containing numeric values
#' @param expand if TRUE, the network is expanded with additional known interactors
#' @param network_type either "functional" or "physical". If "physical", only genes with a physical interaction are connected (e.g. complexes)
#' @param score minimal score for a drawn interaction.
#' @param common_legend If TRUE, minimal and maximal values for the legend are calculated across all columns
#' @param node_deg_above Number specifiying which genes should be deleted from the network (e.g. if 0, all genes with zero interaction partners are removed)
#'
#' @importFrom dplyr mutate select distinct rename 
#' @importFrom rbioapi rba_string_map_ids rba_string_interaction_partners rba_string_interactions_network
#' @importFrom igraph graph_from_data_frame delete.vertices degree
#' @importFrom cowplot plot_grid
#' @import ggraph
#'
#' @return string network with integrated numeric values (fold-changes, lfq-values, ...)
#' @export
#'
get_network <- function(genes, species = 9606, expression_data = NA, expand = FALSE, add_nodes = NA, network_type = "functional", score = 900, common_legend = FALSE, node_deg_above = NA){
  
  # calculate min/max
  if(common_legend){
    min = min(expression_data[, -!is.numeric(expression_data)], na.rm = TRUE)
    max = max(expression_data[, -!is.numeric(expression_data)], na.rm = TRUE)
  } 
  
  # map gene names to stringIDs 
  prots <- genes %>% rbioapi::rba_string_map_ids(species=species)
  
  # If expand = TRUE, all interactions between input proteins and every other STRING protein are returned.
  # If expand = FALSE, interactions among the input set are retrieved. Additional proteins can be added via the add_nodes argument
  if(expand){
    int_net <- rbioapi::rba_string_interaction_partners(prots$stringId, 
                                                        species = 9606, 
                                                        required_score = score,
                                                        network_type = network_type)
  }else{
    int_net <- rbioapi::rba_string_interactions_network(prots$stringId, 
                                                        species = 9606, 
                                                        required_score = score, 
                                                        network_type = network_type, 
                                                        add_nodes = ifelse(is.na(add_nodes),
                                                                           0,
                                                                           add_nodes))
  } 
  
  if(ncol(int_net > 0)){
  
  # get stringIDs for all proteins of the expression table
  all_prots <- expression_data$gene_names %>% unlist() %>%
    rbioapi::rba_string_map_ids(species=species) %>% 
    merge(expression_data, 
          by.x = "queryItem", 
          by.y = "gene_names")
  
  # Create node table
  node_tbl <- int_net[,c("stringId_A", "preferredName_A")] %>% 
    dplyr::rename("stringId_B" = "stringId_A", "preferredName_B" = "preferredName_A") %>%
    rbind(int_net[,c("stringId_B", "preferredName_B")]) %>% 
    dplyr::distinct() %>% 
    dplyr::rename("clean_name" = "preferredName_B", "name" = "stringId_B") %>%
    merge(all_prots, by.x = "name", by.y = "stringId", all.x = TRUE) %>%
    dplyr::rename("gene_names" = "queryItem") %>%
    dplyr::select(name, clean_name, colnames(expression_data))
  
  # Create graph object
  g_obj <- igraph::graph_from_data_frame(int_net, node_tbl, directed = TRUE)
  
  if(!is.na(node_deg_above)){
    g_obj <- igraph::delete.vertices(g_obj , which(igraph::degree(g_obj)<=node_deg_above))
  }
  
  #Plot networks with expression data
  ps <- list()
  
  for(i in  colnames(expression_data[-grep("gene_names", colnames(expression_data))])){
    
    if(!common_legend){
      min = min(expression_data[,i], na.rm = TRUE)
      max = max(expression_data[,i], na.rm = TRUE)
    }
    
    ps[[i]] <- ggraph(g_obj, layout = "stress")+
      geom_edge_link0(aes(edge_width = score), 
                      edge_colour="grey",
                      alpha=0.7) +
      geom_node_point(aes_string(fill = i), 
                      shape = 21, size=12) +
      geom_node_text(aes(label = clean_name))+
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           limits = c(min,max))+
      scale_edge_width_continuous(range = c(0,1))+
      theme_void()+
      labs(title = i)
  }
  
  p <- cowplot::plot_grid(plotlist=ps)
  
  result = list("node_tbl" = node_tbl,
                "graph_obj" = g_obj,
                "plotlist" = ps, 
                "plot" = p,
                "prots" = prots,
                "params" = list("species" = species, 
                                "network_type" = network_type, 
                                "score" = score,
                                "common_legend" = common_legend, 
                                "node_deg_above" = node_deg_above,
                                "add_nodes" = add_nodes))
  
  return(result)
  }else{
    return(NULL)
  }
}


#' Spectronaut to se
#' 
#' Reads in spectronaut or other output and creates a summarized experiment.
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
#' sample names and groups are read automatically (start with "log2quantity_"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP2 make_unique
#' @examples
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#' @export

spectronaut_read_in <- function(file, gene_column = "genes", protein_column = "uni_prot_ids", sep="_rep_", keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  #####read in data (spectronaut output from Fatih Demir)
  df <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
    rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
    janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
    dplyr::rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
    dplyr::rename("condition" = "r_condition", "file_name" = "r_file_name", "replicate" = "r_replicate") #remove the strange r_ from those three column names
  
  #turn dataframe (df) to wide format, which is required for summarised experiment (or perseus)
  data <- df %>% 
    mutate(cond_rep = paste0(condition, "_rep", replicate )) %>% # make yourself a new column with the condition and replicate pasted together
    DEP2::make_unique("genes", "uni_prot_ids", delim = ";") %>%#this is a function of DEP2 to make unique names of genes
    select(cond_rep, ID,genes,uni_prot_ids,  protein_descriptions,log2quantity, ibaq) %>% #I just select these columns because otherwise the 'pivot_wider' function is not working because each row contains non-unique info so everything is matched to everything
    tidyr::pivot_wider(names_from = cond_rep, values_from = c(log2quantity, ibaq)) %>% #make wide
    select(-ID) %>% #now I remove the ID col that was added when I ran make_unique, because I need to make unique again, because later I will need both ID and Name column
    DEP2::make_unique("genes", "uni_prot_ids", delim = ";") %>% #re-run
    janitor::clean_names() %>% #clean_up col-names
    dplyr::rename(ID = id) %>% #clean up col names further
    dplyr::rename_all(list(~as.character(gsub("0_2", "0p2", .)))) %>% #clean up further, as this underscore separation is annoying later on 
    mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x)) %>% # turns Infinite to NA
    mutate_if(is.numeric, function(x) ifelse(is.na(x), NA, x)) %>% #turns NaN to NA
    mutate_if(is.numeric, function(x) 2^x) %>% #inverse log2 of all numerics (careful if you change to include anything else, like the ibaqs that have not been inverse log2)
    mutate_if(is.numeric, function(x) ifelse(is.na(x), NA, x))  %>% #turns NaN to NA
    dplyr::rename("gene_names" = gene_column,
                  "protein_ids" = protein_column)
  
    #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_ids,
           orig_gene_names = gene_names) %>%
    split_genes(colname = "protein_ids", keep_all = keep_all_proteins) %>%
    split_genes(colname = "gene_names", keep_all = keep_all_genes)
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    filename <- df %>% 
      mutate(cond_rep = paste0(tolower(condition), "_rep", replicate)) %>%
      select(file_name, cond_rep) %>%
      unique() %>%
      mutate(cond_rep = janitor::make_clean_names(cond_rep))
    
    LFQ_labels <- colnames(data)[grep("^log2quantity_", colnames(data))]
    
    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("^log2quantity_", "", LFQ_labels),
                                    condition=gsub(paste0("^log2quantity_|",sep, "[0-9].*"), "", LFQ_labels),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE)) %>%
      merge(filename, by.x = "sample", by.y = "cond_rep", all = T)
  }else{
    LFQ_labels<-experimental_design$label
  }
  
  data_se <- make_se(data, which(colnames(data) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}




#' fragpipe to se
#' 
#' Reads in fragpipe or other output and creates a summarized experiment.
#'
#' @param file path to proteinGroups or similar file
#' @param gene_column name of gene_name column after janitor
#' @param protein_column name of protein column after janitor
#' @param sep character describing the separator between sample name and replicate number (e.g. "_rep_", "_r")
#' @param experimental_design dataframe with information regarding samples. If not specified,
#' sample names and groups are read automatically (start with "log2quantity_"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP2 make_unique
#' @examples
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#' @export

fragpipe_read_in <- function(file, gene_column = "gene", protein_column = "protein_id", sep="_rep_", experimental_design = ""){
  
  #####read in data (spectronaut output from Fatih Demir)
  df <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
    rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
    janitor::clean_names() #make column names clean and unique (makes later coding easier)
    
  #turn dataframe (df) to wide format, which is required for summarised experiment (or perseus)
  data <- df %>% 
    dplyr::rename("gene_names" = gene_column,
                  "protein_ids" = protein_column) %>%
    make_unique("gene_names", "protein_ids")
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    LFQ_labels <- colnames(data)[grep(".*_max_lfq_intensity", colnames(data))]
    
    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("_max_lfq_intensity", "", LFQ_labels),
                                    condition=gsub(paste0("_max_lfq_intensity|",sep, "[0-9]*"), "", LFQ_labels),
                                    replicate=gsub(paste0("_max_lfq_intensity|^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE)) 
    }else{
    LFQ_labels<-experimental_design$label
  }
  
  data_se <- make_se(data, which(colnames(data) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}



#' Create a 2D scatterPlot
#' 
#' Plots two columns against each other and marks entries that differ more than standard_dev SD. This calculation is performed as a moving window.
#' @param df data for plotting
#' @param col_x string column name of x-axis data 
#' @param col_y string column name of y-axis data
#' @param col_label string column name of labels
#' @param show_labels Boolean show labels on plot
#' @param title string plot title
#' @param standard_dev int difference to sd that specifies highlighted entries 
#' @param window rolling window for which sd is calculated
#' @return list containing the scatterplot and dfs for highlighted entries
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom stats cor.test
#' @importFrom ggrepel geom_text_repel
#' @importFrom DescTools Closest
#' @importFrom stringr str_wrap
#' @export

scatterPlot <- function(df, col_x, col_y, col_label = "gene_names", show_labels = TRUE, title = "", standard_dev=2, window=1) {
  x <- df[,col_x]
  y <- df[,col_y]
  col_label <- df[,col_label]
  mod = lm(unlist(y)~unlist(x))
  slope=mod$coefficients["x"]
  intercept=if (is.na(mod$coefficients["(Intercept)"])) 0 else mod$coefficients["(Intercept)"]
  angle = atan(slope)
  
  rotate = function (x,y, theta, intercept){
    c(cos(theta) * x + sin(theta) * (y-intercept),
      -sin(theta) * x + cos(theta) * (y-intercept))
  }
  
  ## Rotate so that slope is 0
  rot = mapply(rotate, x, y, angle, intercept)
  xp = rot[1,]
  yp = rot[2,]
  
  ## Compute standard deviation along the x axis at position pos,
  ## for data in a given window (a lower window will be less smooth)
  sdWindow = function (pos, window, x, y) {
    sel = which (x >= pos-window/2 & x <= pos+window/2) 
    sd(y[sel], na.rm = TRUE)
  }
  
  xsd = seq(min(xp, na.rm = TRUE),max(xp, na.rm = TRUE), length=100)
  ysd = sapply(xsd, sdWindow, window, xp, yp)
  
  ## Transform back the data, and two curves for the + and - standard dev
  xpp = mapply(rotate, xp, yp, -angle, -intercept)
  sdp = mapply(rotate, xsd, standard_dev*ysd, -angle, -intercept)
  sdm = mapply(rotate, xsd, -standard_dev*ysd, -angle, -intercept)
  
  sdpdf = data.frame(x=sdp[1,], y=sdp[2,])
  sdmdf = data.frame(x=sdm[1,], y=sdm[2,])
  
  cbPalette <- c("#dddddd","#0000ff", "#00ff00")
  df = data.frame(x=x, y=y)
  
  up<-data.frame("x"=lapply(df$x, function(i) min(Closest(x=sdpdf$x,a=i)))%>%unlist(),"x_old"=df$x,"y"=df$y, "col_label"=col_label)
  up_col_labels<-merge(up,sdpdf, by="x")%>%filter(y.x>=y.y)
  
  down<-data.frame("x"=lapply(df$x, function(i) min(Closest(x=sdmdf$x,a=i)))%>%unlist(),"x_old"=df$x,"y"=df$y, "col_label"=col_label)
  down_col_labels<-merge(down,sdmdf, by="x")%>%filter(y.x<=y.y)
  
  pearson<-cor.test(x,y)
  
  p = ggplot(df, aes(x=x, y=y))
  
  if(show_labels == TRUE){
    p <- p + geom_point (shape=20, color="black", size=3, alpha = 0.5) +
      geom_abline (intercept=intercept, slope=mod$coefficients["x"], color="darkgrey", linetype=2) +
      geom_line (data=sdpdf, aes(x=x, y=y), color="red", alpha=0.5,linetype=2, size=1) +
      geom_line (data=sdmdf, aes(x=x, y=y), color="blue", alpha=0.5, linetype=2, size=1) +
      
      #scale_colour_manual(values=cbPalette) +
      scale_size_manual(values=c(2,3)) +
      geom_point(data=up_col_labels, aes(x=x_old, y=y.x), colour="red")+
      geom_text_repel(data=up_col_labels, aes(x=x_old, y=y.x,label=col_label), nudge_y = 0.1,fill = alpha(c("white"),0.3), label.padding = 0.5, size=3, max.overlaps = Inf)+
      geom_point(data=down_col_labels, aes(x=x_old, y=y.x), colour="blue")+
      geom_text_repel(data=down_col_labels, aes(x=x_old, y=y.x,label=str_wrap(col_label,30)),nudge_y=-0.2,fill = alpha(c("white"),0.3), label.padding = 0.5, size=3, max.overlaps = Inf)+
      geom_label(aes(x=ifelse(any(up_col_labels$x_old<= -0.5), 0.7,-1),y=0.8),label=paste("r: ",as.character(signif(pearson$estimate,digits=3)) , "\np-value: ",as.character(signif(pearson$p.value, digits=3))), hjust=0, fontface="bold", max.overlaps = Inf)+
      theme_bw () + 
      theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5), legend.position = "none", 
            axis.title = element_text(face="bold", lineheight=0.6))+
      labs(x=col_x, y=col_y, title=title)
  }else{
    p <- p + geom_point (shape=20, color="black", size=3, alpha = 0.5) +
      geom_abline (intercept=intercept, slope=mod$coefficients["x"], color="darkgrey", linetype=2) +
      geom_line (data=sdpdf, aes(x=x, y=y), color="red", alpha=0.5,linetype=2, size=1) +
      geom_line (data=sdmdf, aes(x=x, y=y), color="blue", alpha=0.5, linetype=2, size=1) +
      
      #scale_colour_manual(values=cbPalette) +
      scale_size_manual(values=c(2,3)) +
      geom_point(data=up_col_labels, aes(x=x_old, y=y.x), colour="red")+
      geom_point(data=down_col_labels, aes(x=x_old, y=y.x), colour="blue")+
      geom_label(aes(x=ifelse(any(up_col_labels$x_old<= -0.5), 0.7,-1),y=0.8),label=paste("r: ",as.character(signif(pearson$estimate,digits=3)) , "\np-value: ",as.character(signif(pearson$p.value, digits=3))), hjust=0, fontface="bold", max.overlaps = Inf)+
      theme_bw () + 
      theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5), legend.position = "none", 
            axis.title = element_text(face="bold", lineheight=0.6))+
      labs(x=col_x, y=col_y, title=title)    
  }
  return(list(plot = p, up_col_labels = up_col_labels, down_col_labels = down_col_labels))
}


#' Read in Phosphoproteomics Data and Perform Initial Processing
#'
#' This function reads in phosphoproteomics data from the "Phospho (STY)Sites.txt" MaxQuant output, processes it by splitting protein groups,
#' filtering low-quality hits, and expanding the site table. It also creates site IDs and makes gene names unique. In the end, it returns a 
#' SummarizedExperiment.
#' Experimental design can be provided or generated automatically if not provided.
#'
#' @param file A string representing the file path of the tab-delimited phosphoproteomics data.
#' @param gene_column A string representing the name of the gene column in the input data. Default is "gene_names".
#' @param protein_column A string representing the name of the protein column in the input data. Default is "proteins".
#' @param sep A string used as the separator for intensity columns in the input data. Default is "_rep_".
#' @param filt A character vector containing the filters to apply for false and low-quality hits. Default is c("reverse", "potential_contaminant").
#' @param keep_all_proteins A boolean indicating whether to keep all proteins when splitting protein groups. Default is FALSE.
#' @param keep_all_genes A boolean indicating whether to keep all genes when splitting protein groups. Default is FALSE.
#' @param experimental_design A data.frame containing the experimental design information. Default is an empty string.
#' 
#' @importFrom vroom vroom
#' @importFrom janitor make_clean_names
#' @importFrom dplyr select filter
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom DEP2 make_se make_unique
#' 
#' @return A summarizedExperiment containing the phosphoproteomics data.
#' @examples
#' # Assuming "Phospho (STY)Sites.txt" is a valid MaxQuant phosphoproteomics data file
#' processed_data <- phos_read_in_int("Phospho (STY)Sites.txt")
#' @export
phos_read_in_int <- function(file, gene_column = "gene_names", protein_column = "proteins", sep="_rep_", filt = c("reverse", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  data <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names)
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = .[[protein_column]],
           orig_gene_names =  .[[gene_column]]) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes) %>%
    dplyr::select(-contains("score"))
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ is.na(.x)))
  
  # intensity columns
  cols <- grep(paste0(".*",sep,"[0-9]*_[0-9]*$"), colnames(data), value = TRUE)
  
  #expand site_table
  data <-  data %>% tidyr::pivot_longer(cols = cols, names_to = c("set", "multiplicity"), values_to = "int", names_pattern = "(.*)_([0-9]*)$") %>% mutate(set = gsub("intensity", "value", set)) %>% tidyr::pivot_wider(values_from = "int", names_from = "set")
  
  #create site id
  data <- data %>% mutate(site_id = paste0(gene_names, "_", amino_acid, position), 
                          site_id_mult = paste0(gene_names, "_", amino_acid, position, "_", multiplicity))
  
  #Make gene_names unique
  data_unique <- make_unique(data, "site_id_mult", protein_column, delim=";")
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    
    labels <- grep(paste0("intensity_.*", sep, "[0-9]*$"), colnames(data), value = TRUE)
    
    experimental_design <- data.frame(label=labels, 
                                      sample=gsub("intensity_", "", labels),
                                      ID = gsub("intensity_", "", labels),
                                      condition=gsub(paste0("intensity_|",sep, "[0-9]*"), "", labels),
                                      replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", labels, perl = TRUE))
    
  }else{
    experimental_design <- data.frame(experimental_design)
  }
  
  data_se <- DEP2::make_se(data_unique, which(colnames(data_unique) %in% gsub("intensity", "value",labels)), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "intensity_raw"
  return(data_se)
}

#' Read in phosphoproteomics occupancy data and process it
#'
#' This function reads in phosphoproteomics occupancy data from the "Phospho (STY)Sites.txt" MaxQuant output,
#' processes the data, filters false and low-quality hits, and returns a SummarizedExperiment object.
#'
#' @param file A character string specifying the input file path.
#' @param gene_column A character string specifying the column name for gene names in the input file (default: "gene_names").
#' @param protein_column A character string specifying the column name for protein names in the input file (default: "proteins").
#' @param sep A character string specifying the separator for the replicate information in the column names (default: "_rep_").
#' @param filt A character vector specifying the columns to filter for NA values (default: c("reverse", "potential_contaminant")).
#' @param keep_all_proteins A logical value indicating whether to keep all proteins when splitting protein groups (default: FALSE).
#' @param keep_all_genes A logical value indicating whether to keep all genes when splitting gene names (default: FALSE).
#' @param experimental_design A data frame specifying the experimental design. If not provided, it will be generated based on the input data.
#'
#' @importFrom vroom vroom
#' @importFrom janitor make_clean_names
#' @importFrom dplyr select filter
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom DEP2 make_se make_unique
#' 
#' @return A SummarizedExperiment object containing the processed phosphoproteomics occupancy data.
#' @examples
#' # Example input file path (replace with your own file path)
#' input_file <- "path/to/your/input_file.tsv"
#'
#' # Read and process the data
#' data_se <- phos_read_in_occ(input_file)
#'@export
phos_read_in_occ <- function(file, gene_column = "gene_names", protein_column = "proteins", sep="_rep_", filt = c("reverse", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  data <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names)
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = .[[protein_column]],
           orig_gene_names =  .[[gene_column]]) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes)
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ is.na(.x)))
  
  # intensity columns
  labels <- colnames(data)[grepl("occupancy_", colnames(data)) &
                             !grepl("occupancy_(ratio|error)", colnames(data))]
  
  #create site id
  data <- data %>% mutate(site_id = paste0(gene_names, "_", amino_acid, position))
  
  #Make gene_names unique
  data_unique <- make_unique(data, "site_id", protein_column, delim=";")
  
  #2^x because it is later log2 transformed
  data_unique <- data_unique %>% mutate(across(all_of(labels), ~2^.x))
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    experimental_design <- data.frame(label=labels, 
                                      sample=gsub("occupancy_", "", labels),
                                      ID = gsub("occupancy_", "", labels),
                                      condition=gsub(paste0("occupancy_|",sep, "[0-9]*"), "", labels),
                                      replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", labels, perl = TRUE))
    
  }else{
    experimental_design <- data.frame(experimental_design)
  }
  
  data_se <- make_se(data_unique, which(colnames(data_unique) %in% labels), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "occupancy"
  
  assay(data_se)[assay(data_se) == "NaN"] <- NA
  return(data_se)
}

#' Convert a GCT object to a long format data frame
#'
#' This function takes a GCT object (or a file containing GCT data) and converts it to a long format data frame.
#'
#' @param gct A GCT object. If the file parameter is provided, this parameter will be ignored.
#' @param file A file path to a GCT file. If provided, the function will parse the GCT data from the file. Default is an empty string.
#' @return A long format data frame with columns for id.y, id.x, and metadata.
#' 
#' @importFrom cmapR parse_gctx melt_gct
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom plyr colwise
#' @importFrom dplyr filter
#' 
#' @examples
#' # Use an example GCT object from cmapR package
#' example_gct <- cmapR::example_gct
#' long_format_data <- gct_to_long(example_gct)
#' head(long_format_data)
#' @export
gct_to_long <- function(gct, file = ""){
  
  if(!file == ""){
    gct <- cmapR::parse_gctx(file)
  }
  
  m <- cmapR::melt_gct(gct)
  
  m <- m %>% 
    tidyr::pivot_longer(grep(paste(unique(.$id.y), 
                                   collapse = "|"), 
                             colnames(.)), 
                        names_to = c("type", "set"), 
                        values_to = "val", 
                        names_pattern = paste0("^(.*)\\.(.*_[0-9]*)"), 
                        values_transform = list(val = as.character)) %>% 
    filter(set == id.y) %>% 
    tidyr::pivot_wider(names_from = "type", 
                       values_from = "val")
  
  m <- plyr::colwise(type.convert)(m)
  
}


#' Prepare Data for Single Sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' This function prepares the data from a SummarizedExperiment object for ssGSEA analysis. It removes duplicated rows based on the rowData sequence_window, creates a GCT object, and optionally writes the GCT object to a file.
#'
#' @param se A SummarizedExperiment object.
#' @param file A character string indicating the file name to save the GCT object. Defaults to an empty string, which means the GCT object will not be saved to a file.
#' @return A list with the following elements:
#' \itemize{
#'   \item{"data"}{A GCT object with the unique rows from the input SummarizedExperiment object.}
#'   \item{"path"}{A character string indicating the path of the saved GCT file, if applicable.}
#'   \item{"duplicated"}{A SummarizedExperiment object containing the duplicated rows removed from the input object.}
#' }
#' @importFrom cmapR write_gct
#' @examples
#' # Create a SummarizedExperiment object (se) here
#' # ...
#' result <- prep_ssgsea2(se)
#' result <- prep_ssgsea2(se, file = "output")
#' @export
prep_ssgsea2 <- function(se, file = ""){
  
  
  rid <- gsub("^.{8}(.*).{8}$", 
              "\\1-p", 
              rowData(se)$sequence_window)
  
  dupl <- se[duplicated(rid),]
  uni <- se[!duplicated(rid),]
  
  temp <- new("GCT", 
              mat=assay(uni), 
              rdesc=as.data.frame(rowData(uni)), 
              cdesc=as.data.frame(colData(uni)), 
              rid = gsub("^.{8}(.*).{8}$", 
                         "\\1-p", 
                         rowData(uni)$sequence_window))
  
  if(!file == ""){
    cmapR::write_gct(temp, file)
    return(list("data" = temp, "path" = paste0(file, "_n", ncol(se), "x", nrow(uni), ".gct"), "duplicated" = dupl))
  } else{
    return(list("data" = temp, "path" = paste0(file, "_n", ncol(se), "x", nrow(uni), ".gct"), "duplicated" = dupl))
  }
  
}

#' Convert a SummarizedExperiment object to a long format data frame
#'
#' This function takes a SummarizedExperiment object and converts it to a long format data frame.
#' It allows the user to specify the assays to include in the output data frame.
#'
#' @param se A SummarizedExperiment object.
#' @param assays_ A character vector of assay names to include in the output data frame. Default is an empty character vector, which selects all assays.
#' @return A data frame in long format with assay values, rowData, and colData combined.
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object
#' result <- to_long(se, assays_ = c("assay1", "assay2"))
#' 
#' # To select all assays
#' result_all <- to_long(se)
#' @export
to_long <- function(se, assays_ = c("")){
  
  if(length(assays_) == 0){
    assays_ <- names(assays(se)) 
  }
  
  df <- lapply(assays(se)[assays_], "as.data.frame")
  
  df <- do.call("cbind", df) %>% cbind(as.data.frame(rowData(se)))
  
  assays <- grep(paste(assays_, "\\.", sep = ""), colnames(df), value = TRUE) %>% gsub("\\..*","",.) %>% unique()
  
  df_ <- df %>% dplyr::select(-ID) %>% tidyr::pivot_longer(grep(paste(assays_,"\\.", sep = ""), colnames(df), value = TRUE), names_to = c("assays", "ID"), names_pattern = "(.*)\\.(.*)")
  
  df_ <- merge(df_, as.data.frame(colData(se)), by.x = "ID")
  #%>% tidyr::pivot_longer(paste0("assay_",names(assays(se))))
  
  return(df_)
}


#' Get column data from a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @return A data.frame containing the column data.
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object
#' get_coldata(se)
#' @export
get_coldata <- function(se) {
  return(as.data.frame(colData(se)))
}

#' Get row data from a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @return A data.frame containing the row data.
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object
#' get_rowdata(se)
#' @export
get_rowdata <- function(se) {
  return(as.data.frame(rowData(se)))
}

#' Custom theme for ggplot2
#'
#' This function applies a custom theme to a plot by combining the DEP2::theme_DEP1()
#' function with additional theme elements.
#'
#' @return A ggplot2 theme object.
#' @importFrom DEP2 theme_DEP1
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
#' p + my_theme()
#' @export
my_theme <- function() {
  DEP2::theme_DEP1() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Add Significance Information to a SummarizedExperiment Object
#'
#' This function takes a SummarizedExperiment object, calculates significance information,
#' and adds it to the rowData.
#'
#' @param se A SummarizedExperiment object.
#' @param p_thr A numeric value for the p-value threshold (default: 0.05).
#' @param diff_thr A numeric value for the difference threshold (default: 1).
#' @return A SummarizedExperiment object with added significance information.
#' @importFrom tidyr pivot_longer pivot_wider
#' @import dplyr
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object
#' se_with_sign <- add_sign(se)
#' @export
add_sign <- function(se, p_thr = 0.05, diff_thr = 1){
  se_long <- get_rowdata(se) %>% 
  dplyr::select(c(site_id_mult, 
                  ends_with("_diff"), 
                  ends_with("p.adj"), 
                  ends_with("p.val"), 
                  ends_with("CI.L"), 
                  ends_with("CI.R"))) %>% 
    tidyr::pivot_longer(cols = -site_id_mult, 
                        names_to = c("label", "type"), 
                        values_to = "value", 
                        names_pattern = "(.*)_(.*)") %>% 
    filter(!label == "score") %>%
    tidyr::pivot_wider(names_from = "type", 
                       values_from = "value") %>% 
    mutate(significant = ifelse(p.adj < p_thr & (diff > diff_thr | diff < -diff_thr), TRUE, FALSE)) %>% 
    tidyr::pivot_wider(names_from = "label", 
                       values_from = c("p.val","p.adj", "significant", "diff", "CI.L", "CI.R"), 
                       names_glue = "{label}_{.value}") %>% 
    mutate(significant = ifelse(if_any(ends_with("significant"), ~ .x == TRUE), TRUE, FALSE)) %>%
    dplyr::select(-site_id_mult)
  
  rowData(se)[, colnames(se_long)] <- NULL
  
  rowData(se) <- cbind(rowData(se), se_long)
  
  return(se)

}

#' Long Test
#'
#' This function returns test results for an se object in long format.
#'
#' @param se An se object.
#' @return A data frame.
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect ends_with
#' @import dplyr
#' @examples
#' # Example usage:
#' # result <- long_test(se)
#' @export
long_test <- function(se){
  
  res <- se %>% 
    get_rowdata() %>% 
    select(-significant) %>% 
    tidyr::pivot_longer(c(ends_with("_p.adj"), 
                          ends_with("_p.val"), 
                          ends_with("_diff"), 
                          ends_with("_significant"), 
                          ends_with("_CI.L"), 
                          ends_with("_CI.R")), 
                        values_to = "val", 
                        names_to = c("contrast", "pn"), 
                        names_pattern = "(.*)_(.*)$") %>%
    dplyr::select(-ends_with("p.val"), -ends_with("p.adj")) %>%
    tidyr::pivot_wider(names_from = "pn", 
                       values_from = "val") %>% 
    mutate(significant = ifelse(significant == 0, 
                                FALSE, 
                                TRUE))
  return(res)
}

#' Create a clustered heatmap from a given SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @param indicate A character string indicating the variable to be used for annotation.
#' @param type A character string indicating the type of heatmap.
#' @param k An integer specifying the number of clusters to be used in k-means clustering.
#' @param ... Additional arguments to be passed to the underlying complexHeatmap::Heatmap function.
#' 
#' @return A list containing the heatmap plot, a data frame with the heatmap data, and a data frame with the protein clusters.
#'
#' @importFrom DEP2 plot_heatmap
#' @importFrom gdata cbindX
#' @import dplyr
#' @export
#' @examples
#' # Assuming se_diff is a SummarizedExperiment object with appropriate data
#' result <- clustered_heatmap(se_diff, indicate = "condition", type = "centered", k = 5)
#' result$plot
#' result$df
#' result$clusters
clustered_heatmap <- function(se, indicate = "condition", type = "centered", k = 3, ...){
  p_heatmap <- se %>% DEP2::plot_heatmap(indicate = indicate, type = type, kmeans = TRUE, k = k)  

  p_heatmap_data <- se %>% DEP2::plot_heatmap(indicate = indicate, type = type, kmeans = TRUE, k = k, plot = FALSE)  
  
  split <- p_heatmap_data %>% select(k, protein) %>% split.data.frame(.$k)
  split <- lapply(split, function(x) select(x,protein))
  
  split_df <- do.call(gdata::cbindX, split)
  colnames(split_df) <- paste0("cluster_", 1:length(split))
  rownames(split_df) <- NULL
  
  return(list("plot" = p_heatmap, "df" = p_heatmap_data, "clusters" = split_df))
}

#' Perform Over-representation Analysis (ORA) on Phosphoproteomics Data
#'
#' This function takes a SummarizedExperiment object and conducts Over-representation Analysis (ORA) on proteins with at least one significant phosphosite, using the clusterProfiler package.
#'
#' @param se A SummarizedExperiment object containing phosphoproteomics data.
#' @param contr A character vector specifying the contrasts to be analyzed. Default is "all".
#' @param OrgDb A character string specifying the organism database to be used. Default is "org.Hs.eg.db".
#' @param pvalueCutoff Numeric value for the p-value cutoff. Default is 0.4.
#' @param qvalueCutoff Numeric value for the q-value cutoff. Default is 0.8.
#' @param ont A character vector specifying the Gene Ontology (GO) categories to be analyzed. Default is c("BP", "MF", "CC").
#' @return A list containing the ora object, a combined dataframe and a ggplot DEP2icting terms with q < 0.2
#' @import dplyr
#' @import ggplot2
#' @importFrom clusterProfiler enrichGO
#' @importFrom DOSE parse_ratio
#' @export
phospho_ora <- function(se, contr = "all", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.4, qvalueCutoff = 0.8, ont = c("BP", "MF", "CC")){
  
  se_diff_test <- se %>% long_test()
  if(!contr == "all"){
    se_diff_test <- se_diff_test %>% filter(contrast %in% contr)
  }
  
  # Get significant genes
  sign_genes <- se_diff_test %>% filter(significant == TRUE) %>% select(gene_names, contrast, diff)  %>% distinct(gene_names, .keep_all = TRUE) %>% group_by(contrast) %>% summarize(n = n(), gene_names = list(gene_names), diff = list(diff))
  
  # Get background
  background <- se_diff_test$gene_names %>% unique()
  
  #Create list for ora
  names(sign_genes$gene_names) <- sign_genes$contrast
  
  # Perform ORA via clusterProfiler for each contrast and ontology
  ora_res <- list()
  
  for(ont_ in ont){
    ora_res[[ont_]] <- lapply(sign_genes$gene_names, 
                              clusterProfiler::enrichGO, 
                              OrgDb = OrgDb, 
                              keyType = "SYMBOL", 
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff,
                              ont = ont_, 
                              universe = background)
    
    for(sets_ in names(ora_res[[ont_]])){
      ora_res[[ont_]][[sets_]] <- ora_res[[ont_]][[sets_]]@result %>% mutate(contrast = sets_, ont = ont_)
    }
  }
  
  # Create one large table storing the results
  ora_res_df <- do.call("rbind", 
                        lapply(ora_res, 
                               function(x) do.call("rbind", x))
  ) %>%
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
  
  # Test if faceting variables have at least one value
  if(nrow(ora_res_df %>% filter(qvalue < 0.2)) > 0){
    #Keep only Terms that are significant at least once
    ora_plot <- ora_res_df %>% filter(qvalue < 0.2) %>%
      ggplot(aes(x = FoldEnrichment, 
                 y = reorder(Description, FoldEnrichment),
                 color = ifelse(pvalue < 0.05 & qvalue < 0.5, "Significant", "Non-significant"),
                 fill = contrast,
                 size = p.adjust))+
      geom_point(shape = 21,
                 alpha = 0.5) +
      scale_color_manual(values = c("grey", "black"))+
      scale_size_continuous(range = c(5,2)) +
      labs(title = "q < 0.5", color = "p < 0.05 & q < 0.5", fill = "contrast")+
      facet_grid(ont ~., scales = "free_y", space = "free_y") +
      # scale_color_manual(values = c("grey", "firebrick"))+
      geom_segment(aes(xend=0, yend = Description), color = "black", size = 0.1) +
      theme_bw()
    } else{
      ora_plot <- ggplot() + theme_void()
      }
  
  return(list("res" = ora_res, "df" = ora_res_df, "plot" = ora_plot))
}

#' Write phosphoproteomics data to a file
#'
#' This function takes a SummarizedExperiment object, extracts the relevant
#' data, and writes it to a file.
#'
#' @param se A SummarizedExperiment object containing phosphoproteomics data.
#' @param file The filename where the resulting data should be saved (default: current directory).
#' @return relevant data
#' @importFrom DEP2 get_df_wide
#' @import dplyr
#' @export
write_phos <- function(se, file = ""){
  exp <- se %>% 
    get_df_wide() %>%
    select(name,
           gene_names,
           protein,
           protein_names,
           amino_acid, 
           position,
           multiplicity,
           sequence_window,
           ends_with(c("_diff","_p.val","p.adj","significant")),
           2:(ncol(se))) 
  
  colnames(exp) <- gsub("_diff", "_log2FC", colnames(exp))
  
  if(!file == ""){
    openxlsx::write.xlsx(exp, file = file)
  }
  
  return(exp)
  
}

#' Write proteomics data to a file
#'
#' This function takes a SummarizedExperiment object, extracts the relevant
#' data, and writes it to a file.
#'
#' @param se A SummarizedExperiment object containing proteomics data.
#' @param file The filename where the resulting data should be saved (default: current directory).
#' @return relevant data
#' @importFrom DEP2 get_df_wide
#' @import dplyr
#' @export
write_prot <- function(se, file = ""){
  exp <- se %>% 
    get_df_wide() %>%
    select(name,
           gene_names,
           protein_ids,
           protein_descriptions,
           orig_prot_ids,
           orig_gene_names,
           ends_with(c("_diff","_p.val","p.adj","significant")),
           2:(ncol(se))) 
  
  colnames(exp) <- gsub("_diff", "_log2FC", colnames(exp))
  
  if(!file == ""){
    openxlsx::write.xlsx(exp, file = file)
  }
  
  return(exp)
  
}

#' Prepare KSEA input data from SummarizedExperiment object
#'
#' This function processes a SummarizedExperiment object and extracts relevant information
#' for the given contrast. The output is a data frame containing Protein, Gene, Peptide, Residue.Both,
#' p-value, and fold change (FC) columns, filtered to exclude rows with missing values.
#'
#' @param se A SummarizedExperiment object containing the input data.
#' @param contrast A character string specifying the contrast of interest.
#' @return A data frame with columns: Protein, Gene, Peptide, Residue.Both, p, and FC.
#' @importFrom rlang !! 
#' @importFrom ggplot2 sym
#' @import dplyr
#' @export
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'contrast' is a character string
#' # representing the contrast of interest:
#' prep_ksea_data <- prep_ksea(se, contrast)
prep_ksea <- function(se, contrast){
  diff_ <- paste0(contrast, "_diff")
  p_ <- paste0(contrast, "_p.val")
  
  rowData(se) %>% 
    as.data.frame() %>%
    select(proteins, 
           gene_names, 
           amino_acid, 
           position,
           contains(!! p_),
           contains(!! diff_)) %>%
    mutate(Residue.Both = paste0(amino_acid, position), 
           Peptide = "NULL",
           FC = 2^!!sym(diff_)) %>% 
    dplyr::rename("Protein" = "proteins", 
           "Gene" = "gene_names", 
           "p" = p_) %>% 
    select(Protein, 
           Gene, 
           Peptide, 
           Residue.Both, 
           p, 
           FC)  %>% 
    filter(!is.na(Gene)) %>%
    return()
}

#' Extract the centered substring of a given length from an input string
#'
#' This function extracts a substring of length `n` from the center of the input string. If `n` is greater than the length of the input string,
#' the entire input string is returned.
#'
#' @param input_string A character string from which the centered substring will be extracted.
#' @param n An integer specifying the length of the substring to be extracted.
#' @return A character string representing the centered substring of the input string.
#' @export
#' @examples
#' center_substring("hello", 3) # returns "ell"
#' center_substring("world", 2) # returns "or"
#' center_substring("example", 7) # returns "example"
center_substring <- function(input_string, n) {
  string_length <- nchar(input_string)
  
  start_pos <- ifelse(n > string_length, 1, ceiling((string_length - n) / 2) + 1)
  end_pos <- start_pos + n - 1
  
  return(substr(input_string, start_pos, end_pos))
}


#' Prepares Phosphorylation Data for Analysis
#'
#' This function takes a SummarizedExperiment object and prepares the data for
#' phosphorylation analysis using the PhosR package. It also generates a color palette
#' for plotting.
#'
#' @param se A SummarizedExperiment object.
#' @param species A character string specifying the species (default is "human").
#'        Must be one of "human", "mouse", or "rat".
#' @param numMotif An integer specifying the number of motifs (default is 5).
#' @param numSub An integer specifying the number of substrates (default is 1).
#' @param top An integer specifying the top number of substrates to select (default is 30).
#' @return A list containing the following elements:
#'         - kinaseSubstrateScore: A list of kinase substrate scores.
#'         - kinaseSubstratePred: A matrix of kinase substrate predictions.
#'         - mat: A matrix of collapsed phosphorylation data.
#'         - kinase_all_color: A color palette for plotting.
#'         - kinase_signalome_color" = A color palette for plotting.
#'         - seq: A vector of sequences.
#' @importFrom PhosR phosCollapse kinaseSubstrateScore kinaseSubstratePred
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
prep_phosR <- function(se, species = "human", numMotif = 5, numSub = 1, top = 30){
  
  data('PhosphoSitePlus', package = "PhosR")
  data("KinaseFamily", package = "PhosR")
  data("KinaseMotifs", package = "PhosR")
  
  mat <- assay(se)
  rownames(mat) <- paste0(gsub("_", ";", rownames(mat)), ";",gsub("(.*?);.*", "\\1", rowData(se)$sequence_window))
  
  mat <- phosCollapse(assay(se), id=gsub("(.*;.*;)[1-9]*;(.*)", "\\1\\2", rownames(mat), perl = TRUE), 
                      stat=apply(abs(assay(se)), 1, max), by = "max")
  
  seq <- gsub(".*;.*;(.*)", "\\1", rownames(mat))
  
  rownames(mat) <- gsub("(.*;.*;).*","\\1",rownames(mat))
  colnames(mat) <- se$ID
    
  if(species == "human"){
    psite = PhosphoSite.human
  } else if(species == "mouse"){
    psite = PhosphoSite.mouse
  } else if(species == "rat"){
    psite = PhosphoSite.mouse
  } else{
    return("Species must be one of human, mouse or rat")
  }
  
  kss <- kinaseSubstrateScore(psite, mat, seq,numMotif = numMotif, numSub = numSub, species = "human")
  set.seed(1)
  ksp <- kinaseSubstratePred(kss, top=top)
  
  
  # Color palette
  my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  kinase_all_color <- my_color_palette(ncol(kss$combinedScoreMatrix))
  names(kinase_all_color) <- colnames(kss$combinedScoreMatrix)
  kinase_signalome_color <- kinase_all_color[colnames(ksp)]
  
  
  return(list("kinaseSubstrateScore" = kss, "kinaseSubstratePred" = ksp, "mat" = mat, "kinase_all_color" = kinase_all_color, "kinase_signalome_color" = kinase_signalome_color, "seq" = seq))
}



#' Plot Signalome Map
#'
#' This function creates a Signalome Map plot using the PhosR::plotSignalomeMap function
#' and customizes the appearance with additional ggplot2 layers and a custom theme.
#'
#' @param signalome_res A data frame containing the Signalome results.
#' @param kinase_signalome_color A vector of colors for the kinases in the Signalome Map.
#'
#' @return A ggplot object representing the Signalome Map plot.
#' @import ggplot2
#' @importFrom PhosR plotSignalomeMap
#' 
#' @export
# Signalome Map
plot_signalome_map <- function(signalome_res, kinase_signalome_color){
  
  p <- PhosR::plotSignalomeMap(signalome_res, kinase_signalome_color) +
    scale_size_continuous(range = c(1,5)) + 
    scale_y_continuous(labels = paste0("Cluster ", levels(as.factor(signalome_res$proteinModules))), breaks = 1:length(unique(signalome_res$proteinModules))) +
    my_theme()
  
  return(p)
}


#' Create a bar plot with error bars for the given data
#'
#' This function takes a matrix and a signalome result, processes the data, and creates a bar plot with error bars.
#' It returns a list of ggplot objects, one for each module in the signalome result.
#'
#' @param mat A matrix containing the data to be plotted
#' @param signalome_res A list containing the results of the signalome analysis
#' @return A list of ggplot objects, one for each module in the signalome result
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @export

module_barplot <- function(mat, signalome_res){
  
  
  # Get data
  df <- mat %>% as.data.frame() %>% mutate(gene_names = gsub("(.*);(.*);", "\\1", rownames(mat))) %>% tidyr::pivot_longer(-gene_names, names_to = "ID", values_to = "scaled_intensity") %>%
    mutate(condition = gsub("(.*)_.*", "\\1", ID))
  
  mod_list <- split(names(signalome_res$proteinModules), signalome_res$proteinModules)
  names(mod_list) <- paste0("Cluster_", names(mod_list))
  
  mod_list <- lapply(mod_list, function(x) filter(df, gene_names %in% x))
  
  l <- do.call(cbind, lapply(mod_list, function(x) group_by(x, condition) %>% summarize(m = mean(scaled_intensity))))
  
  # Function to calculate mean and standard deviation
  mean_sd <- function(x) {
    y = mean(x, na.rm = TRUE)
    return(c(y = y,
             ymin = y - sd(x, na.rm = TRUE),
             ymax = y + sd(x, na.rm = TRUE)))
  }
  
  # Create the bar plot with error bars
  res <- lapply(seq_along(mod_list), function(z) ggplot(mod_list[[z]], aes(x = factor(condition), y = scaled_intensity, fill = factor(condition))) +
                  geom_boxplot()+
                  xlab("Condition") +
                  ylab("Scaled intensity") +
                  labs(title = names(mod_list)[[z]]) +
                  my_theme())
  
  return(res)
}


#' Prepares phosphoproteomics data for PhosR analysis
#'
#' This function processes phosphoproteomics data, calculates kinase-substrate scores, and generates predictions for the given data.
#' It also generates color palettes for visualization.
#'
#' @param ppe Phosphoproteomics data object
#' @param species Species of the data, one of "human", "mouse", or "rat" (default: "human")
#' @param numMotif Number of motifs to use (default: 5)
#' @param numSub Number of substrates to use (default: 1)
#' @param top Number of top kinases to consider (default: 30)
#' @param assay Type of assay data to use (default: "z_scored")
#'
#' @return A list containing kinaseSubstrateScore, kinaseSubstratePred, mat, kinase_all_color, kinase_signalome_color, and seq
#' @importFrom PhosR PhosphoExperiment kinaseSubstrateScore kinaseSubstratePred
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
prep_phosR_from_ppe <- function(ppe, species = "human", numMotif = 5, numSub = 1, top = 30, assay = "z_scored"){
  
  # Load datasets from PhosR
  data('PhosphoSitePlus', package = "PhosR")
  data("KinaseFamily", package = "PhosR")
  data("KinaseMotifs", package = "PhosR")
  
  # Get data from provided object
  mat <- ppe@assays@data[[assay]]
  
  seq <- gsub(".*;.*;(.*)", "\\1", rownames(mat))
  
  rownames(mat) <- gsub(".*;(.*;.*;).*","\\1",rownames(mat))
  colnames(mat) <- ppe$ID
  
  if(species == "human"){
    psite = PhosphoSite.human
  } else if(species == "mouse"){
    psite = PhosphoSite.mouse
  } else if(species == "rat"){
    psite = PhosphoSite.rat
  } else{
    return("Species must be one of human, mouse or rat")
  }
  
  kss <- kinaseSubstrateScore(psite, mat, seq, numMotif = numMotif, numSub = numSub, species = "human")
  set.seed(1)
  ksp <- kinaseSubstratePred(kss, top=top)
  
  
  # Color palette
  my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  kinase_all_color <- my_color_palette(ncol(kss$combinedScoreMatrix))
  names(kinase_all_color) <- colnames(kss$combinedScoreMatrix)
  kinase_signalome_color <- kinase_all_color[colnames(ksp)]
  
  
  return(list("kinaseSubstrateScore" = kss, "kinaseSubstratePred" = ksp, "mat" = mat, "kinase_all_color" = kinase_all_color, "kinase_signalome_color" = kinase_signalome_color, "seq" = seq))
}

#' AOV Scale PhosR Function
#'
#' This function takes a ppe object and applies a series of transformations,
#' including mean abundance calculation, ANOVA filtering, subsetting, standardization,
#' and z-scoring. The transformed data is then added as a new assay to the ppe object.
#'
#' @param ppe A ppe object.
#' @param p_cut A numeric threshold for selecting rows based on ANOVA p-value (default: 0.05).
#' @param fc_cut A numeric threshold for selecting rows based on fold change (default: 0.5).
#' @param assay A character specifying the assay to use from ppe (default: "normalised").
#' @return A ppe object with the transformed data added as a new assay named "z_scored".
#' @importFrom PhosR PhosphoExperiment meanAbundance matANOVA
#' @export
aov_scale_phosR <- function(ppe, p_cut = 0.05, fc_cut = 0.5, assay = "normalised"){
  mat <- ppe@assays@data[[assay]]
  mat.mean <- meanAbundance(mat, grps = colData(ppe)$condition)
  
  aov <- matANOVA(mat=mat, grps = colData(ppe)$condition)
  idx <- (aov < p_cut) & (rowSums(mat.mean > fc_cut) > 0)
  
  mat.reg <- mat[idx, ,drop = FALSE]
  
  mat.std <- standardise(mat.reg)
  colnames(mat.std) <- ppe$ID
  
  ppe <- ppe[idx,]
  ppe <- sev::add_assay(ppe, mat.std, "z_scored")
  
  return(ppe)
}

#' Test Differential Phosphorylation
#'
#' This function tests differential phosphorylation using a linear model fit and an empirical Bayes approach.
#'
#' @param ppe A PhosR object containing the phosphorylation data.
#' @param contrast A character vector specifying the contrasts to be tested. Default is "treatment_vs_control".
#' @param test_all A logical value indicating whether to test all possible contrasts or just the specified ones. Default is FALSE.
#' @param assay A character string specifying the assay to use for the analysis. Default is "scaled".
#'
#' @return A modified PhosR object with the differential phosphorylation test results added to the rowData.
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom S4Vectors expand.grid
#' @importFrom BiocGenerics cbind
#' @import dplyr
#' @export
test_diff_phosR <- function(ppe, contrast = c("treatment_vs_control"), test_all = FALSE, assay = "scaled"){
  
  # Create design matrix
  design <- model.matrix(~ ppe$condition - 1)
  colnames(design) <- gsub("ppe\\$condition", "", colnames(design))
  
  fit <- lmFit(ppe@assays@data[[assay]], design)
  
  # Create contrasts
  if(test_all == FALSE){
    contrast_ <- gsub("_vs_", "-", contrast)
  } else{
    contrast_ <- expand.grid(colnames(design), colnames(design)) %>% 
      filter(!Var1 == Var2) %>% mutate(contrast_ = paste0(Var1, "-", Var2)) %>%
      select(contrast_) %>%
      unlist()
  }
  
  # Create contrast matrix
  contrast.matrix <- makeContrasts(contrasts = contrast_, levels=design)
  
  # fit model
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  fit2 <- eBayes(fit2)
  
  # Create a named list of contrasts
  names(contrast_) <- contrast_
  
  # Get a list of all tested contrasts and rename the columns accordingly
  l <- lapply(contrast_, function(x) topTable(fit2, coef = x, number = Inf, sort.by = "none", confint = TRUE))
  l <- lapply(seq_along(l), function(i) dplyr::rename_with(l[[i]], ~ paste0(names(l)[i], "_", .x)))
  
  # Merge the results and add the F statistic 
  l_com <- cbind(do.call(cbind, l), topTable(fit2, number = Inf, sort.by = "none", confint = TRUE))
  colnames(l_com) <- gsub("-", "_vs_", 
                          gsub("logFC", "diff", 
                               gsub("adj.P.Val", "p.adj", 
                                    gsub("P.Value", "p.val", colnames(l_com)))))
  
  
  rowData(ppe) <- rowData(ppe)[, !colnames(rowData(ppe)) %in% c(colnames(l_com), paste0(contrast, "_significant"), "significant")]
  rowData(ppe) <- cbind(rowData(ppe), l_com)
  
  return(ppe)
}


#' Convert MaxQuant output to PhosphoExperiment object
#'
#' This function reads a MaxQuant output file, filters the data, and converts it into a PhosphoExperiment object.
#'
#' @param file The input file in tab-delimited format.
#' @param sep The separator used in column names for replicates (default: "_rep_").
#' @param filt A vector of filters to remove unwanted rows (default: c("reverse", "potential_contaminant")).
#' @param experimental_design An optional data frame to provide custom experimental design information.
#' @return A PhosphoExperiment object containing the processed data.
#' @importFrom janitor make_clean_names
#' @importFrom DEP2 make_unique
#' @importFrom BiocGenerics colnames rownames
#' @importFrom PhosR PhosphoExperiment
#' @import dplyr
#' @export
maxq_to_ppe <- function(file, sep="_rep_",
                        filt = c("reverse", "potential_contaminant"), experimental_design = NA){
  
  data <- read.delim(file, sep="\t")
  
  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = proteins,
           orig_gene_names = gene_names, 
           orig_positions_within_proteins = positions_within_proteins,
           orig_sequence_window = sequence_window,
           protein_ids = gsub("(.*?);.*", "\\1", proteins),
           gene_names = gsub("(.*?);.*", "\\1", gene_names),
           positions_within_proteins = gsub("(.*?);.*", "\\1", positions_within_proteins),
           sequence_window = gsub("(.*?);.*", "\\1", sequence_window))
  
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ .x == ""))
  

  # Create assays for ppe
  int <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]*$"), colnames(data))]) %>% log2() 
  mult_1 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_1"), colnames(data))]) %>% log2()
  mult_2 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_2"), colnames(data))]) %>% log2()
  mult_3 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_3"), colnames(data))]) %>% log2()
  
  # replace missing values with NA
  int[is.infinite(int)] <- NA 
  mult_1[is.infinite(mult_1)] <- NA
  mult_2[is.infinite(mult_2)] <- NA
  mult_3[is.infinite(mult_3)] <- NA
  
  # Change colnames
  label<- colnames(int)
  
  ID <- gsub(paste0("intensity_(.*)", sep, "(.*)"), "\\1_\\2", label)
  colnames(int) <- ID
  colnames(mult_1) <- ID
  colnames(mult_2) <- ID
  colnames(mult_3) <- ID
  
  # Create colData if not provided
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    experimental_design<-data.frame(label=label,
                                    sample=gsub("intensity_", "", label),
                                    ID = ID,
                                    condition=gsub(paste0("intensity_|",sep, "[0-9].*"), "", label),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", label, perl = TRUE))
  }
  
  #Create rowData
  rowdata <- data %>% select(-grep("intensity", colnames(data)), -starts_with("score"), -contains("_diff")) %>% mutate(gene_names = ifelse(gene_names == "", "NA", gene_names)) %>% mutate(psite = paste0(gene_names, ";", amino_acid, ";", positions_within_proteins)) %>% make_unique(names = "psite", ids = "psite", delim = "/")
  
  # Create PPE
  ppe <- PhosphoExperiment(assays = list(Quantification = int, 
                                         multiplicity_1 = mult_1,
                                         multiplicity_2 = mult_2,
                                         multiplicity_3 = mult_3),
                           rowData = rowdata,
                           colData = experimental_design,
                           Site = as.numeric(data$positions_within_proteins), 
                           GeneSymbol = data$gene_names, 
                           Residue = data$amino_acid, 
                           Sequence = data$sequence_window, 
                           UniprotID = data$protein_ids,
                           Localisation = as.numeric(data$localization_prob))
  
  rownames(ppe) <- paste(rowdata$protein_ids, rowdata$gene_names, paste0(rowdata$amino_acid, rowdata$positions_within_proteins), rowdata$sequence_window, sep = ";")
  return(ppe)
}


#' Read in Spectronaut Output
#'
#' This function reads in data from Spectronaut output, cleans and formats the data and finally returns the data as a summarized experiment.
#'
#' @param candidates Character string representing the path to the candidates file. Alternatively, a data frame can be passed.
#' @param report Character string representing the path to the report file. Alternatively, a data frame can be passed.
#' @param contrasts Character vector representing the contrasts.
#' @param conditionSetup Character string representing the path to the condition setup file. Alternatively, a data frame can be passed.
#' @param quant_col Character string representing the column name for the quantification data (default: "log2quantity"). If other than "log2quantiy", data is log2 transformed in the process!
#' 
#' @return A data frame with the processed and merged data.
#'
#' @export
#'
#' @importFrom vroom vroom
#' @importFrom janitor clean_names
#' @importFrom tidyr pivot_wider
#' @importFrom DEP2 make_se
#' @import dplyr

spectronaut_to_se <- function(candidates = NULL, report = NULL, contrasts = NULL, conditionSetup = NULL, quant_col = "log2quantity"){
    #####read in data (spectronaut output from Fatih Demir)
    if(typeof(candidates) == "character"){
      df_candidates <- vroom::vroom(candidates, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("r_", "", .)))
    }else{
      df_candidates <- candidates
    }
    #####read in data (spectronaut output from Fatih Demir)  
    if(typeof(report) == "character"){
      df_wide_report <- vroom(report, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("r_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("x[0-9]*_", "", .)))
    }else{
      df_wide_report <- report
    }
    #####read in colData
    if(typeof(conditionSetup) == "character"){
      coldata <- vroom(conditionSetup, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        mutate(label = paste0(condition, "_", replicate)) %>% 
        dplyr::select(run_label, condition, replicate, file_name) %>%
        dplyr::rename(label = run_label) %>% 
        mutate(label = paste0(janitor::make_clean_names(label),"_", quant_col))
    }else{
      coldata <- conditionSetup
    }
    
    ##### get contrasts in correct "orientation"
    if(!is.null(contrasts)){
      df_candidates_new <- data.frame()
      for(con in contrasts){
        condition_num <- strsplit(con, "_vs_")[[1]][1]
        condition_den <- strsplit(con, "_vs_")[[1]][2]
        temp <- df_candidates %>% filter(condition_numerator == condition_num & condition_denominator == condition_den) %>%
          mutate(diff = avg_log2_ratio,
                 contrast = con)
        temp2 <- df_candidates %>% filter(condition_numerator == condition_den & condition_denominator == condition_num) %>%
          mutate(diff = -1 * avg_log2_ratio,
                 contrast = con)
        df_candidates_new <- rbind(df_candidates_new, temp, temp2)
      }
    }else{
      df_candidates_new <- df_candidates %>% 
        mutate(diff = avg_log2_ratio, 
               contrast = paste0(condition_numerator, "_vs_", condition_denominator))
    }
    ##### rename columns
    df_candidates_new <- df_candidates_new %>% dplyr::rename(p.val = pvalue, p.adj = qvalue) %>% 
      dplyr::select(diff, contrast, p.val, p.adj, protein_groups)
    ##### pivot
    df_candidates_new <- df_candidates_new %>% tidyr::pivot_wider(names_from = "contrast", values_from = c("diff", "p.val", "p.adj"), names_glue = "{contrast}_{.value}")
    ##### merge with report
    combined <- merge(df_candidates_new, df_wide_report, by = "protein_groups", all.y = T) %>%
      make_unique("genes", "protein_groups", delim = ";") %>%
      mutate(gene_names = genes,
             protein_ids = protein_groups)
    
    
    return(DEP2::make_se(combined, grep(paste0(".*_", quant_col), colnames(combined)), coldata, log2transform = ifelse(quant_col == "log2quantity", F, T))) 
  }

#' Create a correlation plot with different conditions and labels
#'
#' This function takes a dataframe and creates a scatter plot with correlation
#' analysis, highlighting points based on certain criteria such as being in a
#' gene list or having a rank among the top N differences.
#'
#' @param df A data frame with columns for x, y, and label names.
#' @param x_column Character string specifying the column name for the x-axis values.
#' @param y_column Character string specifying the column name for the y-axis values.
#' @param label_names Character string specifying the column name for the labels.
#' @param gene_list A character vector of gene names to be highlighted. Default is an empty vector.
#' @param top_genes Integer specifying the number of top differing proteins to consider. Default is  20.
#' @param max.overlaps Integer specifying the maximum number of overlapping labels. Default is Inf.
#' 
#' @return A list containing the plot and subsets of the data used for each part of the plot.
#' @export
#'
#' @importFrom patchwork plot_layout
#' @importFrom ggpubr stat_cor
#' @importFrom ggrepel geom_label_repel
#' @import ggplot2
#' @import dplyr

corr_plot <- function(df, x_column, y_column, label_names, gene_list = c(""), top_genes = 20, max.overlaps = Inf) {
  x_char <- df %>% dplyr::select({ {x_column} }) %>% colnames()
  y_char <- df %>% dplyr::select({ {y_column} }) %>% colnames()
  
  df <- df %>% dplyr::select({ {x_column} }, { {y_column} }, { {label_names} }) %>%
    filter(!is.na({ {x_column} }) | !is.na({ {y_column} })) %>%
    mutate(abs_diff = ifelse((is.na({{x_column}}) | is.na({{y_column}})), NA, abs({{x_column}} - {{y_column}})),
           in_gene_list = {{label_names}} %in% gene_list)
  
  
  sub_1 <- df %>% filter(!is.na({ {x_column} }) & !is.na({ {y_column} })) %>% 
    mutate(top_genes_common = rank(-abs_diff, ties.method ="first") < top_genes,
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_common, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_2 <- sub_1 %>% filter(top_genes_common & !in_gene_list)
  sub_3 <- sub_1 %>% filter(in_gene_list)
  
  x_min = min(dplyr::select(sub_1, c({ {x_column} }, { {y_column} })))
  x_max = max(dplyr::select(sub_1, c({ {x_column} }, { {y_column} })))
  
  p1 <- ggplot(sub_1, aes(x = { {x_column} }, y = { {y_column} }, label = { {label_names} })) +
    geom_point(aes(fill = target, shape = factor(top_genes_common, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    scale_shape_manual(values = c("FALSE"= 21, "TRUE" = 23), drop = FALSE) +
    geom_label_repel(data = filter(sub_3, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = max.overlaps, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_2, top_genes_common), color = "firebrick", show.legend = FALSE, max.overlaps = max.overlaps, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    geom_abline(intercept = 0, linetype = "solid", color = "black") +
    ggpubr::stat_cor(data = sub_1, aes(label = ..r.label..)) +
    labs(title = "Common proteins", fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(x_min, x_max), y = c(x_min, x_max)) +
    theme_bw() 
  
  sub_4 <- df %>% filter(!is.na({{x_column}}) & is.na({{y_column}})) %>% 
    mutate(`:=`({{y_column}}, "NA")) %>% 
    mutate(top_genes_x = rank(-{{x_column}}, ties.method ="first") < top_genes,
           jittered = runif(nrow(.), min = -0.4, max = 0.4),
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_x, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_5 <- sub_4 %>% filter(!in_gene_list)
  sub_6 <- sub_4 %>% filter(in_gene_list)
  
  p2 <- ggplot(sub_4, aes(x = jittered, y = {{x_column}}, label = {{label_names}})) +
    geom_point(aes(fill = target, shape = factor(top_genes_x, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    scale_shape_manual(values = c(`FALSE`= 21, `TRUE` = 23), drop = FALSE) +
    
    geom_label_repel(data = filter(sub_6, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = Inf,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_5, top_genes_x), color = "firebrick", show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    labs(title = paste0("Only in ", x_char), fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(-0.4, 0.4)) +
    theme_bw()+
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()    
    )
  
  sub_7 <- df %>% 
    filter(is.na({{x_column}}) & !is.na({{y_column}})) %>% 
    mutate(`:=`({{x_column}}, "NA")) %>% 
    mutate(top_genes_y = rank(-{{y_column}}, ties.method ="first") < top_genes,
           jittered = runif(nrow(.), min = -0.4, max = 0.4),
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_y, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_8 <- sub_7 %>% filter(!in_gene_list)
  sub_9 <- sub_7 %>% filter(in_gene_list)
  
  p3 <- ggplot(sub_7, aes(x = jittered, y = {{y_column}}, label = {{label_names}})) +
    geom_point(aes(fill = target, shape = factor(top_genes_y, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_shape_manual(values = c(`FALSE`= 21, `TRUE` = 23), drop = FALSE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    geom_label_repel(data = filter(sub_9, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = Inf,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_8, top_genes_y), color = "firebrick", show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    
    labs(title = paste0("Only in ", y_char), fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(-0.4, 0.4)) +
    theme_bw()+
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  panel_list <- list()
  panel_list[["common"]] <- p1
  panel_list[[x_char]] <- p2
  panel_list[[y_char]] <- p3
  
  l <- list()
  l[["plot"]] <- (p1 + p2 + p3) + patchwork::plot_layout(guides = "collect", widths = c(0.5, 0.3, 0.3))
  l[["common"]] <- sub_1
  l[[x_char]] <- sub_4
  l[[y_char]] <- sub_7
  l[["top_genes"]] <- sub_2
  l[["panels"]] <- panel_list
  
  return(l)
}

#' advanced_test
#'
#' @param se SummarizedExperiment with groups as condition column of colData
#' @param design_formula formula object specifying the design of the linear model. Best left unchanged
#' @param advanced_contrast character vector specifying the contrasts to be tested. e.g. c("condA_vs_condB" = "condA - condB", "condB_vs_sctrl" = "condB - (condA + condC)/2")
#' @param fdr.type character string specifying the method to use for multiple testing correction. Default is "Strimmer's qvalue(t)".
#'
#' @return SummarizedExperiment with the results of the test added to the rowData
#' @export
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom tidyr gather unite pivot_wider
#' @importFrom tibble rownames_to_column
#' @importFrom qvalue qvalue
#' @importFrom purrr map_df
#' @importFrom stats p.adjust
#' @importFrom fdrtool fdrtool
#' @importFrom SummarizedExperiment rowData colData assay 
#' @importFrom stats model.matrix 
#' @import dplyr

advanced_test <- function (se, design_formula = formula(~0 + condition), advanced_contrast = NULL, 
                       fdr.type = c("Strimmer's qvalue(t)", "Strimmer's qvalue(p)", 
                                    "BH", "Storey's qvalue")) 
{
  
  col_data <- colData(se)
  raw <- assay(se)
  
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>% 
    as.character() %>% .[-1]
  
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  
  cntrst = advanced_contrast
  
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), 
                                      collapse = ", "))
  
  fit <- lmFit(raw, design = design)
  
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if (any(is.na(raw))) {
    for (i in cntrst) {
      
      covariates <- strsplit(i, " - |\\(|\\)| \\+ ") %>% unlist
      covariates <- covariates[!covariates == "" & !grepl("/", covariates)]
      
      single_contrast <- makeContrasts(contrasts = i, levels = covariates)
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 
                                                                             1]
      
    }
  }
  eB_fit <- eBayes(contrast_fit, trend = FALSE)
  retrieve_fun <- function(comp, fit = eB_fit, fdr.type) {
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf, 
                    confint = TRUE)
    res <- res[!is.na(res$t), ]
    if (fdr.type == "Strimmer's qvalue(t)") {
      fdr_res <- fdrtool(res$t, plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if (fdr.type == "Strimmer's qvalue(p)") {
      fdr_res <- fdrtool(res$P.Value, statistic = "pvalue", 
                         plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if (fdr.type == "BH") {
      padj <- p.adjust(res$P.Value, method = "BH")
      res$qval = padj
    }
    if (fdr.type == "Storey's qvalue") {
      qval_res = qvalue::qvalue(res$P.Value)
      res$qval = qval_res$qvalues
    }
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  message(fdr.type)
  limma_res <- purrr::map_df(cntrst, retrieve_fun, fdr.type = fdr.type)
  
  limma_res$comparison = names(cntrst)[match(limma_res$comparison, 
                                             cntrst)]
  
  table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, 
                                       CI.R, t, P.Value, qval, comparison) %>% tidyr::gather(variable, 
                                                                                             value, -c(rowname, comparison)) %>% mutate(variable = recode(variable, 
                                                                                                                                                          logFC = "diff", t = "t.stastic", P.Value = "p.val", qval = "p.adj")) %>% 
    tidyr::unite(temp, comparison, variable) %>% tidyr::pivot_wider(names_from = temp, values_from = value)
  
  rowData(se) <- left_join(as.data.frame(rowData(se)), table, 
                       by = c("name" = "rowname"))
  return(se)
}


#' merge_se
#'
#' @param se A list of SummarizedExperiment objects
#' @param keep_all Logical indicating whether to keep all rows or only those present in all objects. Default is FALSE.
#'
#' @return A SummarizedExperiment object with the merged data.
#' @export
#'
#' @importFrom DEP2 make_se
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom tibble column_to_rownames
#' @import dplyr
merge_se <- function(se = list(), keep_all = FALSE){
  
  old_names <- names(se)
  new_names <- paste0("nof_stripped_sequences_identified_experiment_wide_", names(se))
  mw_names <- paste0("moleculaweight_", names(se))
  
  se <- lapply(seq_along(se), function(x){
    se_rd <- rowData(se[[x]]) %>% as.data.frame() 
    
    colnames(se_rd) <- gsub("nof_stripped_sequences_identified_experiment_wide", new_names[[x]], colnames(se_rd))
    colnames(se_rd) <- gsub("moleculaweight", mw_names[[x]], colnames(se_rd))
    
    se_cd <- colData(se[[x]]) %>% as.data.frame() %>% mutate(set = names(se)[[x]])
    
    se_as <- assay(se[[x]]) %>% as.data.frame()
    
    se_sa <- gsub("(.*)_(quantity|log2quantity)", "\\1", colnames(dplyr::select(se_rd, ends_with("_quantity"), ends_with("_log2quantity"))))
    
    se_rd <- se_rd %>% dplyr::select(name, matches(paste0(se_sa, collapse = "|")), contains("_vs_"),
                                     contains("nof_stripped_sequences_identified_experiment_wide"),
                                     contains("moleculaweight"))
    
    list(cd = se_cd, rd = se_rd, as = se_as, sa = se_sa)
  })
  
  names(se) <- old_names
  
  #Check for duplicated IDs and rename if necessary
  temp <- c()
  for(se_ in se){
    temp <- c(se_$cd$ID, temp)
  }
  
  if(any(duplicated(temp))){
    warning("Duplicated IDs! Renamed to avoid conflicts")
    
    for(i in seq_along(se)){
      se[[i]]$cd$ID <- paste0(names(se)[[i]], "_",se[[i]]$cd$ID)
      se[[i]]$cd$condition <- paste0(names(se)[[i]], "_",se[[i]]$cd$condition)
      colnames(se[[i]]$as) <- paste0(names(se)[[i]], "_", colnames(se[[i]]$as))
      
      
      colnames(se[[i]]$rd) <- gsub("(.*)_vs_(.*)", paste0(names(se)[[i]], "_", "\\1_vs_", names(se)[[i]], "_","\\2"), colnames(se[[i]]$rd))
      
    }
  }
  
  # Merge
  for(i in seq_along(se)){
    if(i == 1){
      rd <- se[[i]]$rd %>% dplyr::select(-name)
      cd <- se[[i]]$cd
      as <- se[[i]]$as
    }else{
      rd <- merge(rd, dplyr::select(se[[i]]$rd, -name), by = "row.names", all = keep_all, sort = FALSE) %>% tibble::column_to_rownames("Row.names") %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      cd <- rbind(cd, se[[i]]$cd) %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      as <- merge(as, se[[i]]$as, by = "row.names", all = TRUE, sort = FALSE, ) %>% tibble::column_to_rownames("Row.names") %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      colnames(as) <- cd$ID  
      rownames(cd) <- cd$ID
    }
  }
  
  
  se <- SummarizedExperiment::SummarizedExperiment(assay = as.matrix(as), rowData = rd, colData = cd)
  rowData(se)$name <- rownames(se)
  
  return(se)
}


#' Plot Antigen Intensities and Differences
#'
#' This function creates patient-wise plots of protein intensity differences for a given contrast between conditions in a `SummarizedExperiment` object. Optionally, intensities can be z-scaled. For z-score calculation, additional conditions can be included. The function also generates an overlap plot that depicts in how many patients the given proteins were above the thresholds.
#'
#' @param se A `SummarizedExperiment` object containing the experimental data.
#' @param contrast A string specifying the comparison between conditions in the format "test_vs_control".
#' @param additional_sets A character vector specifying additional condition sets to include in the analysis. Defaults to "none", meaning no additional sets are included. Set to "all" to include all conditions.
#' @param scale A logical value indicating whether to scale the data by z-transformation. Defaults to `FALSE`.
#' @param min_diff A numeric value specifying the minimum difference required for data points to be considered. Defaults to `1`.
#' @param min_intensity A numeric value specifying the minimum intensity for data points to be considered. Defaults to `10`. If scaled is set to TRUE, this needs to be changed to a lower value. 
#' @param max.overlaps A numeric value indicating the maximum number of label overlaps allowed in the plot. Defaults to `40`.
#' @param targets A data frame containing `name` and `target` columns, specifying target proteins or features of interest. Defaults to an empty data frame.
#'
#' @return A combined plot displaying individual and overlapping data points for the specified contrast.
#' @details The function processes the `SummarizedExperiment` data to calculate the mean of control samples and the differences between test samples and control means. It generates individual plots for test samples and an overlap plot for the specified conditions.
#' 
#' If the mean for the control condition is not already present in the data, it is calculated using `sev::add_stats()`.
#'
#' @examples
#' # Assuming `se` is a valid SummarizedExperiment object and `contrast` is specified:
#' plot_antigen(se, contrast = "treatment_vs_control")
#'
#' @importFrom DEP2 get_df_wide
#' @importFrom dplyr select mutate across full_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @export
plot_antigen <- function(se, contrast, additional_sets = "none", scale = FALSE, min_diff = 1, min_intensity = 10, max.overlaps = 40, targets = data.frame(name = character(0), target = character(0))){
  
  ctrl <- gsub(".*_vs_(.*)$","\\1", contrast)
  ctrl_mean <- paste0("mean_", ctrl)
  
  test_condition <- gsub("(.*)_vs_.*","\\1", contrast)
  
  if(!paste0("mean_", ctrl) %in% colnames(rowData(se))){
    se <- sev::add_stats(se, type = "mean")
    message(paste0("No mean for ctrl condition found. Mean for\"", ctrl, "\" was calculated via sev::add_stats(). This value is not retained in the se object."))
  }
  
  if(additional_sets == "all"){
    additional_sets <- se$condition
  } else if(additional_sets == "none"){
    additional_sets <- NULL
  }
  
  sets <- strsplit(contrast, "_vs_") %>% unlist() %>% c(., additional_sets)  %>% unique()
  
  df_diff <- se %>% DEP2::get_df_wide() 
  
  if(scale){
    temp_df <-  df_diff %>% select(matches(paste0(sets, "_[0-9]*$"))) %>% t() %>% scale() %>% t()%>% cbind(df_diff %>% select(name, contains("_diff"))) %>% mutate(across(everything() , ~ifelse(is.na(.), 0, .)))
  } else {
    temp_df <- df_diff %>% select(matches(paste0(sets, "_[0-9]*$"))) %>% cbind(df_diff %>% select(name, paste0(contrast, "_diff"))) %>% mutate(across(everything() , ~ifelse(is.na(.), 0, .)))
  }
  
  temp_df <- df_diff %>% select(name, matches(paste0(sets, "_[0-9]*$")), !!sym(ctrl_mean)) %>% 
    tidyr::pivot_longer(-c(name, !!sym(ctrl_mean)), names_to = "sample", values_to = "intensity") %>% 
    mutate(diff = intensity - !!sym(ctrl_mean)) %>% 
    tidyr::pivot_wider(names_from = "sample", values_from = c(intensity, diff)) %>% 
    select(name, starts_with("diff")) %>% 
    dplyr::full_join(temp_df, by = c("name" = "name"))
  
  
  test_samples <- grep(paste0("^", test_condition, "_[0-9]*$"), colnames(temp_df), value = TRUE)
  
  p_individuals <- lapply(test_samples, function(x) plot_indiviuals(temp_df, paste0("diff_", x) ,x,  min_diff, min_intensity, max.overlaps = max.overlaps, ctrl = ctrl, scaled = scale)) 
  
  p_overlaps <- plot_overlap(temp_df, samples = grep(paste0("^", sets, "_[0-9]*$", collapse = "|"), colnames(temp_df), value = TRUE), 
                             ctrl, targets = targets, scaled = scale)
  
  p <- patchwork::wrap_plots(p_individuals, ncol = 2) / p_overlaps + 
    patchwork::plot_annotation(title = contrast) + patchwork::plot_layout(heights = c(1, 0.8))
}


#' Plot Individual Antigen Differences
#'
#' This function generates a scatter plot visualizing the difference between the intensity of a sample and the control mean for individual data points in a data frame. Points are colored based on whether they meet specified threshold criteria for both axes.
#'
#' @param df A data frame containing the data to be plotted.
#' @param x A string representing the column name for the x-axis variable (the difference between sample intensity and control mean).
#' @param y A string representing the column name for the y-axis variable (sample intensity).
#' @param cut_x A numeric value specifying the threshold for the x-axis variable, above which points are highlighted.
#' @param cut_y A numeric value specifying the threshold for the y-axis variable, above which points are highlighted.
#' @param max.overlaps A numeric value indicating the maximum number of overlapping labels for the text. Defaults to `10`.
#' @param ctrl A string specifying the control condition name, used for labeling the plot. Defaults to `"ctrl"`.
#' @param scaled A logical value indicating whether the y-axis should be labeled as scaled. Defaults to `FALSE`.
#'
#' @return A `ggplot2` object representing the scatter plot with points color-coded based on whether they meet the threshold criteria. Points above the threshold are labeled with their `name`.
#' @details The function uses `ggplot2` to create a scatter plot where points are highlighted and labeled if they meet or exceed the specified `cut_x` and `cut_y` values. Labels for points are added using `geom_text_repel` to reduce overlap.
#'
#' @examples
#' # Assuming `df` is a valid data frame with relevant columns:
#' plot_indiviuals(df, x = "diff_sample1", y = "sample1", cut_x = 1, cut_y = 10)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs
#' @importFrom dplyr filter
#' @importFrom ggrepel geom_text_repel
plot_indiviuals <- function(df, x, y, cut_x, cut_y, max.overlaps = 10, ctrl = "ctrl", scaled = FALSE){
  
  df <- df %>% filter(!is.na(!!sym(x)) & !is.na(!!sym(y)))
  p <- df %>% ggplot(aes(x = !!sym(x), y = !!sym(y), color = ifelse(!!sym(x) >= cut_x & !!sym(y) > cut_y, "TRUE", "FALSE"))) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "darkgrey")) +
    labs(color = "Above Thresholds", y = paste0(ifelse(scaled, "Scaled ", ""), "Intensity (", y, ")"), x = paste0(y, " - mean(", ctrl, ")")) +
    geom_text_repel(data = filter(df, !!sym(x) >= cut_x & !!sym(y) >= cut_y), aes(label = name), max.overlaps = max.overlaps) +
    my_theme()
  return(p)
}


#' Plot Overlap of Antigen Intensity Differences
#'
#' This function creates a summary plot that visualizes the overlap of proteins based on the mean difference and mean intensity across multiple samples. It highlights the distribution of data points that meet certain criteria and optionally labels specific targets.
#'
#' @param df A data frame containing the processed data, with columns representing sample intensities and differences. Defaults to `temp_df`.
#' @param samples A character vector specifying which sample columns to include. If `NULL`, samples matching the condition pattern are selected automatically.
#' @param condition A string specifying the condition name used for filtering and plotting. Defaults to `"MGN_ag_neg"`.
#' @param ctrl A string representing the control condition name, used for labeling the plot. Defaults to `"ctrl"`.
#' @param targets A data frame with `name` and `target` columns for highlighting specific targets. Defaults to an empty data frame.
#' @param scaled A logical value indicating whether the y-axis should be labeled as scaled. Defaults to `FALSE`.
#'
#' @return A `ggplot2` object representing a scatter plot that shows antigen overlap based on mean differences and intensities, with options for customized labeling and color coding.
#' @details The function processes input data to compute summaries for each antigen and plots them based on the number of samples in which they are detected. Points can be color-coded by target if provided, and labeled to show individual antigens and patients.
#'
#' @examples
#' # Assuming `df` is a valid data frame with appropriate columns:
#' plot_overlap(df, condition = "condition1", ctrl = "control1", targets = targets_df)
#'
#' @importFrom dplyr rename_with select filter group_by summarize left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_discrete labs facet_wrap scale_fill_viridis_d
#' @importFrom ggrepel geom_text_repel
plot_overlap  <- function(df = temp_df, samples = NULL, condition = "MGN_ag_neg", ctrl = "ctrl", targets = data.frame(name = character(0), target = character(0)), scaled = FALSE){ 
  if(is.null(samples)){
    samples <- grep(paste0("^", condition, "_[1-9]*$"), colnames(df), value = TRUE)
  }
  
  
  df_sum <- df %>% dplyr::rename_with(.cols = all_of(samples), ~ gsub("(.*)", "intensity_\\1", .x)) %>% 
    select(name, grep(paste0(samples, collapse = "|"), colnames(.), value = TRUE)) %>%
    tidyr::pivot_longer(-c(name), names_to = c("type", "patient"), values_to = "intensity", names_pattern = "(.*?)_(.*)") %>% 
    tidyr::pivot_wider(names_from = "type", values_from = "intensity") %>%
    filter(diff >= 0 & intensity >= 0.5) %>% 
    group_by(name) %>% 
    summarize(n = factor(n(), levels = 1:length(test_samples)), 
              diff = mean(diff, na.rm = TRUE), 
              patient = paste0(gsub(".*_(.*)", "\\1",patient), collapse = ", "), 
              mean_intensity = mean(intensity, na.rm = TRUE)) %>% filter(n != 0)
  
  sample = paste0(".*", condition, "_[1-9]*$")
  
  df_sum <- df_sum %>% dplyr::left_join(targets, by = "name") 
  
  if(targets %>% nrow() > 0){
    p_summary <- df_sum %>% ggplot(aes(x = diff, y = mean_intensity, fill = target))+
      geom_point(shape = 21, size = 3, alpha = 0.8)+
      scale_fill_discrete(na.value = "grey90") +
      ggrepel::geom_text_repel(data = filter(df_sum, !n ==1 & !n == 2), aes(label = name), min.segment.length = 0) +
      ggrepel::geom_text_repel(data = filter(df_sum, n == 1 | n == 2), aes(label = paste0(name, " | ", patient)), min.segment.length = 0) +
      facet_wrap(.~n, scales = "free") +
      labs(y = paste0("Mean(", ifelse(scaled, "Scaled ", ""), "Intensity)"), x = "Mean(Diff)", title = paste0(condition, "vs", ctrl)) +
      my_theme() 
  } else {
    p_summary <- df_sum %>% ggplot(aes(x = diff, y = mean_intensity, fill = n))+
      geom_point(shape = 21, size = 3, alpha = 0.8)+
      scale_fill_viridis_d(option = "inferno", end = 0.8) +
      ggrepel::geom_text_repel(data = filter(df_sum, !n ==1 & !n == 2), aes(label = name), min.segment.length = 0) +
      ggrepel::geom_text_repel(data = filter(df_sum, n == 1 | n == 2), aes(label = paste0(name, " | ", patient)), min.segment.length = 0) +
      facet_wrap(.~n, scales = "free") +
      labs(y = paste0("Mean(", ifelse(scaled, "Scaled ", ""), "Intensity)"), x = "Mean(Diff)", title = paste0(condition, "vs", ctrl)) +
      my_theme()   
  }
  
  return(p_summary)
}


#' Plot Antigen Missing Data Visualization
#'
#' This function creates a detailed plot to visualize missing data patterns, highlighting the percentage of missing data in test and control conditions. It also allows for the labeling of specific targets.
#'
#' @param se A `SummarizedExperiment` object containing the data to be analyzed.
#' @param test_condition A string specifying the test condition to analyze. Defaults to `"lcm_igan"`.
#' @param ctrl_condition A string specifying the control condition to compare against. Defaults to `"lcm_ctrl"`.
#' @param quantile A numeric value indicating the quantile cutoff for low intensity in the control condition. Defaults to `0.1`.
#' @param perc_low_ctrl A numeric threshold for filtering antigens based on the percentage of low-intensity measurements in the control condition. Defaults to `80`.
#' @param targets A data frame with `name` and `target` columns to highlight specific targets. Defaults to an empty data frame.
#'
#' @return A `ggplot2` object showing a plot with points representing antigens, color-coded based on their missing data percentages, and optionally labeled.
#' @details The function calculates missing and measured data percentages, applies quantile-based cutoffs for low intensity, and generates a scatter plot with facets showing different categories of missing data in the test and control conditions. The visualization helps identify patterns of data availability.
#'
#' @examples
#' # Assuming `se` is a `SummarizedExperiment` object with relevant data:
#' plot_antigen_missing(se, test_condition = "test_cond", ctrl_condition = "control_cond", quantile = 0.05, perc_low_ctrl = 70)
#'
#' @importFrom dplyr filter mutate group_by ungroup summarize_all summarize left_join
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_brewer facet_grid labs theme element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_rev
#' @export
plot_antigen_missing <- function(se, test_condition = "lcm_igan", ctrl_condition = "lcm_ctrl", 
                                 quantile = 0.1, perc_low_ctrl = 80, targets = data.frame(name = character(0), target = character(0))){
  
  df_long <- get_df_long(se[, se$condition %in% c(test_condition, ctrl_condition)]) %>% select(name, condition, intensity, label, replicate)
  n_ctrl <- sum(se$condition == ctrl_condition)
  n_test <- sum(se$condition == test_condition)
  
  df_pat <- df_long %>% 
    filter(!is.na(intensity), condition == test_condition) %>% 
    mutate(condition = "test_condition") %>%
    group_by(name) %>% 
    mutate(n = n()) %>% 
    ungroup() %>%  
    mutate(perc_miss = ((n_test-n)/n_test)*100, 
           perc_meas = ((n/n_test)*100)) %>% 
    group_by(condition, name) %>% 
    mutate(mean = mean(intensity, na.rm = TRUE)) %>% 
    ungroup() %>%   
    tidyr::pivot_wider(names_from = condition, 
                       values_from = c(intensity, perc_miss, perc_meas, mean, n))
  
  df_ctrl <- df_long %>% 
    filter(condition == ctrl_condition) %>% 
    mutate(condition = "ctrl_condition") %>%
    group_by(replicate, condition) %>% 
    mutate(cut_off = quantile(intensity, probs = quantile, na.rm = TRUE)) %>% 
    ungroup() %>%
    group_by(name, condition) %>% 
    summarize(n = sum(!is.na(intensity)), 
              n_low = sum(intensity < cut_off |is.na(intensity)), 
              n_miss = sum(is.na(intensity)),
              perc_miss = ((n_ctrl-sum(!is.na(intensity)))/n_ctrl)*100, 
              perc_low =  ((sum(intensity < cut_off |is.na(intensity)))/n_ctrl)*100, 
              mean = mean(intensity, na.rm = TRUE)) %>% 
    ungroup() %>% 
    tidyr::pivot_wider(names_from = condition, 
                       values_from = c(perc_miss, perc_low, mean, n_low, n))
  
  df_final <- inner_join(df_pat, df_ctrl, by = "name") %>% 
    #filter(perc_low_ctrl_condition >= perc_low_ctrl) %>% 
    group_by(name) %>% 
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(perc_miss_ctrl_cut = cut(perc_miss_ctrl_condition, 
                                    breaks = c(0, 70, 80, 99, 100), 
                                    labels = c("0-70 %","70-80 %", "80-99 %", "100 %"), 
                                    include.lowest = TRUE)) %>% 
    mutate(x_coord = runif(nrow(.))) %>% 
    mutate(perc_low_ctrl_cut = cut(perc_low_ctrl_condition, 
                                   breaks = c(0, 80, 99, 100), 
                                   labels = c("0-80 %", "80-99 %", "100 %"), 
                                   include.lowest = TRUE, include.highest = TRUE)) %>%
    mutate(perc_miss_igan_cut = cut(perc_miss_test_condition, 
                                    breaks = c(0, 10, 30, 50, 100), 
                                    labels = c("0-10 %", "10-30 %", "30-50 %", "50-100 %"), 
                                    include.lowest = TRUE))%>%
    mutate(perc_meas_igan_cut = cut(perc_meas_test_condition, 
                                    breaks = c(0, 30, 50, 70, 90, 100), 
                                    labels = c("0-30 %", "30-50 %", "50-70 %", "70-90 %", "90-100 %"), 
                                    include.lowest = TRUE)) %>%
    mutate(perc_miss_ctrl_cut = forcats::fct_rev(factor(perc_miss_ctrl_cut, levels = c("0-70 %","70-80 %", "80-99 %", "100 %"))),
           perc_low_ctrl_cut = forcats::fct_rev(factor(perc_low_ctrl_cut, levels = c("0-80 %", "80-99 %", "100 %"))),
           perc_miss_igan_cut = factor(perc_miss_igan_cut, levels = c("0-10 %", "10-30 %", "30-50 %", "50-100 %")),
           perc_meas_igan_cut = factor(perc_meas_igan_cut, levels = c("0-30 %", "30-50 %", "50-70 %", "70-90 %", "90-100 %"))) %>% 
    dplyr::left_join(targets, by = "name") 
  
  set.seed(1)
  
  if(nrow(targets) == 0){
    p_final <- df_final  %>% 
      #filter(perc_meas_igan > 30) %>%
      ggplot(aes(x = x_coord, y = mean_test_condition, fill = perc_miss_ctrl_cut)) +
      geom_point(shape = 21, size = 3)+
      scale_fill_brewer(palette = "OrRd", direction = -1) +
      ggrepel::geom_text_repel(aes(label = name), min.segment.length = 0, max.overlaps = 10) +
      facet_grid(perc_low_ctrl_cut ~ perc_meas_igan_cut, scales = "free") +
      labs(fill = paste0("Percent missing in ", ctrl_condition), y = paste0("Mean(Intensity) in ", test_condition), x = "") +
      my_theme() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  } else{
    p_final <- df_final  %>% 
      #filter(perc_meas_igan > 30) %>%
      ggplot(aes(x = x_coord, y = mean_test_condition, fill = target)) +
      geom_point(aes(shape = perc_miss_ctrl_cut), size = 3)+
      scale_fill_brewer(palette = "OrRd", direction = -1) +
      scale_shape_manual(values = c(21, 22, 24, 23, 25)) +
      ggrepel::geom_text_repel(aes(label = name), min.segment.length = 0, max.overlaps = 10) +
      facet_grid(perc_low_ctrl_cut ~ perc_meas_igan_cut, scales = "free") +
      labs(fill = paste0("Percent missing in ", ctrl_condition), y = paste0("Mean(Intensity) in ", test_condition), x = "") +
      my_theme() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  return(p_final)
}


############################################################################
############################################################################
############################################################################
#' Retrieve Gene Data from a Specified KEGG Pathway
#'
#' This function uses the R-package KEGGREST to retrieve data for a specified KEGG pathway and then returns a tidy data frame
#' that includes: "kegg_gene_id", "description", "gene_symbol", "kegg_ontology_number", "enzyme_commission_number", "kegg_path_number" and "kegg_path_name" 
#'
#' @param keggpath_id A character string specifying the KEGG pathway ID.
#' Default is "hsa00190" for oxidative phosphorylation in Homo sapiens.
#'
#' @return A data frame containing Gene IDs, Descriptions, Gene symbols, KEGG ontology numbers,
#' Enzyme Commission numbers, the pathway ID queried, and the pathway name.
#'
#' @examples
#' genes_from_kegg("hsa00190")
#' oxphos <- genes_from_kegg(keggpath_id = "hsa00190")  # default, or specify another pathway ID https://www.kegg.jp/entry/map00190
#' cellcycle <- genes_from_kegg(keggpath_id = "hsa04110") #https://www.genome.jp/pathway/map04110
#' mtor <- genes_from_kegg(keggpath_id = "hsa04150") #https://www.kegg.jp/pathway/map=map04150
#' ampk <- genes_from_kegg(keggpath_id = "hsa04152") #https://www.kegg.jp/pathway/map=map04152
#' print(cellcycle)
#' 
#' 
#' # Or mutliple keggpaths, assuming you have a vector of 4 keggpath_ids
#' keggpath_ids <- c("hsa00190", "hsa04110", "hsa04150", "hsa04152")
#' 
#' # Initialize an empty list to store the data frames
#' list_of_dataframes <- list()
#' 
#' # Loop over each keggpath_id
#' for (keggpath_id in keggpath_ids) {
#'   # Call the genes_from_kegg function for each keggpath_id
#'   df <- genes_from_kegg(keggpath_id)
#' 
#'   # Name the data frame according to the keggpath_id and add it to the list
#'   list_of_dataframes[[keggpath_id]] <- df
#' }
#' 
#' # Now, list_of_dataframes is a list of data frames, each named after a keggpath_id
#' # Use bind_rows to combine all data frames into a single data frame
#' combined_df <- bind_rows(list_of_dataframes, .id = "keggpath_id")
#' @importFrom dplyr mutate
#' @export
genes_from_kegg <- function(keggpath_id = "hsa00190") {
  # Helper to manage library loading
  ensurePackage <- function(pkg, bioconductor = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioconductor) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  
  # Ensure required libraries are loaded
  ensurePackage("KEGGREST", bioconductor = TRUE)
  ensurePackage("tidyverse")
  ensurePackage("magrittr")
  
  
  # Fetching KEGG pathway data only once
  pathway_data <- KEGGREST::keggGet(keggpath_id)[[1]]
  
  # Get human genes
  genes <- pathway_data$GENE
  
  # Separate IDs and Descriptions
  kegg_gene_id <- genes[seq(1, length(genes), by = 2)]
  description <- genes[seq(2, length(genes), by = 2)]
  
  # Combine into a data frame
  gene_data <- data.frame(kegg_gene_id, description, stringsAsFactors = FALSE)
  
  # Extract additional details
  gene_data <- gene_data %>%
    dplyr::mutate(
      gene_symbol = sapply(description, function(x) {
        strsplit(x, ";")[[1]][1]
      }),
      kegg_ontology_number = sapply(description, function(x) {
        matches <- regmatches(x, regexpr("(?<=\\[KO:)[^\\]]+", x, perl = TRUE))
        if (length(matches) > 0) matches[1] else NA
      }),
      enzyme_commission_number = sapply(description, function(x) {
        matches <- regmatches(x, regexpr("(?<=\\[EC:)[^\\]]+", x, perl = TRUE))
        if (length(matches) > 0) matches[1] else NA
      }),
      kegg_path_number = keggpath_id, # Add the new column 'kegg'
      kegg_path_name = pathway_data$NAME # Directly add pathway name from fetched data
    )
  
  return(gene_data)
}


################################################
################################################
################################################
#' Get Data in Wide Format with Specified/Available Assays Information from summarized experiment
#'
#' This function converts a SummarizedExperiment object into a wide-format data frame.
#' It appends all or selected assays data to the regular DEP2 table created with DEP2::get_df_wide()
#'
#' @param se A SummarizedExperiment object containing the assay(s).
#' @param assays Optional character vector specifying the assay(s) to include.
#'        If NULL or not provided, all available assays from the SummarizedExperiment
#'        are used.
#'
#' @return A data frame in wide format with selected/all assays data joined via dplyr::full_join().
#' 
#' @examples
#' # Assuming se_diff is an existing SummarizedExperiment object
#' df_out <- get_df_wide_append_assay(se_diff, assays = c("lfq_raw", "imputed_DEP2"))
#' df_out_all <- get_df_wide_append_assay(se_diff)  # Uses all available assays
#'
#' @export

get_df_wide_append_assay <- function(se, assays = NULL) {
  # If assays is NULL, retrieve all assay names from SummarizedExperiment
  if (is.null(assays)) {
    assays <- names(assays(se))
    print("Using all available assays.")
  } else {
    print("Using specified assays.")
  }
  
  # Print the assays being processed
  message("Assay names: ", paste(assays, collapse = ", "))
  
  if (length(assays) < 1) {
    stop("No valid assays provided or available in the SummarizedExperiment object.")
  }
  
  # Proceed with the original dataframe and merging process using selected assays
  # Initial df using DEP2::get_df_wide(), and some sorting
  df <- get_df_wide(se) %>% 
    dplyr::select(name,gene_names,protein_ids, protein_descriptions = description, ID, orig_prot_ids = protein, contains('_vs_'),1:ncol(.))
  
  # Initialize an empty tibble for the eventual joins
  result <- tibble(name = df$name)
  
  # Loop through each specified assay, process data, and join
  for (assay_name in assays) {
    if (!assay_name %in% names(assays(se))) {
      warning(paste("Assay", assay_name, "not found in SummarizedExperiment; it will be skipped."))
      next
    }
    
    assay_data <- as.data.frame(assay(se, assay_name)) %>%
      rownames_to_column("name") %>%
      as_tibble() %>%
      rename_with(~paste0(.x, "_", assay_name), -name)
    
    # Join with the main result DataFrame
    result <- result %>%
      full_join(assay_data, by = "name")
  }
  
  # Join all together with initial df
  result2 <- df %>%
    full_join(result, by = "name")
  
  return(result2)
}




############################################################################
############################################################################
############################################################################
#' Evenly spread/randomize samples into manageable blocks for batch-processing
#'
#' The function considers batches, conditions and bioreplicates, and evenly spreads/randomizes samples into manageable experimental blocks.
#' For example, you have an experiment with 200 samples (4 conditions,n=50 per condition).
#' But you cannot process all of these in a single block, because you can only process a maximum of 30 samples at one time.
#' Therefore you need to split your lab work into blocks of a certain size..
#' This function assists in the blocking, so that conditions (and any other prior batching or cohorts) are spread EVENLY across your blocks.
#' This ensures that 
#' i) the effects you are measuring between your conditions are not clouded/overshadowed by blocking-effects and 
#' ii) the effects of the blocking are not accidentally interpreted by you as 'real' effects caused by your conditions
#'
#' @param df Input Data frame that contains your blocked.
#' @param batch_col Column name representing batches (in case the samples come from different batches, e.g. cohorts or prior batching during processing), defaults to "batch".
#' @param condition_col Column name representing your test conditions, defaults to a column called "condition" in your dataframe.
#' @param bio_replicate_col Column name representing biological replicates, defaults to "bio_replicate".
#' @param scramble Logical flag (TRUE or FALSE) to determine whether rows should be shuffled, defaults to TRUE.
#' @param seed Optional numeric seed for reproducibility, defaults to NULL (no seed). But the recommendation is to set the seed to make reproducible.
#' @return A data frame with rows rearranged into blocks. 
#' The main outcome is column called "block". 
#' The "block" column contains a number for your experimental blocks. Importantly, the "block_size" column contains the minimum of samples that you have to run together in one block.
#' If this number exceeds what you can process in one block, then you need to strike some compromises. Compromise options include i) leaving out batches all togehter or,
#' ii) in-silico reducing of any prior batching/cohorting, e.g. if you have 4 prior batches, turn them into 2 evenly mixed in-silico batches).
#' 
#' @examples
#' set.seed(42)
#' 
#' #' Define the conditions, batches, and maximum bioreplicates
#' conditions <- c("KO_ctrl", "WT_ctrl", "KO_treated", "WT_treated")
#' num_batches <- 3
#' max_bioreplicates <- 5
#' 
#' #' Initialize lists to store data
#' conditions_list <- c()
#' batches_list <- c()
#' bioreplicates_list <- c()
#' 
#' #' Generate random data for the sample data frame
#' for (condition in conditions) {
#'   for (batch in 1:num_batches) {
#'     num_bioreplicates <- sample(1:max_bioreplicates, 1)
#' 
#'     conditions_list <- c(conditions_list, rep(condition, num_bioreplicates))
#'     batches_list <- c(batches_list, rep(paste0("Batch", batch), num_bioreplicates))
#'     bioreplicates_list <- c(bioreplicates_list, seq(1, num_bioreplicates))
#'   }
#' }
#' 
#' #' Create a data frame
#' sample_df <- data.frame(
#'   condition = factor(conditions_list),
#'   batch = factor(batches_list),
#'   bio_replicate = bioreplicates_list
#' )
#' 
#' #' Display the sample data frame
#' print(sample_df)
#' 
#' 
#' 
#' result <- block_randomize(df = sample_df,
#'                                 batch_col = "batch",
#'                                 condition_col = "condition",
#'                                 bio_replicate_col = "bio_replicate",
#'                                 scramble = T,
#'                                 seed = 1000)
#'                                 
#' @export

block_randomize <- function(df, batch_col = "batch", condition_col = "condition", bio_replicate_col = "bio_replicate", scramble = TRUE, seed = NULL) {
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("janitor", quietly = TRUE)) {
    install.packages("janitor")
  }
  if (!requireNamespace("anticlust", quietly = TRUE)) {
    install.packages("anticlust")
  }
  
  library(dplyr)
  library(janitor)
  library(anticlust)
  
  # If batch_col is left empty, create a new 'batch' column with all 1s as values to indicate it is a single batch.
  # If batch_col is left empty, any pre-existing column called 'batch' is renamed to 'user_provided_batch'. (This is to ensure it is not overwritten or discarded)
  if (missing(batch_col) || batch_col == "") {
    batch_col <- "batch"
    if (exists(batch_col, where = df)) {
      df <- df %>% rename(NOTUSED_user_defined_batch_NOTUSED = !!sym(batch_col))
      batch_col <- "batch"
    }
    df <- df %>% mutate(!!batch_col := 1)
  }
  
  
  # Set seed for reproducibility if scramble is TRUE
  if (scramble) {
    if (!is.null(seed)) {
      set.seed(as.numeric(seed))
    } else {
      set.seed(12345)
    }
  }
  
  # Clean and the dataframe up and then scramble if desired
  tryCatch({
    df_cleaned <- df %>%
      dplyr::rename(batch = {{batch_col}},
                    condition = {{condition_col}},
                    bio_replicate = {{bio_replicate_col}}) %>%
      dplyr::filter(if_any(everything(), ~ !is.na(.))) %>%
      dplyr::filter(if_any(condition, ~ !is.na(.))) %>%
      dplyr::arrange(if (scramble) sample(nrow(.)) else NULL) %>% #scramble if desired
      dplyr::arrange(batch, condition) #order rows by batch, condition
    
    # Calculate parameters that are important for the anticlust package
    min_set <- (length(unique(df_cleaned$batch)) * length(unique(df_cleaned$condition))) #this is the number of minimum samples to process simultaneously. The code multiplies the number of unique batches by the number of unique conditions. This operation calculates the total number of unique combinations of batches and conditions
    number_sets <- ceiling(nrow(df_cleaned) / min_set) #this number indicates how many min sets you have. performs a calculation to determine how many minimum sets (min_set) are needed to cover the total number of rows in a data frame (df_cleaned)
    
    # Create a new dataframe with categorical sampling using the anticlust package
    df_new <- df_cleaned %>%
      mutate(block = anticlust::categorical_sampling(., number_sets)) %>% #if number of sets =1 it doesn't work. not sure how it handles not even numbers.
      dplyr::arrange(block) %>%
      mutate(new_sample_number = row_number()) %>%
      dplyr::group_by(block) %>% 
      mutate(block_size = n()) %>% 
      ungroup() %>% 
      dplyr::group_by(block, condition) %>% 
      mutate(block_condition_size = n()) %>% 
      ungroup() %>% 
      mutate(bio_reps_scrambled = scramble) %>% 
      dplyr::select(new_sample_number, block, block_size,dplyr::contains("batch"), condition, block_condition_size, bio_replicate, bio_reps_scrambled, 1:ncol(.) )
    
    # Add seed_used column only if scramble is TRUE
    if (scramble) {
      df_new <- df_new %>%
        mutate(seed_used = ifelse(!is.null(seed), as.character(seed), "12345"))
    }
    
    return(df_new)
  }, error = function(e) {
    cat("\n")
    cat("The following Error occurred in the block_randomiz function: \n")
    cat("\n")
    cat(conditionMessage(e), "\n")
    cat("\n")
    cat("Most likely, you have mis-spelled a column name in the arguments. \n")
    cat("\n")
    cat("In case of error around argument K needing to be bigger than 1, you will have to process all samples in one go. \n")
    cat("If you don't want to do that, leave out the 'batch_col'argument.\n")
    cat("However, doing so jeopordizes your entire future pipeline, because any prior batching of samples will not be accounted for. \n")
    
    return(NULL)
  })
}



############################################################################
############################################################################
############################################################################
#' Construct a PubMed Query String
#'
#' This function constructs a PubMed query string based on the provided
#' publication date range and combinations of two sets of keywords.
#'
#' @param date_range A character vector of length 2 indicating the start and end publication years.
#' @param keywords1 A character vector of the first set of keywords.
#' @param keywords2 A character vector of the second set of keywords.
#' @return A character string representing the PubMed query.
#' @examples
#' date_range <- c("2015", "2025")
#' keywords1 <- c("kidney", "renal")
#' keywords2 <- c("organoid", "spheroid", "tubuloid", "organ-on-a-chip", "microfluidic")
#' query <- pubmed_query(date_range, keywords1, keywords2)
#' cat(query)
#' @export
pubmed_query <- function(date_range, keywords1, keywords2) {
  # Validate input
  if (length(date_range) != 2) {
    stop("date_range should be a vector of two elements (start and end year).")
  }
  if (!all(sapply(date_range, is.character))) {
    stop("date_range should be a vector of characters.")
  }
  
  # Construct date part
  date_part <- paste0('("', date_range[1], '"[Publication Date] : "', date_range[2], '"[Publication Date])')
  
  # Construct keyword combinations
  keyword_combinations <- unlist(lapply(keywords1, function(k1) {
    sapply(keywords2, function(k2) {
      paste0('("', k1, ' ', k2, '")')
    })
  }))
  
  # Combine keyword combinations with OR
  keyword_part <- paste(keyword_combinations, collapse = " OR ")
  
  # Combine date part and keyword part
  query <- paste0(date_part, " AND (", keyword_part, ")")
  
  return(query)
}



############################################################################
############################################################################
############################################################################
#' Construct a Scopus Query String
#'
#' This function constructs a Scopus query string based on the provided
#' publication date range and combinations of two sets of keywords.
#'
#' @param date_range A character vector of length 2 indicating the start and end publication years.
#' @param keywords1 A character vector of the first set of keywords.
#' @param keywords2 A character vector of the second set of keywords.
#' @return A character string representing the Scopus query.
#' @examples
#' date_range <- c("2015", "2025")
#' keywords1 <- c("kidney", "renal")
#' keywords2 <- c("organoid", "spheroid", "tubuloid", "organ-on-chip", "microfluidic")
#' query <- scopus_query(date_range, keywords1, keywords2)
#' cat(query)
#' @export
scopus_query <- function(date_range, keywords1, keywords2) {
  # Validate input
  if (length(date_range) != 2) {
    stop("date_range should be a vector of two elements (start and end year).")
  }
  if (!all(sapply(date_range, is.character))) {
    stop("date_range should be a vector of characters.")
  }
  
  # Construct date part
  date_part <- paste0('PUBYEAR AFT ', date_range[1], ' AND PUBYEAR BEF ', date_range[2])
  
  # Construct keyword combinations
  keyword_combinations <- unlist(lapply(keywords1, function(k1) {
    sapply(keywords2, function(k2) {
      paste0('TITLE-ABS-KEY("', k1, ' ', k2, '")')
    })
  }))
  
  # Combine keyword combinations with OR
  keyword_part <- paste(keyword_combinations, collapse = " OR ")
  
  # Combine date part and keyword part
  query <- paste0(date_part, " AND (", keyword_part, ")")
  
  return(query)
}


############################################################################
############################################################################
############################################################################
#' Add statistical summaries to SummarizedExperiment
#'
#' This function calculates various statistical summaries for protein intensities
#' in a SummarizedExperiment object. The user can specify which statistics to compute
#' via the `type` argument, with the default being "all" to calculate all available statistics.
#' The computed statistics are merged back into the `rowData` of the SummarizedExperiment object.
#'
#' @param se A `SummarizedExperiment` object containing the protein data to be summarized.
#' @param type A character vector specifying which statistics to calculate. Options include:
#'   \describe{
#'     \item{"all"}{(default) calculates all available statistics (mean, sd, min, max, percentiles, IQR, counts).}
#'     \item{Custom set of statistics, e.g. `c("mean", "sd", "min")`} selects only the specified statistics. Available options are:
#'         \itemize{
#'           \item "mean": Mean intensity
#'           \item "sd": Standard deviation
#'           \item "min": Minimum intensity
#'           \item "max": Maximum intensity
#'           \item "p25": 25th percentile intensity
#'           \item "p50": 50th percentile (median) intensity
#'           \item "p75": 75th percentile intensity
#'           \item "iqr": Interquartile range
#'           \item "missing_count": Count of missing (NA) values
#'           \item "measured_count": Count of non-NA values
#'         }
#'   }
#'
#' @return The updated `SummarizedExperiment` object with the selected statistics added to the `rowData`.
#'
#' @examples
#' # Example usage with default (all statistics):
#' updated_se <- add_stats(se)
#'
#' # Example usage to return only mean, sd, and min statistics:
#' updated_se_selected <- add_stats(se, type = c("mean", "sd", "min"))
#'
#' # Example usage to calculate mean and missing_count only:
#' updated_se_mean_missing <- add_stats(se, type = c("mean", "missing_count"))
#'
#' @export
add_stats <- function(se, type = "all") {
  # Extract data in long format from the SummarizedExperiment object
  df_long <- get_df_long(se)
  
  # Convert all non-finite values (NaN, Inf, -Inf) to NA
  df_long <- df_long %>%
    mutate(intensity = ifelse(is.finite(intensity), intensity, NA))
  
  # Define the list of possible statistics and their respective calculation methods with the is.na check
  stats_to_calculate <- list(
    mean_intensity = ~ ifelse(all(is.na(.x)), NA, mean(.x, na.rm = TRUE)),
    sd_intensity = ~ ifelse(all(is.na(.x)), NA, sd(.x, na.rm = TRUE)),
    min_intensity = ~ ifelse(all(is.na(.x)), NA, min(.x, na.rm = TRUE)),
    max_intensity = ~ ifelse(all(is.na(.x)), NA, max(.x, na.rm = TRUE)),
    p25_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.25, na.rm = TRUE)),
    p50_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.50, na.rm = TRUE)),  # Median
    p75_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.75, na.rm = TRUE)),
    iqr_intensity = ~ ifelse(all(is.na(.x)), NA, IQR(.x, na.rm = TRUE)),  # Interquartile range
    missing_count = ~ sum(is.na(.x)),  # Count of NA values
    measured_count = ~ sum(!is.na(.x))  # Count of non-NA values
  )
  
  # If type is "all", calculate all stats; otherwise, filter the stats to calculate
  if (type == "all") {
    selected_stats <- names(stats_to_calculate)
  } else {
    selected_stats <- intersect(names(stats_to_calculate), paste0(type, "_intensity"))
    if ("missing_count" %in% type) selected_stats <- c(selected_stats, "missing_count")
    if ("measured_count" %in% type) selected_stats <- c(selected_stats, "measured_count")
  }
  
  # Summarize the selected statistics
  df_summaries <- df_long %>%
    group_by(name, condition) %>%
    summarize(across(
      .cols = intensity,
      .fns = stats_to_calculate[selected_stats],
      .names = "{.fn}"
    )) %>%
    ungroup()
  
  # Pivot the data to a wide format, appending suffixes to columns
  df_summaries_wide <- df_summaries %>%
    tidyr::pivot_wider(
      names_from = condition, 
      values_from = selected_stats,
      names_glue = "{condition}_{.value}"  # Automatically create column names with suffixes
    )
  
  # Arrange columns alphabetically (with 'name' as the first column)
  df_summaries_wide <- df_summaries_wide %>%
    dplyr::select(name, sort(colnames(df_summaries_wide)[-1]))  # Keep 'name' first
  
  # Extract the existing rowData from the SummarizedExperiment object
  row_data <- rowData(se)
  
  # Merge the new summary data into rowData
  row_data <- as.data.frame(row_data) %>%
    left_join(df_summaries_wide, by = "name")  # Merge by 'name' which is the protein identifier
  
  # Update the rowData of the SummarizedExperiment object
  rowData(se) <- as(row_data, "DataFrame")  # Convert it back to a DataFrame
  
  # Return the updated SummarizedExperiment object
  return(se)
}

############################################################################
############################################################################
############################################################################
#' Create SummarizedExperiment object from Spectronaut reports
#'
#' This function processes Spectronaut reports and returns a SummarizedExperiment object. 
#' It allows for optional transformations like collapsing extra underscores in the 
#' condition column and converting condition values to lowercase.
#'
#' @param candidates Character or dataframe. Path to the candidates file or a dataframe.
#' @param report Character or dataframe. Path to the Spectronaut report file or a dataframe.
#' @param contrasts Character vector. Optional. Specify contrasts for differential analysis. 
#' If NULL, default diff calculations will be used.
#' @param conditionSetup Character or dataframe. Path to the condition setup file or a dataframe.
#' @param quant_col Character. Column name indicating the quantity column. Default is "log2quantity".
#' @param suggest_contrasts Logical. If TRUE, the function will suggest possible contrasts based on conditions. Default is TRUE.
#' @param collapse_conditions Logical. If TRUE, extra underscores in the condition column are collapsed to the minimum number of underscores. Default is FALSE.
#' @param cond_tolower Logical. If TRUE, the condition column is converted to lowercase. Default is FALSE.
#'
#' @return A SummarizedExperiment object containing the processed data.
#'
#' @details 
#' The function processes Spectronaut reports, condition setups, and contrasts (if provided) 
#' to create a SummarizedExperiment object for downstream analysis. It allows for condition-specific 
#' transformations such as collapsing underscores and converting conditions to lowercase.
#' 
#' If `collapse_conditions` is TRUE, any extra underscores in the condition column will be collapsed 
#' to the minimum number of underscores found in the column.
#' 
#' If `cond_tolower` is TRUE, the values in the condition column will be converted to lowercase.
#'
#' @examples
#' # Example with custom arguments
#' se <- spectronaut_to_se(candidates = "path_to_candidates.txt", 
#'                         report = "path_to_report.txt", 
#'                         conditionSetup = "path_to_conditions.txt", 
#'                         quant_col = "log2quantity", 
#'                         collapse_conditions = TRUE, 
#'                         cond_tolower = TRUE)
#'                         
#' # Example with defaults
#' se <- spectronaut_to_se()
#'
#' @export
optimized_spectronaut_to_se <- function(candidates = NULL, report = NULL, contrasts = NULL, 
                              conditionSetup = NULL, quant_col = "log2quantity", 
                              suggest_contrasts = TRUE, collapse_conditions = FALSE, 
                              conditions_tolower = FALSE) {
  
  # Check if candidates data is provided and load if necessary
  if (!is.null(candidates)) {
    if (typeof(candidates) == "character") {
      df_candidates <- vroom::vroom(candidates, delim = "\t", col_names = TRUE, 
                                    guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
        dplyr::rename_all(tolower) %>% 
        janitor::clean_names() %>% 
        dplyr::rename_all(~ gsub("pg_", "", .)) %>% 
        dplyr::rename_all(~ gsub("^r_", "", .)) %>%  # Remove 'r_' only at the start of string
        dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA)))  # Replace NaN, -Inf, Inf with NA
    } else {
      df_candidates <- candidates
    }
  } else {
    df_candidates <- NULL
  }
  
  # Load and clean report data if it's a file path
  if (typeof(report) == "character") {
    df_wide_report <- vroom(report, delim = "\t", col_names = TRUE, 
                            guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
      dplyr::rename_all(tolower) %>% 
      janitor::clean_names() %>% 
      dplyr::rename_all(~ gsub("pg_", "", .)) %>% 
      dplyr::rename_all(~ gsub("^r_", "", .)) %>%  # Remove 'r_' only at the start of string
      dplyr::rename_all(~ gsub("x[0-9]*_", "", .)) %>%  # Remove "x" followed by numbers
      dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA)))  # Replace NaN, -Inf, Inf with NA
  } else {
    df_wide_report <- report
  }
  
  # Load and clean condition setup data if it's a file path
  if (typeof(conditionSetup) == "character") {
    coldata <- vroom(conditionSetup, delim = "\t", col_names = TRUE, 
                     guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
      dplyr::rename_all(tolower) %>% 
      janitor::clean_names() %>% 
      dplyr::select(run_label, condition, replicate, file_name) %>%
      mutate(condition = gsub("[^[:alnum:]]", "_", condition)) %>% 
      dplyr::rename(label = run_label) %>% 
      dplyr::mutate(label = paste0(janitor::make_clean_names(label), "_", quant_col)) %>% 
      dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA))) %>%  # Replace NaN, -Inf, Inf with NA
      dplyr::mutate(renamingcol = sub("^x", "", label)) %>% 
      dplyr::mutate(renamingcol = gsub("_raw_log2quantity$", "", renamingcol)) %>% 
      dplyr::mutate(condition_rep = paste0(condition, "_rep", replicate))
  } else {
    coldata <- conditionSetup
  }
  
  # If conditions_tolower is TRUE, convert the condition column to lowercase
  if (conditions_tolower) {
    coldata <- coldata %>%
      dplyr::mutate(condition = tolower(condition))
  }
  
  # If collapse_conditions is TRUE, collapse underscores in the condition column
  if (collapse_conditions) {
    # Find the minimum number of underscores in the condition column
    min_underscores <- min(str_count(coldata$condition, "_"))
    
    # Collapse extra underscores to match the minimum number
    coldata <- coldata %>%
      dplyr::mutate(
        condition = ifelse(
          stringr::str_count(condition, "_") > min_underscores,
          stringr::str_replace(condition, "(.*)_(.*)", "\\1\\2"), 
          condition
        )
      )
  }
  
  # Suggest possible contrasts using expand.grid based on unique conditions
  if (suggest_contrasts) {
    unique_conditions <- coldata %>% dplyr::distinct(condition) %>% dplyr::pull(condition)
    possible_contrasts <- expand.grid(condition_num = unique_conditions, 
                                      condition_den = unique_conditions) %>%
      dplyr::filter(condition_num != condition_den) %>%
      dplyr::mutate(contrast = paste0(condition_num, "_vs_", condition_den)) %>%
      dplyr::arrange(contrast) %>% 
      dplyr::pull(contrast)
    
    # Print possible contrasts to the console
    message("Possible contrasts based on conditions in the condition setup:")
    message(paste(possible_contrasts, collapse = "\n"))
  }
  
  # If candidates is not NULL, process contrasts or use default diff calculations
  if (!is.null(df_candidates)) {
    if (!is.null(contrasts)) {
      df_candidates_new <- data.frame()
      for (con in contrasts) {
        condition_num <- strsplit(con, "_vs_")[[1]][1]
        condition_den <- strsplit(con, "_vs_")[[1]][2]
        temp <- df_candidates %>% 
          dplyr::filter(condition_numerator == condition_num & condition_denominator == condition_den) %>% 
          dplyr::mutate(diff = avg_log2_ratio, contrast = con)
        temp2 <- df_candidates %>% 
          dplyr::filter(condition_numerator == condition_den & condition_denominator == condition_num) %>% 
          dplyr::mutate(diff = -1 * avg_log2_ratio, contrast = con)
        df_candidates_new <- dplyr::bind_rows(df_candidates_new, temp, temp2)
      }
    } else {
      df_candidates_new <- df_candidates %>% 
        mutate(diff = avg_log2_ratio, contrast = paste0(condition_numerator, "_vs_", condition_denominator))
    }
    
    # Prepare the candidates data
    df_candidates_new <- df_candidates_new %>% 
      dplyr::rename(p.val = pvalue, p.adj = qvalue) %>% 
      dplyr::select(diff, contrast, p.val, p.adj, protein_groups) %>% 
      tidyr::pivot_wider(names_from = "contrast", values_from = c("diff", "p.val", "p.adj"), names_glue = "{contrast}_{.value}")
  } else {
    # If no candidates are provided, create an empty placeholder
    df_candidates_new <- data.frame(protein_groups = unique(df_wide_report$protein_groups))
  }
  
  # Merge candidates with the report data
  combined <- df_candidates_new %>%
    dplyr::right_join(df_wide_report, by = "protein_groups") %>% 
    DEP2::make_unique("genes", "protein_groups", delim = ";") %>% 
    dplyr::mutate(gene_names = genes, protein_ids = protein_groups)
  
  # Update column names based on condition setup
  for (i in 1:nrow(coldata)) {
    current_prefix <- coldata$renamingcol[i]  # The prefix you want to replace
    new_prefix <- coldata$condition_rep[i]    # The new prefix
    
    # Replace the exact prefix only in the columns that match it fully
    colnames(combined) <- str_replace(colnames(combined), 
                                      paste0("^", current_prefix, "(?=_|$)"), 
                                      new_prefix)
  }
  
  # Adjust labels in coldata
  coldata <- coldata %>% 
    dplyr::mutate(old_label = label) %>% 
    dplyr::mutate(label = paste0(condition_rep, "_raw_", quant_col))
  
  # Create the SummarizedExperiment object
  se <- DEP2::make_se(combined, grep(paste0(".*_", quant_col), colnames(combined)), coldata, log2transform = ifelse(quant_col == "log2quantity", FALSE, TRUE))
  
  # Add 'sample' column for compatibility with SEV filtering functions
  colData(se)$sample <- colData(se)$ID
  
  # Return the SE object
  return(se)
}



#########
#' Subset a Dataframe Based on Conditions, Contrasts, and Additional Columns
#'
#' This function subsets a dataframe by selecting specified condition columns, contrast columns, 
#' and optionally extra columns. It also supports selecting summary intensity columns for conditions, 
#' and ensures that the intensity summary columns are selected in a specific order.
#'
#' @param df A dataframe from which to subset.
#' @param id_col A string specifying the name of the ID column. Default is "name".
#' @param condition_list A character vector of conditions. Columns matching these conditions followed by numbers 
#'   will be selected.
#' @param contrast_list A character vector of contrasts. Columns starting with these contrasts will be selected.
#' @param get_summaries A logical value indicating whether to include summary intensity columns for the specified conditions.
#'   The intensity summary columns include mean, standard deviation, min, max, quartiles, IQR, missing counts, and measured counts.
#'   Default is TRUE.
#' @param extra_cols A character vector of additional column names to include in the subset. The columns must exist in the dataframe.
#'   Default is NULL, meaning no extra columns will be selected.
#'
#' @return A dataframe containing the selected columns in the following order: ID column, extra columns (if specified), 
#'   condition columns (including intensity summary columns if `get_summaries` is TRUE in the specified order), and contrast columns.
#'
#' @examples
#' # Example dataframe
#' df <- data.frame(
#'   name = c("A", "B", "C"),
#'   cond1_1 = c(1, 2, 3),
#'   cond1_mean_intensity = c(10, 20, 30),
#'   contrastA_1 = c(100, 200, 300),
#'   col2 = c("X", "Y", "Z")
#' )
#'
#' # Subset the dataframe with conditions and contrasts
#' colab_subset_df(df, id_col = "name", condition_list = c("cond1"), contrast_list = c("contrastA"), extra_cols = c("col2"))
#'
#' @export
colab_subset_df <- function(df, id_col = "name", condition_list, contrast_list, get_summaries = TRUE, extra_cols = NULL) {
  
  # Define intensity suffixes to look for in condition columns if get_summaries is TRUE
  intensity_suffixes <- c("mean_intensity", "sd_intensity", "min_intensity", 
                          "max_intensity", "p25_intensity", "p50_intensity", 
                          "p75_intensity", "iqr_intensity", "missing_count", 
                          "measured_count")
  
  # Create a pattern for condition columns: match the condition followed by any number at the end
  condition_cols <- unlist(lapply(condition_list, function(cond) {
    grep(paste0("^", cond, "_[0-9]+$"), names(df), value = TRUE)
  }))
  
  # Optionally add any columns that start with the condition and end with one of the intensity suffixes
  if (get_summaries) {
    intensity_cols <- unlist(lapply(condition_list, function(cond) {
      unlist(lapply(intensity_suffixes, function(suffix) {
        grep(paste0("^", cond, "_", suffix, "$"), names(df), value = TRUE)
      }))
    }))
    all_condition_cols <- c(condition_cols, intensity_cols)
  } else {
    all_condition_cols <- condition_cols
  }
  
  # Create a pattern for contrast columns: match columns starting with contrast strings
  contrast_cols <- unlist(lapply(contrast_list, function(contrast) {
    grep(paste0("^", contrast), names(df), value = TRUE)
  }))
  
  # Check if extra_cols is provided and exists in the dataframe
  if (!is.null(extra_cols)) {
    extra_cols <- extra_cols[extra_cols %in% names(df)]  # Only keep valid column names
  } else {
    extra_cols <- character(0)  # Empty vector if extra_cols is NULL
  }
  
  # Select columns in the desired order: first 'id_col', then extra_cols, then condition columns (including intensity if selected in the correct order), and finally contrast columns
  df_subset <- df %>%
    dplyr::select(all_of(id_col), all_of(extra_cols), all_of(all_condition_cols), all_of(contrast_cols))
  
  # Return the subsetted dataframe
  return(df_subset)
}


############
#' Volcano plot with targets and significant genes
#'
#' This function creates a volcano plot based on the input summarized experiment (SE) object, highlighting significant genes and specific target genes.
#'
#' @param se A SummarizedExperiment object containing the gene expression data.
#' @param contrast A string specifying the contrast or comparison of interest. This should match the column prefixes for p-values, adjusted p-values, and fold changes in `rowData(se)`.
#' @param id_col A string specifying the column in `rowData(se)` that contains gene IDs or names. Defaults to "gene_names".
#' @param target_names A character vector of target gene names to be highlighted. Defaults to an empty vector.
#' @param label_sign A logical value indicating whether significant genes should be labeled on the plot. Defaults to TRUE.
#' @param label_targets A logical value indicating whether target genes should be labeled on the plot. Defaults to TRUE.
#' @param max.overlaps Maximum number of label overlaps allowed for `ggrepel`. Defaults to 15.
#' @param labelsize Numeric value for the text size of labels. Defaults to 2.
#'
#' @return A ggplot2 object representing the volcano plot.
#' 
#' @details
#' The volcano plot shows the log2 fold changes on the x-axis and the -log10 p-values on the y-axis. Significant genes and target genes are highlighted with different colors and shapes. Optionally, labels can be added for significant and target genes using `ggrepel`.
#'
#' @examples
#' # Example usage:
#' # volcano_plot <- lars_volcano(se, contrast = "condition_A_vs_B", target_names = c("gene1", "gene2"))
#' # print(volcano_plot)
#'
#' @export
lars_volcano <- function(se, contrast, id_col = "gene_names", target_names = c(""), 
                         label_sign = TRUE, label_targets = TRUE, max.overlaps = 15, labelsize = 2) 
{
  pval <- paste0(contrast, "_p.val")
  padj <- paste0(contrast, "_p.adj")
  diff <- paste0(contrast, "_diff")
  sign <- paste0(contrast, "_significant")
  plot_data <- rowData(se) %>%
    as.data.frame() %>%
    dplyr::select(all_of(c(id_col, "name")), starts_with(contrast)) %>%
    mutate(target_gene = factor(ifelse(!!sym(id_col) %in% target_names, "Target", "Non-Target"),
                                levels = c("Non-Target", "Target"))) %>%
    mutate(sign = !!sym(sign)) %>%
    mutate(target_sign = factor(ifelse(!!sym(id_col) %in% target_names, "Target", 
                                       ifelse(!!sym(sign) == TRUE, "Significant", "None")),
                                levels = c("None", "Significant", "Target"))) %>%
    mutate(alpha = ifelse(target_sign == "None", "None", "goi")) %>%
    arrange(target_gene, sign)
  p1 <- ggplot(plot_data, aes_string(x = diff, y = paste0("-log10(", pval, ")"))) +
    geom_point(aes(fill = target_sign, alpha = alpha, shape = sign), 
               size = 1.5, stroke = 0.1, show.legend = TRUE) + 
    scale_fill_manual(values = c(Significant = "#F8756D", 
                                 Target = "#7570b3", None = "grey90"), drop = FALSE) + 
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), drop = FALSE) + 
    scale_alpha_manual(values = c(None = 0.2, goi = 0.7), drop = FALSE) + 
    guides(alpha = FALSE, fill = guide_legend(override.aes = list(shape = 21, size = 2)), 
           shape = guide_legend(override.aes = list(size = 2))) + 
    labs(title = contrast, fill = "Targets", shape = "Significant", 
         x = "Log2(Fold Change)", y = "-Log10(P-value)") + 
    theme_bw()
  # Labels for significant genes
  if (label_sign) {
    p1 <- p1 + geom_text_repel(data = subset(plot_data, sign), aes(label = name), 
                               size = labelsize, color = "black", 
                               min.segment.length = 0, max.overlaps = max.overlaps,
                               segment.size = 0.2, 
                               box.padding = 0.3,  # Adds padding between points and labels
                               point.padding = 0.3,  # Padding between point and label
                               bg.color = "white", # White background for text
                               bg.r = 0.15)  # Border radius around the label
  }
  # Labels for target genes
  if (label_targets) {
    p1 <- p1 + geom_label_repel(data = subset(plot_data, target_gene == "Target"), 
                                aes(label = name), size = labelsize, color = "#7570b3", 
                                min.segment.length = 0, max.overlaps = max.overlaps,
                                segment.size = 0.2, 
                                box.padding = 0.3,  # Adds padding between points and labels
                                point.padding = 0.3,  # Padding between point and label
                                bg.color = "#7570b3",  # White background for label
                                bg.r = 0.15)  # Border radius around the label
  }
  return(p1)
}


######
#' Convert Log2 Differences to Fold Change in SummarizedExperiment Object
#'
#' This function takes a SummarizedExperiment object and adds new columns containing fold-changes. 
#' The function looks for columns ending with "_diff" in the `rowData` slot and uses dplyr::mutate to 'unlog' them from log2 scale. 
#' It preserves `NA` values in the input and returns a modified SummarizedExperiment object.
#' The function also checks to ensure that the rowData of the summarizedExperiment object remains
#' identical before and after the transformation, except for the newly created fold change columns.
#' Upon success, it prints a confirmation message and lists the names of the newly added columns.
#'
#' @param se A SummarizedExperiment object. This object should contain 
#'   rowData with columns that end with "_diff", representing log2 differences.
#' @return A SummarizedExperiment object with the same structure, but with 
#'   additional columns where log2 values are converted to fold change.
#' @examples
#' # Assuming `se` is a SummarizedExperiment object:
#' # se <- log2fc_to_fc(se)
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr mutate across
#' @importFrom stats ifelse
#' @importFrom digest digest
#' @export
log2fc_to_fc <- function(se) {
  # Check if input is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("The input must be a SummarizedExperiment object.")
  }
  
  # Extract rowData
  rowdata_df <- as.data.frame(rowData(se))
  
  # Check if rowData is empty
  if (nrow(rowdata_df) == 0) {
    stop("rowData is empty. The SummarizedExperiment object must contain data in rowData.")
  }
  
  # Check if any columns end with '_diff'
  diff_cols <- grep('_diff$', colnames(rowdata_df))
  if (length(diff_cols) == 0) {
    stop("No columns ending with '_diff' found in rowData.")
  }
  
  # Store a digest (hash) of the original rowData before modifications
  original_rowdata_hash <- digest::digest(rowdata_df, algo = "sha256")
  
  # Perform mutation: convert log2 to fold change, while handling NA values
  rowdata_df <- rowdata_df %>%
    mutate(across(ends_with('_diff'), 
                  ~ ifelse(is.na(.x), NA, 2^.x),  # Convert log2 to fold change and keep NA
                  .names = "{col}_foldchange"))
  
  # Capture the names of the newly created columns
  new_foldchange_cols <- grep('_foldchange$', colnames(rowdata_df), value = TRUE)
  
  # Remove '_diff' from the newly created columns' names
  colnames(rowdata_df) <- sub('_diff_foldchange$', '_foldchange', colnames(rowdata_df))
  
  # Reassign the mutated rowData back to the SummarizedExperiment object
  rowData(se) <- rowdata_df
  
  # Check that the object is identical apart from the newly created columns
  modified_rowdata_hash <- digest::digest(rowdata_df[ , !grepl('_foldchange$', colnames(rowdata_df))], algo = "sha256")
  
  if (original_rowdata_hash != modified_rowdata_hash) {
    stop("The SummarizedExperiment object has been modified in a way other than adding fold change columns.")
  }
  
  # If successful, print confirmation and list newly added columns
  message("The RowData of the SummarizedExperiment object is identical apart from the newly added columns.")
  message("Newly added columns: ", paste(new_foldchange_cols, collapse = ", "))
  
  # Return the updated SummarizedExperiment object
  return(se)
}

#####
library(SummarizedExperiment)
library(dplyr)

#' Modify rowData or colData of a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @param datatype A character string specifying whether to modify "rowData" or "colData". Must be one of "rowData" or "colData".
#' @param modify_function A function to apply to the selected data (e.g., dplyr functions like mutate or left_join).
#' @param ... Additional arguments passed to the modify_function.
#'
#' @return A modified SummarizedExperiment object with updated rowData or colData.
#'
#' @examples
#' se <- SummarizedExperiment(assays = list(counts = matrix(1:12, nrow = 3)))
#' colData(se) <- DataFrame(sample = c("S1", "S2", "S3"))
#' rowData(se) <- DataFrame(name = c("gene1", "gene2", "gene3"))
#'
#' # Modify colData by adding a new column
#' se <- se_modify(se, "colData", mutate, new_col = c("A", "B", "C"))
#'
#' # Modify rowData by adding a new column
#' se <- se_modify(se, "rowData", mutate, new_col = c(10, 20, 30))
#'
#' # Example that would trigger an error by scrambling 'sample' column in colData
#' try(se <- se_modify(se, "colData", mutate, sample = sample(sample)))
#'
#' # Example that would trigger an error by scrambling 'name' column in rowData
#' try(se <- se_modify(se, "rowData", mutate, name = sample(name)))
#'
#' @import SummarizedExperiment
#' @import dplyr
#' @export
se_modify <- function(se, datatype = c("rowData", "colData"), modify_function, ...) {
  # Ensure datatype is either 'rowData' or 'colData'
  datatype <- match.arg(datatype)
  
  # Extract the relevant data (either rowData or colData)
  data <- if (datatype == "rowData") {
    as.data.frame(rowData(se))
  } else {
    as.data.frame(colData(se))
  }
  
  # Store original rownames for safety checks
  original_rownames <- rownames(data)
  
  # Store additional column check depending on datatype
  if (datatype == "rowData") {
    original_name_column <- data$name
  } else {
    original_sample_column <- data$sample
  }
  
  # Apply the specified function to the data
  modified_data <- modify_function(data, ...)
  
  # Safety checks for colData
  if (datatype == "colData") {
    # Check if rownames are still in the same order
    if (!identical(original_rownames, rownames(modified_data))) {
      stop("Error: The rownames of colData have changed. Ensure that the operation preserves rownames.")
    }
    
    # Check if rownames match the 'sample' column
    if (!all(original_rownames == modified_data$sample)) {
      stop("Error: The rownames of colData do not match the 'sample' column after modification.")
    }
  }
  
  # Safety checks for rowData
  if (datatype == "rowData") {
    # Check if rownames are still in the same order
    if (!identical(original_rownames, rownames(modified_data))) {
      stop("Error: The rownames of rowData have changed. Ensure that the operation preserves rownames.")
    }
    
    # Check if rownames match the 'name' column
    if (!all(original_rownames == modified_data$name)) {
      stop("Error: The rownames of rowData do not match the 'name' column after modification.")
    }
  }
  
  # Assign the modified data back to the appropriate slot in SummarizedExperiment
  if (datatype == "rowData") {
    rowData(se) <- DataFrame(modified_data)
  } else {
    colData(se) <- DataFrame(modified_data)
  }
  
  return(se)
}

