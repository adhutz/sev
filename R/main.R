
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
#' @param id_col id column
#'
#' @return volcano plot
#' @importFrom ggplot2 ggplot geom_point scale_color_manual theme_bw labs lims
#' @importFrom ggrepel geom_text_repel
#'

se_volcano<-function(se, contrast_, id_col = "gene_names"){
  LFC <- rowData(se)[,paste0(contrast_, "_diff")]
  p <- rowData(se)[,paste0(contrast_, "_p.val")]%>%-log2(.)
  p.adj <- rowData(se)[,paste0(contrast_, "_p.adj")]%>%-log2(.)
  sig <- rowData(se)[,paste0(contrast_, "_significant")]
  
  change<-ifelse(sig & (LFC>0), "Increase",
                 ifelse(sig & (LFC < 0), "Decrease","None"))
  
  df <- data.frame(LFC,p,sig,change, "SYMBOL" = rowData(se)[[id_col]])
  
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

#' Handling of columns with multiple entries
#'
#' Split columns with multiple entries (need to be separated by a semicolon) into multiple rows with one entry each.
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
  
  if(fun == "mixed"){
    
    temp <- DEP2::impute(se, fun = fun, randna = rowData(se)$randna, ...)
    
  }else{
    
    temp <- DEP2::impute(se, fun = fun, ...)
    
  }
  
  se <- add_assay(se, assay(temp), "imputed_DEP2")
  
  assays(se, withDimnames = FALSE)$imputed <- assay(se) %>% as.data.frame() %>% dplyr::mutate_all(~ ifelse(is.na(.x), 1,0)) %>% as.matrix()
  
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
  mod = lm(y~x)
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
#' @importFrom dplyr select filter mutate
#' @importFrom tidyr pivot_longer pivot_wider
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
#' @importFrom dplyr select mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect ends_with
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
#' @importFrom dplyr select
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
#' @importFrom dplyr select ends_with
#' @importFrom DEP2 get_df_wide
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
#' @importFrom dplyr select ends_with
#' @importFrom DEP2 get_df_wide
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
#' @importFrom dplyr select filter mutate contains rename
#' @importFrom rlang !! 
#' @importFrom ggplot2 sym
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
#' @importFrom dplyr mutate filter group_by summarize
#' @importFrom tidyr pivot_longer
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
#' @importFrom dplyr rename_with mutate select filter
#' @importFrom BiocGenerics cbind
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
#' @importFrom dplyr mutate filter starts_with contains select 
#' @importFrom DEP2 make_unique
#' @importFrom BiocGenerics colnames rownames
#' @importFrom PhosR PhosphoExperiment
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
#'
#' @return A data frame with the processed and merged data.
#'
#' @export
#'
#' @importFrom vroom vroom
#' @importFrom janitor clean_names
#' @importFrom dplyr rename_all, select, mutate, filter
#' @importFrom tidyr pivot_wider
#' @importFrom DEP2 make_se

spectronaut_to_se <- function(candidates = NULL, report = NULL, contrasts = NULL, conditionSetup = NULL){
  
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
      dplyr::rename(label = run_label)
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
  
  
  return(DEP2::make_se(combined, grep(".*_quantity$", colnames(combined)), coldata)) 
}
