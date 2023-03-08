
gat<-GeneAnnoTable(PanelWidth=8L)
got<-GOTable(PanelWidth=8L)
rst <- RowDataTable(PanelWidth = 12L)


#' isee_mini()
#' 
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

#' filter_perseus()
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


#'Impute_perseus()
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

#' se_volcano()
#' 
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

#' Split_genes()
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

#' se_read_in()
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
  data <- read.delim(file, sep="\t")

  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()

  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_ids,
           orig_gene_names = gene_names) %>%
    split_genes(colname = "protein_ids", keep_all = keep_all_proteins) %>%
    split_genes(colname = "gene_names", keep_all = keep_all_genes) %>%
    dplyr::rename("perseus_intensity" = "intensity")

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
#' 
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
#' 
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


#' impute_DEP()
#' 
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



#' Fix_maxq_pig()
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


#' get_network()
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
  
}


#' spectronaut_read_in()
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
#' sample names and groups are read automatically (start with "log2quantity_"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP make_unique
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
    DEP::make_unique("genes", "uni_prot_ids", delim = ";") %>%#this is a function of DEP to make unique names of genes
    select(cond_rep, ID,genes,uni_prot_ids,  protein_descriptions,log2quantity, ibaq) %>% #I just select these columns because otherwise the 'pivot_wider' function is not working because each row contains non-unique info so everything is matched to everything
    tidyr::pivot_wider(names_from = cond_rep, values_from = c(log2quantity, ibaq)) %>% #make wide
    select(-ID) %>% #now I remove the ID col that was added when I ran make_unique, because I need to make unique again, because later I will need both ID and Name column
    DEP::make_unique("genes", "uni_prot_ids", delim = ";") %>% #re-run
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
      mutate(cond_rep = paste0(tolower(condition), "_rep", replicate )) %>%
      select(file_name, cond_rep) %>%
      unique() 
    
    LFQ_labels <- colnames(data)[grep("^log2quantity_", colnames(data))]
    
    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("^log2quantity_", "", LFQ_labels),
                                    condition=gsub(paste0("^log2quantity_|",sep, "[0-9].*"), "", LFQ_labels),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE)) %>%
      merge(filename, by.x = "sample", by.y = "cond_rep")
  }else{
    LFQ_labels<-experimental_design$label
  }
  
  data_se <- make_se(data, which(colnames(data) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}

