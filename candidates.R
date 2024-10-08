

###this is for plotting Venn and Euler diagrams for gene_ids enriched in multiple comparisons

#load libraries
library(ggvenn)
library(tidyverse)
library(readr)
library(stringr)

#Human
  #upload output tables from Deseq2
  #this script is currently set up to be performed on the human comparisons
  #replace h for m to perform mouse comparison
    #human Sm-RIP vs human polyA-RNA
    hrvhx <- read_delim("/filepath to Deseq2 ouptus.tsv", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

    #Modified Sm-RIP (mouse extract + human polyA-RNA + ATP) vs human Sm-RIP
    hrvmxhra <- read_delim("/filepath to Deseq2 ouptus.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

    #Modified Sm-RIP (mouse extract + human polyA-RNA + ATP) vs mouse extract aligned to human genome to remove genomically ambiguous counts
    mxvmxhra <- read_delim("/filepath to Deseq2 ouptus.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

    #Modified Sm-RIP ATP vs no ATP (mouse extract + human polyA-RNA + ATP vs mouse extract + human polyA-RNA)
    mxhrvmxhra <- read_delim("/filepath to Deseq2 ouptus.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

  #make new column that prints "UP", "DOWN", or NA for diffexpressed genes
    #add a column of NOs
      hrvhx$hrvhx <- "NO"
      hrvmxhra$hrvmxhra <- "NO"
      mxvmxhra$mxvmxhra <- "NO"
      mxhrvmxhra$mxhrvmxhra <- "NO"

  #if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    hrvhx$hrvhx[hrvhx$log2FoldChange > 0.6 & hrvhx$padj < 0.05] <- "UP"
    hrvmxhra$hrvmxhra[hrvmxhra$log2FoldChange > 0.6 & hrvmxhra$padj < 0.05] <- "UP"
    mxvmxhra$mxvmxhra[mxvmxhra$log2FoldChange > 0.6 & mxvmxhra$padj < 0.05] <- "UP"
    mxhrvmxhra$mxhrvmxhra[mxhrvmxhra$log2FoldChange > 0.6 & mxhrvmxhra$padj < 0.05] <- "UP"

  #if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    hrvhx$hrvhx[hrvhx$log2FoldChange < -0.6 & hrvhx$padj < 0.05] <- "DOWN"
    hrvmxhra$hrvmxhra[hrvmxhra$log2FoldChange < -0.6 & hrvmxhra$padj < 0.05] <- "DOWN"
    mxvmxhra$mxvmxhra[mxvmxhra$log2FoldChange < -0.6 & mxvmxhra$padj < 0.05] <- "DOWN"
    mxhrvmxhra$mxhrvmxhra[mxhrvmxhra$log2FoldChange < -0.6 & mxhrvmxhra$padj < 0.05] <- "DOWN"

  #remove unnecessary columns
    hrvhx_removed <- hrvhx[,c(1,16)]
    hrvmxhra_removed <- hrvmxhra[,c(1,16)]
    mxvmxhra_removed <- mxvmxhra[,c(1,16)]
    mxhrvmxhra_removed <- mxhrvmxhra[,c(1,16)]

  #merge using "Human_gene_stable_ID" thus yielding columns with differential expression status for each comparison in a column
    updown <- inner_join(hrvhx_removed, hrvmxhra_removed, by = "gene_id")
    updown <- inner_join(updown, mxvmxhra_removed, by = "gene_id")
    updown <- inner_join(updown, mxhrvmxhra_removed, by = "gene_id")

  #function to count occurrences of "UP" in a row
    count_UP <- function(row) {
      sum(row == "UP")}
  #apply the count_UP function to each row of the data frame
    updown$count <- apply(updown, 1, count_UP)

  #make tables that includes only "UP" "gene_id"
    hrvhx_up <- updown%>%filter(hrvhx =="UP")
      hrvhx_up <- hrvhx_up[,c(1)] # 5073 obs
    hrvmxhra_up <- updown%>%filter(hrvmxhra =="UP")
      hrvmxhra_up <- hrvmxhra_up[,c(1)] # 3317 obs
    mxvmxhra_up <- updown%>%filter(mxvmxhra =="UP")
      mxvmxhra_up <- mxvmxhra_up[,c(1)] # 8968 obs
    mxhrvmxhra_up <- updown%>%filter(mxhrvmxhra =="UP")
      mxhrvmxhra_up <- mxhrvmxhra_up[,c(1)] # 485 obs

  #make a list containing each of the columns above
    up_genes <- list('Sm-RIP vs polyA-RNA' = hrvhx_up$gene_id,
                       'ATP vs polyA-RNA' = hrvmxhra_up$gene_id,
                       'ATP vs background' = mxvmxhra_up$gene_id,
                       'ATP vs no ATP' = mxhrvmxhra_up$gene_id)
      #plot venn diagram
      venn <- ggvenn(up_genes,
                       show_percentage = FALSE,
                       fill_color = c("#AADC3275", "#20A48675", "#365D8D75", "#44015475"),
                       stroke_size = 0,
                       set_name_color =c("#AADC3200", "#20A48600", "#365D8D00", "#44015400"),
                       set_name_size = 6,
                       text_color = "white",
                       text_size = 6
      )
      ggsave(filename = "venn.png", plot = venn, device = "png", dpi = 300, height = 4, width = 4, path = "/filepath/")
      
      #make list containing only those gene_id's enriched in all 4 comparisons
      all4 <- updown%>%filter(count == "4") 
      
      # to give all pertinent information, merge with masterfile from atp_vs_noatp.R 
      deseq_atp <- read_delim("filepath to master file.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
      
      #want the ATP-dependent ids that are in the candidate list
      candidates <- deseq_atp %>% filter(ensembl_gene_id %in% all4$gene_id) 
      write_tsv(as.data.frame(candidates), file = "/filepath/candidates.tsv", quote = "none")
   
#make tables for biotype and sm-sites
  #biotype
    h_candidates_biotype <- as.data.frame(table(h_candidates$gene_biotype))    
  #smsite
    h_candidates_sm <- as.data.frame(table(h_candidates$smsite))
  #smtype
    h_candidates_st <- as.data.frame(table(h_candidates$smtype))
   
