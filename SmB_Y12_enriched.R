

#------------------------
# Identify the overlap between transcripts enriched between Y12 and SmB/B'/N RIP and distribution of Sm-sites, Splice-sites, and Branch-points
#------------------------

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(gtable)
library(patchwork)

#------------------------
# Import datasets
#------------------------
  h_cyto_vs_SmB_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/h_cyto_vs_SmB_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  h_cyto_vs_Y12_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/h_cyto_vs_Y12_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # simplify datasets
    hcB <- h_cyto_vs_SmB_deseq_results[,c(1,3,7)] 
    hcB <- na.omit(hcB) 
    
    hcY <- h_cyto_vs_Y12_deseq_results[,c(1,3,7)] 
    hcY <- na.omit(hcY) 

  # import smsite master table
    human_smsite_master <- read_delim("PATH TO/human_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # inner_join with human_smsite_master
    hcB <- inner_join(hcB, human_smsite_master, by = "ensembl_gene_id")
    hcY <- inner_join(hcY, human_smsite_master, by = "ensembl_gene_id")
    
    
#------------------------   
# Make Venn Diagrams for RNAs enriched in both
#------------------------  
# PART 1: make Venn Diagram
  # make data.frames including only significantly enriched genes
    hcB_up <- hcB %>% filter(log2FoldChange > 0.6 & padj < 0.05) 
    hcY_up <- hcY %>% filter(log2FoldChange > 0.6 & padj < 0.05) 
    
  # protein_coding only
    hcB_up_prot <- hcB_up %>% filter(gene_biotype == "protein_coding") 
    hcY_up_prot <- hcY_up %>% filter(gene_biotype == "protein_coding") 
    
    # convert to list for venn plotting
    hBY_vn <- list("SmB" = hcB_up$ensembl_gene_id,
                  "Y12" = hcY_up$ensembl_gene_id)
    
    hBY_prot_vn <- list("SmB" = hcB_up_prot$ensembl_gene_id,
                        "Y12" = hcY_up_prot$ensembl_gene_id)
    
    # plot venn diagram
      #hBY_venn <- ggvenn(hBY_vn,
      hBY_prot_venn <- ggvenn(hBY_prot_vn,
                        show_percentage = FALSE,
                        fill_color = c("#20A48675", "#44015475"),
                        stroke_size = 0,
                        set_name_color = c("#20A486", "#440154"),
                        set_name_size = 6,
                        text_color = "white",
                        text_size = 6)
      
      ggsave(filename = "h_SmB_Y12_enriched_venn.png", plot = hBY_venn, device = "png", dpi = 300, height = 4, width = 4, path = "/Users/dapperwhiskyman/Desktop/deseq2/tables_figures/reviewer_response/SmB_Y12/human/")
      ggsave(filename = "h_SmB_Y12_prot_enriched_venn.png", plot = hBY_prot_venn, device = "png", dpi = 300, height = 4, width = 4, path = "/Users/dapperwhiskyman/Desktop/deseq2/tables_figures/reviewer_response/SmB_Y12/human/")
      
    
    # Define gene sets for Fisher's Exact Test
      h_smb_genes <- hcB_up$ensembl_gene_id
      h_y12_genes <- hcY_up$ensembl_gene_id
      
      h_smb_prot_genes <- hcB_up_prot$ensembl_gene_id
      h_y12_prot_genes <- hcY_up_prot$ensembl_gene_id
    
    # Define the universe of genes in the DESeq2-tested analysis
      h_gene_universe <- union(hcB$ensembl_gene_id, hcY$ensembl_gene_id)

      # for only protein_coding genes
        hcB_prot <- hcB %>% filter(gene_biotype == "protein_coding") 
        hcY_prot <- hcY %>% filter(gene_biotype == "protein_coding") 
        h_gene_prot_universe <- union(hcB_prot$ensembl_gene_id, hcY_prot$ensembl_gene_id)
    
    # Get counts for Fisher's test
      a <- length(intersect(h_smb_genes, h_y12_genes))                  # In both
      b <- length(setdiff(h_smb_genes, h_y12_genes))                    # Only in SmB
      c <- length(setdiff(h_y12_genes, h_smb_genes))                    # Only in Y12
      d <- length(setdiff(h_gene_universe, union(h_smb_genes, h_y12_genes)))  # In neither
      
      e <- length(intersect(h_smb_prot_genes, h_y12_prot_genes))                  # In both
      f <- length(setdiff(h_smb_prot_genes, h_y12_prot_genes))                    # Only in SmB
      g <- length(setdiff(h_y12_prot_genes, h_smb_prot_genes))                    # Only in Y12
      h <- length(setdiff(h_gene_prot_universe, union(h_smb_prot_genes, h_y12_prot_genes)))  # In neither
    
    # Create the contingency table
      h_contingency_table <- matrix(c(a, b, c, d), nrow = 2)
      h_prot_contingency_table <- matrix(c(e, f, g, h), nrow = 2)
      
    # Run Fisher's exact test
      h_fisher_result <- fisher.test(h_contingency_table)
      h_prot_fisher_result <- fisher.test(h_prot_contingency_table)
      # View results
        h_fisher_result$p.value
        h_fisher_result 
        
        h_prot_fisher_result$p.value 
        h_prot_fisher_result 

    
# PART 2: get breakdown of biotype and Sm-sites in dualy enriched genes
  
  # SUB-PART: Import intron counts for cytoplasmic RNA to annotate smsite_master
    # Human  
      h_2025_intron_counts <- read_delim("PATH TO featureCounts OUTPUT FOR COUNTING INTRONIC SEQUENCE/h_2025_intron_counts.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      # remove unneeded columns and simplify sample headers
        h_2025_intron_counts <- h_2025_intron_counts[,c(1,7:9)]
        names(h_2025_intron_counts) <- c("ensembl_gene_id", "M1", "M2", "M3") 
    
        # Compute mean counts by condition
          mean_intron <- h_2025_intron_counts %>%
            rowwise() %>%
            mutate(
              h_cyto_intron    = mean(c_across(starts_with("M")))
              ) %>%
            ungroup() %>%
            dplyr::select(ensembl_gene_id, ends_with("_intron")) 
    
        # Process and only keep gene_ids that have >0 counts
          mean_intron$ensembl_gene_id <- str_split(mean_intron$ensembl_gene_id, "\\.", simplify = T)[,1] 
          h_2025_intronic <- mean_intron %>% filter(h_cyto_intron > 0) 
          
  # Annnotate smsite-master with cyto_intron
    human_smsite_master <- human_smsite_master %>% mutate(intron_cyto = case_when(ensembl_gene_id %in% h_2025_intronic$ensembl_gene_id ~"present", TRUE ~ "absent"))
        
    write_tsv(as.data.frame(human_smsite_master), file = "PATH TO SAVING DIRECTORY/human_smsite_master.tsv", quote = "none")
   
  # select for retaining only those gene_ids pertaining to enrichment with both antibodies
    h_master_both_enriched <- human_smsite_master %>% filter(ensembl_gene_id %in% hcB_up$ensembl_gene_id, ensembl_gene_id %in% hcY_up$ensembl_gene_id) 
    
  # convert to list for venn plotting by canon, 5'ss, intron, bpt
    h_enriched_csm <- h_master_both_enriched %>% filter(smsite == "canon") 
    h_enriched_5ss <- h_master_both_enriched %>% filter(five_ss == "5'ss") 
    h_enriched_intron <- h_master_both_enriched %>% filter(intron_cyto == "present") 
    h_enriched_bpt <- h_master_both_enriched %>% filter(bpt == "stringent") 
    
    h_enriched_vn <- list("Canonical" = h_enriched_csm$ensembl_gene_id,
                          "5'ss" = h_enriched_5ss$ensembl_gene_id,
                          "Bpt" = h_enriched_bpt$ensembl_gene_id,
                          "Intron" = h_enriched_intron$ensembl_gene_id)
    
    # plot venn diagram
      hBY_enriched_venn <- ggvenn(h_enriched_vn,
                         show_percentage = FALSE,
                         fill_color = c("#44015475", "#F9741575", "#63B8FF75", "#ecba0075"),
                         stroke_size = 0,
                         set_name_color =c("#440154", "#F97415", "#63B8FF", "#ecba00"),
                         set_name_size = 6,
                         text_color = "white",
                         text_size = 6)
    
      ggsave(filename = "h_SmB_Y12_enriched_site_venn.png", plot = hBY_enriched_venn, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")


#------------------------   
# Make CDF plots to look at enrichment for different sites
#------------------------  
# PART 1: Inner_join RIP-seq datasets with smsite master
  hB <- inner_join(hcB, human_smsite_master, by = "ensembl_gene_id")
  hY <- inner_join(hcY, human_smsite_master, by = "ensembl_gene_id") 
  
# PART 2: Add columns for separating out different sites
  # add columns for graphing comparisons of five_ss w/ canon Sm-site or intron w/ canon Sm-site
    # to make it easier, for 'obs", order is hB, hY, mB, mY
      x <- human_smsite_master
      x <- hB 
      x <- hY 
      # for bpt
        x_bpt_w_c <- x %>% filter(bpt == "stringent" & smsite == "canon") 
        x_bpt_wo_c <- x %>% filter(bpt == "stringent" & smsite != "canon") 
        x_c_wo_bpt <- x %>% filter(bpt != "stringent" & smsite == "canon")
        x_no_b_no_c <- x %>% filter(bpt != "stringent" & smsite != "canon")
        x <- x %>% mutate(bpt_canon = case_when(ensembl_gene_id %in% x_bpt_w_c$ensembl_gene_id ~ "Bpt_Canon", ensembl_gene_id %in% x_bpt_wo_c$ensembl_gene_id ~ "Bpt_only", ensembl_gene_id %in% x_c_wo_bpt$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))
  
      # for five_ss
        x_5ss_w_c <- x %>% filter(five_ss == "5'ss" & smsite == "canon") 
        x_5ss_wo_c <- x %>% filter(five_ss == "5'ss" & smsite != "canon") 
        x_c_wo_5ss <- x %>% filter(five_ss != "5'ss" & smsite == "canon") 
        x_no_5_no_c <- x %>% filter(five_ss != "5'ss" & smsite != "canon") 
        x <- x %>% mutate(five_ss_canon = case_when(ensembl_gene_id %in% x_5ss_w_c$ensembl_gene_id ~ "5'ss_Canon", ensembl_gene_id %in% x_5ss_wo_c$ensembl_gene_id ~ "5'ss_only", ensembl_gene_id %in% x_c_wo_5ss$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))
  
      # for intron  
        x_i_w_c <- x %>% filter(intron_cyto == "present" & smsite == "canon") 
        x_i_wo_c <- x %>% filter(intron_cyto == "present" & smsite != "canon") 
        x_c_wo_i <- x %>% filter(intron_cyto != "present" & smsite == "canon")
        x_no_i_no_c <- x %>% filter(intron_cyto != "present" & smsite != "canon") 
        x <- x %>% mutate(cyto_intron_canon = case_when(ensembl_gene_id %in% x_i_w_c$ensembl_gene_id ~ "Intron_Canon", ensembl_gene_id %in% x_i_wo_c$ensembl_gene_id ~ "Intron_only", ensembl_gene_id %in% x_c_wo_i$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))
  
    # to make it easier
      h <- x
      hB <- x
      hY <- x
      
  # save master table
    write_tsv(as.data.frame(h), file = "PATH TO SAVING DIRECTORY/human_smsite_master.tsv", quote = "none")
    write_tsv(as.data.frame(hB), file = "PATH TO SAVING DIRECTORY/SmB_Y12/h_SmB.tsv", quote = "none")
    write_tsv(as.data.frame(hY), file = "PATH TO SAVING DIRECTORY/SmB_Y12/h_Y12.tsv", quote = "none")
   
    
    
# PART 3: Import deseq_atp datasets, and make a column for ATP-enriched and Canonical Sm-site containing gene_ids
    
    hdeseq_atp <- read_delim("PATH TO/hdeseq_atp.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # filter for ATP-enriched
      h_atp_up <- hdeseq_atp %>% filter(log2FoldChange >= 0.6 & padj <= 0.05 & smsite == "canon" & gene_biotype == "protein_coding") 
      
    # Annotate columns for ATP-enriched protein_coding genes that have Canonical Sm-sites
      hB <- hB %>% mutate(candidate = case_when(ensembl_gene_id %in% h_atp_up$ensembl_gene_id ~"present", TRUE ~"absent"))
      hY <- hY %>% mutate(candidate = case_when(ensembl_gene_id %in% h_atp_up$ensembl_gene_id ~"present", TRUE ~"absent"))
      
    # save master table
      write_tsv(as.data.frame(hB), file = "PATH TO SAVING DIRECTORY/h_SmB.tsv", quote = "none")
      write_tsv(as.data.frame(hY), file = "PATH TO SAVING DIRECTORY/h_Y12.tsv", quote = "none")
      
      
    
# PART 4: Prepare data for plotting  
  # isolate only protein_coding
    hB_prot <- hB %>% filter(gene_biotype == "protein_coding") 
    hY_prot <- hY %>% filter(gene_biotype == "protein_coding") 
   
  # remove U7
    hB_prot <- hB_prot %>% filter(smsite != "U7") 
    hY_prot <- hY_prot %>% filter(smsite != "U7") 
    
    
# PART 4: Plot CDFs by sites
  # to make it easier
    x <- hB_prot
    x <- hY_prot
    
  # CDF plot
    #bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt)) +
    #fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss)) +
    #int_plot <- ggplot(x, aes(log2FoldChange, color=intron_cyto)) +
    #sm_plot <- ggplot(x, aes(log2FoldChange, color=smsite)) +
    #sm_bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt_canon)) +
    #sm_fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss_canon)) +
    #sm_int_plot <- ggplot(x, aes(log2FoldChange, color=cyto_intron_canon)) +
    cand_plot <- ggplot(x, aes(log2FoldChange, color=candidate)) +
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) anti-Sm RIP vs Cytoplasmic RNA", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      #coord_cartesian(xlim=c(-1.5,1.5)) + # for B_prot
      coord_cartesian(xlim=c(-4,4)) + # for Y_prot
        #scale_color_manual("Bpt:", values = c("#2828cc", "#36648B", "#63B8FF"), breaks = c("No Bpt", "consensus", "stringent"), labels = c("No Bpt", "Consensus", "Stringent"))+ 
        #scale_color_manual("5'ss:", values = c("#2828cc", "#B80069", "#F97415"), breaks = c("No 5'ss", "5'ss_bulge", "5'ss"), labels = c("No 5'ss", "5'ss_bulge", "5'ss"))+ 
        #scale_color_manual("Intron:", values = c("#a9a2a2", "#ecba00"), breaks = c("absent", "present"), labels = c("No Intron", "Intron"))+
        #scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#440154"), breaks = c("absent", "noncanon", "canon" ), labels = c("Absent", "Noncanonical", "Canonical")) + 
        #scale_color_manual("Site", values = c("#7c7887", "#63B8FF", "#440154", "#f0199e"), breaks = c("Absent", "Bpt_only", "Canon_only", "Bpt_Canon" ), labels = c("Absent", "Bpt_only", "Canonical-Sm only", "Bpt & Canonical-Sm")) + 
        #scale_color_manual("Site", values = c("#4C4545", "#F97415", "#440154", "#CE2727"), breaks = c("Absent", "5'ss_only", "Canon_only", "5'ss_Canon" ), labels = c("Absent", "5'ss only", "Canonical-Sm only", "5'ss & Canonical-Sm")) + 
        #scale_color_manual("Site:", values = c("#a9a2a2", "#ecba00","#440154", "#3bbf5e"), breaks = c("Absent", "Intron_only", "Canon_only", "Intron_Canon"), labels = c("Absent", "Intron_only", "Canon_only", "Intron & Canon"))+
        scale_color_manual("Candidate", values = c("#FDE725", "#440154"), breaks = c("absent", "present" ), labels = c("Absent", "Candidate")) + 
        theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            panel.grid = element_blank(),
            axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none",
            #legend.box.background = element_rect("white"), 
            #legend.key = element_blank(), 
            #legend.direction = "horizontal",
            #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
            #legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
      
      # for hB_prot
        ggsave(filename = "h_SmB_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_SmB_Candidates_CDF.png", plot = cand_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        
      # for hY_prot
        ggsave(filename = "h_Y12_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_Y12_Candidates_CDF.png", plot = cand_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        
  
# PART 5: Perform Wilcoxon Rank Sum Test with Continuity Correction and Benjamini-Hochberg False Discovery Rate corrections for multiple testing
  # to make it easier
    x <- hB_prot
    x <- hY_prot
    
      table(x$bpt)
      
      n <- subset(x, bpt == "No Bpt")$log2FoldChange
      c <- subset(x, bpt == "consensus")$log2FoldChange
      s <- subset(x, bpt == "stringent")$log2FoldChange
      cn <- wilcox.test(c, n, exact = FALSE, alternative = 'greater')$p.value 
      sn <- wilcox.test(s, n, exact = FALSE, alternative = 'greater')$p.value 
      sc <- wilcox.test(s, c, exact = FALSE, alternative = 'greater')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(cn, sn, sc)
        p.adjust(raw_pvals, method = "BH") 
      
      table(x$five_ss)
      
      n <- subset(x, five_ss == "No 5'ss")$log2FoldChange
      b <- subset(x, five_ss == "5'ss_bulge")$log2FoldChange
      s <- subset(x, five_ss == "5'ss")$log2FoldChange
      bn <- wilcox.test(b, n, exact = FALSE, alternative = 'greater')$p.value 
      sn <- wilcox.test(s, n, exact = FALSE, alternative = 'greater')$p.value 
      sb <- wilcox.test(s, b, exact = FALSE, alternative = 'greater')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(bn, sn, sb)
        p.adjust(raw_pvals, method = "BH") 
      
      table(x$intron_cyto)
      
      a <- subset(x, intron_cyto == "absent")$log2FoldChange
      p <- subset(x, intron_cyto == "present")$log2FoldChange
        wilcox.test(p, a, exact = FALSE, alternative = 'greater')$p.value 
      
      table(x$smsite)
     
      a <- subset(x, smsite == "absent")$log2FoldChange
      n <- subset(x, smsite == "noncanon")$log2FoldChange
      c <- subset(x, smsite == "canon")$log2FoldChange
      na <- wilcox.test(n, a, exact = FALSE, alternative = 'greater')$p.value 
      ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
      cn <- wilcox.test(c, n, exact = FALSE, alternative = 'greater')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(na, ca, cn)
        p.adjust(raw_pvals, method = "BH") 
      
      table(x$bpt_canon)
      
      a <- subset(x, bpt_canon == "Absent")$log2FoldChange
      b <- subset(x, bpt_canon == "Bpt_only")$log2FoldChange
      c <- subset(x, bpt_canon == "Canon_only")$log2FoldChange
      bc <- subset(x, bpt_canon == "Bpt_Canon")$log2FoldChange
      ba <- wilcox.test(b, a, exact = FALSE, alternative = 'greater')$p.value 
      ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
      bca <- wilcox.test(bc, a, exact = FALSE, alternative = 'greater')$p.value 
      cb <- wilcox.test(c, b, exact = FALSE, alternative = 'greater')$p.value
      bcb <- wilcox.test(bc, b, exact = FALSE, alternative = 'greater')$p.value
      bcc <- wilcox.test(bc, c, exact = FALSE, alternative = 'greater')$p.value
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(ba, ca, bca, cb, bcb, bcc)
        p.adjust(raw_pvals, method = "BH")
      
      table(x$five_ss_canon)
     
      a <- subset(x, five_ss_canon == "Absent")$log2FoldChange
      s <- subset(x, five_ss_canon == "5'ss_only")$log2FoldChange
      c <- subset(x, five_ss_canon == "Canon_only")$log2FoldChange
      sc <- subset(x, five_ss_canon == "5'ss_Canon")$log2FoldChange
      sa <- wilcox.test(s, a, exact = FALSE, alternative = 'greater')$p.value 
      ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
      sca <- wilcox.test(sc, a, exact = FALSE, alternative = 'greater')$p.value 
      cs <- wilcox.test(c, s, exact = FALSE, alternative = 'greater')$p.value
      scs <- wilcox.test(sc, s, exact = FALSE, alternative = 'greater')$p.value
      scc <- wilcox.test(sc, c, exact = FALSE, alternative = 'greater')$p.value
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(sa, ca, sca, cs, scs, scc)
        p.adjust(raw_pvals, method = "BH")
          
      table(x$cyto_intron_canon)
      
      a <- subset(x, cyto_intron_canon == "Absent")$log2FoldChange
      i <- subset(x, cyto_intron_canon == "Intron_only")$log2FoldChange
      c <- subset(x, cyto_intron_canon == "Canon_only")$log2FoldChange
      ic <- subset(x, cyto_intron_canon == "Intron_Canon")$log2FoldChange
      ia <- wilcox.test(i, a, exact = FALSE, alternative = 'greater')$p.value 
      ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
      ica <- wilcox.test(ic, a, exact = FALSE, alternative = 'greater')$p.value 
      ci <- wilcox.test(c, i, exact = FALSE, alternative = 'greater')$p.value
      ici <- wilcox.test(ic, i, exact = FALSE, alternative = 'greater')$p.value
      icc <- wilcox.test(ic, c, exact = FALSE, alternative = 'greater')$p.value
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(ia, ca, ica, ci, ici, icc)
        p.adjust(raw_pvals, method = "BH") 
              
      table(x$candidate)
      
      a <- subset(x, candidate == "absent")$log2FoldChange
      p <- subset(x, candidate == "present")$log2FoldChange
        wilcox.test(p, a, exact = FALSE, alternative = 'greater')$p.value 
      
              
#------------------------   
# Collect the distribution of total, enriched, and depleted gene_ids with site information
#------------------------   
  
# PART 1: Filter for only enriched protein_coding genes
  # to make it easier
    x <- hB_prot
    x <- hY_prot
   
    # filter for significantly enriched or depleted gene_ids       
      xu <- x %>% filter(log2FoldChange >= 0.6 & padj <= 0.05) 
      xd <- x %>% filter(log2FoldChange <= -0.6 & padj <= 0.05)
      
# PART 2: Create data.frames for each site-type      
  # for bpt
    b <- as.data.frame(table(x$bpt))
      names(b) <- c("Site", "All")
    bu <- as.data.frame(table(xu$bpt))
      names(bu) <- c("Site", "Sm-enriched")
    bd <- as.data.frame(table(xd$bpt))
      names(bd) <- c("Site", "Sm-depleted")
      # full_join into one table
        b <- full_join(b, bu, by = "Site")
        b <- full_join(b, bd, by = "Site")
    
  # for five_ss
    f <- as.data.frame(table(x$five_ss))
      names(f) <- c("Site", "All")
    fu <- as.data.frame(table(xu$five_ss))
      names(fu) <- c("Site", "Sm-enriched")
    fd <- as.data.frame(table(xd$five_ss))
      names(fd) <- c("Site", "Sm-depleted")
      # full_join into one table
        f <- full_join(f, fu, by = "Site")
        f <- full_join(f, fd, by = "Site")
      
  # for Sm-site
    s <- as.data.frame(table(x$smsite))
      names(s) <- c("Site", "All")
    su <- as.data.frame(table(xu$smsite))
      names(su) <- c("Site", "Sm-enriched")
    sd <- as.data.frame(table(xd$smsite))
      names(sd) <- c("Site", "Sm-depleted")
      # full_join into one table
        s <- full_join(s, su, by = "Site")
        s <- full_join(s, sd, by = "Site")
          
  # for Sm-site
    i <- as.data.frame(table(x$intron_cyto))
      names(i) <- c("Site", "All")
    iu <- as.data.frame(table(xu$intron_cyto))
      names(iu) <- c("Site", "Sm-enriched")
    id <- as.data.frame(table(xd$intron_cyto))
      names(id) <- c("Site", "Sm-depleted")
      # full_join into one table
      i <- full_join(i, iu, by = "Site")
      i <- full_join(i, id, by = "Site")
      
  # Part 3: Bind all together and write tsv
      x <- bind_rows(b, f, s, i)
      
      hSmB <- x
      hY12 <- x
      
      write_tsv(as.data.frame(hSmB), file = "PATH TO SAVING DIRECTORY/h_SmB_site_breakdown.tsv", quote = "none")
      write_tsv(as.data.frame(hY12), file = "PATH TO SAVING DIRECTORY/h_Y12_site_breakdown.tsv", quote = "none")

    
#------------------------   
# Create volcano, density, and CDF plots for RNA species
#------------------------     
      
# PART 1: Make a column grouping different RNA species
  # to make it easier
    x <- hB
    x <- hY
   
  # filter by biotype, then create new column
      sn <- x %>% filter(gene_biotype == "snRNA") 
      m <- x %>% filter(gene_biotype == "protein_coding") 
      lnc <- x %>% filter(gene_biotype == "lncRNA") 
      sca <- x %>% filter(gene_biotype == "scaRNA") 
      sno <- x %>% filter(gene_biotype == "snoRNA") 
      snc <- bind_rows(sca, sno) 
      x <- x%>%mutate(rna_species = case_when(ensembl_gene_id %in% sn$ensembl_gene_id ~ "snRNA",
                                                                  ensembl_gene_id %in% upsmc$ensembl_gene_id ~ "Candidate",
                                                                  ensembl_gene_id %in% m$ensembl_gene_id ~ "mRNA",
                                                                  ensembl_gene_id %in% lnc$ensembl_gene_id ~ "lncRNA",
                                                                  ensembl_gene_id %in% snc$ensembl_gene_id ~ "sncRNA",
                                                                  TRUE ~ "Other"))   
  
  
# PART 2: Volcano Plot, Density Plot, & CDF Plot
  # Volcano
    volcano <- ggplot(data = x, aes(x = log2FoldChange, y = -log10(padj), col = rna_species)) + 
      geom_point(alpha=0.35) + 
      labs(x = "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "-log10(p-adj)") +
      theme_minimal() +
      scale_color_manual("Biotype:", values = c("#FDE725", "#5DC863", "#21908c", "#3B528B", "#440154"), breaks = c("Other", "sncRNA", "snRNA", "lncRNA", "mRNA"), labels = c("Other", "sno/scaRNA", "snRNA", "lncRNA", "mRNA"))+
      coord_cartesian(xlim = c(-10,15)) +
      theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
            panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.x = element_blank(),
            axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none")+ 
      #legend.box.background = element_blank(), 
      #legend.key = element_blank(), 
      #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
      #legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))+
      geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")
      
  # Density for volcano plot
    density <- ggplot(data = x, aes(x = log2FoldChange, color = rna_species)) + 
      geom_density(linewidth = 1) +
      labs(x= "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "Density") +
      coord_cartesian(xlim = c(-10,15)) +
      scale_color_manual(" ", values = c("#FDE725", "#5DC863", "#21908c", "#3B528B", "#440154"), breaks = c("Other", "sncRNA", "snRNA", "lncRNA", "mRNA"), labels = c("Other", "sno/scaRNA", "snRNA", "lncRNA", "mRNA"))+
      theme(legend.position = "top",
            legend.box.background = element_blank(), 
            legend.key = element_blank(), 
            legend.box = "horizontal",
            legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.25),
            panel.grid.minor = element_blank()) +
      geom_vline(xintercept = c(-0.6, 0.6), col="black", linetype = "dashed") +
      guides(color = guide_legend(nrow = 2))
      
  # CDF (don't know if this will be used)
    cdf <- ggplot(x, aes(log2FoldChange, color=rna_species)) +
      stat_ecdf(geom = "step", linewidth = 1.25) +
      labs(x = "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      coord_cartesian(xlim=c(-10,15)) +
      scale_color_manual("Biotype:", values = c("#FDE725", "#5DC863", "#21908c", "#3B528B", "#440154"), breaks = c("Other", "sncRNA", "snRNA", "lncRNA", "mRNA"), labels = c("Other", "sno/scaRNA", "snRNA", "lncRNA", "mRNA"))+
      theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            panel.grid = element_blank(),
            axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none",
            #legend.box.background = element_rect("white"), 
            #legend.key = element_blank(), 
            #legend.direction = "horizontal",
            #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
            #legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
  
  # Store plots for assembling figures
    Bv <- volcano
    Bd <- density
    Bc <- cdf
    
    Yv <- volcano
    Yd <- density
    Yc <- cdf
    
    design <- "
      11
      23
      45
      "
    p <- guide_area() + Bv + Yv + Bd + Yd + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect", axis_titles = "collect")
    ggsave(filename = "h_SmB_Y12_volcano_density.png", p = p, device = "png", dpi = 300, height = 6, width = 8, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "m_SmB_Y12_volcano_density.png", p = p, device = "png", dpi = 300, height = 6, width = 8, path = "PATH TO SAVING DIRECTORY")
    
  
#------------------------
# Plots for exon number as a proxy for probability of intron retention
#------------------------  
    
# PART 1: import exon number tables Dr. Manu Sanjeev. These are only protein_coding. Annotate Y12 and SmB/B'/N datasets
    
  hdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/hdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # select only enseml_gene_id and ExonNumber columns
    en <- hdeseq_equil_PC_ExonNumbers %>% dplyr::select(ensembl_gene_id, ExonNumber) 
    
    hB_prot <- inner_join(hB_prot, en, by = "ensembl_gene_id") 
    hY_prot <- inner_join(hY_prot, en, by = "ensembl_gene_id") 
  
# PART 2: CDF Plots
    
  x <- hB_prot
  x <- hY_prot

  # CDF plot
    x_plot <- ggplot(x, aes(log2FoldChange, color=ExonNumber)) +
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) anti-Sm RIP vs Cytoplasmic RNA", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      #coord_cartesian(xlim=c(-1.5,1.5)) + # for B_prot
      coord_cartesian(xlim=c(-4,4)) + # for Y_prot
      scale_color_manual("Exon Number:", values = c("#413d3d", "#7e514a","#b46a48", "#e18c39", "#feb701"), breaks = c("1 to 4", "5 to 7","8 to 11", "12 to 17", "17+"), labels = c("1 to 4", "5 to 7","8 to 11", "12 to 17", "17+"))+
      theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            panel.grid = element_blank(),
            axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none",
            #legend.box.background = element_rect("white"), 
            #legend.key = element_blank(), 
            #legend.direction = "horizontal",
            #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
            #legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
    
    hBxn_plot <- x_plot
    hYxn_plot <- x_plot

    # save plots
    ggsave(filename = "h_SmB_ExonNumber_CDF.png", plot = hBxn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "h_Y12_ExonNumber_CDF.png", plot = hYxn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")

# PART 3: Wilcoxon Rank Sum Test and Benjamini-Hochberg Fasle discovery Rate correction for multiple testing
    
    x <- hB_prot
    x <- hY_prot
    
    table(x$ExonNumber)

    one <- subset(x, ExonNumber == "1 to 4")$log2FoldChange
    two <- subset(x, ExonNumber == "5 to 7")$log2FoldChange
    three <- subset(x, ExonNumber == "8 to 11")$log2FoldChange
    four <- subset(x, ExonNumber == "12 to 17")$log2FoldChange
    five <- subset(x, ExonNumber == "17+")$log2FoldChange
    w1 <- wilcox.test(two, one, exact = FALSE, alternative = 'greater')$p.value 
    w2 <- wilcox.test(three, two, exact = FALSE, alternative = 'greater')$p.value 
    w3 <- wilcox.test(four, three, exact = FALSE, alternative = 'greater')$p.value
    w4 <- wilcox.test(five, four, exact = FALSE, alternative = 'greater')$p.value
    # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
    raw_pvals <- c(w1, w2, w3, w4)
    p.adjust(raw_pvals, method = "BH")
    
#------------------------
# Plots for intronic reads as quintile
#------------------------  
    
# PART 1: import intronic counts
    
  h_2025_intronic <- read_delim("Desktop/deseq2/tables_figures/reviewer_response/intron/h_2025_intronic.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  x <- hB_prot
  x <- hY_prot
  
  x <- inner_join(x, h_2025_intronic, by = "ensembl_gene_id") # 11404 obs
    
  x <- x %>% mutate(h_cyto_intron_quintile = ntile(h_cyto_intron, 5))  
    
# PART 2: CDF Plot
    
  plot <- ggplot(x, aes(log2FoldChange, color= factor(h_cyto_intron_quintile))) +
    stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
    labs(x = "log2(f∆) anti-Sm RIP vs Human Cytoplasmic RNA", y = "Fraction of Data") +
    scale_y_continuous(position = "right") +
    coord_cartesian(xlim=c(-1.5,1.5)) + # for B_prot
    #coord_cartesian(xlim=c(-4,4)) + # for Y_prot
    scale_color_manual("Intron Quintile:", values = c("#a9a2a2", "#d09f94","#ea9f78", "#f4a750", "#ebb800"), breaks = c("1", "2", "3", "4", "5"), labels = c("1st", "2nd", "3rd", "4th", "5th"))+
    theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
          panel.grid = element_blank(),
          axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
          axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
          axis.title.x = element_blank(),
          legend.position = "none",
          #legend.box.background = element_rect("white"), 
          #legend.key = element_blank(), 
          #legend.direction = "horizontal",
          #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
          #legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
          panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
    
  # save plots
    ggsave(filename = "h_SmB_IntronQuintile_CDF.png", plot = plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "h_Y12_IntronQuintile_CDF.png", plot = plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
  
# PART 3: Wilcoxon Rank Sum Test and Benjamini-Hochberg False discovery Rate correction for multiple testing
    
    table(x$h_cyto_intron_quintile)
   
    one <- subset(x, h_cyto_intron_quintile == "1")$log2FoldChange
    two <- subset(x, h_cyto_intron_quintile == "2")$log2FoldChange
    three <- subset(x, h_cyto_intron_quintile == "3")$log2FoldChange
    four <- subset(x, h_cyto_intron_quintile == "4")$log2FoldChange
    five <- subset(x, h_cyto_intron_quintile == "5")$log2FoldChange
    w1 <- wilcox.test(two, one, exact = FALSE, alternative = 'greater')$p.value 
    w2 <- wilcox.test(three, two, exact = FALSE, alternative = 'greater')$p.value 
    w3 <- wilcox.test(four, three, exact = FALSE, alternative = 'greater')$p.value
    w4 <- wilcox.test(five, four, exact = FALSE, alternative = 'greater')$p.value
    # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
    raw_pvals <- c(w1, w2, w3, w4)
    p.adjust(raw_pvals, method = "BH")
  
  
  
    
