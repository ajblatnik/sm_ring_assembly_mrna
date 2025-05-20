
#------------------------
# Identify the effect Splice site, branch-point, and intron has in Sm-protein enrichment from 2023 RSB-500 + 0.1% NP-40 RIP-seq datasets
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
  hr_vs_hx_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/hr_vs_hx_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  mxhr_vs_mxhra_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/mx+hr_vs_mx+hr+a_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

  mr_vs_mx_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/mr_vs_mx_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  hx_mr_vs_hx_mr_a_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/hx+mr_vs_hx+mr+a_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

  # simplify datasets
    hrx <- hr_vs_hx_deseq_results[,c(1,3,7)] 
    hrx <- na.omit(hrx) 
    names(hrx)[1] <- "ensembl_gene_id"

    hxa <- mxhr_vs_mxhra_deseq_results[,c(1,3,7)] #
    hxa <- na.omit(hxa) 
    names(hxa)[1] <- "ensembl_gene_id"

    mrx <- mr_vs_mx_deseq_results[,c(1,3,7)] 
    mrx <- na.omit(mrx) 
    names(mrx)[1] <- "ensembl_gene_id"

    mxa <- hx_mr_vs_hx_mr_a_deseq_results[,c(1,3,7)] 
    mxa <- na.omit(mxa) 
    names(mxa)[1] <- "ensembl_gene_id"
    

#------------------------   
# Make CDF plots to look at enrichment for different sites
#------------------------  
  # PART 1: Inner_join RIP-seq datasets with smsite master
          
    # Import smsite_master sheets
      human_smsite_master <- human_smsite_ss_bpt_intron_master <- read_delim("PATH TO/human_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      mouse_smsite_master <- mouse_smsite_ss_bpt_intron_master <- read_delim("PATH TO/mouse_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
    # Inner_join smsite_master with deseq comparison results        
      hrx <- inner_join(hrx, human_smsite_master, by = "ensembl_gene_id")
      hxa <- inner_join(hxa, human_smsite_master, by = "ensembl_gene_id")

      mrx <- inner_join(mrx, mouse_smsite_master, by = "ensembl_gene_id")
      mxa <- inner_join(mxa, mouse_smsite_master, by = "ensembl_gene_id")

  # PART 2: Add columns for separating out different sites
    # add columns for graphing comparisons of five_ss w/ canon Sm-site or intron w/ canon Sm-site
    # to make it easier, for 'obs", order is hrx, hxa, mrx, mxa
      x <- hrx 
      x <- hxa 
      x <- mrx 
      x <- mxa 
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
        x_i_w_c <- x %>% filter(intron_equil == "present" & smsite == "canon")
        x_i_wo_c <- x %>% filter(intron_equil == "present" & smsite != "canon") 
        x_c_wo_i <- x %>% filter(intron_equil != "present" & smsite == "canon") 
        x_no_i_no_c <- x %>% filter(intron_equil != "present" & smsite != "canon") 
        x <- x %>% mutate(equil_intron_canon = case_when(ensembl_gene_id %in% x_i_w_c$ensembl_gene_id ~ "Intron_Canon", ensembl_gene_id %in% x_i_wo_c$ensembl_gene_id ~ "Intron_only", ensembl_gene_id %in% x_c_wo_i$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))

        x_i_w_c <- x %>% filter(intron_polyA == "present" & smsite == "canon") 
        x_i_wo_c <- x %>% filter(intron_polyA == "present" & smsite != "canon") 
        x_c_wo_i <- x %>% filter(intron_polyA != "present" & smsite == "canon") 
        x_no_i_no_c <- x %>% filter(intron_polyA != "present" & smsite != "canon")
        x <- x %>% mutate(polyA_intron_canon = case_when(ensembl_gene_id %in% x_i_w_c$ensembl_gene_id ~ "Intron_Canon", ensembl_gene_id %in% x_i_wo_c$ensembl_gene_id ~ "Intron_only", ensembl_gene_id %in% x_c_wo_i$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))
        
    # to make it easier
      hrx <- x
      hxa <- x
      mrx <- x
      mxa <- x

  # save master table
    write_tsv(as.data.frame(hrx), file = "PATH TO SAVING DIRECTORY/h_equil.tsv", quote = "none")
    write_tsv(as.data.frame(hxa), file = "PATH TO SAVING DIRECTORY/h_atp.tsv", quote = "none")
    write_tsv(as.data.frame(mrx), file = "PATH TO SAVING DIRECTORY/m_equil.tsv", quote = "none")
    write_tsv(as.data.frame(mxa), file = "PATH TO SAVING DIRECTORY/m_atp.tsv", quote = "none")
    
  # PART 4: import exon number plots
    
    hdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/hdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    hdeseq_atp_PC_ExonNumbers <- read_delim("PATH TO/hdeseq_atp_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    mdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/mdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    mdeseq_atp_PC_ExonNumbers <- read_delim("PATH TO/mdeseq_atp_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      
  # PART 3: Prepare data for plotting  
    # isolate only protein_coding
      hrx_prot <- hrx %>% filter(gene_biotype == "protein_coding") 
      hxa_prot <- hxa %>% filter(gene_biotype == "protein_coding") 
      mrx_prot <- mrx %>% filter(gene_biotype == "protein_coding") 
      mxa_prot <- mxa %>% filter(gene_biotype == "protein_coding") 
    # remove U7
      hrx_prot <- hrx_prot %>% filter(smsite != "U7") 
      hxa_prot <- hxa_prot %>% filter(smsite != "U7") 
      mrx_prot <- mrx_prot %>% filter(smsite != "U7") 
      mxa_prot <- mxa_prot %>% filter(smsite != "U7") 

  # PART 4: Plot CDFs by sites
    # to make it easier
      x <- hrx_prot
      x <- hxa_prot
      x <- mrx_prot
      x <- mxa_prot

    # CDF plot
      #bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt)) +
      #fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss)) +
      #int_plot <- ggplot(x, aes(log2FoldChange, color=intron_equil)) +
      #int_plot <- ggplot(x, aes(log2FoldChange, color=intron_polyA)) +
      #sm_plot <- ggplot(x, aes(log2FoldChange, color=smsite)) +
      #sm_bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt_canon)) +
      #sm_fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss_canon)) +
      #sm_int_plot <- ggplot(x, aes(log2FoldChange, color=equil_intron_canon)) +
      sm_int_plot <- ggplot(x, aes(log2FoldChange, color=polyA_intron_canon)) +
        stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
        labs(x = "log2(f∆) anti-Sm RIP vs Cytoplasmic RNA", y = "Fraction of Data") +
        scale_y_continuous(position = "right") +
        coord_cartesian(xlim=c(-1,1)) + # for atp
        #coord_cartesian(xlim=c(-4,4)) + # for equil
        #scale_color_manual("Bpt:", values = c("#2828cc", "#36648B", "#63B8FF"), breaks = c("No Bpt", "consensus", "stringent"), labels = c("No Bpt", "Consensus", "Stringent"))+ 
        #scale_color_manual("5'ss:", values = c("#2828cc", "#B80069", "#F97415"), breaks = c("No 5'ss", "5'ss_bulge", "5'ss"), labels = c("No 5'ss", "5'ss_bulge", "5'ss"))+ 
        #scale_color_manual("Intron:", values = c("#a9a2a2", "#ecba00"), breaks = c("absent", "present"), labels = c("No Intron", "Intron"))+
        #scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#440154"), breaks = c("absent", "noncanon", "canon" ), labels = c("Absent", "Noncanonical", "Canonical")) + 
        #scale_color_manual("Site", values = c("#7c7887", "#63B8FF", "#440154", "#f0199e"), breaks = c("Absent", "Bpt_only", "Canon_only", "Bpt_Canon" ), labels = c("Absent", "Bpt_only", "Canonical-Sm only", "Bpt & Canonical-Sm")) + 
        #scale_color_manual("Site", values = c("#4C4545", "#F97415", "#440154", "#CE2727"), breaks = c("Absent", "5'ss_only", "Canon_only", "5'ss_Canon" ), labels = c("Absent", "5'ss only", "Canonical-Sm only", "5'ss & Canonical-Sm")) + 
        scale_color_manual("Site:", values = c("#a9a2a2", "#ecba00","#440154", "#3bbf5e"), breaks = c("Absent", "Intron_only", "Canon_only", "Intron_Canon"), labels = c("Absent", "Intron_only", "Canon_only", "Intron & Canon"))+
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

      # for hrx_prot (equil)
        ggsave(filename = "h_equil_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_equil_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")

      # for hxa_prot (atp)
        ggsave(filename = "h_atp_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "h_atp_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")

      # for mrx_prot (equil)
        ggsave(filename = "m_equil_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_equil_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")

      # for mxa_prot (atp)
        ggsave(filename = "m_atp_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
        ggsave(filename = "m_atp_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")


  # PART 5: Perform Wilcoxon Rank Sum Test with Continuity Correction and Benjamini-Hochberg False Discovery Rate corrections for multiple testing
    # to make it easier
      x <- hrx_prot
      x <- hxa_prot
      x <- mrx_prot
      x <- mxa_prot

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
        d <- subset(x, five_ss == "5'ss_bulge")$log2FoldChange
        s <- subset(x, five_ss == "5'ss")$log2FoldChange
        dn <- wilcox.test(d, n, exact = FALSE, alternative = 'greater')$p.value 
        sn <- wilcox.test(s, n, exact = FALSE, alternative = 'greater')$p.value 
        sd <- wilcox.test(s, d, exact = FALSE, alternative = 'greater')$p.value 
        # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
          raw_pvals <- c(dn, sn, sd)
          p.adjust(raw_pvals, method = "BH") 

        table(x$intron_equil)
       
        a <- subset(x, intron_equil == "absent")$log2FoldChange
        p <- subset(x, intron_equil == "present")$log2FoldChange
          wilcox.test(p, a, exact = FALSE, alternative = 'greater')$p.value 
          
        table(x$intron_polyA)
      
        a <- subset(x, intron_polyA == "absent")$log2FoldChange
        p <- subset(x, intron_polyA == "present")$log2FoldChange
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

        table(x$equil_intron_canon)
       
        a <- subset(x, equil_intron_canon == "Absent")$log2FoldChange
        i <- subset(x, equil_intron_canon == "Intron_only")$log2FoldChange
        c <- subset(x, equil_intron_canon == "Canon_only")$log2FoldChange
        ic <- subset(x, equil_intron_canon == "Intron_Canon")$log2FoldChange
        ia <- wilcox.test(i, a, exact = FALSE, alternative = 'greater')$p.value 
        ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
        ica <- wilcox.test(ic, a, exact = FALSE, alternative = 'greater')$p.value 
        ci <- wilcox.test(c, i, exact = FALSE, alternative = 'greater')$p.value
        ici <- wilcox.test(ic, i, exact = FALSE, alternative = 'greater')$p.value
        icc <- wilcox.test(ic, c, exact = FALSE, alternative = 'greater')$p.value
        # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
          raw_pvals <- c(ia, ca, ica, ci, ici, icc)
          p.adjust(raw_pvals, method = "BH") 
          
        table(x$polyA_intron_canon)
     
        a <- subset(x, polyA_intron_canon == "Absent")$log2FoldChange
        i <- subset(x, polyA_intron_canon == "Intron_only")$log2FoldChange
        c <- subset(x, polyA_intron_canon == "Canon_only")$log2FoldChange
        ic <- subset(x, polyA_intron_canon == "Intron_Canon")$log2FoldChange
        ia <- wilcox.test(i, a, exact = FALSE, alternative = 'greater')$p.value 
        ca <- wilcox.test(c, a, exact = FALSE, alternative = 'greater')$p.value 
        ica <- wilcox.test(ic, a, exact = FALSE, alternative = 'greater')$p.value 
        ci <- wilcox.test(c, i, exact = FALSE, alternative = 'greater')$p.value
        ici <- wilcox.test(ic, i, exact = FALSE, alternative = 'greater')$p.value
        icc <- wilcox.test(ic, c, exact = FALSE, alternative = 'greater')$p.value
        # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
          raw_pvals <- c(ia, ca, ica, ci, ici, icc)
          p.adjust(raw_pvals, method = "BH")


#------------------------   
# Collect the distribution of total, enriched, and depleted gene_ids with site information
#------------------------   

  # PART 1: Filter for only enriched protein_coding genes
    # to make it easier
      x <- hrx_prot
      x <- hxa_prot
      x <- mrx_prot
      x <- mxa_prot        
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

    # for intron_equil
      i <- as.data.frame(table(x$intron_equil))
        names(i) <- c("Site", "All")
      iu <- as.data.frame(table(xu$intron_equil))
        names(iu) <- c("Site", "Sm-enriched")
      id <- as.data.frame(table(xd$intron_equil))
        names(id) <- c("Site", "Sm-depleted")
        # full_join into one table
          i <- full_join(i, iu, by = "Site")
          i <- full_join(i, id, by = "Site")
          
    # for intron_polyA
      j <- as.data.frame(table(x$intron_polyA))
        names(j) <- c("Site", "All")
      ju <- as.data.frame(table(xu$intron_polyA))
        names(ju) <- c("Site", "Sm-enriched")
      jd <- as.data.frame(table(xd$intron_polyA))
        names(jd) <- c("Site", "Sm-depleted")
        # full_join into one table
          j <- full_join(j, ju, by = "Site")
          j <- full_join(j, jd, by = "Site")

    # PART 3: Bind all together and write tsv
      x <- bind_rows(b, f, s, i)
      x <- bind_rows(b, f, s, j)
      
      h_rx <- x
      h_xa <- x
      m_rx <- x
      m_xa <- x

      write_tsv(as.data.frame(h_rx), file = "PATH TO SAVING DIRECTORY/h_equil_site_breakdown.tsv", quote = "none")
      write_tsv(as.data.frame(h_xa), file = "PATH TO SAVING DIRECTORY/h_atp_site_breakdown.tsv", quote = "none")
      write_tsv(as.data.frame(m_rx), file = "PATH TO SAVING DIRECTORY/m_equil_site_breakdown.tsv", quote = "none")
      write_tsv(as.data.frame(m_xa), file = "PATH TO SAVING DIRECTORY/m_atp_site_breakdown.tsv", quote = "none")
      

#------------------------
# Plots for exon number as a proxy for probability of intron retention
#------------------------  
      
# PART 1: import exon number tables Dr. Manu Sanjeev. These are only protein_coding.
      
    hdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/hdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    hdeseq_atp_PC_ExonNumbers <- read_delim("PATH TO/hdeseq_atp_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    mdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/mdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    mdeseq_atp_PC_ExonNumbers <- read_delim("PATH TO/mdeseq_atp_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      
# PART 2: CDF Plots
    
    x <- hdeseq_equil_PC_ExonNumbers
    x <- hdeseq_atp_PC_ExonNumbers
    x <- mdeseq_equil_PC_ExonNumbers
    x <- mdeseq_atp_PC_ExonNumbers
      
  # CDF plot
    x_plot <- ggplot(x, aes(log2FoldChange, color=ExonNumber)) +
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) anti-Sm RIP vs Cytoplasmic RNA", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      #coord_cartesian(xlim=c(-1,1)) + # for atp
      coord_cartesian(xlim=c(-4,4)) + # for equil
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
    
    hexn_plot <- x_plot
    haxn_plot <- x_plot
    mexn_plot <- x_plot
    maxn_plot <- x_plot
    
  # save plots
    ggsave(filename = "h_equil_exon_number_CDF.png", plot = hexn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "h_atp_exon_number_CDF.png", plot = haxn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "m_equil_exon_number_CDF.png", plot = mexn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "m_atp_exon_number_CDF.png", plot = maxn_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    
# PART 2: Wilcoxon Rank Sum Test and Benjamini-Hochberg Fasle discovery Rate correction for multiple testing
    
    x <- hdeseq_equil_PC_ExonNumbers
    x <- hdeseq_atp_PC_ExonNumbers
    x <- mdeseq_equil_PC_ExonNumbers
    x <- mdeseq_atp_PC_ExonNumbers
    
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
        
  h_2023_intronic <- read_delim("PATH TO featureCounts INTRONIC OUTPUT/h_2023_intronic.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  m_2023_intronic <- read_delim("PATH TO featureCounts INTRONIC OUTPUT/m_2023_intronic.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  x <- hrx_prot
  x <- hxa_prot
  x <- mrx_prot
  x <- mxa_prot
        
  x <- inner_join(x, h_2023_intronic, by = "ensembl_gene_id") 
  x <- inner_join(x, m_2023_intronic, by = "ensembl_gene_id") 
  
  x <- x %>% mutate(hx_intron_quintile = ntile(hx_intron, 5))  
  x <- x %>% mutate(hr_intron_quintile = ntile(hr_intron, 5))  
  x <- x %>% mutate(mx_intron_quintile = ntile(mx_intron, 5))  
  x <- x %>% mutate(mr_intron_quintile = ntile(mr_intron, 5))  
  
# PART 2: CDF Plot
        
  plot <- ggplot(x, aes(log2FoldChange, color= factor(h_cyto_intron_quintile))) +
    stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
    labs(x = "log2(f∆) anti-Sm RIP vs Human Cytoplasmic RNA", y = "Fraction of Data") +
    scale_y_continuous(position = "right") +
    #coord_cartesian(xlim=c(-1,1)) + # for atp
    coord_cartesian(xlim=c(-4,4)) + # for equil
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
      
