
#--------------------------
# Process ATP vs no ATP conditions for hybrid species Sm-ring assembly reactions
# this corresponds to 2025 samples treated with 2M urea, 5 mg/mL heparin and washed in RSB-1000 + 0.4% NP-40
#--------------------------

library(tidyverse)
library(pheatmap)
library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(gtable)
library(patchwork)
library(ggvenn)

#--------------------------
# Import ATP vs no ATP condition, treated with 2M urea, 5mg/mL heparin, and 1M NaCl washes, captured with Y12 antibody
#--------------------------

# PART 1: Combine deseq2 comparison output with mouse_smsite_master
hxmr_vs_hxmra_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/hxmr_vs_hxmra_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # remove unnecessary columns and remove 'NA'
    m_atp <- hxmr_vs_hxmra_deseq_results[,c(1:3,7)]
    m_atp <- na.omit(m_atp) 
  
    # inner_join with mouse_smsite_master
    mouse_smsite_master <- read_delim("PATH TO/mouse_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
      m_atp <- inner_join(m_atp, mouse_smsite_master, by = "ensembl_gene_id") 

      
# PART 2: Add columns for separating out different sites
      
  x <- m_atp 
      
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
    x_i_w_c <- x %>% filter(intron_polyA == "present" & smsite == "canon") 
    x_i_wo_c <- x %>% filter(intron_polyA == "present" & smsite != "canon") 
    x_c_wo_i <- x %>% filter(intron_polyA != "present" & smsite == "canon") 
    x_no_i_no_c <- x %>% filter(intron_polyA != "present" & smsite != "canon")
    x <- x %>% mutate(polyA_intron_canon = case_when(ensembl_gene_id %in% x_i_w_c$ensembl_gene_id ~ "Intron_Canon", ensembl_gene_id %in% x_i_wo_c$ensembl_gene_id ~ "Intron_only", ensembl_gene_id %in% x_c_wo_i$ensembl_gene_id ~ "Canon_only", TRUE ~ "Absent"))
      
  m_atp <- x
      
  # save master table
    write_tsv(as.data.frame(m_atp), file = "PATH TO SAVING DIRECTORY/m_atp.tsv", quote = "none")
    
# PART 3: Import mdeseq_atp datasets for RSB-500 + 0.1% NP-40 data, and make a column for ATP-enriched and Canonical Sm-site containing gene_ids
    
    mdeseq_atp <- read_delim("PATH TO/mdeseq_atp.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
    
    # filter for ATP-enriched
      m_atp_up <- mdeseq_atp %>% filter(log2FoldChange >= 0.6 & padj <= 0.05 & smsite == "canon" & gene_biotype == "protein_coding") 
    
    # Annotate columns for ATP-enriched protein_coding genes that have Canonical Sm-sites
      m_atp <- m_atp %>% mutate(candidate = case_when(ensembl_gene_id %in% m_atp_up$ensembl_gene_id ~"present", TRUE ~"absent"))
    
    # save master table
      write_tsv(as.data.frame(hB), file = "PATH TO SAVING DIRECTORY/m_atp.tsv", quote = "none")
    
# PART 4: Filter to keep protein_coding and remove smsite U7
  
  m_atp_prot <- m_atp %>% filter(gene_biotype == "protein_coding") 
  m_atp_prot <- m_atp_prot %>% filter(smsite != "U7")
      
#--------------------------
# CDF plot
#--------------------------      
# PART 1: Plots  
  x <- m_atp_prot
      
  #bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt)) +
  #fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss)) +
  #int_plot <- ggplot(x, aes(log2FoldChange, color=intron_polyA)) +
  #sm_plot <- ggplot(x, aes(log2FoldChange, color=smsite)) +
  #sm_bpt_plot <- ggplot(x, aes(log2FoldChange, color=bpt_canon)) +
  #sm_fss_plot <- ggplot(x, aes(log2FoldChange, color=five_ss_canon)) +
  sm_int_plot <- ggplot(x, aes(log2FoldChange, color=polyA_intron_canon)) +
  #cand_plot <- ggplot(x, aes(log2FoldChange, color=candidate)) +
  stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
    labs(x = "log2(f∆) anti-Sm RIP: ATP vs No ATP", y = "Fraction of Data") +
    scale_y_continuous(position = "right") +
    #coord_cartesian(xlim=c(-2,2)) + 
    coord_cartesian(xlim=c(-1,1)) + # for canon vs ___
    #scale_color_manual("Bpt:", values = c("#2828cc", "#36648B", "#63B8FF"), breaks = c("No Bpt", "consensus", "stringent"), labels = c("No Bpt", "Consensus", "Stringent"))+ 
    #scale_color_manual("5'ss:", values = c("#2828cc", "#B80069", "#F97415"), breaks = c("No 5'ss", "5'ss_bulge", "5'ss"), labels = c("No 5'ss", "5'ss_bulge", "5'ss"))+ 
    #scale_color_manual("Intron:", values = c("#a9a2a2", "#ecba00"), breaks = c("absent", "present"), labels = c("No Intron", "Intron"))+
    #scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#440154"), breaks = c("absent", "noncanon", "canon" ), labels = c("Absent", "Noncanonical", "Canonical")) + 
    #scale_color_manual("Site", values = c("#7c7887", "#63B8FF", "#440154", "#f0199e"), breaks = c("Absent", "Bpt_only", "Canon_only", "Bpt_Canon" ), labels = c("Absent", "Bpt_only", "Canonical-Sm only", "Bpt & Canonical-Sm")) + 
    #scale_color_manual("Site", values = c("#4C4545", "#F97415", "#440154", "#CE2727"), breaks = c("Absent", "5'ss_only", "Canon_only", "5'ss_Canon" ), labels = c("Absent", "5'ss only", "Canonical-Sm only", "5'ss & Canonical-Sm")) + 
    scale_color_manual("Site:", values = c("#a9a2a2", "#ecba00","#440154", "#3bbf5e"), breaks = c("Absent", "Intron_only", "Canon_only", "Intron_Canon"), labels = c("Absent", "Intron_only", "Canon_only", "Intron & Canon"))+
    #scale_color_manual("Candidate", values = c("#FDE725", "#440154"), breaks = c("absent", "present" ), labels = c("Absent", "Candidate")) + 
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

      
  # for saving plots
    ggsave(filename = "matp_Y12_Bpt_CDF.png", plot = bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_5ss_CDF.png", plot = fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Int_CDF.png", plot = int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_CDF.png", plot = sm_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_Bpt_CDF.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_5ss_CDF.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_Int_CDF.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Candidate_CDF.png", plot = cand_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_Bpt_CDF_coor1.png", plot = sm_bpt_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_5ss_CDF_coor1.png", plot = sm_fss_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "matp_Y12_Sm_Int_CDF_coor1.png", plot = sm_int_plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    
    
# PART 2: Perform Wilcoxon Rank Sum Test with Continuity Correction and Benjamini-Hochberg False Discovery Rate corrections for multiple testing
    
  x <- m_atp_prot
    
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
  ba <- wilcox.test(s, a, exact = FALSE, alternative = 'greater')$p.value 
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
    
  table(x$candidate)
  
  a <- subset(x, candidate == "absent")$log2FoldChange
  p <- subset(x, candidate == "present")$log2FoldChange
    wilcox.test(p, a, exact = FALSE, alternative = 'greater')$p.value 
    
    
#------------------------   
# Collect the distribution of total, enriched, and depleted gene_ids with site information
#------------------------   
    
# PART 1: Filter for only enriched protein_coding genes

  x <- m_atp_prot
         
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
    i <- as.data.frame(table(x$intron_polyA))
      names(i) <- c("Site", "All")
    iu <- as.data.frame(table(xu$intron_polyA))
      names(iu) <- c("Site", "Sm-enriched")
    id <- as.data.frame(table(xd$intron_polyA))
      names(id) <- c("Site", "Sm-depleted")
      # full_join into one table
        i <- full_join(i, iu, by = "Site")
        i <- full_join(i, id, by = "Site")
    
# Part 3: Bind all together and write tsv
    x <- bind_rows(b, f, s, i)
    
    matp <- x
    
    write_tsv(as.data.frame(matp), file = "PATH TO SAVING DIRECTORY/matp_Y12_enriched_site_breakdown.tsv", quote = "none")

    
#------------------------   
# Create volcano, density, and CDF plots for RNA species
#------------------------     
    
# PART 1: Make a column grouping different RNA species
  
  x <- m_atp
    
  # filter by biotype, then create new column
    sn <- x %>% filter(gene_biotype == "snRNA") 
    m <- x %>% filter(gene_biotype == "protein_coding") 
    lnc <- x %>% filter(gene_biotype == "lncRNA") 
    sca <- x %>% filter(gene_biotype == "scaRNA") 
    sno <- x %>% filter(gene_biotype == "snoRNA") 
    snc <- bind_rows(sca, sno) 
    x <- x%>%mutate(rna_species = case_when(ensembl_gene_id %in% sn$ensembl_gene_id ~ "snRNA",
                                            ensembl_gene_id %in% m$ensembl_gene_id ~ "mRNA",
                                            ensembl_gene_id %in% lnc$ensembl_gene_id ~ "lncRNA",
                                            ensembl_gene_id %in% snc$ensembl_gene_id ~ "sncRNA",
                                            TRUE ~ "Other"))       
    
    
# PART 2: Volcano Plot, Density Plot, & CDF Plot
  # Volcano
    volcano <- ggplot(data = x, aes(x = log2FoldChange, y = -log10(padj), col = rna_species)) + 
      geom_point(alpha=0.35) + 
      labs(x = "log2(f∆) anti-Sm RIP: ATP vs No ATP", y = "-log10(p-adj)") +
      theme_minimal() +
      scale_color_manual("Biotype:", values = c("#FDE725", "#5DC863", "#21908c", "#3B528B", "#440154"), breaks = c("Other", "sncRNA", "snRNA", "lncRNA", "mRNA"), labels = c("Other", "sno/scaRNA", "snRNA", "lncRNA", "mRNA"))+
      coord_cartesian(xlim = c(-4,4)) +
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
      labs(x= "log2(f∆) anti-Sm RIP: ATP vs Npo ATP", y = "Density") +
      coord_cartesian(xlim = c(-4,4)) +
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
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) anti-Sm RIP: vs ATP vs No ATP", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      coord_cartesian(xlim=c(-2,2)) +
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
    v <- volcano
    d <- density
    c <- cdf
    
    design <- "
      11
      24
      35
      "
    p <- guide_area() + v + d + c + plot_spacer() + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect", axis_titles = "collect")
    ggsave(filename = "matp_Y12_volcano_density.png", p = p, device = "png", dpi = 300, height = 6, width = 8, path = "PATH TO SAVING DIRECTORY")
    
    
#------------------------   
# Overlap between ATP-enriched RNAs between 0.5M NaCl and 2M urea, 5mg/mL heparin, 1M NaCl with Y12 antibody
#------------------------       
  
# PART 1: define datasets to filter for only ATP-enriched RNAs
    
  # from 2M urea, 5mg/mL heparin, 1M NaCl
    hxmr_vs_hxmra_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/hxmr_vs_hxmra_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
  # from 0.5M NaCl
    hx_mr_vs_hx_mr_a_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT/hx+mr_vs_hx+mr+a_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
    
  # remove unnecessary columns and remove 'NA'
    m_atp <- hx_mr_vs_hx_mr_a_deseq_results[,c(1:3,7)] 
    m_atp <- na.omit(m_atp)
    names(m_atp)[1] <- "ensembl_gene_id"
    
    m_atp_1 <- hxmr_vs_hxmra_deseq_results[,c(1:3,7)] 
    m_atp_1 <- na.omit(m_atp_1) 
    
  # filter for only enriched gene_ids
    m_atp_up <- m_atp %>% filter(log2FoldChange >= 0.6 & padj <= 0.05) 
    m_atp_1_up <- m_atp_1 %>% filter(log2FoldChange >= 0.06 & padj <= 0.05) 
    
  # make a list of Canonical Sm-site gene_ids from mouse_smsite_master
    m_canon <- mouse_smsite_master %>% filter(smsite == "canon") 
    
  # inner_join deseq analyses with mouse_smsite_master
    m_atp_up <- inner_join(m_atp_up, mouse_smsite_master, by = "ensembl_gene_id")
    m_atp_1_up <- inner_join(m_atp_1_up, mouse_smsite_master, by = "ensembl_gene_id")
    
  # for protein_coding genes only
    m_atp_up_prot <- m_atp_up %>% filter(gene_biotype == "protein_coding") 
    m_atp_1_up_prot <- m_atp_1_up %>% filter(gene_biotype == "protein_coding")
    m_canon_prot <- m_canon %>% filter(gene_biotype == "protein_coding")
    
# PART 2: Venn Plotting
    
  # convert to list for venn plotting
    vn <- list("Standard" = m_atp_up$ensembl_gene_id,
               "Urea/Heparin/1M NaCl" = m_atp_1_up$ensembl_gene_id,
               "Canonical-Sm" = m_canon$ensembl_gene_id)
    
    prot_vn <- list("Standard" = m_atp_up_prot$ensembl_gene_id,
                    "Urea/Heparin/1M NaCl" = m_atp_1_up_prot$ensembl_gene_id,
                    "Canonical-Sm" = m_canon_prot$ensembl_gene_id)
    
  # plot venn diagram
    #venn <- ggvenn(vn,
    venn <- ggvenn(prot_vn,        
                   show_percentage = FALSE,
                   fill_color = c("#20A48675","#3B528B75", "#44015475"),
                   stroke_size = 0,
                   set_name_color = c("#20A486", "#3B528B", "#440154"),
                   set_name_size = 6,
                   text_color = "white",
                   text_size = 6)
    
    ggsave(filename = "all_enriched_rna_venn.png", plot = venn, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    ggsave(filename = "mRNA_enriched_venn.png", plot = venn, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    

#------------------------
# Plots for exon number as a proxy for probability of intron retention
#------------------------  
    
# PART 1: import exon number tables Dr. Manu Sanjeev. These are only protein_coding. Annotate Y12 and SmB/B'/N datasets
    
  mdeseq_equil_PC_ExonNumbers <- read_delim("PATH TO/mdeseq_equil_PC_ExonNumbers.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # select only enseml_gene_id and ExonNumber columns
    en <- mdeseq_equil_PC_ExonNumbers %>% dplyr::select(ensembl_gene_id, ExonNumber) 
    
    matp_prot <- inner_join(m_atp_prot, en, by = "ensembl_gene_id")

# PART 2: CDF Plot
    
    plot <- ggplot(matp_prot, aes(log2FoldChange, color=ExonNumber)) +
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) anti-Sm RIP vs Cytoplasmic RNA", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      coord_cartesian(xlim=c(-2,2)) + 
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
    
    # save plots
      ggsave(filename = "matp_Y12_ExonNumber_CDF.png", plot = plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")

# PART 3: Wilcoxon Rank Sum Test and Benjamini-Hochberg Fasle discovery Rate correction for multiple testing
    
    x <- matp_prot
    
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
    
  m_2023_intronic <- read_delim("PATH TO featureCounts INTRONIC OUTPUT/m_2023_intronic.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  m_atp_prot <- inner_join(m_atp_prot, m_2023_intronic, by = "ensembl_gene_id") 
  
  m_atp_prot <- m_atp_prot %>% mutate(mr_intron_quintile = ntile(mr_intron, 5))  
    
# PART 2: CDF Plot
    
  plot <- ggplot(m_atp_prot, aes(log2FoldChange, color= factor(mr_intron_quintile))) +
    stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
    labs(x = "log2(f∆) anti-Sm RIP: ATP vs no ATP", y = "Fraction of Data") +
    scale_y_continuous(position = "right") +
    coord_cartesian(xlim=c(-2,2)) + 
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
    ggsave(filename = "matp_Y12_IntronQuintile_CDF.png", plot = plot, device = "png", dpi = 300, height = 4, width = 4, path = "PATH TO SAVING DIRECTORY")
    
# PART 3: Wilcoxon Rank Sum Test and Benjamini-Hochberg False discovery Rate correction for multiple testing
    
  table(m_atp_prot$mr_intron_quintile)
    
    one <- subset(m_atp_prot, mr_intron_quintile == "1")$log2FoldChange
    two <- subset(m_atp_prot, mr_intron_quintile == "2")$log2FoldChange
    three <- subset(m_atp_prot, mr_intron_quintile == "3")$log2FoldChange
    four <- subset(m_atp_prot, mr_intron_quintile == "4")$log2FoldChange
    five <- subset(m_atp_prot, mr_intron_quintile == "5")$log2FoldChange
    w1 <- wilcox.test(two, one, exact = FALSE, alternative = 'greater')$p.value 
    w2 <- wilcox.test(three, two, exact = FALSE, alternative = 'greater')$p.value 
    w3 <- wilcox.test(four, three, exact = FALSE, alternative = 'greater')$p.value
    w4 <- wilcox.test(five, four, exact = FALSE, alternative = 'greater')$p.value
    # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
    raw_pvals <- c(w1, w2, w3, w4)
    p.adjust(raw_pvals, method = "BH")    
    
    
    
      
