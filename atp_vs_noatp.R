


###this is a script for processing the ATP vs No ATP data 
# for mouse B(hxmr) vs C(hxmra)
# for human E(mxhr) vs F(mxhra)

#load libraries
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(gridExtra)
library(gtable)
library(patchwork)

#Human
#import deseq analysis for Sm-RIP: ATP vs noATP comparison
  filename <- read_delim("filepath to filename.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  #prepare tables for joining with smsite tables and remove raw data values, keeping only gene_id, baseMean, log2FoldChange, and padj
    deseq_atp <-filename[,c(1:3,7)] 
      names(deseq_atp)[1] <- "ensembl_gene_id"
      deseq_atp <- na.omit(deseq_atp) 
#import smsite master table
  smsite_master <- read_delim("/filepath to smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  smsite_master <- smsite_master %>% mutate(across(everything(), ~replace_na(., 'NA')))
  #inner_join deseq analysis with smsite master tables
    deseq_atp <- inner_join(deseq_atp, human_smsite_master, by = "ensembl_gene_id") 
   
  #lets get the descriptive analysis done here for biotypes and sm_sites of total and of those enriched
    #make a table for the biotype of rnas for types of smsites
      #for human total classification of RNAs 
        table(deseq_atp$smsite) 
    #now for rnas enriched >= 0.6 log2FoldChange
      deseq_atp_UP <- deseq_atp%>%filter(log2FoldChange >= 0.6 & padj <= 0.05) 
      deseq_atp_DOWN <- deseq_atp%>%filter(log2FoldChange <= -0.6 & padj <= 0.05) 
        table(deseq_atp_UP$smsite)
        table(deseq_atp_DOWN$smsite)
    #make master table with each
      ba <- as.data.frame(table(deseq_atp$gene_biotype)) 
        names(ba) <- c("gene_biotype", "All")
      baU <- as.data.frame(table(deseq_atp_UP$gene_biotype))
        names(hbaU) <- c("gene_biotype", "All_UP")
      baD <- as.data.frame(table(deseq_atp_DOWN$gene_biotype))
        names(baD) <- c("gene_biotype", "All_DOWN")

      baa <- deseq_atp%>%filter(smsite == "absent")    
        baa <- as.data.frame(table(baa$gene_biotype))
        names(baa) <- c("gene_biotype","Absent")
      baaU <- deseq_atp_UP%>%filter(smsite == "absent")    
        baaU <- as.data.frame(table(baaU$gene_biotype))
        names(baaU) <- c("gene_biotype","Absent_UP")
      baaD <- deseq_atp_DOWN%>%filter(smsite == "absent")   
        baaD <- as.data.frame(table(baaD$gene_biotype))
        names(baaD) <- c("gene_biotype","Absent_DOWN")

      ban <- deseq_atp%>%filter(smsite == "noncanon") 
        ban <- as.data.frame(table(ban$gene_biotype))
        names(ban) <- c("gene_biotype", "Noncanonical")
      banU <- deseq_atp_UP%>%filter(smsite == "noncanon") 
        banU <- as.data.frame(table(banU$gene_biotype))
        names(banU) <- c("gene_biotype", "Noncanonical_UP")
      banD <- deseq_atp_DOWN%>%filter(smsite == "noncanon") 
        banD <- as.data.frame(table(banD$gene_biotype))
        names(banD) <- c("gene_biotype", "Noncanonical_DOWN")

      bac <- deseq_atp%>%filter(smsite == "canon") 
        bac <- as.data.frame(table(bac$gene_biotype))
        names(bac) <- c("gene_biotype", "Canonical")
      bacU <- deseq_atp_UP%>%filter(smsite == "canon") 
        bacU <- as.data.frame(table(bacU$gene_biotype))
        names(bacU) <- c("gene_biotype", "Canonical_UP")
      bacD <- deseq_atp_DOWN%>%filter(smsite == "canon") 
        bacD <- as.data.frame(table(bacD$gene_biotype))
        names(bacD) <- c("gene_biotype", "Canonical_DOWN")

      #combine into a single table
        atpbio <- full_join(ba, baa, by = "gene_biotype")
          atpbio <- full_join(atpbio, ban, by = "gene_biotype")
          atpbio <- full_join(atpbio, bac, by = "gene_biotype")
          atpbio <- full_join(atpbio, baU, by = "gene_biotype")
          atpbio <- full_join(atpbio, baaU, by = "gene_biotype")
          atpbio <- full_join(atpbio, banU, by = "gene_biotype")
          atpbio <- full_join(atpbio, bacU, by = "gene_biotype")
          atpbio <- full_join(atpbio, baD, by = "gene_biotype")
          atpbio <- full_join(atpbio, baaD, by = "gene_biotype")
          atpbio <- full_join(atpbio, banD, by = "gene_biotype")
          atpbio <- full_join(atpbio, bacD, by = "gene_biotype")
        write_tsv(as.data.frame(atpbio), file = "/filepath/atpbio.tsv", quote = "none")
   
# reannotate with another column for gene_biotype that combines sno/scaRNA and can plot from single column
  snRNA <- smsite_master%>%filter(gene_biotype == "snRNA") 
  mRNA <- smsite_master%>%filter(gene_biotype == "protein_coding") 
  lncRNA <- smsite_master%>%filter(gene_biotype == "lncRNA") 
  scaRNA <- smsite_master%>%filter(gene_biotype == "scaRNA") 
  snoRNA <- smsite_master%>%filter(gene_biotype == "snoRNA") 
  sncRNA <- bind_rows(scaRNA, snoRNA) 
  deseq_atp <- deseq_atp%>%mutate(rna_class = case_when(ensembl_gene_id %in% snRNA$ensembl_gene_id ~ "snRNA",
                                                        ensembl_gene_id %in% mRNA$ensembl_gene_id ~ "mRNA",
                                                        ensembl_gene_id %in% lncRNA$ensembl_gene_id ~ "lncRNA",
                                                        ensembl_gene_id %in% sncRNA$ensembl_gene_id ~ "sncRNA",
                                                        TRUE ~ "Other",))

# cdf plots of anti-Sm-RIP: ATP vs no ATP, color coding by Sm-site
  #for all RNAs
  #exclude U7 sm-sites
    deseq_atp_noU7 <- deseq_atp%>%filter(smsite != "U7")
    #cdf 
      c_atp_sm <- ggplot(deseq_atp_noU7, aes(log2FoldChange, color = smsite)) +
        stat_ecdf(geom = "step", linewidth = 1.25) +
        labs(x = "log2(F∆) Sm-RIP: ATP vs no ATP", y = "Fraction of Data") +
        scale_y_continuous(position = "left") +
        coord_cartesian(xlim=c(-1,1)) +
        # use either line for plotting by Sm-site type or combined Sm-site
        #scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#3B528B", "#44015450", "#44015475","#440154"), breaks = c("NA", "NC", "U7", "U5", "U2", "U1" ), labels = c("Absent", "Noncanonical", "U7", "U5", "U2", "U1")) + 
        scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#440154"), breaks = c("absent", "noncanon", "canon" ), labels = c("Absent", "Noncanonical", "Canonical")) + 
        theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
              panel.grid = element_blank(),
              axis.text = element_text( size = rel(1.25), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              axis.title.x = element_blank(),
              legend.position = "top",
              legend.box.background = element_blank(), 
              legend.key = element_blank(), 
              legend.direction = "horizontal",
              legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
              legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
      # gather number of events for smsite
      table(deseq_atp_noU7$smsite)
        #wilcox test 
          absent <- subset(deseq_atp_noU7, smsite == "absent")$log2FoldChange
          noncanon <- subset(deseq_atp_noU7, smsite == "noncanon")$log2FoldChange
          canon <- subset(deseq_atp_noU7, smsite == "canon")$log2FoldChange
            w1 <- wilcox.test(noncanon, absent, exact = FALSE, alternative = 'greater') 
            w2 <- wilcox.test(canon, absent, exact = FALSE, alternative = 'greater') 
            w3 <- wilcox.test(canon, noncanon, exact = FALSE, alternative = 'greater') 
      # gather number of events for smsite type  
        table(hdeseq_atp$smtype)
        # wilcox test
          na <- subset(hdeseq_atp, smtype == "NA")$log2FoldChange
          nc <- subset(hdeseq_atp, smtype == "NC")$log2FoldChange
          u7 <- subset(hdeseq_atp, smtype == "U7")$log2FoldChange
          u5 <- subset(hdeseq_atp, smtype == "U5")$log2FoldChange
          u2 <- subset(hdeseq_atp, smtype == "U2")$log2FoldChange
          u1 <- subset(hdeseq_atp, smtype == "U1")$log2FoldChange
            w1 <- wilcox.test(nc, na, exact = FALSE, alternative = 'greater') 
            w2 <- wilcox.test(u7, na, exact = FALSE, alternative = 'greater') 
            w3 <- wilcox.test(u5, na, exact = FALSE, alternative = 'greater')
            w4 <- wilcox.test(u2, na, exact = FALSE, alternative = 'greater') 
            w5 <- wilcox.test(u1, na, exact = FALSE, alternative = 'greater') 

  # repeat for:
    #for protein_coding genes
      deseq_atp_prot <- deseq_atp%>%filter(gene_biotype == "protein_coding") 
      deseq_atp_prot_noU7 <- deseq_atp_prot%>%filter(smsite != "U7") 
    
    #for lncRNAs
      deseq_atp_lncRNA <- deseq_atp%>%filter(gene_biotype == "lncRNA") 
      deseq_atp_lncRNA_noU7 <- deseq_atp_lncRNA%>%filter(smsite != "U7")
   
    #for sncRNAs (sno/scaRNA)
      deseq_atp_sncRNA <- deseq_atp%>%filter(rna_class == "sncRNA")
      deseq_atp_sncRNA_noU7 <- deseq_atp_sncRNA%>%filter(smsite != "U7") 
    
    #for snRNAs
      deseq_atp_snRNA <- deseq_atp%>%filter(rna_class == "snRNA") 
      deseq_atp_snRNA_noU7 <- deseq_atp_snRNA%>%filter(smsite != "U7") 
      
  #to plot cdfs for frequency of sm-sites vs anti-Sm-RIP log2(F∆): ATP vs no ATP
    #Only polyA-RNAs with smsites, frequency of canonical smsites
    # obtain breakdown of number of genes that contain a specified amount of canonical smsites
      table(deseq_atp_noU7$canon_freq)
      #group by 0, 1, 2, 3+
        deseq_atp_polyA <- deseq_atp_noU7%>%filter(gene_biotype %in% c("protein_coding", "lncRNA"))
        deseq_atp_noA <- deseq_atp_polyA%>%filter(smsite != "absent") 
        deseq_atp_noA$canon_freq[deseq_atp_noA$canon_freq %in% c(3:80)] <- '3+'
      #cdf
        c_atp_ccf <- ggplot(deseq_atp_noA, aes(log2FoldChange, color = canon_freq)) +
          stat_ecdf(geom = "step", linewidth = 1.25) +
          labs(x = "log2(F∆) Sm-RIP: ATP vs no ATP", y = "Fraction of Data") +
          scale_y_continuous(position = "left") +
          coord_cartesian(xlim=c(-1,1)) +
          scale_color_manual("Canonical Frequency", values = c("#21908C", "#31699E", "#443A83", "#440154"), breaks = c("0", "1", "2", "3+"), labels = c("0", "1", "2", "3+"))+
          theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                panel.grid = element_blank(),
                axis.text = element_text( size = rel(1.25), family = "Arial", face = "bold", color = "black"),
                axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                legend.box.background = element_blank(), 
                legend.key = element_blank(), 
                legend.direction = "vertical",
                legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
                legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                legend.position = "top",
                panel.grid.major = element_line(colour = "gray90", linetype = "dashed")) +
          guides(color=guide_legend(nrow=1, byrow=TRUE))
        # gather number of events for smsite type
          table(hdeseq_atp_noA$canon_freq)
          #wilcox test
            zero <- subset(deseq_atp_noA, canon_freq == "0")$log2FoldChange
            one <- subset(deseq_atp_noA, canon_freq == "1")$log2FoldChange
            two <- subset(deseq_atp_noA, canon_freq == "2")$log2FoldChange
            three <- subset(deseq_atp_noA, canon_freq == "3+")$log2FoldChange
              w1 <- wilcox.test(three, two, exact = FALSE, alternative = 'greater') 
              w2 <- wilcox.test(two, one, exact = FALSE, alternative = 'greater') 
              w3 <- wilcox.test(one, zero, exact = FALSE, alternative = 'greater') 
    # repeat for:
      # only polyA-RNAs with canonical Sm-sites, frequency of noncanonical Sm-sites
        deseq_atp_canon <- deseq_atp_noA%>%filter(smsite == "canon") 
        table(deseq_atp_canon$nc_freq)
      
      # only RNA, frequency of noncanonical Sm-sites
        deseq_atp_noncanon <- deseq_atp_polyA%>%filter(smsite != "canon") 
        table(deseq_atp_noncanon$nc_freq)
      
  # plot for transcript and UTR length
  # scatter plot for log10(transcript_length) vs log2(f∆)
  as_tl <- ggplot(data=deseq_atp, aes(x=log2FoldChange, y=log10(transcript_length))) + 
    #geom_point(alpha=0.25, color = "black") + 
    geom_pointdensity(alpha = 0.5) +
    scale_color_viridis_c(name = "Transcripts") +
    theme_minimal() +
    coord_cartesian(xlim=c(-2,3), ylim =c(1,5)) +
    labs(y = "log10(transcript length)", x = "log2(f∆) anti-Sm-RIP: ATP vs no ATP") +
    #stat_cor(method = "pearson", label.x = -1.95, label.y = 4.5, size = 5, label.sep="\n")+
    theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
          panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
          axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
          legend.position = c(0.825, 0.775), 
          legend.box.background = element_rect(fill = "white"), 
          legend.key = element_blank(), 
          legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
          legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"))
  # Pearson correlation test
  r <-cor.test(log10(deseq_atp$transcript_length), deseq_atp$log2FoldChange, method="pearson") 

  #repeat for 3utr_length in protein_coding genes
  deseq_atp_prot <- deseq_atp%>%filter(gene_biotype == "protein_coding") 
    
#plotting themes for figures
  #plotting theme for Figure 5
    design1 <- "
          12
          "
    f5 <- c_atp_sm_lncrna + c_atp_sm_prot + plot_layout(design = design1, axis_titles = "collect")
    ggsave(filename = "deseq_atp_sm_lnc_m_cdf.png", plot = f5, device = "png", dpi = 300, height = 6, width = 10, path = "/filepath")
    
  #plotting theme for Supplement Figure 5 
    design2 <- "
          123
          "
    sf5 <- c_atp_sm + c_atp_sm_sncrna + c_atp_sm_snrna + plot_layout(design = design2, axis_titles = "collect")
    ggsave(filename = "deseq_atp_sm_cdf_sf.png", plot = sf5, device = "png", dpi = 300, height = 6, width = 15, path = "/filepath")
   
    # for smtype
      design_type <- "
        123
        "
        at <- (c_atp_sm + c_atp_sm_lncrna + c_atp_sm_prot) + plot_layout(design = design_type, axis_titles = "collect")
        ggsave(filename = "deseq_atp_smtype_cdf.png", plot = at, device = "png", dpi = 300, height = 6, width = 14.5, path = "/filepath")
    
  #plotting theme for Supplement Figure 6
    design3 <- "
        123
        "
    sf6 <- c_atp_ccf + c_atp_cnf + c_atp_ncf + plot_layout(design = design3, axis_titles = "collect")
    ggsave(filename = "deseq_atp_freq_cdf.png", plot = sf6, device = "png", dpi = 300, height = 6, width = 14, path = "/filepath")
     
  #plot theme for Supplement Figure 7
    sf7 <- as_tl + plot_spacer() + as_u3l + plot_layout(widths = c(2.25,0.25,2.25), ncol = 3, nrow = 1)
    ggsave(filename = "atp_lengths.png", plot = sf7, device = "png", dpi = 300, height = 5, width = 10, path = "/filepath")
