

###Analysis of anti-Sm-RIP vs polyA-RNA 
# for mouse: M(mr) vs D(mx), both aligned to mouse genome
# for human: H(hr) vs A(hx), both aligned to human genome

#load R libraries
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

#import smsite master sheet
  smsite_master <- read_delim("filepath to smsite_master.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
  # ensure NA is 'NA' and not an empty cell
    smsite_master <- smsite_master %>% mutate(across(everything(), ~replace_na(., 'NA')))

#import Deseq2 output for Sm-RIP vs PolyA-RNA comparison
  deseq_results <- read_delim("filepath to deseq_results.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
  # remove unnecessary columns
    deseq_equil <- deseq_results[,c(1:3,7)] 
  # rename gene_id column
    names(deseq_equil)[1] <- "ensembl_gene_id"
  # omit any rows with NA
    deseq_equil <- na.omit(deseq_equil) 
   
  #inner_join deseq analysis with smsite master, becomes master file used to generate the following figures
  deseq_equil <- inner_join(deseq_equil, smsite_master, by = "ensembl_gene_id") 
    write_tsv(as.data.frame(deseq_equil), file = "/filepath/deseq_equil.tsv", quote = "none")

#lets get the descriptive analysis done here for biotypes and sm_sites of total and of those enriched
  #make a table for the biotype of rnas for types of smsites
    #for human total classification of RNAs 
      table(deseq_equil$smsite) 
      table(hdeseq_equil$smtype)
      
    #now for rnas enriched >= 0.6 log2FoldChange
      deseq_equil_UP <- deseq_equil%>%filter(log2FoldChange >= 0.6 & padj <= 0.05) 
      deseq_equil_DOWN <- deseq_equil%>%filter(log2FoldChange <= -0.6 & padj <= 0.05) 
      table(deseq_equil_UP$smsite)
      table(deseq_equil_DOWN$smsite)
        
  #make master table for biotypes and smsites
    be <- as.data.frame(table(deseq_equil$gene_biotype)) 
      names(be) <- c("gene_biotype", "All")
    beU <- as.data.frame(table(deseq_equil_UP$gene_biotype)) 
      names(beU) <- c("gene_biotype", "All_UP")
    beD <- as.data.frame(table(deseq_equil_DOWN$gene_biotype)) 
      names(beD) <- c("gene_biotype", "All_DOWN")

    bea <- deseq_equil%>%filter(smsite == "absent")    
      bea <- as.data.frame(table(bea$gene_biotype))
      names(bea) <- c("gene_biotype","absent")
    beaU <- deseq_equil_UP%>%filter(smsite == "absent")   
      beaU <- as.data.frame(table(beaU$gene_biotype))
      names(beaU) <- c("gene_biotype","Absent_UP")
    beaD <- deseq_equil_DOWN%>%filter(smsite == "absent")    
      beaD <- as.data.frame(table(beaD$gene_biotype))
      names(beaD) <- c("gene_biotype","Absent_DOWN")

    ben <- deseq_equil%>%filter(smsite == "noncanon")
      ben <- as.data.frame(table(ben$gene_biotype))
      names(ben) <- c("gene_biotype", "Noncanonical")
    benU <- deseq_equil_UP%>%filter(smsite == "noncanon") 
      benU <- as.data.frame(table(benU$gene_biotype))
      names(benU) <- c("gene_biotype", "Noncanonical_UP")
    benD <- deseq_equil_DOWN%>%filter(smsite == "noncanon") 
      benD <- as.data.frame(table(benD$gene_biotype))
      names(benD) <- c("gene_biotype", "Noncanonical_DOWN")

    bec <- deseq_equil%>%filter(smsite == "canon") 
      bec <- as.data.frame(table(bec$gene_biotype))
      names(bec) <- c("gene_biotype", "Canonical")
    becU <- deseq_equil_UP%>%filter(smsite == "canon") 
      becU <- as.data.frame(table(becU$gene_biotype))
      names(becU) <- c("gene_biotype", "Canonical_UP")
    becD <- deseq_equil_DOWN%>%filter(smsite == "canon")
      becD <- as.data.frame(table(becD$gene_biotype))
      names(becD) <- c("gene_biotype", "Canonical_DOWN")
      
  #combine into a single table
    equilbio <- full_join(be, bea, by = "gene_biotype")
      equilbio <- full_join(equilbio, ben, by = "gene_biotype")
      equilbio <- full_join(equilbio, bec, by = "gene_biotype")
      equilbio <- full_join(equilbio, beU, by = "gene_biotype")
      equilbio <- full_join(equilbio, beaU, by = "gene_biotype")
      equilbio <- full_join(equilbio, benU, by = "gene_biotype")
      equilbio <- full_join(equilbio, becU, by = "gene_biotype")
      equilbio <- full_join(equilbio, beD, by = "gene_biotype")
      equilbio <- full_join(equilbio, beaD, by = "gene_biotype")
      equilbio <- full_join(equilbio, benD, by = "gene_biotype")
      equilbio <- full_join(equilbio, becD, by = "gene_biotype")
    write_tsv(as.data.frame(equilbio), file = "/filepath/equil_biotype.tsv", quote = "none")


#annotate to assign a new column with gene_biotype names that combines sca/sno RNAs into a single group
  snRNA <- smsite_master%>%filter(gene_biotype == "snRNA") 
  mRNA <- smsite_master%>%filter(gene_biotype == "protein_coding") 
  lncRNA <- smsite_master%>%filter(gene_biotype == "lncRNA") 
  scaRNA <- smsite_master%>%filter(gene_biotype == "scaRNA") 
  snoRNA <- smsite_master%>%filter(gene_biotype == "snoRNA") 
  sncRNA <- bind_rows(scaRNA, snoRNA) # 1071 obs
  deseq_equil <- deseq_equil%>%mutate(rna_class = case_when(ensembl_gene_id %in% snRNA$ensembl_gene_id ~ "snRNA",
                                                              ensembl_gene_id %in% mRNA$ensembl_gene_id ~ "mRNA",
                                                              ensembl_gene_id %in% lncRNA$ensembl_gene_id ~ "lncRNA",
                                                              ensembl_gene_id %in% sncRNA$ensembl_gene_id ~ "sncRNA",
                                                              TRUE ~ "Other"))

  #volcano and density plots of anti-Sm-RIP vs polyA-RNA, color coding gene_biotype--renamed rna_class
  #volcano plot 
  v_equil_rclass <- ggplot(data = deseq_equil, aes(x = log2FoldChange, y = -log10(padj), col = rna_class)) + 
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

  #density for volcano plot
  d_equil_rclass <- ggplot(data = deseq_equil, aes(x = log2FoldChange, color = rna_class)) + 
    geom_density(linewidth = 1.25) +
    labs(x= "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "Density") +
    coord_cartesian(xlim = c(-10,15), ylim = c(0,0.4)) +
    scale_color_manual(" ", values = c("#FDE725", "#5DC863", "#21908c", "#3B528B", "#440154"), breaks = c("Other", "sncRNA", "snRNA", "lncRNA", "mRNA"), labels = c("Other", "sno/scaRNA", "snRNA", "lncRNA", "mRNA"))+
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.key = element_blank(), 
          #legend.box = "horizontal",
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
  
  # cdf plots of anti-Sm-RIP vs polyA-RNA,, color coding by Sm-site
    # for all RNAs
      # cdf *make sure to switch between smsite and smtype
      # remove U7 smsite containing rnas
      deseq_equil_noU7 <- deseq_equil %>% filter(smsite != "U7")
      c_equil_sm <- ggplot(deseq_equil_noU7, aes(log2FoldChange, color = smtype)) +
        stat_ecdf(geom = "step", linewidth = 1.25) +
        labs(x = "log2(f∆) Sm-RIP/polyA-RNA", y = "Fraction of Data") +
        scale_y_continuous(position = "left") +
        coord_cartesian(xlim=c(-4,4)) +
        #choose to plot by smsite or smsite type
        #scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#440154"), breaks = c("absent", "noncanon", "canon" ), labels = c("Absent", "Noncanonical", "Canonical")) + 
        scale_color_manual("Sm-site", values = c("#FDE725", "#21908C", "#3B528B", "#44015450", "#44015475","#440154"), breaks = c("NA", "NC", "U7", "U5", "U2", "U1" ), labels = c("Absent", "Noncanical", "U7", "U5", "U2", "U1")) + 
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
      #obtain contribution to plot
      table(deseq_equil_noU7$smsite)
        # absent    canonical   noncanonical 
      # wilcox test 
        absent <- subset(deseq_equil_noU7, smsite == "absent")$log2FoldChange
        noncanon <- subset(deseq_equil_noU7, smsite == "noncanon")$log2FoldChange
        canon <- subset(deseq_equil_noU7, smsite == "canon")$log2FoldChange
          wilcox.test(noncanon, absent, exact = FALSE, alternative = 'greater') 
          wilcox.test(canon, absent, exact = FALSE, alternative = 'greater') 
          wilcox.test(canon, noncanon, exact = FALSE, alternative = 'greater')
      #repeat for smtype
      table(hdeseq_equil$smtype)
      # wilcox test
        na <- subset(deseq_equil, smtype == "NA")$log2FoldChange
        nc <- subset(deseq_equil, smtype == "NC")$log2FoldChange
        u7 <- subset(deseq_equil, smtype == "U7")$log2FoldChange
        u5 <- subset(deseq_equil, smtype == "U5")$log2FoldChange
        u2 <- subset(deseq_equil, smtype == "U2")$log2FoldChange
        u1 <- subset(deseq_equil, smtype == "U1")$log2FoldChange
          wilcox.test(nc, na, exact = FALSE, alternative = 'greater') 
          wilcox.test(u7, na, exact = FALSE, alternative = 'greater') 
          wilcox.test(u5, na, exact = FALSE, alternative = 'greater') 
          wilcox.test(u2, na, exact = FALSE, alternative = 'greater') 
          wilcox.test(u1, na, exact = FALSE, alternative = 'greater') 
          
          
    # repeat the above cdf plot and statistical tests for each of the following, inserting the dataframe into the data variable of the plot function
    #for protein_coding genes
      deseq_equil_prot <- deseq_equil %>% filter(gene_biotype == "protein_coding") 
      deseq_equil_prot_noU7 <- deseq_equil_prot %>% filter(smsite != "U7")
     
    #for lncRNAs
      deseq_equil_lncRNA <- deseq_equil%>%filter(gene_biotype == "lncRNA") 
      deseq_equil_lncRNA_noU7 <- deseq_equil_lncRNA %>% filter(smsite != "U7")
     
    #for sncRNAs (sno/scaRNA)
      deseq_equil_sncRNA <- deseq_equil %>% filter(rna_class == "sncRNA")
      deseq_equil_sncRNA_noU7 <- deseq_equil_sncRNA %>% filter(smsite != "U7") 
      
    #for snRNAs
      deseq_equil_snRNA <- deseq_equil%>%filter(rna_class == "snRNA") 
      deseq_equil_snRNA_noU7 <- deseq_equil_snRNA%>%filter(smsite != "U7") 
        
  #cdf plots for frequency of sm-sites vs anti-Sm-RIP vs polyA-RNA log2(F∆)
  #Only RNAs with smsites, frequency of canonical smsites
    #obtain number of genes that have a specific number of smsites
    table(deseq_equil_noU7$canon_freq)
    #group by 0, 1, 2, 3+
      deseq_equil_noA <- deseq_equil_noU7%>%filter(smsite != "absent") 
      deseq_equil_noA$canon_freq[hdeseq_equil_noA$canon_freq %in% c(3:80)] <- '3+'
      #cdf
        c_equil_ccf <- ggplot(deseq_equil_noA, aes(log2FoldChange, color = canon_freq)) +
          stat_ecdf(geom = "step", linewidth = 1.25) +
          labs(x = "log2(f∆) anti-Sm-RIP vs polyA-RNA", y = "Fraction of Data") +
          scale_y_continuous(position = "left") +
          coord_cartesian(xlim=c(-3,3)) +
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
        #obtain number contributing to plots
        table(deseq_equil_noA$canon_freq)
          #wilcox test
            zero <- subset(deseq_equil_noA, canon_freq == "0")$log2FoldChange
            one <- subset(deseq_equil_noA, canon_freq == "1")$log2FoldChange
            two <- subset(deseq_equil_noA, canon_freq == "2")$log2FoldChange
            three <- subset(deseq_equil_noA, canon_freq == "3+")$log2FoldChange
              w1 <- wilcox.test(three, two, exact = FALSE, alternative = 'greater') 
              w2 <- wilcox.test(two, one, exact = FALSE, alternative = 'greater') 
              w3 <- wilcox.test(one, zero, exact = FALSE, alternative = 'greater') 
    #repeat the above cdf plot for the following:
      #only RNAs with canonical Sm-sites, frequency of noncanonical Sm-sites
        deseq_equil_canon <- deseq_equil_noU7%>%filter(smsite == "canon") # 8322
        table(deseq_equil_canon$nc_freq)
        #group by 0, 1-2, 3-6, 7-13, 14+
          deseq_equil_canon$nc_freq[deseq_equil_canon$nc_freq %in% c(1:5)] <- '1-5'
          deseq_equil_canon$nc_freq[deseq_equil_canon$nc_freq %in% c(6:10)] <- '6-10'
          deseq_equil_canon$nc_freq[deseq_equil_canon$nc_freq %in% c(11:15)] <- '11-15'
          deseq_equil_canon$nc_freq[deseq_equil_canon$nc_freq %in% c(16:764)] <- '16+'
        
      #only RNA, frequency of noncanonical Sm-sites
        deseq_equil_noncanon <- deseq_equil_noU7 %>% filter(smsite != "canon") 
        table(deseq_equil_noncanon$nc_freq)
        #group by 0, 1, 2-3, 4-8, 9+
          deseq_equil_noncanon$nc_freq[deseq_equil_noncanon$nc_freq %in% c(2:3)] <- '2-3'
          deseq_equil_noncanon$nc_freq[deseq_equil_noncanon$nc_freq %in% c(4:8)] <- '4-8'
          deseq_equil_noncanon$nc_freq[deseq_equil_noncanon$nc_freq %in% c(9:77)] <- '9+'
          
   #scatter plot for log10(transcript_length) vs log2(f∆)
    es_tl <- ggplot(data=deseq_equil, aes(x=log2FoldChange, y=log10(transcript_length))) + 
      #geom_point(alpha=0.25, color = "black") + 
      geom_pointdensity(alpha = 0.5) +
      scale_color_viridis_c(name = "Transcripts") +
      theme_minimal() +
      coord_cartesian(xlim=c(-5,10), ylim =c(1,6)) +
      labs(y = "log10(transcript length)", x = "log2(f∆) anti-Sm-RIP vs polyA-RNA") +
      #stat_cor(method = "pearson", label.x = -4, label.y = 5.5, size = 5, label.sep="\n")+
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
    #pearson correlation test
    r <-cor.test(log10(deseq_equil$transcript_length), deseq_equil$log2FoldChange, method="pearson") 
    
    #repeat for 3utr_length in protein_coding genes
    es_u3l <- ggplot(data=deseq_equil_prot, aes(x=log2FoldChange, y= log10(utr3_length))) + 
      #geom_point(alpha=0.25, color = "black") + 
      geom_pointdensity(alpha = 0.5) +
      scale_color_viridis_c(name = "Transcripts") +
      theme_minimal() +
      coord_cartesian(xlim=c(-5,10), ylim =c(0.5,5)) +
      labs(y = "log10(3'UTR length)", x = "log2(f∆) anti-Sm-RIP vs polyA-RNA") +
      #stat_cor(method = "pearson", label.x = -4, label.y = 4.5, size = 5, label.sep="\n")+
      theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
            panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.x = element_blank(),
            axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            #legend.position.inside = c(0.17, 0.78), 
            legend.position = c(0.825, 0.775), 
            legend.box.background = element_rect(fill = "white"), 
            legend.key = element_blank(), 
            legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"))
    # remove rows with '0' utr3_length
      esu3l <- deseq_equil_prot %>% filter(utr3_length != 0)
    # pearson correlation test
    r <-cor.test(log10(esu3l$utr3_length), esu3l$log2FoldChange, method="pearson") 



#plotting themes for figures
    # plotting theme for Figure 4
    design <- "
      3
      1
      2
      "
    f4 <- v_equil_rclass + d_equil_rclass + guide_area() + plot_layout(heights = c(1,6,2), design = design, guides = "collect", axis_titles = "collect")
    ggsave(filename = "deseq_equil_biotype_volcano_density_cdf.png", plot = f4, device = "png", dpi = 300, height = 6, width = 5, path = "/filepath")
    
  #plotting theme for Figure 5
    design1 <- "
      12345
      "
    f5 <- (c_equil_sm + c_equil_sm_sncrna + c_equil_sm_snrna + c_equil_sm_lncrna + c_equil_sm_prot) + plot_layout(design = design1, axis_titles = "collect")
    ggsave(filename = "deseq_equil_sm_cdf.png", plot = f5, device = "png", dpi = 300, height = 6, width = 22, path = "/filepath")
     
  #plotting theme for Supplement Figure 3
    design2 <- "
        123
        " 
    sf3 <- (c_equil_ccf + c_equil_cnf + c_equil_ncf) + plot_layout(design = design2, axis_titles = "collect")
    ggsave(filename = "deseq_equil_freq_cdf.png", plot = sf3, device = "png", dpi = 300, height = 6, width = 14, path = "/filepath")
    
  #plotting theme for Supplement Figure 4
    sf4 <- es_tl + plot_spacer() + es_u3l + plot_layout(widths = c(2.25,0.25,2.25), ncol = 3, nrow = 1)
    ggsave(filename = "equil_lengths.png", plot = sf4, device = "png", dpi = 300, height = 5, width = 10, path = "/filepath")

      
      
      
      
      
      
      
