
## scripts to characterize smsite data from the idenfiication of smsites. 

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(gridExtra)
library(gtable)
library(ggpubr)
library(patchwork)

#import human_smsite_master and mouse_smsite_master
human_smsite_master <- read_delim("/filepath/human_smsite_master.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
# NA sometimes populates as an empty cell, NA is currently used to determine genes in which smsites are absent
human_smsite_master <- human_smsite_master %>% mutate(across(everything(), ~replace_na(., 'NA')))
  # numbers of smsites and smsite types
    table(human_smsite_master$smsite)
    table(human_smsite_master$smtype)
    
mouse_smsite_master <- read_delim("/filepath/mouse_smsite_master.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)
# NA sometimes populates as an empty cell, NA is currently used to determine genes in which smsites are absent
mouse_smsite_master <- mouse_smsite_master %>% mutate(across(everything(), ~replace_na(., 'NA')))
  # numbers of smsites and smsite types
    table(mouse_smsite_master$smsite)
    table(mouse_smsite_master$smtype)
  
# make the tables for bar graphs
  # switch between either master as the human data or mouse data to reuse the same code for analysis of both
    master <- human_smsite_master
    # master <- mouse_smsite_master
  # biotype for All, N/A, Absent, Noncanonical, and Canonical smsites
    h <- as.data.frame(table(master$gene_biotype))
      names(h) <- c("gene_biotype", "All")
    hna <-  master%>%filter(smtype == "NA")
      hna <- as.data.frame(table(hna$gene_biotype))
      names(hna) <- c("gene_biotype", "NA")
    hnc <-  master%>%filter(smtype == "NC")
      hnc <- as.data.frame(table(hnc$gene_biotype))
      names(hnc) <- c("gene_biotype", "NC")
    hu1 <-  master%>%filter(smtype == "U1")
      hu1 <- as.data.frame(table(hu1$gene_biotype))
      names(hu1) <- c("gene_biotype", "U1")
    hu2 <-  master%>%filter(smtype == "U2")
      hu2 <- as.data.frame(table(hu2$gene_biotype))
      names(hu2) <- c("gene_biotype", "U2")
    hu5 <-  master%>%filter(smtype == "U5")
      hu5 <- as.data.frame(table(hu5$gene_biotype))
      names(hu5) <- c("gene_biotype", "U5")
    hu7 <-  master%>%filter(smtype == "U7")
      hu7 <- as.data.frame(table(hu7$gene_biotype))
      names(hu7) <- c("gene_biotype", "U7")
    # now combine into a single table
      h_insilicobio <- full_join(h, hna, by = "gene_biotype")
      h_insilicobio <- full_join(h_insilicobio, hnc, by = "gene_biotype")
      h_insilicobio <- full_join(h_insilicobio, hu1, by = "gene_biotype")
      h_insilicobio <- full_join(h_insilicobio, hu2, by = "gene_biotype")
      h_insilicobio <- full_join(h_insilicobio, hu5, by = "gene_biotype")
      h_insilicobio <- full_join(h_insilicobio, hu7, by = "gene_biotype")
        write_tsv(as.data.frame(h_insilicobio), file = "/filepath/h_insilico_biotype.tsv", quote = "none")
        # write_tsv(as.data.frame(h_insilicobio), file = "/filepath/m_insilico_biotype.tsv", quote = "none")
  
  # smtype numbers within specific biotypes
    h_all <- as.data.frame((table(master$smtype)))
      names(h_all) <- c("smtype", "All")
    h_mrna <- master%>%filter(gene_biotype == "protein_coding") 
      h_mrna <- as.data.frame(table(h_mrna$smtype))
      names(h_mrna) <- c("smtype", "mRNA")
    h_lncrna <- master%>%filter(gene_biotype == "lncRNA") 
      h_lncrna <- as.data.frame(table(h_lncrna$smtype)) 
      names(h_lncrna) <- c("smtype", "lncRNA")
    h_scarna <- master%>%filter(gene_biotype == "scaRNA") 
      h_scarna <- as.data.frame(table(h_scarna$smtype)) 
      names(h_scarna) <- c("smtype", "scaRNA")
    h_snorna <- master%>%filter(gene_biotype == "snoRNA") 
      h_snorna <- as.data.frame(table(h_snorna$smtype)) 
      names(h_snorna) <- c("smtype", "snoRNA")
    h_snrna <- master%>%filter(gene_biotype == "snRNA") 
      h_snrna <- as.data.frame(table(h_snrna$smtype)) 
      names(h_snrna) <- c("smtype", "snRNA")
    # now combine into a single table
      h_smbio <- full_join(h_all, h_mrna, by = "smtype")
      h_smbio <- full_join(h_smbio, h_lncrna, by = "smtype")
      h_smbio <- full_join(h_smbio, h_scarna, by = "smtype")
      h_smbio <- full_join(h_smbio, h_snorna, by = "smtype")
      h_smbio <- full_join(h_smbio, h_snrna, by = "smtype")
        write_tsv(as.data.frame(h_smbio), file = "/filepath/h_smsite_biotype.tsv", quote = "none")
        # write_tsv(as.data.frame(h_smbio), file = "/filepath/m_smsite_biotype.tsv", quote = "none")
        
  
  # smtype frequency for regions
    # import output files
    # *make sure to swap human and mouse
      u1 <- read_delim("/filepath to u1_species.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      u2 <- read_delim("/filepath to u2_species.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      u5 <- read_delim("/filepath to u5_species.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      u7 <- read_delim("/filepath to u7_species.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      pe <- read_delim("/filepath to noncanonical_species.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      # only keep those transcript_ids in the master file
        u1f <- u1 %>% filter(ensembl_transcript_id %in% master$ensembl_transcript_id)
        u2f <- u2 %>% filter(ensembl_transcript_id %in% master$ensembl_transcript_id)
        u5f <- u5 %>% filter(ensembl_transcript_id %in% master$ensembl_transcript_id)
        u7f <- u7 %>% filter(ensembl_transcript_id %in% master$ensembl_transcript_id)
        pef <- pe %>% filter(ensembl_transcript_id %in% master$ensembl_transcript_id)
        # only keep protein coding genes
          u1fp <- u1f %>% filter(gene_biotype == "protein_coding")
          u2fp <- u2f %>% filter(gene_biotype == "protein_coding")
          u5fp <- u5f %>% filter(gene_biotype == "protein_coding")
          u7fp <- u7f %>% filter(gene_biotype == "protein_coding")
          pefp <- pef %>% filter(gene_biotype == "protein_coding")
        
    # now filter by smtype and table by region, only protein coding genes
    u1 <- as.data.frame(table(u1fp$region))
      names(u1) <- c("region", "U1" )
    u2 <- as.data.frame(table(u2fp$region))
      names(u2) <- c("region", "U2" )
    u5 <- as.data.frame(table(u5fp$region))
      names(u5) <- c("region", "U5" )
    u7 <- as.data.frame(table(u7fp$region))
      names(u7) <- c("region", "U7" )
    pe <- as.data.frame(table(pefp$region))
      names(pe) <- c("region", "NC" )
    # full_join by region
      region <- full_join(u1, u2, by = "region")
      region <- full_join(region, u5, by = "region")
      region <- full_join(region, u7, by = "region")
      region <- full_join(region, pe, by = "region")
        write_tsv(as.data.frame(region), file = "/filepath/h_smsite_region.tsv", quote = "none")
        # write_tsv(as.data.frame(region), file = "/filepath/m_smsite_region.tsv", quote = "none")
        
  
  
  # make density plots of transcript length (y) vs density (x), fill based on smtype or smsite designation
  # density plot showing distribution of rna lengths and whether they are different between canonical, noncanonical, or absent sm-sites.
    # smsite_master is either human_smsite_master or mouse_smsite_master
      NA <- smsite_master %>% filter(smsite == "NA")   
      NC <- smsite_master %>% filter(smsite == "NC")    
      U1 <- smsite_master %>% filter(smtype == "U1")      
      U2 <- smsite_master %>% filter(smtype == "U2") 
      U5 <- smsite_master %>% filter(smtype == "U5")
      U7 <- smsite_master %>% filter(smtype == "U7")
      c <- smsite_master %>% filter(smsite == "canon")
        NA_med <- log10(median(NA$transcript_length)) 
        NC_med <- log10(median(NC$transcript_length)) 
        U1_med <- log10(median(U1$transcript_length)) 
        U2_med <- log10(median(U2$transcript_length)) 
        U5_med <- log10(median(U5$transcript_length)) 
        U7_med <- log10(median(U7$transcript_length)) 
        c_med <- log10(median(c$transcript_length)) 
      
          tl_med <- data.frame(Sm_site = c("NA", "NC", "U1", "U2", "U5", "U7"),
                               median_transcript_length = c("NA_med", "NC_med", "U1_med", "U2_med", "U5_med", "U7_med"))
          tl_med1 <- data.frame(Sm_site = c("absent", "noncanon", "canon"),
                                median_transcript_length = c("NA_med", "NC_med", "c_med"))
      
    
    # remove U7
      no_u7 <- smsite_master %>% filter(smsite != "U7")
      
    # density plot for transcript length
      d_tl <- ggplot(data = no_u7, aes(y = log10(transcript_length), color = smsite)) + 
        geom_density(linewidth = 1.5) + 
        labs(y = "log10(transcript length)", x = "Density") +
        coord_cartesian(xlim = c(0,1.5), ylim = c(1,4.5)) +
        geom_hline(data = h_tl_med1, aes(yintercept = as.numeric(median_transcript_length), color = Sm_site), alpha = 1, linewidth = 1.25, linetype ="dotted")+
        # choose to plot either smsite or smtype
        #scale_color_manual("Sm-site", values = c("#440154","#3B528B", "#21908C", "#FDE725"), breaks = c("canon","U7", "noncanon", "absent"), labels = c("Canonical", "U7", "Noncanonical", "Absent")) +
        scale_color_manual("Sm-site", values = c("#440154","#21908C", "#FDE725"), breaks = c("canon", "noncanon", "absent"), labels = c("Canonical", "Noncanonical", "Absent")) +
        theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),
              legend.position = "bottom", 
              legend.direction = "vertical",
              legend.box.background = element_blank(), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))

    # frequency vs log10(transcript_length) scatter plot
      cs_tl <- ggplot(data = no_u7, aes(x=log10(transcript_length), y=canon_freq)) + 
        geom_point(alpha=0.35, color = "#440154") + 
        theme_minimal() +
        coord_cartesian(ylim=c(0,15)) +
        labs(y = "Canonical Frequency", x = "log10(transcript length)") +
        stat_cor(method = "kendall", label.x = 0.25, label.y = 12.5, size = 5, label.sep="\n")+
        theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
              panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.position.inside = c(0.17, 0.78), 
              legend.box.background = element_rect("white"), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))
  
      ns_tl <- ggplot(data = no_u7, aes(x=log10(transcript_length), y=nc_freq)) + 
        geom_point(alpha=0.35, color = "#21908C") + 
        theme_minimal() +
        coord_cartesian(ylim=c(0,100)) +
        labs(y = "Noncanonical Frequency", x = "log10(transcript length)") +
        stat_cor(method = "kendall", label.x = 0.25, label.y = 84, size = 5, label.sep="\n")+
        theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
              panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.position.inside = c(0.17, 0.78), 
              legend.box.background = element_rect("white"), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))
  
    # the same graphs for 3utr_length
      # only look at protein_coding
        p <- smsite_master%>%filter(gene_biotype == "protein_coding") 
          pNA <- p%>%filter(smtype == "NA")           
          pNC <- p%>%filter(smtype == "NC")     
          pU1 <- p%>%filter(smtype == "U1")        
          pU2 <- p%>%filter(smtype == "U2")
          pU5 <- p%>%filter(smtype == "U5")
          pU7 <- p%>%filter(smtype == "U7")
          pc <- p%>%filter(smsite == "canon")
            pNA_med <- log10(median(pNA$utr3_length))      
            pNC_med <- log10(median(pNC$utr3_length))      
            pU1_med <- log10(median(pU1$utr3_length))      
            pU2_med <- log10(median(pU2$utr3_length))      
            pU5_med <- log10(median(pU5$utr3_length))     
            pU7_med <- log10(median(pU7$utr3_length))      
            pc_med <- log10(median(pc$utr3_length))       
            
            u3l_med <- data.frame(Sm_site = c("NA", "NC", "U1", "U2", "U5", "U7"),
                                  median_utr3_length = c("pNA_med", "pNC_med", "pU1_med", "pU2_med", "pU5_med", "pU7_med"))
            u3l_med1 <- data.frame(Sm_site = c("absent", "noncanon", "canon"),
                                   median_utr3_length = c("pNA_med", "pNC_med", "pc_med"))
            
    # remove U7
      p_no_u7 <- p %>% filter(smsite != "U7")
            
    # Density plot, for utr3_length
      d_u3l <- ggplot(data = p_no_u7, aes(y = log10(utr3_length), color = smsite)) + 
        geom_density(linewidth = 1.5) + 
        coord_cartesian(xlim=c(0,1.25)) +
        labs(y = "log10(3'UTR length)", x = "Density") +
        geom_hline(data = h_u3l_med1, aes(yintercept = as.numeric(median_utr3_length), color = Sm_site), alpha = 1, linewidth = 1.25, linetype ="dotted")+
        # choose either to plot smsite or smtype
        #scale_color_manual("Sm-type", values = c("#440154", "#3B528B", "#21908C","#5DC863", "#FDE725", "grey"), breaks = c("U1", "U2", "U5", "U7", "NC", "NA"), labels = c("U1", "U2", "U5", "U7", "NC", "NA")) +
        scale_color_manual("Sm-site", values = c("#440154", "#21908C", "#FDE725"), breaks = c("canon", "noncanon", "absent"), labels = c("Canonical", "Noncanonical", "Absent")) +
        theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),
              legend.position = "bottom", 
              legend.direction = "vertical",
              legend.box.background = element_blank(), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))

    # Scatter plot for frequency of canonical sm-sites in 3'UTR vs log10(utr3_length)
      cs_u3l <- ggplot(data = p_no_u7, aes(x=log10(utr3_length), y= as.numeric(canon_3utr))) + 
        geom_point(alpha=0.35, color = "#440154") + 
        coord_cartesian(xlim=c(0,5), ylim=c(0,15)) +
        theme_minimal() +
        labs(y = "Canonical Frequency", x = "log10(3'UTR length)") +
        stat_cor(method = "kendall", label.x = 0.25, label.y = 12.5, size = 5, label.sep="\n")+
        theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
              panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.position.inside = c(0.17, 0.78), 
              legend.box.background = element_rect("white"), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))
      
      ns_u3l <- ggplot(data = p_no_u7, aes(x=log10(utr3_length), y= as.numeric(nc_3utr))) + 
        geom_point(alpha=0.35, color = "#21908C") + 
        coord_cartesian(xlim=c(0,5), ylim=c(0,100)) +
        theme_minimal() +
        labs(y = "Noncanonical Frequency", x = "log10(3'UTR length)") +
        stat_cor(method = "kendall", label.x = 0.25, label.y = 84, size = 5, label.sep="\n")+
        theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
              panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.position.inside = c(0.17, 0.78), 
              legend.box.background = element_rect("white"), 
              legend.key = element_blank(), 
              legend.text = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))
  

#compiling graphs for smsite frequency vs length
  f1 <- d_tl + plot_spacer() + cs_tl + plot_spacer() + ns_tl + plot_layout(widths = c(1.25,0.25,2.25,0.25,2.25), ncol = 5, nrow = 1)
  ggsave(filename = "insilico_transcript_length.png", plot = f1, device = "png", dpi = 300, height = 5, width = 9, path = "/filepath")
 
#compiling graphs for smsite frequency vs 3'UTR length
  sf1 <- d_u3l + plot_spacer() + cs_u3l + plot_spacer() + ns_u3l + plot_layout(widths = c(1,0.25,2.25,0.25,2.25), ncol = 5, nrow = 1)
  ggsave(filename = "insilico_3utr_length.png", plot = sf1, device = "png", dpi = 300, height = 5, width = 9, path = "/filepath")


# want to get the breakdown of Sm-sites found within each named snRNA group
  # filter for snRNA biotype
    snRNA <- smsite_master %>% filter(gene_biotype == "snRNA") 
      u1 <- snRNA %>% filter(!str_detect(external_gene_name, "U11")) 
        u1 <- u1 %>% filter(str_detect(external_gene_name, "U1")) 
        
      u2 <- snRNA %>% filter(str_detect(external_gene_name, "U2"))
      
      u4 <- snRNA %>% filter(!str_detect(external_gene_name, "U4ATAC")) 
        u4 <- u4 %>% filter(str_detect(external_gene_name, "U4"))
        
      u4atac <- snRNA %>% filter(str_detect(external_gene_name, "U4ATAC")) 
      
      u5 <- snRNA %>% filter(str_detect(external_gene_name, "U5"))
      
      u6 <- snRNA %>% filter(!str_detect(external_gene_name, "U6ATAC")) 
        u6 <- u6 %>% filter(str_detect(external_gene_name, "U6"))
      
      u6atac <- snRNA %>% filter(str_detect(external_gene_name, "U6ATAC")) 
      
      u7 <- snRNA %>% filter(str_detect(external_gene_name, "U7"))
      
      u11 <- snRNA %>% filter(str_detect(external_gene_name, "U11")) 
      
      u12 <- snRNA %>% filter(str_detect(external_gene_name, "U12")) 
    
 
