
### this is script to analyze sequencing performed in sma models and asking whether there is a change in abundance for mRNAs with sm-sites in sma
# cdf plots for smsite 
# filtering of downregulated transcripts shared in each and abundance of sm-sites.

# import libraries
  library(readr)
  library(stringr)
  library(tidyr)
  library(dplyr)
  library(tidyverse)
  library(biomaRt)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(gtable)
  library(patchwork)
  library(ggvenn)

# import datasets
  # doktor 
    # spinal cord 16691 obs PC1 51% disease, PC2 21%
      DESeqDoktorMice <- read_delim("/filepath/DESeqDoktorMice.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    # hela cells 12733 obs PC1 47% ?, PC2 45% ?
      DESeqDoktorHela <- read_delim("/filepath/DESeqDoktorHela.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    # liver 19496 obs PC1 39% ?, PC2 30% ?
      DESeqDoktorMiceLiver <- read_delim("/filepath/DESeqDoktorMiceLiver.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    # muscle 16853 obs PC1 70% disease, PC2 12%
      DESeqDoktorMiceMuscle <- read_delim("/filepath/DESeqDoktorMiceMuscle.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # maeda
    # 25218 obs PC1 99% disease, PC2 0%
      DESeqMaedaMice <- read_delim("/filepath/DESeqMaedaMice.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # nichterwitz
    # 19271 obs PC1 28% sample, PC2 17% disease
      NichterwitzDESeqSMA_MNs <- read_delim("/filepath/NichterwitzDESeqSMA_MNs.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # zhang
    # 19174 obs PC1 62% sample, PC2 25% disease
      DESeqSMA_MNs <- read_delim("/filepath/DESeqSMA_MNs.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # human/mouse_smsite_master 
    # 70611 obs
      human_smsite_master <- read_delim("/filepath/human_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      # remove U7 from smsite
        human_smsite_master_noU7 <- human_smsite_master %>% filter(smsite != "U7") # 70547 obs
    # 57186 obs
      mouse_smsite_master <- read_delim("/filepath/mouse_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      # remove U7 from smsite
        mouse_smsite_master_noU7 <- mouse_smsite_master %>% filter(smsite != "U7") # 57159 obs
      
# process tables
  # rename last column to enembl_gene_id
    names(DESeqDoktorMice)[8] <- "ensembl_gene_id"
    names(DESeqDoktorHela)[8] <- "ensembl_gene_id"
    names(DESeqDoktorMiceLiver)[7] <- "ensembl_gene_id"
    names(DESeqDoktorMiceMuscle)[7] <- "ensembl_gene_id"
    names(DESeqMaedaMice)[8] <- "ensembl_gene_id"
    names(NichterwitzDESeqSMA_MNs)[8] <- "ensembl_gene_id"
    names(DESeqSMA_MNs)[7] <- "ensembl_gene_id"
    
  # inner join with smsite master tables
    dsc <- inner_join(DESeqDoktorMice, mouse_smsite_master, by = "ensembl_gene_id") # 16689
      dsc <- dsc %>% filter(!(smsite == "U7"))
    dhe <- inner_join(DESeqDoktorHela, human_smsite_master, by = "ensembl_gene_id") # 12692
      dhe <- dhe %>% filter(!(smsite == "U7"))
    dlv <- inner_join(DESeqDoktorMiceLiver, mouse_smsite_master, by = "ensembl_gene_id") # 19496
      dlv <- dlv %>% filter(!(smsite == "U7"))
    dml <- inner_join(DESeqDoktorMiceMuscle, mouse_smsite_master, by = "ensembl_gene_id") # 16853
      dml <- dml %>% filter(!(smsite == "U7"))
    mmn <- inner_join(DESeqMaedaMice, mouse_smsite_master, by = "ensembl_gene_id") # 25216
      mmn <- mmn %>% filter(!(smsite == "U7"))
    nmn <- inner_join(NichterwitzDESeqSMA_MNs, mouse_smsite_master, by = "ensembl_gene_id") # 19269
      nmn <- nmn %>% filter(!(smsite == "U7"))
    zmn <- inmmnzmn <- inner_join(DESeqSMA_MNs, mouse_smsite_master, by = "ensembl_gene_id") # 19172
      zmn <- zmn %>% filter(!(smsite == "U7"))
      write_tsv(as.data.frame(dsc), file = "/filepath/DoktorSC/Doktor_SC_master.tsv", quote = "none")
      write_tsv(as.data.frame(dhe), file = "/filepath/Doktor_HeLa_master.tsv", quote = "none")
      write_tsv(as.data.frame(dlv), file = "/filepath/Doktor_LV_master.tsv", quote = "none")
      write_tsv(as.data.frame(dml), file = "/filepath/Doktor_ML_master.tsv", quote = "none")
      write_tsv(as.data.frame(mmn), file = "/filepath/Maeda_MN_master.tsv", quote = "none")
      write_tsv(as.data.frame(nmn), file = "/filepath/Nichterwitz_MN_master.tsv", quote = "none")
      write_tsv(as.data.frame(zmn), file = "/filepath/Zhang_MN_master.tsv", quote = "none")
      
    # cdf for all rnas 
      dsc_p <- ggplot(dsc, aes(log2FoldChange, color = smsite)) +
        stat_ecdf(geom = "step", linewidth = 1) +
        labs(x = "log2(f∆) SC: Taiwanese SMA vs Het", y = "Fraction of Data") +
        scale_y_continuous(position = "left") +
        coord_cartesian(xlim=c(-1,1)) +
        scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanonical", "Canonical"))+
        theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
              panel.grid = element_blank(),
              axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              #axis.title.x = element_blank(),
              legend.position = "top",
              legend.box.background = element_blank(), 
              legend.key = element_blank(), 
              legend.direction = "horizontal",
              legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
      table(dsc$smsite)
        # wilcox test 
          absent <- subset(dsc, smsite == "absent")$log2FoldChange
          noncanon <- subset(dsc, smsite == "noncanon")$log2FoldChange
          canon <- subset(dsc, smsite == "canon")$log2FoldChange
            w1 <- wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less')
            w2 <- wilcox.test(canon, absent, exact = FALSE, alternative = 'less')
            w3 <- wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
            
        dml_p <- ggplot(dml, aes(log2FoldChange, color = smsite)) +
              stat_ecdf(geom = "step", linewidth = 1) +
              labs(x = "log2(f∆) Muscle: Taiwanese SMA vs Het", y = "Fraction of Data") +
              scale_y_continuous(position = "left") +
              coord_cartesian(xlim=c(-1,1)) +
              scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
              theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                    panel.grid = element_blank(),
                    axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                    axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                    #axis.title.x = element_blank(),
                    legend.position = "top",
                    legend.box.background = element_blank(), 
                    legend.key = element_blank(), 
                    legend.direction = "horizontal",
                    legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                    legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
            table(dml$smsite)
            # wilcox test 
              absent <- subset(dml, smsite == "absent")$log2FoldChange
              noncanon <- subset(dml, smsite == "noncanon")$log2FoldChange
              canon <- subset(dml, smsite == "canon")$log2FoldChange
                wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                wilcox.test(canon, absent, exact = FALSE, alternative = 'less') 
                wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
          
          dlv_p <- ggplot(dlv, aes(log2FoldChange, color = smsite)) +
                  stat_ecdf(geom = "step", linewidth = 1) +
                  labs(x = "log2(f∆) Liver: Taiwanese SMA vs Het", y = "Fraction of Data") +
                  scale_y_continuous(position = "left") +
                  coord_cartesian(xlim=c(-1.5,1.5)) +
                  scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
                  theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                        panel.grid = element_blank(),
                        axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                        #axis.title.x = element_blank(),
                        legend.position = "top",
                        legend.box.background = element_blank(), 
                        legend.key = element_blank(), 
                        legend.direction = "horizontal",
                        legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
                table(dlv$smsite)
                # wilcox test 
                  absent <- subset(dlv, smsite == "absent")$log2FoldChange
                  noncanon <- subset(dlv, smsite == "noncanon")$log2FoldChange
                  canon <- subset(dlv, smsite == "canon")$log2FoldChange
                    wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
                    
            dhe_p <- ggplot(dhe, aes(log2FoldChange, color = smsite)) +
                  stat_ecdf(geom = "step", linewidth = 1) +
                  labs(x = "log2(f∆) Liver: Taiwanese SMA vs Het", y = "Fraction of Data") +
                  scale_y_continuous(position = "left") +
                  coord_cartesian(xlim=c(-1,1)) +
                  scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
                  theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                        panel.grid = element_blank(),
                        axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                        #axis.title.x = element_blank(),
                        legend.position = "top",
                        legend.box.background = element_blank(), 
                        legend.key = element_blank(), 
                        legend.direction = "horizontal",
                        legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                        panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
                table(dhe$smsite)
                  # wilcox test 
                    absent <- subset(dhe, smsite == "absent")$log2FoldChange
                    noncanon <- subset(dhe, smsite == "noncanon")$log2FoldChange
                    canon <- subset(dhe, smsite == "canon")$log2FoldChange
                        wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                        wilcox.test(canon, absent, exact = FALSE, alternative = 'less')
                        wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
                    
                
          mmn_p <- ggplot(mmn, aes(log2FoldChange, color = smsite)) +
                stat_ecdf(geom = "step", linewidth = 1) +
                labs(x = "log2(f∆) mESC MN: SMA vs Het", y = "Fraction of Data") +
                scale_y_continuous(position = "left") +
                coord_cartesian(xlim=c(-3,3)) +
                scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
                theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                      panel.grid = element_blank(),
                      axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                      #axis.title.x = element_blank(),
                      legend.position = "top",
                      legend.box.background = element_blank(), 
                      legend.key = element_blank(), 
                      legend.direction = "horizontal",
                      legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
                table(mmn$smsite)
                # wilcox test 
                  absent <- subset(mmn, smsite == "absent")$log2FoldChange
                  noncanon <- subset(mmn, smsite == "noncanon")$log2FoldChange
                  canon <- subset(mmn, smsite == "canon")$log2FoldChange
                    wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
                    
          nmn_p <- ggplot(nmn, aes(log2FoldChange, color = smsite)) +
                stat_ecdf(geom = "step", linewidth = 1) +
                labs(x = "log2(f∆) MN: ∆7SMA vs ∆7Het", y = "Fraction of Data") +
                scale_y_continuous(position = "left") +
                coord_cartesian(xlim=c(-2,2)) +
                scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
                theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                      panel.grid = element_blank(),
                      axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                      #axis.title.x = element_blank(),
                      legend.position = "top",
                      legend.box.background = element_blank(), 
                      legend.key = element_blank(), 
                      legend.direction = "horizontal",
                      legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
                table(nmn$smsite)
                # wilcox test 
                  absent <- subset(nmn, smsite == "absent")$log2FoldChange
                  noncanon <- subset(nmn, smsite == "noncanon")$log2FoldChange
                  canon <- subset(nmn, smsite == "canon")$log2FoldChange
                    wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, absent, exact = FALSE, alternative = 'less') 
                    wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
                    
          zmn_p <- ggplot(zmn, aes(log2FoldChange, color = smsite)) +
                stat_ecdf(geom = "step", linewidth = 1) +
                labs(x = "log2(f∆) MN: ∆7SMA vs ∆7Het", y = "Fraction of Data") +
                scale_y_continuous(position = "left") +
                coord_cartesian(xlim=c(-4,4)) +
                scale_color_manual("Sm-site:", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanon", "Canon"))+
                theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                      panel.grid = element_blank(),
                      axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                      #axis.title.x = element_blank(),
                      legend.position = "top",
                      legend.box.background = element_blank(), 
                      legend.key = element_blank(), 
                      legend.direction = "horizontal",
                      legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                      panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
                table(zmn$smsite)
                    # wilcox test 
                      absent <- subset(zmn, smsite == "absent")$log2FoldChange
                      noncanon <- subset(zmn, smsite == "noncanon")$log2FoldChange
                      canon <- subset(zmn, smsite == "canon")$log2FoldChange
                        wilcox.test(noncanon, absent, exact = FALSE, alternative = 'less') 
                        wilcox.test(canon, absent, exact = FALSE, alternative = 'less') 
                        wilcox.test(canon, noncanon, exact = FALSE, alternative = 'less') 
                    

    # repeat cdfs for protein_coding genes
      dscp <- dsc %>% filter(gene_biotype == "protein_coding") 
      dhep <- dhe %>% filter(gene_biotype == "protein_coding") 
      dmlp <- dml %>% filter(gene_biotype == "protein_coding") 
      dlvp <- dlv %>% filter(gene_biotype == "protein_coding") 
      mmnp <- mmn %>% filter(gene_biotype == "protein_coding") 
      nmnp <- nmn %>% filter(gene_biotype == "protein_coding") 
      zmnp <- zmn %>% filter(gene_biotype == "protein_coding")
      
    #repeat cdfs for lncRNAs
      dscl <- dsc %>% filter(gene_biotype == "lncRNA") 
      dhel <- dhe %>% filter(gene_biotype == "lncRNA") 
      dmll <- dml %>% filter(gene_biotype == "lncRNA") 
      dlvl <- dlv %>% filter(gene_biotype == "lncRNA") 
      mmnl <- mmn %>% filter(gene_biotype == "lncRNA") 
      nmnl <- nmn %>% filter(gene_biotype == "lncRNA") 
      zmnl <- zmn %>% filter(gene_biotype == "lncRNA") 
      
  
#plotting themes          
  #sma model cdf delineating by smsite
    design6 <- "
    11
    22
    33
    44
    55
    66
    "
  sma_cdf <- guide_area() + (mmnl_p + mmnp_p + plot_layout(axis_titles = "collect")) + plot_spacer() + (dscl_p + dscp_p + plot_layout(axis_titles = "collect")) + plot_spacer() + (dlvl_p + dlvp_p + plot_layout(axis_titles = "collect")) + plot_layout(design = design6, axis_titles = "collect", guides = "collect", heights = c(0.4,4,0.2,4,0.2,4))
  ggsave(filename = "sma_cdf.png", plot = sma_cdf, device = "png", dpi = 300, height = 15, width = 10, path = "/filepath")

# if we make the hypothesis that the affected RNAs should be reduced in SMA, lets make a list of gene_ids which are all reduced in SMA, then find overlap with smsite datasets.
  # first make a table containing only canon smsites
    canon <- mouse_smsite_master %>% filter(smsite == "canon") # 11129 obs
  # or should be physiologically associated with Sm-proteins
    # import mdeseq_equil
      mdeseq_equil <- read_delim("/filepath/mdeseq_equil.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
       # keep only things enriched
        physm <- mdeseq_equil %>% filter(log2FoldChange >= 0.6 & padj <= 0.5) # 5641 obs
  
  # focus only on Maeda, Doktor SC, and Doktor LV datasets
  #keep only those gene_ids corresponding to a log2F∆ < -0.6 & padj < 0.05
    dsc_DOWN <- dsc %>% filter(log2FoldChange <= 0 & padj <= 0.05) # 774 obs
    dlv_DOWN <- dlv %>% filter(log2FoldChange <= 0 & padj <= 0.05) # 2365 obs
    mmn_DOWN <- mmn %>% filter(log2FoldChange <= 0 & padj <= 0.05) # 6281 obs

  #make a list containing each of the columns above
    down_genes <- list('Doktor SC' = dsc_DOWN$ensembl_gene_id,
                     'Doktor LV' = dlv_DOWN$ensembl_gene_id,
                     'Maeda MN' = mmn_DOWN$ensembl_gene_id)

  #plot venn diagram
    sma_venn <- ggvenn(down_genes,
                     show_percentage = FALSE,
                     fill_color = c("#5dc86375", "#21908C75", "#44015475"),
                     stroke_size = 0,
                     set_name_color =c("#5dc863", "#21908C", "#440154"),
                     set_name_size = 6,
                     text_color = "white",
                     text_size = 6)
                        ggsave(filename = "sma_down_venn.png", plot = sma_venn, device = "png", dpi = 300, height = 4, width = 4, path = "/filepath")
    
   
    # make tables of the gene_ids shared by each comparison in mdeseq_equil
      # keep only gene_ids for each dataset
        dsc_id <- dsc_DOWN[,c(8)] # 774 obs
        dlv_id <- dlv_DOWN[,c(7)] # 2365 obs
        mmn_id <- mmn_DOWN[,c(8)] # 6281 obs
      sma_down <- inner_join(dsc_id, dlv_id, by = "ensembl_gene_id")
        sma_down <- inner_join(sma_down, mmn_id, by = "ensembl_gene_id")
        # inner_join sma_down with mouse_smsite_master
          sma_down_master <- inner_join(sma_down, mouse_smsite_master, by = "ensembl_gene_id")
          sma_down_equil <- inner_join(sma_down, mdeseq_equil, by = "ensembl_gene_id") 
            write_tsv(as.data.frame(sma_down), file = "/filepath/sma_down_SC_MN_LV.tsv", quote = "none")
            write_tsv(as.data.frame(sma_down_equil), file = "/filepath/sma_down_equil.tsv", quote = "none")
         
    # make volcano plots for sma_down
      # median log2FoldChange
        lg2a <- sma_down_equil %>% filter(smsite == "absent") 
        lg2n <- sma_down_equil %>% filter(smsite == "noncanon")
        lg2c <- sma_down_equil %>% filter(smsite == "canon")
          a_med <- median(lg2a$log2FoldChange)
          n_med <- median(lg2n$log2FoldChange) 
          c_med <- median(lg2c$log2FoldChange)
            lg2 <- data.frame(smsite = c("absent", "noncanon", "canon"),median_lg2fc = c("a_med", "n_med", "c_med"))
       
      # sma_down_equil
        sma_down_equil_volplot <- ggplot(data = sma_down_equil, aes(x = log2FoldChange, y = -log10(padj), col = smsite)) + 
          geom_point(alpha=1) + 
          labs(x = "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "-log10(p-adj)") +
          theme_minimal() +
          scale_color_manual(" ", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanonical", "Canonical"))+
          coord_cartesian(xlim = c(-3,3)) +
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
        plot(sma_down_equil_volplot)
        
        sma_down_equil_density <- ggplot(data = sma_down_equil, aes(x = log2FoldChange, color = smsite)) + 
          geom_density(linewidth = 1.25) +
          labs(x= "log2(f∆) anti-Sm RIP vs polyA-RNA", y = "Density") +
          coord_cartesian(xlim = c(-3,3), ylim = c(0,1)) +
          geom_vline(data = lg2, aes(xintercept = as.numeric(median_lg2fc), color = smsite), alpha = 1, linewidth = 1.25, linetype ="dotted")+
          scale_color_manual(" ", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanonical", "Canonical"))+
          theme(legend.position = "top",
                legend.box.background = element_blank(), 
                legend.key = element_blank(), 
                legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
                axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
                axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
                axis.title.x = element_blank(),
                panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),
                panel.grid.minor = element_blank()) +
          geom_vline(xintercept = c(-0.6, 0.6), col="black", linetype = "dashed")
        plot(sma_down_equil_density)
        
        #plotting theme for volcano and density plot for down-regulated genes in SMA
        design <- "
      3
      1
      2
      "
        sde <- ((sma_down_equil_volplot/sma_down_equil_density) + plot_layout(heights = c(2,2))) + guide_area() + plot_layout(heights = c(0.5,4,1.5), design = design, guides = "collect", axis_titles = "collect")
        ggsave(filename = "sma_down_equil_smsite_volcano_density_cdf.png", plot = sde, device = "png", dpi = 300, height = 6, width = 5, path = "/filepath")
        
