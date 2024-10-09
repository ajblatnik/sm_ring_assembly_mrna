

## this is a script used to analize the SmB-CLIP-Seq data, Nijssen axon/soma seq data, and Lu 2014 Y12 RIP data used for supplementary figures.

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
  # SmB-CLIP-Seq
    SmBclipMediumvsStringentDESeq <- read_delim("/filepath/SmBclipMediumvsStringentDESeq.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # nijssen 2018
    njiissen_2018_soma_vs_axon <- read_delim("/filepath/njiissen_2018_soma_vs_axon.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # smsite_master
    mouse_smsite_master <- read_delim("/filepath/mouse_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    human_smsite_master <- read_delim("/filepath/human_smsite_master.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# analysis for SmB-CLIP-Seq data
  # rename Human_gene_stable_id to ensembl_gene_id
    SmBclipMediumvsStringentDESeq <- SmBclipMediumvsStringentDESeq[,c(1:7)]
    names(SmBclipMediumvsStringentDESeq)[7] <- "ensembl_gene_id"
  # inner join SmB-CLIP-Seq with mouse_smsite_master
    smb <- inner_join(SmBclipMediumvsStringentDESeq, human_smsite_master, by = "ensembl_gene_id") # 11283 obs
    write_tsv(as.data.frame(smb), file = "/filepath/SmB-CLIP-Seq_master.tsv", quote = "none")
    
    # remove U7
      smb <- smb %>% filter(smsite != "U7")
    # make lists to analyze for different biotypes
      smb_snRNA <- smb %>% filter(gene_biotype == "snRNA") 
      smb_mRNA <- smb %>% filter(gene_biotype == "protein_coding") 
      smb_lncRNA <- smb %>% filter(gene_biotype == "lncRNA") 
      
    # plot cdfs 
      # all rnas
        smb_p <- ggplot(smb, aes(log2FoldChange, color = smsite)) +
          stat_ecdf(geom = "step", linewidth = 1) +
          labs(x = "log2(f∆) SmB-CLIP: Medium vs Stringent", y = "Fraction of Data") +
          scale_y_continuous(position = "left") +
          coord_cartesian(xlim=c(-10,10)) +
          scale_color_manual(" ", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanonical", "Canonical"))+
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
        # obtain numbers contributing to each function
        table(smb$smsite)
        # wilcox test 
          absent <- subset(smb, smsite == "absent")$log2FoldChange
          noncanon <- subset(smb, smsite == "noncanon")$log2FoldChange
          canon <- subset(smb, smsite == "canon")$log2FoldChange
            w1 <- wilcox.test(noncanon, absent, exact = FALSE, alternative = 'greater')
            w2 <- wilcox.test(canon, absent, exact = FALSE, alternative = 'greater') 
            w3 <- wilcox.test(canon, noncanon, exact = FALSE, alternative = 'greater') 
        # repeat for snRNA, mRNA, and lncRNA
     
      # compare rna's enriched in SmB-CLIP-Seq and Sm-RIP 
          smb_3rdquart_enriched <- smb %>% filter(log2FoldChange >= 2.8467 & padj <= 0.05) 
        # keep only those gene_ids corresponding to a log2F∆ >= 0.6 & padj < 0.05 enriched with sm-proteins
          hequil_enriched <- hdeseq_equil %>% filter(log2FoldChange >= 0.6 & padj <= 0.05) 
        # inner_join tables
          smb_3rd_equil <- inner_join(smb_3rdquart_enriched, hee, by = "ensembl_gene_id") # 549 obs, 578 obs
            write_tsv(as.data.frame(smb_3rd_equil), file = "/filepath/smb_3rd_equil.tsv", quote = "none")

          # make a list containing each of the columns above
            third <- list('SmB-CLIP' = smb_3rdquart_enriched$ensembl_gene_id,
                          'Sm-enriched' = hequil_enriched$ensembl_gene_id) 
                
          # plot venn diagram
            third_venn <- ggvenn(third,
                                  show_percentage = FALSE,
                                  fill_color = c("#5dc86375", "#44015475"),
                                  stroke_size = 0,
                                  set_name_color =c("#5dc863", "#440154"),
                                  set_name_size = 6,
                                  text_color = "white",
                                  text_size = 6)
            ggsave(filename = "smb_3rdquart_noATP_venn.png", plot = third_venn, device = "png", dpi = 300, height = 4, width = 4, path = "/filepath/")
            
              
# analysis for nijssen axon vs soma seq
  # rename Mouse_gene_name to external_gene_name
    names(njiissen_2018_soma_vs_axon)[1] <- "external_gene_name"
  # inner_join with other datasets
    # mouse_smsite_master
      nijssen <- inner_join(njiissen_2018_soma_vs_axon, mouse_smsite_master, by = "external_gene_name") # 9848 obs
    # make lists of enriched gene_id from mdeseq_equil 
      equil_enriched <- mdeseq_equil %>% filter(log2FoldChange >= 0.6 & padj <= 0.05)
    # make columns with a Sm-enriched or ATP-enriched yay or nay
      nijssen <- nijssen %>% mutate(Sm = case_when(ensembl_gene_id %in% equil_enriched$ensembl_gene_id ~"Sm-enriched", TRUE ~"Un-enriched"))
      nijssen <- nijssen %>% mutate(ATP = case_when(ensembl_gene_id %in% atp_enriched$ensembl_gene_id ~"ATP-enriched", TRUE ~"Un-enriched"))
      write_tsv(as.data.frame(nijssen), file = "/filepath/nijssen_2018_master.tsv", quote = "none")

    # remove U7
      nijssen <- nijssen %>% filter(smsite != "U7")
    # cdf for smsite
      nijssen_p <- ggplot(nijssen, aes(log2FoldChange, color = smsite)) +
        stat_ecdf(geom = "step", linewidth = 1) +
        labs(x = "log2(f∆) Soma vs Axon", y = "Fraction of Data") +
        scale_y_continuous(position = "left") +
        coord_cartesian(xlim=c(-6,10)) +
        scale_color_manual(" ", values = c("#FDE725", "#21908c", "#440154"), breaks = c("absent", "noncanon", "canon"), labels = c("Absent", "Noncanonical", "Canonical"))+
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
      table(nijssen$smsite)
      # wilcox test 
        absent <- subset(nijssen, smsite == "absent")$log2FoldChange
        noncanon <- subset(nijssen, smsite == "noncanon")$log2FoldChange
        canon <- subset(nijssen, smsite == "canon")$log2FoldChange
          w1 <- wilcox.test(noncanon, absent, exact = FALSE, alternative = 'greater') 
          w2 <- wilcox.test(canon, absent, exact = FALSE, alternative = 'greater') 
          w3 <- wilcox.test(canon, noncanon, exact = FALSE, alternative = 'greater') 
          
    # cdf for equil
      nijssen_equil_p <- ggplot(nijssen, aes(log2FoldChange, color = Sm)) +
            stat_ecdf(geom = "step", linewidth = 1) +
            labs(x = "log2(f∆) Soma vs Axon", y = "Fraction of Data") +
            scale_y_continuous(position = "left") +
            coord_cartesian(xlim=c(-6,10)) +
            scale_color_manual(" ", values = c("#FDE725", "#440154"), breaks = c("Un-enriched", "Sm-enriched"), labels = c("Un-enriched", "Sm-enriched"))+
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
          table(nijssen$Sm)
          # wilcox test 
            en <- subset(nijssen, Sm == "Sm-enriched")$log2FoldChange
            un <- subset(nijssen, Sm == "Un-enriched")$log2FoldChange
              w1 <- wilcox.test(en, un, exact = FALSE, alternative = 'greater') #true location shift is less than 0:  W = 9042258, p-value = 1.393262e-46




  # upload Todd HMG 2013 COP-RIP data
    # had to convert PDF to excel online, then open in Numbers and replace all '--' with '-', then output as tsv
      todd <- read_delim("Desktop/aCOP-SMN/Todd HMG 2013/dds480supp_table2.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
      # this is currently enriched in aCOP-RIP over NSC-34 input, but in a negative binomial...so - is in fact more enriched with aCOP
        # mutate a column to mark which gene_names are in the aCOP-RIP dataset
          mdeseq_equil <- mdeseq_equil %>% mutate(aCOP = case_when(external_gene_name %in% todd$`Gene Symbol` ~"aCOP-enriched", TRUE ~"Un-enriched"))
          # 1315 gene_names are shared in both datasets.
          
          # now plot cdf for log2f∆ Sm-enrichment delineating by aCOP
          acop_p <- ggplot(mdeseq_equil, aes(log2FoldChange, color = aCOP)) +
            stat_ecdf(geom = "step", linewidth = 1) +
            labs(x = "log2(f∆) anti-Sm-RIP vs polyA-RNA", y = "Fraction of Data") +
            scale_y_continuous(position = "left") +
            coord_cartesian(xlim=c(-4,4)) +
            scale_color_manual(" ", values = c("#FDE725", "#440154"), breaks = c("Un-enriched", "aCOP-enriched"), labels = c("Un-enriched", "aCOP-enriched"))+
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
          table(mdeseq_equil$aCOP)
          # wilcox test 
          en <- subset(mdeseq_equil, aCOP == "aCOP-enriched")$log2FoldChange
          un <- subset(mdeseq_equil, aCOP == "Un-enriched")$log2FoldChange
          w1 <- wilcox.test(en, un, exact = FALSE, alternative = 'greater') 
          
          ggsave(filename = "Sm_aCOP_cdf.png", plot = acop_p, device = "png", dpi = 300, height = 6, width = 5, path = "/filepath/")
          
          # output the aCOP-enriched RNAs for biotype and smsite
          acop_enriched <- mdeseq_equil %>% filter(aCOP == "aCOP-enriched")
          acop_enriched <- acop_enriched %>% filter(log2FoldChange >= 0.6 & padj < 0.05)
          write_tsv(as.data.frame(acop_enriched), file = "/filepath/Sm-enriched_aCOP_master.tsv", quote = "none")
          write_tsv(as.data.frame(table(acop_enriched$gene_biotype)), file = "/filepath/mequil_aCOP_biotype.tsv", quote = "none")
          write_tsv(as.data.frame(table(acop_enriched$smtype)), file = "/filepath/mequil_aCOP_smtype.tsv", quote = "none")
          
          
  # to compare with Lu Genome Biology 2014, HeLa Y12 IP vs IgG isotype control
    Lu_HeLa_Y12 <- read_delim("Desktop/deseq2/tables_figures/supplemental_data/Lu_HeLa_Y12/Lu_HeLa_Y12.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
    # make venn diagram for comparison
      # keep human genes enriched
        hee <- hdeseq_equil %>% filter(log2FoldChange >= 0.6 & padj <= 0.05)
        
        # make a list containing each of the columns above
          lu <- list('Lu 2014 ' = Lu_HeLa_Y12$gene,
                      'Sm-enriched' = hee$external_gene_name)
        
        # plot venn diagram
        lu_venn <- ggvenn(lu,
                          show_percentage = FALSE,
                          fill_color = c("#5dc86375", "#44015475"),
                          stroke_size = 0,
                          set_name_color =c("#5dc863", "#440154"),
                          set_name_size = 6,
                          text_color = "white",
                          text_size = 6)
        ggsave(filename = "lu_venn.png", plot = lu_venn, device = "png", dpi = 300, height = 4, width = 4, path = "/filepath/Lu_HeLa_Y12")
        
        # mutate Lu_HeLa_Y12 to find those enriched with sm-proteins
          lu <- Lu_HeLa_Y12 %>% mutate(Sm_enriched = case_when(gene %in% hee$external_gene_name ~ "yes", TRUE ~ "no"))
          
        # inner_join with human_smsite_master
          names(lu)[1] <- "external_gene_name"
          lu <- lu[,c(1,3:6)]
          lu_master <- inner_join(lu, human_smsite_master, by = "external_gene_name")
          write_tsv(as.data.frame(lu_master), file = "/filepath/Lu_HeLa_Y12/lu_master.tsv", quote = "none")
   
# plotting themes          
  # SmB-CLIP
    design_smb <- "
    1111
    2345
    "
    smb_clip <- guide_area() + smb_p + smb_snRNA_p + smb_lncRNA_p + smb_mRNA_p + plot_layout(design = design_smb, axis_titles = "collect", guides = "collect", heights = c(0.4,4,4,4,4))
    ggsave(filename = "smb_clip.png", plot = smb_clip, device = "png", dpi = 300, height = 6, width = 15, path = "/filepath/")
    
  # nijssen
    design_nij <- "
    12
    "
    nij <- nijssen_p + nijssen_equil_p + plot_layout(design = design_nij, axis_titles = "collect")
    ggsave(filename = "nijssen_soma_vs_axon.png", plot = nij, device = "png", dpi = 300, height = 6, width = 10, path = "/filepath/")
    
  
