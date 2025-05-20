#-----------------------------
#-----------------------------
# standard DESeq2 and PCA analyzing new sequencing 
#-----------------------------
#-----------------------------


#-----------------------------
# Load required libraries
#-----------------------------

library(readr)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(gtable)
library(patchwork)
library(ggpointdensity)
library(MASS)

#-----------------------------
# Load featureCounts outputs !! remove header prior to import!!
#-----------------------------

  counts <- read_delim("PATH TO featureCounts OUTPUT.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # remove unneeded columns and simplify sample headers
    counts <- counts[,c(1,7:15)] 
    names(counts) <- c("ensembl_gene_id", "Condition 1 Sample 1", "Condition 1 Sample 2", "Condition 1 Sample 3", "Condition 2 Sample 1", "Condition 2 Sample 2", "Condition 2 Sample 3")
    
#-----------------------------
# Prepare for creation of dds object
#----------------------------- 

# Step 1: remove rows pertaining to genes located in pseudoautosomal regions, which are duplicated.
  counts <- counts %>% filter(str_detect(ensembl_gene_id, "PAR") == F)
  
# split version and stable IDs ('.__'), keeping only the stable ID
  counts$ensembl_gene_id <- str_split(counts$ensembl_gene_id, "\\.", simplify = T)[,1]

# Set rownames to ensembl_gene_id
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$ensembl_gene_id
  counts$ensembl_gene_id <- NULL

# Remove rows that have no counts
  counts <- counts[rowSums(counts) >= 1, ] 
  counts[] <- lapply(counts, as.integer)
  colnames(counts)


# ----------------------------
# Create sample metadata
# ----------------------------
  
# Make sure column names match your BAM file names (from counts)
sample_info <- data.frame(
  sample = colnames(counts),
  condition = c("A", "A", "A", "B", "B", "B")
)
rownames(sample_info) <- sample_info$sample
sample_info$sample <- NULL
# check that row names for sample_info are the same as the column names for counts (and in same order!)
all(rownames(sample_info) == colnames(counts)) # should be TRUE

# ----------------------------
# Create DESeq2 object
# ----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

# Filter out low-count introns
#dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq normalization and analysis
dds <- DESeq(dds)
# pull normalized counts from the dds object
tmp = counts(dds, normalized = TRUE)
# create dataframe with these normalized counts and convert the rownames to the ensembl_gene_id column
dds_counts = as.data.frame(tmp) %>% rownames_to_column(var = "ensembl_gene_id")
colData(dds) #check that group is associating correctly

write_tsv(as.data.frame(dds_counts), file = "~/Desktop/deseq2/tables_figures/reviewer_response/deseq2/human_deseq2_SmB_Y12_counts.tsv", quote = "none")


# ----------------------------
# Perform direct comparisons for future analysis
# ----------------------------
# define direct comparisons to be performed
comparisons = data.frame(
  x = c("A"),
  y = c("B")
)
head(comparisons) # double check data.frame 

# define the gene_id data frame  by select columns that are non-numeric/numeric
gene_ids <- data.frame(ensembl_gene_id = dds_counts$ensembl_gene_id)

# iterate through comparisons
for(i in 1:nrow(comparisons)){
  
  sample_x = as.character(comparisons[i,1])
  sample_y = as.character(comparisons[i,2])
  curr = paste0(sample_x,"_vs_",sample_y)
  print(curr)
  
  xnames <- rownames(subset(sample_info, condition == sample_x))
  ynames <- rownames(subset(sample_info, condition == sample_y))
  
  res = results(dds, contrast = c("condition", sample_y, sample_x))
  
  # Subset normalized counts to samples of interest
  c_sub <- subset(dds_counts, select = c("ensembl_gene_id", xnames, ynames))
  
  # Merge DESeq2 results with normalized counts
  resdata <- merge(as.data.frame(res), c_sub, by.x = "row.names", by.y = "ensembl_gene_id", sort = FALSE)
  names(resdata)[1] <- "ensembl_gene_id"
  
  # Compute average expression across replicates
  resdata[[sample_x]] <- rowMeans(resdata[, xnames])
  resdata[[sample_y]] <- rowMeans(resdata[, ynames])
  
  # Reorder columns if needed
  deseq_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  
  write.table(resdata, file =  paste0("~/PATH TO SAVING DIRECTORY/", curr, "_deseq_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# ----------------------------
# PCA Plot
# ----------------------------
# Perform variance stabilizing transformation
  vsd <- vst(dds, blind=FALSE)

# PCA plot using DESeq2 built-in function
  p <- plotPCA(vsd, intgroup="condition") +
    ggtitle("TITLE") +
    theme_minimal()
  
  ggsave(filename = "NAME.png", plot = p, device = "png", dpi = 300, height = 6, width = 7, path = "PATH TO SAVING DIRECTORY")
 
  # Optional: save PCA data to CSV for external plotting
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
write.csv(pcaData, file = "PATH TO SAVING DIRECTORY/filename.csv", row.names=FALSE)


# ----------------------------
# Scatter Plots for Correlation between SmB and Y12 captured RNAs
# ----------------------------

# import comparison files
  human_Y12_vs_SmB_deseq_results <- read_delim("PATH TO DESEQ2 OUTPUT FILE/human_Y12_vs_SmB_deseq_results.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# inner_join with human_smsite_master to select only protein_coding genes
  protein <- inner_join(human_Y12_vs_SmB_deseq_results, human_smsite_master, by = "ensembl_gene_id") 
  protein <- protein %>% filter(gene_biotype == "protein_coding") 

# scatter plot
  hYB <- ggplot(data = human_Y12_vs_SmB_deseq_results, aes(x = log10(Y12), y = log10(SmB))) + 
  #hYB_prot <- ggplot(data = protein, aes(x = log10(Y12), y = log10(SmB))) + 
    geom_pointdensity(alpha = 0.5) +
    scale_color_viridis_c(name = "Transcripts") +
    geom_smooth(method = "rlm",     # robust linear model
                method.args = list(method = "M"), # Huber M-estimator
                formula = y ~ x,
                se     = TRUE,     # shaded 95% CI
                color  = "black",  # line color
                linetype = "dashed",
                linewidth   = 0.8) +
    theme_minimal() +
    labs(y = "SmB/B'/N", x = "Y12") +
    #stat_cor(method = "pearson", label.x = 0, label.y = 6, size = 5, label.sep="\n")+
    theme(panel.background = element_rect(fill = NULL, color = "black", linewidth = rel(3)), 
          panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
          axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
          #legend.position = "none"
          legend.position = c(0.8, 0.25), 
          legend.box.background = element_rect(fill = "white"), 
          legend.key = element_blank(), 
          legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
          legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black")
          )
  
  ggsave(filename = "h_SmB_Y12_pearson_corr.png", plot = hYB, device = "png", dpi = 300, height = 5, width = 5.625, path = "Desktop/deseq2/tables_figures/reviewer_response/SmB_Y12/human")
  ggsave(filename = "h_SmB_Y12_prot_pearson_corr.png", plot = hYB_prot, device = "png", dpi = 300, height = 5, width = 5.625, path = "Desktop/deseq2/tables_figures/reviewer_response/SmB_Y12/human")
  
