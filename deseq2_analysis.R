#!/apps/R/gnu/9.1/3.6.3/bin/Rscript
setwd("/filepath")

library(ggpubr)
library(tidyverse)
library(ggplot2)
library(lemon)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(ggpmisc)
library(ggseqlogo)
library(rstatix)
library(ggdist)
library(readr)
library(DESeq2)
library(stringr)

##### Deseq2 and PCA analysis starts here

# Run these to store global variables for scatter plot theme
BIGPOINT=1
MEDPOINT=0.6
SMALLPOINT=0.2
ALPHAPOINT=0.7
LINETYPE="solid"
LINECOLOR="grey70"
LINEALPHA=1
LINEWIDTH =0.4
FC=1
PVALUE=0.05
RPM = 0

# Run these to store as function for scatter plot theme
my_theme = function() {
    
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.13, "cm"),
          text = element_text(size=10, color = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, vjust = -0.5),
          axis.text = element_text(size = 11, color = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0) ,color = "black"),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), color = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(size = .5, colour = "black"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1,"line"),
          #legend.justification = c(0, 1),
          #legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank())
}

xy_plot_log2 = function(res, x, y, cols, labs, rho, m, M){

    p = ggplot(data = res %>% arrange(order), aes(x = X, y = Y, color = lab)) + 
    geom_point(aes(size = point_size), alpha = 1) + 
    scale_size_identity() + 
    scale_alpha_identity() + 
    my_theme() + 
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2, shape = 15))) + 
    coord_cartesian() + 
    #coord_capped_cart(bottom = "both", left = "both") + 
    scale_y_continuous(
      limits = c(m, M),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)
    ) + 
    scale_x_continuous(
      limits = c(m, M),
      trans = "log2",
      labels = trans_format("log2", math_format(2^.x)),
      breaks = trans_breaks("log2", function(x) 2^x)
    ) +
    scale_color_manual(values = cols, labels = labs) + 
    labs(x = paste0(x), y = paste0(y), color = "") 
    
    if (rho != "None"){
        p = p + labs(subtitle = paste0("Pearson's rho: ", round(rho, 2)))
    } 
    
    return(p)
}

xy_dge = function(res, x, y, axmin, axmax, deseq = FALSE, fold_change = FALSE) {
     point_size = 0.6
#    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    cor = cor.test(res$X, res$Y, method = "pearson")
    rho = cor(res$X, res$Y, method = "pearson")
    
    if (deseq) {
      res = res %>% select(X, Y, log2FoldChange, padj)
      res = res %>% dplyr::rename(pvalue = padj)
    } else {
      res = res %>% select(X, Y, log2FoldChange, pvalue)
    }
    
    if (fold_change == FALSE){
      fc = FC
    } else {
      fc = fold_change
    }

    res = res %>% 
    mutate(lab = ifelse(log2FoldChange >= fc & pvalue < PVALUE & Y >= RPM, "Up",
                        ifelse(log2FoldChange <= -1*fc & pvalue < PVALUE & X >= RPM , "Down", "None"))) %>% 
    mutate(lab = ifelse(is.na(lab), "None", lab)) %>% 
    mutate(order = ifelse(lab == "None", 0, 1)) %>%
    mutate(size = ifelse(lab == "Up" | lab == "Down", .8, 0.4)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 0.7, 0.3))
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    up = res %>% filter(lab == "Up") 
    down = res %>% filter(lab == "Down") 
    unchanged = res %>% filter(lab == "None")
    
    up_percent = round(100*( nrow(up) / nrow(res) ))
    down_percent = round(100*( nrow(down) / nrow(res) ))
    
    unchanged_percent = 100-(up_percent + down_percent)
    
    d = data.frame(x = c("left", "center", "right"), 
                 y = c("middle", "bottom", "bottom"), 
                 label = c(paste0(nrow(up), " (", up_percent,"%)"), 
                           paste0(nrow(down), " (", down_percent,"%)"), 
                           paste0("Pearson's Rho:",round(rho,3))))
    
    cols = c("Up" = "springgreen3", 
           "Down" = "violetred1", 
           "None" = "grey80")
    
    labs = c(
      paste0("log2FC >= ",1*fc," & p < ",PVALUE," (n =",nrow(up),")"),
      paste0("log2FC <= ",-1*fc," & p < ",PVALUE," (n =",nrow(down),")"),
      "Other")
    
    res['point_size'] = point_size
    p = xy_plot_log2(res, x, y, cols, labs, rho, axmin, axmax) 
    p = p +
      geom_abline(slope = 1, intercept = fc, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
      geom_abline(slope = 1, intercept = -1*fc, lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
    
    return(p)
}


# upload the featureCounts master table for comparisons as "counts". 
counts <- read_delim("/filepath to featurecounts output.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

#### Convert Gene Version IDs to Stable IDs

#remove rows pertaining to genes located in pseudoautosomal regions, which are duplicated.
counts<- counts%>%filter(str_detect(gene_id, "PAR") == F)

#split version and stable IDs ('.__'), keeping only the stable ID
counts$Human_Gene_Stable_ID <- str_split(counts$gene_id, "\\.", simplify = T)[,1]

#replace values in "gene_id" with the stable gene ID in "Human_Gene_Stable_ID"
counts$gene_id <- counts$Human_Gene_Stable_ID

#remove now duplicate column "Human_Gene_Stable_ID"
counts <- select(counts, -"Human_Gene_Stable_ID")

#remove ".featureCounts.tsv" from all column header names
colnames(counts) = gsub(".featureCounts.tsv","",colnames(counts)) 

# coldata
# coldata must be inform 

# sample, Condition
# X_rep1, X
# X_rep2, X
# Y_rep1, Y
# Y_rep2, Y

# upload the sample groupings file, make it in excel, save as .csv
coldata <- read_csv("coldata.csv")
coldata = coldata %>% tibble::column_to_rownames(var = "Sample")

# select cols that are non-numeric/numeric
tmp = data.frame(counts, row.names = 1, check.names=FALSE)
names_df = select_if(tmp, Negate(is.numeric))
names_df['gene_id'] = rownames(names_df)
cts = select_if(tmp, is.numeric)
cts = cts[ , sort(colnames(cts))]

# remove rows that have no counts
x = rowSums(cts) >= 1 
cts_filt = cts[x,]
cts_filt[] <- sapply(cts_filt, as.integer)
list(colnames(cts_filt))[[1]]

# set Condition equal to coldata Condition
Condition = coldata$Condition

# make sure columns in counts are in same order as they appear in coldata
rownames(coldata)
cts_filt = cts_filt[,rownames(coldata)]
all( rownames(coldata) == colnames(cts_filt) ) # IF THIS IS FALSE YOU HAVE BIG ISSUE, MUST MAKE THE ORDER OF THE COLUMNS IN TABLE WITH COUNTS SAME ORDER AS ROWS IN COLDATA!!!!!!!!

# DESeq
deobj <- DESeqDataSetFromMatrix(countData = cts_filt, colData = coldata, design = ~Condition)
dds <- DESeq(deobj) #for Wald significance test
tmp = counts(dds, normalized = TRUE)
dds_counts = as.data.frame(tmp) %>% rownames_to_column(var = "gene_id")
colData(dds) #check that group is associating correctly

# PCA plot to assess variance within sample groups and between sample groups 
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(aspect.ratio = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PCA

# results and plot for each comparison
# here you can set your comparisons, look at mine as an example:
comparisons = data.frame(x = c("mr", "mr", "mr", "mx", "mx", "hx", "hx", "hx", "hx+mr", "hx+mr", "mx+hr"), y = c("mx", "hx+mr", "hx+mr+a", "hx+mr", "hx+mr+a", "mr", "hx+mr", "hx+mr+a", "hx+mr+a", "hx+mr+a", "mx+hr+a"))

head(comparisons)

# iterate through comparisons
for(i in 1:nrow(comparisons)){
  
  sample_x = as.character(comparisons[i,1])
  sample_y = as.character(comparisons[i, 2])
  
  curr = paste0(sample_x,"_vs_",sample_y)
  print(curr)
  
  xdf = subset(coldata, Condition == paste0(sample_x))
  xnames = rownames(xdf)
  
  ydf = subset(coldata, Condition == paste0(sample_y))
  ynames = rownames(ydf)
  
  res = results(dds, contrast = c("Condition",paste0(sample_y),paste0(sample_x)))
  
  c = as.data.frame(counts(dds,normalized = TRUE))
  keep = c(paste0(xnames), paste0(ynames))
  
  c_sub = subset(c, select = keep)
  
  resdata = merge(as.data.frame(res), c_sub, by = 'row.names', sort = FALSE)
  names(resdata)[1] <- "gene_id"
  
  resdata[paste0(sample_x)] = rowMeans( resdata[ , xnames] )
  resdata[paste0(sample_y)] = rowMeans( resdata[ , ynames] )

  resdata = merge(resdata, names_df, by = "gene_id", all.x = TRUE)

  deseq_cols = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  cols = c(names_df, deseq_cols)
  
  #ord = c(cols, deseq_cols, list(xnames)[[1]], list(ynames)[[1]], sample_x, sample_y)
  write.table(resdata, paste0(curr,"_deseq_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  p = xy_dge(resdata, paste0(sample_x), paste0(sample_y), 2^-5, 2^20, deseq = TRUE, fold_change = FALSE) # When fold change is set to FALSE it uses 2 fold as a cut off by default
  ggsave(p, filename = paste0(curr,"_deseq_results.pdf"), dpi = 300, height = 12, width = 12)
  
}
