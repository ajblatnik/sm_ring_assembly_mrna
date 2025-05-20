
# -----------------------------
# this script creates a SAF of intronic regions for featureCounts to assign mapped reads to. 
# -----------------------------


library(GenomicFeatures)
library(GenomicRanges)

# Load txdb and get gene/exon ranges
  txdb <- makeTxDbFromGFF("PATH TO GTF/gencode.v40.chr_patch_hapl_scaff.annotation.gtf")
  txdb <- makeTxDbFromGFF("PATH TO GTF/gencode.vM29.chr_patch_hapl_scaff.annotation.gtf")
  genes_gr <- genes(txdb)
  exons_by_gene <- exonsBy(txdb, by = "gene")

# Reduce overlapping exons
  reduced_exons_by_gene <- reduce(exons_by_gene)

# Make sure both are named GRangesList
  gene_ids <- intersect(names(genes_gr), names(reduced_exons_by_gene))
  genes_gr <- genes_gr[gene_ids]
  reduced_exons_by_gene <- reduced_exons_by_gene[gene_ids]

# Convert genes_gr to a GRangesList by splitting
  genes_by_id <- split(genes_gr, genes_gr$gene_id)

# Filter to matching gene IDs
  gene_ids <- intersect(names(genes_by_id), names(reduced_exons_by_gene))
  genes_by_id <- genes_by_id[gene_ids]
  reduced_exons_by_gene <- reduced_exons_by_gene[gene_ids]

# Subtract reduced exons from full gene range to get introns
  introns_by_gene <- Map(setdiff, genes_by_id, reduced_exons_by_gene)

    # sanity check
      all(names(reduced_exons_by_gene) %in% names(genes_by_id))  # Should be TRUE
      all(names(genes_by_id) %in% names(reduced_exons_by_gene))  # Should also be TRUE
      all(names(introns_by_gene) %in% names(genes_by_id))  # Should also be TRUE
      all(names(introns_by_gene) %in% names(reduced_exons_by_gene))  # Should also be TRUE
    
    # code to double-check everything worked
      # Pick a gene with known exons/introns
        gene_id <- "ENSG00000172062.17"  # SMN1, for example
    
        cat("Gene range:\n")
        genes_gr[gene_id]
    
        cat("\nReduced exon ranges:\n")
        reduced_exons_by_gene[[gene_id]]
    
        cat("\nInferred intron ranges:\n")
        intron_ranges <- introns_by_gene[[gene_id]]
        intron_ranges
    
      # Check for any overlaps
        any_overlap <- any(overlapsAny(intron_ranges, reduced_exons_by_gene[[gene_id]]))
        if (any_overlap) {
          message("⚠️ Introns overlap exons! Something went wrong.")
        } else {
          message("✅ Introns are disjoint from exons.")
        }
    
# Flatten the list to a single GRanges object
  introns_gr <- unlist(GRangesList(introns_by_gene))
    
# Assign gene IDs back based on names of the original list
  introns_gr$gene_id <- rep(names(introns_by_gene), lengths(introns_by_gene))
    
# Convert to SAF data frame
  intron_saf <- data.frame(
    GeneID = introns_gr$gene_id,
    Chr    = as.character(seqnames(introns_gr)),
    Start  = start(introns_gr),
    End    = end(introns_gr),
    Strand = as.character(strand(introns_gr)),
    stringsAsFactors = FALSE
    )
    
# Optional: sort for sanity
  intron_saf <- intron_saf[order(intron_saf$Chr, intron_saf$Start), ]
  
  # check the SAF built correctly
    str(intron_saf)
    head(intron_saf, 10)
    
    anyNA(intron_saf)                     # should be FALSE
    any(grepl("\\s", intron_saf$GeneID))  # should be FALSE
    unique(intron_saf$Strand)             # should be "+" "-" "*" — no surprises
    any(intron_saf$Start >= intron_saf$End)  # should be FALSE
    any(intron_saf$Start == intron_saf$End)  # should be FALSE
    
    # create dataframe of only introns where the Start is >= the End
      bad_introns <- intron_saf[intron_saf$Start == intron_saf$End, ]
      bad_introns # to view the dataframe
      nrow(bad_introns) # only 64 for human, 40 for mouse, all are equal start and end, therefore let's remove
      
      intron_saf <- intron_saf[intron_saf$Start < intron_saf$End, ]
    
# write SAF to file for use with featureCounts
  write.table(intron_saf, file="PATH TO SAVING DIRECTORY/gencode.v40.cphsa.by_gene.saf", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(intron_saf, file="PATH TO SAVING DIRECTORY/gencode.vM29.cphsa.by_gene.saf", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
