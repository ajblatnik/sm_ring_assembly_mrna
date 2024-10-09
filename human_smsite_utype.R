

## this script takes the Sm-site identification pipeline output and assembles a new document based off transcript isoforms for human
## takes find_smsites.py output files (total of 10 files (Refseq and Gencode, U1, U2, U5, U7, Noncanonical) and crossreferences identified 
## Sm-sites between both Refseq and Gencode annotations, yielding a single file for each Sm-site type (U1, U2, U5, U7). Then these are 
## processed to generate a single master file yielding the frequency of Sm-sites identified in the MANE/Canonical transcript of a gene,
## yielding a master table giving Sm-site breakdowns for a single transcript per gene and used as a metagene analysis for the study.

# import libraries
  library(readr)
  library(stringr)
  library(tidyr)
  library(dplyr)
  library(tidyverse)
  library(biomaRt)

# in each step below, use variations of the following script to check the number of gene_ids, transcript_ids, etc within each new iteration of the tables
  table <- as.data.frame(table(TABLE$COLUMN))

# create a table to cross-reference refseq and ensembl ids
  # access BioMart and find attributes of the database
    ensembl <- useMart("ensembl")
    ensembl_h <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
    attributes <- listAttributes(ensembl_h)
  # create table with ensembl_gene_id, ensembl_transcript_id, transcript_is_canonical, refseq_mrna, refseq_ncrna, transcript_mane_select, gene_biotype
    refseq <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical", "refseq_mrna", "refseq_ncrna", "transcript_mane_select", "gene_biotype"),
                    mart = ensembl_h,
                    useCache = F) # 318178 obs, 70611 gene_ids, 278220 transcript_ids, however multiple refseq_mrna and reseq_ncrna ids
    write_tsv(as.data.frame(refseq), file = "/Users/dapperwhiskyman/Desktop/deseq2/from_BioMart/h_refseq.tsv", quote = "none")
  # create table with ensembl_transcript_id, transcript_length, 5_utr_start, 5_utr_end, CDS_length, 3_utr_start, 3_utr_end)
    lengths <- getBM(attributes = c("ensembl_transcript_id", "transcript_length", "5_utr_start", "5_utr_end", "cds_length", "3_utr_start", "3_utr_end"),
                    mart = ensembl_h,
                    useCache = F) # 673496 obs, 278220 transcript_ids
    # subtract 5'utr and 3'utr start and ends to give distances
    # rename header to allow proper column calling
      names(lengths) <- c("ensembl_transcript_id", "transcript_length", "utr5_start", "utr5_end", "cds_length", "utr3_start", "utr3_end")
    # make 'NA' = '0'
      missing_values <- is.na(lengths$cds_length)
      lengths$cds_length[missing_values] <- 0
      # make two new columns equal to utr5_end - utr5_start and utr3_end - utr3_start
        lengths <- lengths %>% mutate(utr5_length = utr5_end - utr5_start)
        lengths <- lengths %>% mutate(utr3_length = utr3_end - utr3_start)
        lengths <- lengths[,c(1:2,5,8:9)]
        # split apart to select the longest lengths corresponding to each transcript
          lengths_tl <- lengths[,c(1:2)]
          lengths_cl <- lengths[,c(1,3)]
          lengths_5l <- lengths[,c(1,4)]
          lengths_3l <- lengths[,c(1,5)]
        # sort by descending length and slice to top (ie longest length)
          lengths_tl <- lengths_tl %>% group_by(ensembl_transcript_id) %>% arrange(desc(transcript_length)) %>% slice(1) %>% ungroup() 
          lengths_cl <- lengths_cl %>% group_by(ensembl_transcript_id) %>% arrange(desc(cds_length)) %>% slice(1) %>% ungroup() 
          lengths_5l <- lengths_5l %>% group_by(ensembl_transcript_id) %>% arrange(desc(utr5_length)) %>% slice(1) %>% ungroup() 
          lengths_3l <- lengths_3l %>% group_by(ensembl_transcript_id) %>% arrange(desc(utr3_length)) %>% slice(1) %>% ungroup() 
      # inner_join by ensembl_transcript_id to provide a table with the longest lengths for each measure for each transcript
        lengths_final <- inner_join(lengths_tl, lengths_5l, by = "ensembl_transcript_id")
        lengths_final <- inner_join(lengths_final, lengths_cl, by = "ensembl_transcript_id")
        lengths_final <- inner_join(lengths_final, lengths_3l, by = "ensembl_transcript_id") , 
      write_tsv(as.data.frame(lengths_final), file = "/filepath/h_lengths.tsv", quote = "none")
      
  # inner_join refseq with lengths_final by ensembl_transcript_id since the same number of transcript_ids
    refleng <- inner_join(refseq, lengths_final, by = "ensembl_transcript_id") 
    # make a table containing only MANE transcripts
      mane <- refleng %>% filter(!is.na(transcript_mane_select) & transcript_mane_select != "") 
      # need to filter this down to a single gene_id to single transcript_id to single refseq_id
      # write a function to remove the suffix on the transcript_mane_select refseq_id and then run it on mane
        remove_suffix <- function(x) {sub("\\.[^.]*$", "", x)}
        mane$transcript_mane_select <- sapply(mane$transcript_mane_select, remove_suffix)
        # now keep rows in which refseq_mrna is equal to transcript_mane_select
          mane <- mane %>% filter(refseq_mrna == transcript_mane_select) 
          # the remaining duplicates have multiple refseq_ncrna, which can be removed as the refseq_mrna will take precedent
            mane <- mane %>% group_by(ensembl_gene_id) %>% slice(1) %>% ungroup()
            # remove uneccessary columns and rename refseq_mrna to refseq_id
              mane <- mane[c(1:2,4,7:11)]
              names(mane)[3] <- "refseq_id"
    # make a table of gene_ids that are not MANE
      tisc <- refleng %>% anti_join(mane, by = "ensembl_gene_id") 
        # keep only those with transcript_is_canonical == 1
          tisc <- tisc %>% filter(transcript_is_canonical == "1") 
          # sort by descending transcript_length and slice to top (ie longest length)
            tisc <- tisc %>% group_by(ensembl_gene_id) %>% arrange(desc(transcript_length)) %>% slice(1) %>% ungroup() 
            # convert empty strings to NA
              tisc <- tisc %>% mutate(refseq_mrna = na_if(refseq_mrna, ""))
              # replace NA values in refseq_mrna with corresponding values from refseq_ncrna, keeping refseq_mrna values if they exist
                tisc <- tisc %>% mutate(refseq_mrna = if_else(is.na(refseq_mrna), refseq_ncrna, refseq_mrna)) 
                # remove uneccessary columns and rename refseq_mrna to refseq_id
                  tisc <- tisc[c(1:2,4,7:11)]
                  names(tisc)[3] <- "refseq_id"
      
    # bind_rows from mane and tisc to build the complete refseq table yielding a single gene_id, transcript_id, and refseq_id, etc.
      refseq_final <- bind_rows(mane, tisc)  
      write_tsv(as.data.frame(refseq_final), file = "/filepath/h_refseq_master.tsv", quote = "none")
 
      
      
      
      
      
# import human refseq and gencode sm-site output files - called u1g for U1 Sm-site from gencode

  input_dir_r <- "/filepath to output files from find_smsites.py"
  input_dir_g <- "/filepath to output files from find_smsites.py"
  output_dir_r <- "/filepath to desired output directory"
  output_dir_g <- "/filepath to desired output directory"

  # get all .txt files in the input directory
  txt_files_r <- list.files(input_dir_r, pattern = "*_stem_loops_nih_human.txt", full.names = TRUE)
  txt_files_g <- list.files(input_dir_g, pattern = "*_stem_loops_gencode_human.txt", full.names = TRUE)
  
  remove_suffix <- function(x) {sub("\\.[^.]*$", "", x)}
  
  # Loop through each refseq output file
  for (file_r in txt_files_r) {
        
    ur <- read_delim(file_r, delim = "|", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
    ur <- ur[,c(1:3,5,7)]
    names(ur) <- c("refseq_id", "position", "region", "gene_name", "gene_biotype")
    # the sequence is not added with '|' but on a separate space, so need to remove and paste into the last column
      for(i in 1:(nrow(ur)-1)){ur$Sequence[i] <- as.character(ur$refseq_id[i+1])}
      # remove second row as these are blank or duplicated
        ur <- ur %>% filter(!is.na(region))
        # remove '>' in front of refseq_id
          ur$refseq_id <- gsub(">", "", ur$refseq_id)
          # remove '.x' from refseq_id (this is the accession information) using remove_suffix function written in line 61
            ur$refseq_id <- sapply(ur$refseq_id, remove_suffix)
    
          # Generate output filename using the first two characters of the basename
            base_name <- tools::file_path_sans_ext(basename(file_r))
            output_filename_r <- paste0(output_dir_r, "/", substr(base_name, 1, 2), "_nih.tsv")  
          
          # write the processed data to a new file
            write_tsv(as.data.frame(ur), file = output_filename_r, quote = "none")
  }
     
  # Loop through each gencode output file
   for (file_g in txt_files_g) {     
    ug <- read_delim(file_g, delim = "|", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    ug <- ug[,c(1:5)]
    names(ug) <- c("ensembl_transcript_id", "position", "region", "gene_biotype", "gene_name")
    # the sequence is not added with '|' but on a separate space, so need to remove and paste into the last column
      for(i in 1:(nrow(ug)-1)){ug$Sequence[i] <- as.character(ug$ensembl_transcript_id[i+1])}
      # remove second row as these are blank or duplicated
        ug <- ug %>% filter(!is.na(region))
        # remove '>' in front of ensembl_transcript_id
          ug$ensembl_transcript_id <- gsub(">", "", ug$ensembl_transcript_id)
          # remove '.x' from refseq_id (this is the accession information) using remove_suffix function written in line 61
            ug$ensembl_transcript_id <- sapply(ug$ensembl_transcript_id, remove_suffix)
            
            # Generate output filename using the first two characters of the basename
            base_name <- tools::file_path_sans_ext(basename(file_g))
            output_filename_g <- paste0(output_dir_g, "/", substr(base_name, 1, 2), "_gencode.tsv")  
            
            # write the processed data to a new file
            write_tsv(as.data.frame(ug), file = output_filename_g, quote = "none")
  }
 
  
# now to combine tables and make one master smsite list that crossrefrences both annotations
  input_dir_r2 <- "/filepath to output from above code processing refseq data"
  input_dir_g2 <- "/filepath to output from above code processing gencode data"
  output_dir <- "/filepath to desired output directory"
  
  # need to remember to import refseq_final (h_refseq_master or m_refseq_master)
  
  # get all .txt files in the input directory
    txt_files_r2 <- list.files(input_dir_r2, pattern = "*_nih.tsv", full.names = TRUE)
    txt_files_g2 <- list.files(input_dir_g2, pattern = "*_gencode.tsv", full.names = TRUE)
  
  # Ensure both directories have the same number of files
    if (length(txt_files_r2) != length(txt_files_g2)) {
      stop("The number of files in the two directories do not match.")
    }
  
  # Loop through each file
    for (i in seq_along(txt_files_r2)) { 
   
    file_r2 <- txt_files_r2[i]
    file_g2 <- txt_files_g2[i]
    
    ur <- read_delim(file_r2, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
    ug <- read_delim(file_g2, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # for each, keep only those refseq_ids or transcript_ids present in the refseq_master table (do XR annotated genes have a corresponding gencode annotation) then full_join by sequence
    ur <- ur %>% filter(refseq_id %in% refseq_final$refseq_id) 
    ug <- ug %>% filter(ensembl_transcript_id %in% refseq_final$ensembl_transcript_id) 
  # full_join by sequence column
    uf <- full_join(ur, ug, by = "Sequence")
    # some rows are duplicate sequences but the annotated positions between refseq and gencode are different, and gene_name.x is different from gene_name.y
      # make a table of only unique sequences
        uf_unique <- uf %>% group_by(Sequence) %>% filter(n() == 1) %>% ungroup()
        # make table including only those rows that either do or do not contain an ensembl_transcript_id
        uf_unique_do <- uf_unique %>% filter(!is.na(ensembl_transcript_id))
        uf_unique_donot <- uf_unique %>% filter(is.na(ensembl_transcript_id)) 
        # remove empty columns and then merge with refseq_final
          uf_unique_donot <- uf_unique_donot[,c(1:6)]  
          uf_unique_donot <- inner_join(uf_unique_donot, refseq_final, by = "refseq_id") 
          # for each unique sequence, slice to one
            uf_unique_donot <- uf_unique_donot %>% group_by(Sequence) %>% slice(1) %>% ungroup() 
            uf_unique_donot <- uf_unique_donot[,c(1,3,4,6,8:9)]
            names(uf_unique_donot) <- c("refseq_id", "region", "gene_name", "Sequence", "ensembl_transcript_id", "gene_biotype")
      # make a table of sequences that are duplicated
        uf_duplicate <- uf %>% filter(!(Sequence %in% uf_unique$Sequence)) 
        # filter for those with duplicate ensembl_transcript_ids
          ensembl_transcript_id <- "ensembl_transcript_id"
          uf_duplicate1 <- uf_duplicate %>% group_by(.data[[ensembl_transcript_id]]) %>% filter(n() > 1) %>% ungroup() 
          # for duplicate sequences, filter out those in which the gene names in x and y aren't the same
            uf_duplicate2 <- uf_duplicate1 %>% group_by(Sequence) %>% filter(n() > 1 & gene_name.x == gene_name.y) %>% ungroup() 
          # table of remaining sequences
            uf_duplicate3 <- uf_duplicate1 %>% filter(!(Sequence %in% uf_duplicate2$Sequence)) 
            # group and slice to 1
              uf_duplicate3 <- uf_duplicate3 %>% group_by(Sequence) %>% slice(1) %>% ungroup() 
        # table of remaining sequences from uf_duplicate
          uf_duplicate4 <- uf_duplicate %>% filter(!(Sequence %in% uf_duplicate1$Sequence)) 
          # for duplicate sequences, filter out those in which the gene names in x and y aren't the same
            uf_duplicate5 <- uf_duplicate4 %>% group_by(Sequence) %>% filter(n() > 1 & gene_name.x == gene_name.y) %>% ungroup() 
            uf_duplicate5 <- uf_duplicate5 %>% group_by(Sequence) %>% slice(1) %>% ungroup() 
          # for remaining sequences from uf_duplicate4
            uf_duplicate6 <- uf_duplicate4 %>% filter(!(Sequence %in% uf_duplicate5$Sequence)) 
            # remove gene_names starting with 'ENSG'
              sequence <- "Sequence"
              uf_duplicate7 <- uf_duplicate6 %>% group_by(.data[[sequence]]) %>% filter(!(duplicated(.data[[sequence]]) | duplicated(.data[[sequence]], fromLast = TRUE)) | !grepl("^ENSG", gene_name.y)) %>% ungroup() # 67 obs, 41 unique sequences
            # for remaining sequences in uf_duplicate6
              uf_duplicate8 <- uf_duplicate6 %>% filter(!(Sequence %in% uf_duplicate7$Sequence)) 
        
      # make a table of duplicate sequences, filtering for same gene_name !! does not pick-up anything without an ensembl_transcript_id or refseq_id
        # arrange by sequence, slice to row containing mane transcript
          uf_duplicate_mane <- uf_duplicate %>% filter(refseq_id %in% mane$refseq_id) 
        # make a table of those sequences that are not mane transcripts 
          uf_duplicate_notmane <- uf_duplicate %>% filter(!(Sequence %in% uf_duplicate_mane$Sequence)) 
      # make table of sequences that do not have a refseq_id
        uf_no_rsi <- uf %>% filter(is.na(refseq_id) | refseq_id == "") 
        # remove those that are in uf_unique or uf_duplicate
          uf_no_rsi <- uf_no_rsi %>% filter(Sequence %in% uf_unique$Sequence | Sequence %in% uf_duplicate) 
      # bind_rows
        uf_final <- bind_rows(uf_unique_do, uf_duplicate2)
        uf_final <- bind_rows(uf_final, uf_duplicate3)
        uf_final <- bind_rows(uf_final, uf_duplicate5) 
        uf_final <- bind_rows(uf_final, uf_duplicate7)
        uf_final <- bind_rows(uf_final, uf_duplicate8) 
        uf_final <- uf_final[,c(1,6:7,9:11)]
        names(uf_final) <- c("refseq_id", "Sequence", "ensembl_transcript_id", "region", "gene_biotype", "gene_name")
        uf_final <- bind_rows(uf_final, uf_unique_donot)
        uf_final <- uf_final[,c("ensembl_transcript_id", "region", "gene_biotype", "gene_name", "refseq_id", "Sequence")]
      # there are some sequences identified with a refseq_id that are not identified in gencode, want to pull the gencode information
        
        # Generate output filename using the first two characters of the basename
          base_name <- tools::file_path_sans_ext(basename(file_g2))
          output_filename <- paste0(output_dir, "/", substr(base_name, 1, 2), "_human.tsv")  
        
        # write the processed data to a new file
        write_tsv(as.data.frame(uf_final), file = output_filename, quote = "none")
  }    
        
      
    

# compile a single master table from refseq_final (h_refseq_master or m_refseq_master) that provides smsite-frequency data
  
  u1 <- read_delim("/filepath to u1_human.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  u2 <- read_delim("/filepath to u2_human.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  u5 <- read_delim("/filepath to u5_human.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  u7 <- read_delim("/filepath to u7_human.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  pe <- read_delim("filepath to noncanonical_human.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # get smsite frequency per ensembl_transcript_id
    u1_freq <- as.data.frame(table(u1$ensembl_transcript_id)) 
    u2_freq <- as.data.frame(table(u2$ensembl_transcript_id)) 
    u5_freq <- as.data.frame(table(u5$ensembl_transcript_id)) 
    u7_freq <- as.data.frame(table(u7$ensembl_transcript_id)) 
    pe_freq <- as.data.frame(table(pe$ensembl_transcript_id)) 
    # rename columns for each
      names(u1_freq) <- c("ensembl_transcript_id", "u1_freq")
      names(u2_freq) <- c("ensembl_transcript_id", "u2_freq")
      names(u5_freq) <- c("ensembl_transcript_id", "u5_freq")
      names(u7_freq) <- c("ensembl_transcript_id", "u7_freq")
      names(pe_freq) <- c("ensembl_transcript_id", "pe_freq")
      # full_join to make a comprehensive table
        smsite <- full_join(u1_freq, u2_freq, by = "ensembl_transcript_id")
        smsite <- full_join(smsite, u5_freq, by = "ensembl_transcript_id")
        smsite <- full_join(smsite, u7_freq, by = "ensembl_transcript_id")
        smsite <- full_join (smsite, pe_freq, by = "ensembl_transcript_id")
        # replace all NA with 0
          replace_na_u1 <- is.na(smsite$u1_freq)
            smsite$u1_freq[replace_na_u1] <- 0
          replace_na_u2 <- is.na(smsite$u2_freq)
            smsite$u2_freq[replace_na_u2] <- 0
          replace_na_u5 <- is.na(smsite$u5_freq)
            smsite$u5_freq[replace_na_u5] <- 0
          replace_na_u7 <- is.na(smsite$u7_freq)
            smsite$u7_freq[replace_na_u7] <- 0
          replace_na_pe <- is.na(smsite$perm_freq)
            smsite$perm_freq[replace_na_pe] <- 0
          # make a new column with the difference between perm_freq and the rest, to make a column of truly permissive smsites
            smsite <- smsite %>% mutate(nc_freq = pe_freq - u1_freq - u2_freq - u5_freq - u7_freq) 
            # remove perm_freq column
              smsite <- smsite[,c(1:5,7)] # there is a single transcript_id represented
              
  # append smsite frequency for each region of rna
      list <- list(u1 = u1, u2 = u2, u5 = u5, u7 = u7, pe = pe)
      
    # for 5'UTR
      process_5utr <- function(data, source) {
        data_5utr <- data %>% filter(region == "5'UTR")
        data_5utr_freq <- as.data.frame(table(data_5utr$ensembl_transcript_id))
        names(data_5utr_freq) <- c("ensembl_transcript_id", paste0(source, "_5utr_freq"))
        data_5utr_freq[, paste0(source, "_5utr_freq")][is.na(data_5utr_freq[, paste0(source, "_5utr_freq")])] <- 0
        return(data_5utr_freq)
      }
      # Apply the function to each dataset in the list
        results <- lapply(names(list), function(name) process_5utr(list[[name]], name))
        utr5 <- Reduce(function(x, y) full_join(x, y, by = "ensembl_transcript_id"), results)
        utr5 <- utr5 %>% mutate(across(ends_with("_5utr_freq"), ~ replace_na(., 0)))
        utr5 <- utr5 %>% mutate(nc_5utr_freq = pe_5utr_freq - u1_5utr_freq - u2_5utr_freq - u5_5utr_freq - u7_5utr_freq)
        utr5 <- utr5[,c(1:5,7)]
        
    # for CDS
      process_cds <- function(data, source) {
        data_cds <- data %>% filter(region == "CDS")
        data_cds_freq <- as.data.frame(table(data_cds$ensembl_transcript_id))
        names(data_cds_freq) <- c("ensembl_transcript_id", paste0(source, "_cds_freq"))
        data_cds_freq[, paste0(source, "_cds_freq")][is.na(data_cds_freq[, paste0(source, "_cds_freq")])] <- 0
        return(data_cds_freq)
      }
        # Apply the function to each dataset in the list
          results <- lapply(names(list), function(name) process_cds(list[[name]], name))
          cds <- Reduce(function(x, y) full_join(x, y, by = "ensembl_transcript_id"), results)
          cds <- cds %>% mutate(across(ends_with("_cds_freq"), ~ replace_na(., 0)))
          cds <- cds %>% mutate(nc_cds_freq = pe_cds_freq - u1_cds_freq - u2_cds_freq - u5_cds_freq - u7_cds_freq)
          cds <- cds[,c(1:5,7)]
              
    # for 3'UTR
      process_3utr <- function(data, source) {
        data_3utr <- data %>% filter(region == "3'UTR")
        data_3utr_freq <- as.data.frame(table(data_3utr$ensembl_transcript_id))
        names(data_3utr_freq) <- c("ensembl_transcript_id", paste0(source, "_3utr_freq"))
        data_3utr_freq[, paste0(source, "_3utr_freq")][is.na(data_3utr_freq[, paste0(source, "_3utr_freq")])] <- 0
        return(data_3utr_freq)
      }
        # Apply the function to each dataset in the list
          results <- lapply(names(list), function(name) process_3utr(list[[name]], name))
          utr3 <- Reduce(function(x, y) full_join(x, y, by = "ensembl_transcript_id"), results)
          utr3 <- utr3 %>% mutate(across(ends_with("_3utr_freq"), ~ replace_na(., 0)))
          utr3 <- utr3 %>% mutate(nc_3utr_freq = pe_3utr_freq - u1_3utr_freq - u2_3utr_freq - u5_3utr_freq - u7_3utr_freq)
          utr3 <- utr3[,c(1:5,7)]
        
    # for Exon_Mismatch (u7 does not have any so need to modify code, otherwise the function errors)
      list1 <- list(u1 = u1, u2 = u2, u5 = u5, pe = pe)
      process_exmis <- function(data, source) {
        data_exmis <- data %>% filter(region == "Exon_Mismatch")
        data_exmis_freq <- as.data.frame(table(data_exmis$ensembl_transcript_id))
        names(data_exmis_freq) <- c("ensembl_transcript_id", paste0(source, "_exon_mismatch_freq"))
        data_exmis_freq[, paste0(source, "_exon_mismatch_freq")][is.na(data_exmis_freq[, paste0(source, "_exon_mismatch_freq")])] <- 0
        return(data_exmis_freq)
        }
        # Apply the function to each dataset in the list
          results <- lapply(names(list1), function(name) process_exmis(list1[[name]], name))
          exmis <- Reduce(function(x, y) full_join(x, y, by = "ensembl_transcript_id"), results)
          exmis <- exmis %>% mutate(across(ends_with("_exon_mismatch_freq"), ~ replace_na(., 0)))
          exmis <- exmis %>% mutate(nc_exon_mismatch_freq = pe_exon_mismatch_freq - u1_exon_mismatch_freq - u2_exon_mismatch_freq - u5_exon_mismatch_freq)
          exmis <- exmis[,c(1:4,6)]
     
    # for N/A
      process_na <- function(data, source) {
        data_na <- data %>% filter(region == "N/A")
        data_na_freq <- as.data.frame(table(data_na$ensembl_transcript_id))
        names(data_na_freq) <- c("ensembl_transcript_id", paste0(source, "_na_freq"))
        data_na_freq[, paste0(source, "_na_freq")][is.na(data_na_freq[, paste0(source, "_na_freq")])] <- 0
        return(data_na_freq)
      }
        # Apply the function to each dataset in the list
          results <- lapply(names(list), function(name) process_na(list[[name]], name))
          na <- Reduce(function(x, y) full_join(x, y, by = "ensembl_transcript_id"), results)
          na <- na %>% mutate(across(ends_with("_na_freq"), ~ replace_na(., 0)))
          na <- na %>% mutate(nc_na_freq = pe_na_freq - u1_na_freq - u2_na_freq - u5_na_freq - u7_na_freq)
          na <- na[,c(1:5,7)]
      
    # full_join into one table 
      region <- full_join(utr5,cds, by = "ensembl_transcript_id")
      region <- full_join(region, utr3, by = "ensembl_transcript_id")
      region <- full_join(region, exmis, by = "ensembl_transcript_id")
      region <- full_join(region, na, by = "ensembl_transcript_id")
      names(region) <-c("ensembl_transcript_id", "u1_5utr", "u2_5utr", "u5_5utr", "u7_5utr", "nc_5utr", "u1_cds", "u2_cds", "u5_cds", "u7_cds", "nc_cds", "u1_3utr", "u2_3utr", "u5_3utr", "u7_3utr", "nc_3utr", "u1_exon_mismatch", "u2_exon_mismatch", "u5_exon_mismatch", "nc_exon_mismatch", "u1_na", "u2_na", "u5_na", "u7_na", "nc_na")
          
    
    # combine with smsite 
      smsite_region <- inner_join(smsite, region, by = "ensembl_transcript_id") 
      #combine with refseq_final
        master <- full_join(refseq_final, smsite_region, by = "ensembl_transcript_id") 
        # replace all NA with 0
          master <- master %>% mutate(across(everything(), ~replace_na(., 0)))
        # mutate to add an smsite column: canon (u1, u2, u5, u7), noncanon (nc), absent, case_when is ordered first to last, therefore if found in u1, will label u1 and remove from the query
          # *keep U7 separate as it recieves a different sm-ring
          master <- master %>% mutate(smsite = case_when(ensembl_transcript_id %in% u1$ensembl_transcript_id ~ "canon",
                                                         ensembl_transcript_id %in% u2$ensembl_transcript_id ~ "canon",
                                                         ensembl_transcript_id %in% u5$ensembl_transcript_id ~ "canon",
                                                         ensembl_transcript_id %in% u7$ensembl_transcript_id ~ "U7",
                                                         ensembl_transcript_id %in% pe$ensembl_transcript_id ~ "noncanon",
                                                         TRUE ~ "absent"))
          master <- master %>% mutate(smtype = case_when(ensembl_transcript_id %in% u2$ensembl_transcript_id ~ "U2",
                                                         ensembl_transcript_id %in% u1$ensembl_transcript_id ~ "U1",
                                                         ensembl_transcript_id %in% u5$ensembl_transcript_id ~ "U5",
                                                         ensembl_transcript_id %in% u7$ensembl_transcript_id ~ "U7",
                                                         ensembl_transcript_id %in% pe$ensembl_transcript_id ~ "NC",
                                                         TRUE ~ "NA"))
          
        # now add gene_name to make the true master table 
          gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = ensembl_h,
                           useCache = F)
          master <- inner_join(master, gene_name, by = "ensembl_gene_id")
          # reorder to double check frequencies
            column_order <- c("ensembl_gene_id", "ensembl_transcript_id", "refseq_id", "external_gene_name", "gene_biotype", "transcript_length", "utr5_length", "cds_length", "utr3_length", "smsite", "smtype", "u1_freq", "u1_5utr", "u1_cds", "u1_3utr", "u1_exon_mismatch", "u1_na", "u2_freq", "u2_5utr", "u2_cds", "u2_3utr", "u2_exon_mismatch", "u2_na", "u5_freq", "u5_5utr", "u5_cds", "u5_3utr", "u5_exon_mismatch", "u5_na", "u7_freq", "u7_5utr", "u7_cds", "u7_3utr", "u7_na", "nc_freq", "nc_5utr", "nc_cds", "nc_3utr", "nc_exon_mismatch", "nc_na")
            master <- master[, column_order]
            # make a column combining U1, U2, and U5 freq to be canon_freq, and also canon_3utr_freq
              master <- master %>% mutate(canon_freq = u1_freq + u2_freq + u5_freq)
              master <- master %>% mutate(canon_3utr = u1_3utr + u2_3utr + u5_3utr)
                write_tsv(as.data.frame(master), file = "/filepath/human_smsite_master.tsv", quote = "none")
            
    

  

 
  
  
  
  
  
