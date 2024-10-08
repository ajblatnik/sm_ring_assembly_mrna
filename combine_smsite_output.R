

## this script is used to combine Sm-site output data into a single file

# load libraries
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)

 # import Sm-site outputs from find_smsites.py
    u1 <- read_delim("filepath to u1.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    u2 <- read_delim("filepath to u2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    u5 <- read_delim("filepath to u5.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    u7 <- read_delim("filepath to u7.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    noncanonical <- read_delim("filepath to noncanonical.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # make a column name providing the smtype
      u1_u <- mutate(u1, smtype = "U1") 
      u2_u <- mutate(u2, smtype = "U2") 
      u5_u <- mutate(u5, smtype = "U5") 
      u7_u <- mutate(u7, smtype = "U7") 

      # use bind_rows() to combine u1, u2, u5, and u7 
        combined <- bind_rows(u1, u2, u5, u7) 
        
        # use setdiff() to only keep rows unique to noncanonical
          nc <- setdiff(noncanonical, combined) 
          nc <- mutate(nc, smtype = "Noncanonical") 
      
      # use bind_rows() to combine u1, u2, u5, u7, and nc
        smsite_output <- bind_rows(u1_u, u2_u, u5_u, u7_u, nc) 
        
      # write output to .tsv
        write_tsv(as.data.frame(smsite_output), file = "/filepath/smsite_output.tsv", quote = "none")
     
        
        
