# Identify transcripts that contain Sm-sites from both RefSeq and Gencode genome and transcriptome references
# Script order: create_database.py -> sm_sequence.py -> rnafold -> stem_loops.py 
# copy and paste each of the following python sripts into their own file to run
# scripts for human_refseq, human_gencode, mouse_refseq, and mouse_gencode are separated by ####
# Organization of the full document:
    #####################

    # Bash setup of python environment

    #####################
    #####################

    # create_database.py

        # human_refseq

    #####################

        # mouse_refseq

    #####################

        # human_gencode

    #####################

        # mouse_gencode

    #####################
    #####################
    #####################

    # sm_sequence.py

        # Refseq

    #####################

        # Gencode

    #####################
    #####################
    #####################

    # rnafold

    #####################
    #####################
    #####################

    # stem_loops.py

####################################################################################################################################################################################### 

# Setup Python environment using terminal. These are commented out as they should be supplied in command line and not run in the .py 
# Install Miniconda3 (if not already installed)
# Download and install from: https://docs.conda.io/en/latest/miniconda.html
    # module load miniconda 3
# Create a new conda environment named 'smsite_env'
    # conda create -n smsite_env python=3.9

# Activate the environment
    # conda activate smsite_env

# Install the required packages
    # conda install -c conda-forge biopython pandas numpy sqlite 

# Verify the installations
    # python -c "import Bio; import pandas; import numpy; import sqlite3; print('All packages are installed correctly')"

# Run your script using the code below
    # python your_script.py

####################################################################################################################################################################################### 
####################################################################################################################################################################################### 

# create_database.py 

# following is the python script used to generate the database to easily retrieve gene annotation information.

# Import packages 
# sqlite3: Provides a lightweight disk-based database, which does not require a separate server process.
# Bio.SeqIO, Bio.Seq, Bio.SeqUtils: From Biopython, used for handling sequence data.
# pandas: Data analysis library, though it's not used in the script.
# os: Used to change the current working directory.
import sqlite3
from sqlite3 import Error
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqUtils
import os
import re

# for human NCBI RefSeq
HUMAN_REFSEQ_DATADIRECTORY = '/filepath'
HUMAN_REFSEQ_GFFFILENAME = '/filepath/GCF_000001405.40_GRCh38_genomic.gff'
human_refseq_db = '/filepath/GRCh38_refseq.db'

os.chdir(HUMAN_REFSEQ_DATADIRECTORY)

conn = sqlite3.connect(human_refseq_db)
cur = conn.cursor()


#create database tables
cur.execute("""CREATE TABLE TRANSCRIPTS (
            TRANSCRIPT_ID text not null,
            TRANSCRIPT_TYPE text not null,
            GENE_NAME text not null,
            PARENT_GENE_NAME text not null,
            CDS_START_INDEX integer,
            CDS_LEN integer,
            EXON_LENS integer,
            primary key(TRANSCRIPT_ID)
)""")

cur.execute("""CREATE TABLE GENES (
            GENE_ID text not null,
            GENE_NAME text not null,
            GENE_TYPE text not null,
            primary key(GENE_ID)
)""")

#set collum numbers for gff file
typeidx = 2
infoidx = 8
start_idx = 3
end_idx = 4
strand_idx = 6

id_info_idx = 0
gene_name_info_idx = 2
parent_gene_info_idx = 1
parent_rna_idx = 1

with open(HUMAN_REFSEQ_GFFFILENAME, 'r') as f:
    for whole_line in f:
        if not whole_line.startswith('#'):
            line = whole_line.strip().split('\t')
            info = line[infoidx].split(';')

            # if entry is of type transcript add it to table
            if "RNA" in line[typeidx].upper() or 'transcript' in line[typeidx]:
                rna_id = info[id_info_idx].replace("ID=rna-", "")
                gene_name = ''
                parent_gene_name = ''
                type = line[typeidx]

                for value in info:
                    if value.startswith("gene="):
                        gene_name = value.replace("gene=", "")
                        break

                if info[parent_gene_info_idx].startswith("Parent=gene-"):
                    parent_gene_name = info[parent_gene_info_idx].replace("Parent=gene-", "")

                cur.execute(
                    "insert into TRANSCRIPTS(TRANSCRIPT_ID,GENE_NAME,PARENT_GENE_NAME,TRANSCRIPT_TYPE,EXON_LENS) values('{}','{}','{}','{}',0)".format(rna_id, gene_name,
                                                                                                                                                       parent_gene_name, type))

                exons = []
                found_cds_start = False


            # if entry is of type exon update exon length
            # trancript entry will appear before any exon entries
            elif line[typeidx] == 'exon' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):
                exons.append([int(line[start_idx]), int(line[end_idx])])
                cur.execute("UPDATE TRANSCRIPTS SET EXON_LENS = EXON_LENS + {} WHERE TRANSCRIPT_ID = '{}'".format(((abs(int(line[end_idx]) - int(line[start_idx]))) + 1), rna_id))

            # if entry is of type cds update cds length
            elif line[typeidx] == 'CDS' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):

                if not found_cds_start:

                    if len(exons) == 0:
                        print(rna_id + ": error no exons")

                    #get cds start postion by counting number of nucletoides in exons before first cds occurs
                    cds_start = 0
                    if line[strand_idx] == '+':
                        for index in range(0, len(exons)):
                            if exons[index][0] <= int(line[start_idx]) <= exons[index][1]:
                                cds_start = int(line[start_idx]) - exons[index][0] + cds_start
                                break
                            else:
                                cds_start += exons[index][1] - exons[index][0] + 1
                    else:
                        for index in range(0, len(exons)):
                            if exons[index][1] >= int(line[end_idx]) >= exons[index][0]:
                                cds_start = exons[index][1] - int(line[end_idx]) + cds_start
                                break
                            else:
                                cds_start += exons[index][0] - exons[index][1] + 1


                    cur.execute("UPDATE TRANSCRIPTS SET CDS_START_INDEX = {} WHERE TRANSCRIPT_ID = '{}'".format(cds_start, rna_id))
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = {} WHERE TRANSCRIPT_ID = '{}'".format(abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

                    exons = []
                    found_cds_start = True

                else:
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = CDS_LEN + {} WHERE TRANSCRIPT_ID = '{}'".format(abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))


            elif "gene" == line[typeidx] or 'pseudogene' == line[typeidx]:
                gene_id = info[id_info_idx].replace("ID=gene-", "")
                gene_name = info[gene_name_info_idx].replace('Name=', '')

                gene_func = ''
                for value in info:
                    if value.startswith("gene_biotype="):
                        gene_func = value.replace("gene_biotype=", "")
                        break

                cur.execute("insert into GENES(GENE_ID,GENE_NAME,GENE_TYPE) values('{}','{}','{}')".format(gene_id, gene_name, gene_func))


    print("done")

f.close()

conn.commit()
conn.close()

#######################################################################################################################################################################################

# for mouse NCBI RefSeq
MOUSE_REFSEQ_DATADIRECTORY = '/filepath'
MOUSE_REFSEQ_GFFFILENAME = '/filepath/GCF_000001635.27_GRCm39_genomic.gff'
mouse_refseq_db = '/filepath/GRCm39_refseq.db'

os.chdir(MOUSE_REFSEQ_DATADIRECTORY)

conn = sqlite3.connect(mouse_refseq_db)
cur = conn.cursor()


#create database tables
cur.execute("""CREATE TABLE TRANSCRIPTS (
            TRANSCRIPT_ID text not null,
            TRANSCRIPT_TYPE text not null,
            GENE_NAME text not null,
            PARENT_GENE_NAME text not null,
            CDS_START_INDEX integer,
            CDS_LEN integer,
            EXON_LENS integer,
            primary key(TRANSCRIPT_ID)
)""")

cur.execute("""CREATE TABLE GENES (
            GENE_ID text not null,
            GENE_NAME text not null,
            GENE_TYPE text not null,
            primary key(GENE_ID)
)""")

#set collum numbers for gff file
typeidx = 2
infoidx = 8
start_idx = 3
end_idx = 4
strand_idx = 6

id_info_idx = 0
gene_name_info_idx = 2
parent_gene_info_idx = 1
parent_rna_idx = 1

with open(MOUSE_REFSEQ_GFFFILENAME, 'r') as f:
    for whole_line in f:
        if not whole_line.startswith('#'):
            line = whole_line.strip().split('\t')
            info = line[infoidx].split(';')

            # if entry is of type transcript add it to table
            if "RNA" in line[typeidx].upper() or 'transcript' in line[typeidx]:
                rna_id = info[id_info_idx].replace("ID=rna-", "")
                gene_name = ''
                parent_gene_name = ''
                type = line[typeidx]

                for value in info:
                    if value.startswith("gene="):
                        gene_name = value.replace("gene=", "")
                        break

                if info[parent_gene_info_idx].startswith("Parent=gene-"):
                    parent_gene_name = info[parent_gene_info_idx].replace("Parent=gene-", "")

                cur.execute(
                    "insert into TRANSCRIPTS(TRANSCRIPT_ID,GENE_NAME,PARENT_GENE_NAME,TRANSCRIPT_TYPE,EXON_LENS) values('{}','{}','{}','{}',0)".format(rna_id, gene_name,
                                                                                                                                                       parent_gene_name, type))

                exons = []
                found_cds_start = False


            # if entry is of type exon update exon length
            # trancript entry will appear before any exon entries
            elif line[typeidx] == 'exon' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):
                exons.append([int(line[start_idx]), int(line[end_idx])])
                cur.execute("UPDATE TRANSCRIPTS SET EXON_LENS = EXON_LENS + {} WHERE TRANSCRIPT_ID = '{}'".format(((abs(int(line[end_idx]) - int(line[start_idx]))) + 1), rna_id))

            # if entry is of type cds update cds length
            elif line[typeidx] == 'CDS' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):

                if not found_cds_start:

                    if len(exons) == 0:
                        print(rna_id + ": error no exons")

                    #get cds start postion by counting number of nucletoides in exons before first cds occurs
                    cds_start = 0
                    if line[strand_idx] == '+':
                        for index in range(0, len(exons)):
                            if exons[index][0] <= int(line[start_idx]) <= exons[index][1]:
                                cds_start = int(line[start_idx]) - exons[index][0] + cds_start
                                break
                            else:
                                cds_start += exons[index][1] - exons[index][0] + 1
                    else:
                        for index in range(0, len(exons)):
                            if exons[index][1] >= int(line[end_idx]) >= exons[index][0]:
                                cds_start = exons[index][1] - int(line[end_idx]) + cds_start
                                break
                            else:
                                cds_start += exons[index][0] - exons[index][1] + 1


                    cur.execute("UPDATE TRANSCRIPTS SET CDS_START_INDEX = {} WHERE TRANSCRIPT_ID = '{}'".format(cds_start, rna_id))
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = {} WHERE TRANSCRIPT_ID = '{}'".format(abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

                    exons = []
                    found_cds_start = True

                else:
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = CDS_LEN + {} WHERE TRANSCRIPT_ID = '{}'".format(abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))


            elif "gene" == line[typeidx] or 'pseudogene' == line[typeidx]:
                gene_id = info[id_info_idx].replace("ID=gene-", "")
                gene_name = info[gene_name_info_idx].replace('Name=', '')

                gene_func = ''
                for value in info:
                    if value.startswith("gene_biotype="):
                        gene_func = value.replace("gene_biotype=", "")
                        break

                cur.execute("insert into GENES(GENE_ID,GENE_NAME,GENE_TYPE) values('{}','{}','{}')".format(gene_id, gene_name, gene_func))


    print("done")

f.close()

conn.commit()
conn.close()

#######################################################################################################################################################################################

# for human Gencode
HUMAN_GENCODE_DATADIRECTORY = '/filepath'
HUMAN_GENCODE_GFFFILENAME = '/human_filepath/gencode.v46.chr_patch_hapl_scaff.annotation.gff3'
human_gencode_db = 'filepath/GRCh38_gencode.db'

# Change the current working directory
os.chdir(HUMAN_GENCODE_DATADIRECTORY)

# Connect to the SQLite database
conn = sqlite3.connect(human_gencode_db)
cur = conn.cursor()

# Create tables for storing transcripts and genes
cur.execute("""CREATE TABLE IF NOT EXISTS TRANSCRIPTS (
            TRANSCRIPT_ID text not null,
            TRANSCRIPT_TYPE text not null,
            GENE_NAME text not null,
            PARENT_GENE_NAME text not null,
            CDS_START_INDEX integer,
            CDS_LEN integer,
            EXON_LENS integer,
            primary key(TRANSCRIPT_ID)
)""")

cur.execute("""CREATE TABLE IF NOT EXISTS GENES (
            GENE_ID text not null,
            GENE_NAME text not null,
            GENE_TYPE text not null,
            primary key(GENE_ID)
)""")

# Index positions for parsing the GFF file
typeidx = 2
infoidx = 8
start_idx = 3
end_idx = 4
strand_idx = 6

gene_name_info_idx = 2
parent_gene_info_idx = 1
parent_rna_idx = 1

# Open and read the GFF file
with open(HUMAN_GENCODE_GFFFILENAME, 'r') as f:
    rna_id = None
    found_cds_start = False
    exons = []

    for whole_line in f:
        if not whole_line.startswith('#'):  # Skip comments
            line = whole_line.strip().split('\t')
            info = line[infoidx].split(';')

            # Process RNA/transcript lines
            if "RNA" in line[typeidx].upper() or 'transcript' in line[typeidx]:
                for value in info:
                    if value.startswith("gene_name="):
                        gene_name = value.replace("gene_name=", "")
                        break
                
                for value in info:
                    if value.startswith("Parent="):
                        parent_gene_name = value.replace("Parent=", "")
                        break
                
                for value in info:
                    if value.startswith("transcript_type="):
                        type = value.replace("transcript_type=", "")
                        break
                
                for value in info:
                    if value.startswith("ID=ENST"):
                        rna_id = value.replace("ID=", "")
                        break

                cur.execute(
                    "INSERT INTO TRANSCRIPTS(TRANSCRIPT_ID, GENE_NAME, PARENT_GENE_NAME, TRANSCRIPT_TYPE, EXON_LENS) VALUES (?, ?, ?, ?, 0)", (rna_id, gene_name, parent_gene_name, type)
                )

                exons = []
                found_cds_start = False

            # Process exon lines
            elif line[typeidx] == 'exon' and rna_id and rna_id == info[parent_rna_idx].replace("Parent=", ""):
                exons.append([int(line[start_idx]), int(line[end_idx])])
                cur.execute("UPDATE TRANSCRIPTS SET EXON_LENS = EXON_LENS + ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

            # Process CDS lines
            elif line[typeidx] == 'CDS' and rna_id and rna_id == info[parent_rna_idx].replace("Parent=", ""):
                if not found_cds_start:
                    if len(exons) == 0:
                        print(rna_id + ": error no exons")

                    cds_start = 0
                    if line[strand_idx] == '+':
                        for index in range(len(exons)):
                            if exons[index][0] <= int(line[start_idx]) <= exons[index][1]:
                                cds_start = int(line[start_idx]) - exons[index][0] + cds_start
                                break
                            else:
                                cds_start += exons[index][1] - exons[index][0] + 1
                    else:
                        for index in range(len(exons)):
                            if exons[index][1] >= int(line[end_idx]) >= exons[index][0]:
                                cds_start = exons[index][1] - int(line[end_idx]) + cds_start
                                break
                            else:
                                cds_start += exons[index][0] - exons[index][1] + 1

                    cur.execute("UPDATE TRANSCRIPTS SET CDS_START_INDEX = ? WHERE TRANSCRIPT_ID = ?", (cds_start, rna_id))
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

                    exons = []
                    found_cds_start = True
                else:
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = CDS_LEN + ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

            # Process gene and pseudogene lines
            elif "gene" == line[typeidx] or 'pseudogene' == line[typeidx]:
                gene_name = info[gene_name_info_idx].replace('Name=', '')
                
                for value in info:
                    if value.startswith("gene_id="):
                        gene_id = value.replace("gene_id=", "")
                        break

                gene_func = ''
                for value in info:
                    if value.startswith("gene_biotype="):
                        gene_func = value.replace("gene_biotype=", "")
                        break

                cur.execute("INSERT INTO GENES(GENE_ID, GENE_NAME, GENE_TYPE) VALUES (?, ?, ?)", (gene_id, gene_name, gene_func))

    print("done")

# Close the file and commit changes to the database
conn.commit()
conn.close()

#######################################################################################################################################################################################

# for mouse Gencode
MOUSE_GENCODE_DATADIRECTORY = '/filepath'
MOUSE_GENCODE_GFFFILENAME = '/filepath/gencode.vM35.chr_patch_hapl_scaff.annotation.gff3'
mouse_gencode_db = '/filepath/GRCm39_gencode.db'

# Change the current working directory
os.chdir(MOUSE_GENCODE_DATADIRECTORY)

# Connect to the SQLite database
conn = sqlite3.connect(mouse_gencode_db)
cur = conn.cursor()

# Create tables for storing transcripts and genes
cur.execute("""CREATE TABLE IF NOT EXISTS TRANSCRIPTS (
            TRANSCRIPT_ID text not null,
            TRANSCRIPT_TYPE text not null,
            GENE_NAME text not null,
            PARENT_GENE_NAME text not null,
            CDS_START_INDEX integer,
            CDS_LEN integer,
            EXON_LENS integer,
            primary key(TRANSCRIPT_ID)
)""")

cur.execute("""CREATE TABLE IF NOT EXISTS GENES (
            GENE_ID text not null,
            GENE_NAME text not null,
            GENE_TYPE text not null,
            primary key(GENE_ID)
)""")

# Index positions for parsing the GFF file
typeidx = 2
infoidx = 8
start_idx = 3
end_idx = 4
strand_idx = 6

gene_name_info_idx = 2
parent_gene_info_idx = 1
parent_rna_idx = 1

# Open and read the GFF file
with open(GFFFILENAME, 'r') as f:
    rna_id = None
    found_cds_start = False
    exons = []

    for whole_line in f:
        if not whole_line.startswith('#'):  # Skip comments
            line = whole_line.strip().split('\t')
            info = line[infoidx].split(';')

            # Process RNA/transcript lines
            if "RNA" in line[typeidx].upper() or 'transcript' in line[typeidx]:
                for value in info:
                    if value.startswith("gene_name="):
                        gene_name = value.replace("gene_name=", "")
                        break
                
                for value in info:
                    if value.startswith("Parent="):
                        parent_gene_name = value.replace("Parent=", "")
                        break
                
                for value in info:
                    if value.startswith("transcript_type="):
                        type = value.replace("transcript_type=", "")
                        break
                
                for value in info:
                    if value.startswith("ID=ENSMUST"):
                        rna_id = value.replace("ID=", "")
                        break

                cur.execute(
                    "INSERT INTO TRANSCRIPTS(TRANSCRIPT_ID, GENE_NAME, PARENT_GENE_NAME, TRANSCRIPT_TYPE, EXON_LENS) VALUES (?, ?, ?, ?, 0)", (rna_id, gene_name, parent_gene_name, type)
                )

                exons = []
                found_cds_start = False

            # Process exon lines
            elif line[typeidx] == 'exon' and rna_id and rna_id == info[parent_rna_idx].replace("Parent=", ""):
                exons.append([int(line[start_idx]), int(line[end_idx])])
                cur.execute("UPDATE TRANSCRIPTS SET EXON_LENS = EXON_LENS + ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

            # Process CDS lines
            elif line[typeidx] == 'CDS' and rna_id and rna_id == info[parent_rna_idx].replace("Parent=", ""):
                if not found_cds_start:
                    if len(exons) == 0:
                        print(rna_id + ": error no exons")

                    cds_start = 0
                    if line[strand_idx] == '+':
                        for index in range(len(exons)):
                            if exons[index][0] <= int(line[start_idx]) <= exons[index][1]:
                                cds_start = int(line[start_idx]) - exons[index][0] + cds_start
                                break
                            else:
                                cds_start += exons[index][1] - exons[index][0] + 1
                    else:
                        for index in range(len(exons)):
                            if exons[index][1] >= int(line[end_idx]) >= exons[index][0]:
                                cds_start = exons[index][1] - int(line[end_idx]) + cds_start
                                break
                            else:
                                cds_start += exons[index][0] - exons[index][1] + 1

                    cur.execute("UPDATE TRANSCRIPTS SET CDS_START_INDEX = ? WHERE TRANSCRIPT_ID = ?", (cds_start, rna_id))
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

                    exons = []
                    found_cds_start = True
                else:
                    cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = CDS_LEN + ? WHERE TRANSCRIPT_ID = ?", (abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))

            # Process gene and pseudogene lines
            elif "gene" == line[typeidx] or 'pseudogene' == line[typeidx]:
                gene_name = info[gene_name_info_idx].replace('Name=', '')
                
                for value in info:
                    if value.startswith("gene_id="):
                        gene_id = value.replace("gene_id=", "")
                        break

                gene_func = ''
                for value in info:
                    if value.startswith("gene_biotype="):
                        gene_func = value.replace("gene_biotype=", "")
                        break

                cur.execute("INSERT INTO GENES(GENE_ID, GENE_NAME, GENE_TYPE) VALUES (?, ?, ?)", (gene_id, gene_name, gene_func))

    print("done")

# Close the file and commit changes to the database
conn.commit()
conn.close()

####################################################################################################################################################################################### 
####################################################################################################################################################################################### 
####################################################################################################################################################################################### 

# sm_sequence.py

# this script requires the databases created in the previous section. 
# this script identifies all transcripts containing the specified Sm-site sequence, and outputs two files for each type:
#   1) essential - sequences 3' of the identified Sm-site
#   2) nonessential - sequences 5' of the identified Sm-site
# outputs are used as input into RNAfold to predict secondary folding around each Sm-site. 
# the script is written such that the databases and transcriptome fasta files can be substituted in order to output sequences surrounding the sm-sites for each species.
# there is one script written for NCBI Refseq .fna and one script writted for Gencode .fa

# Import necessary modules
from Bio import SeqIO  # For parsing sequence files
import pandas as pd  # For data manipulation
from Bio.Seq import Seq  # For handling sequences
from Bio import SeqUtils  # For sequence utility functions
import os  # For operating system dependent functionality
import numpy as np  # For numerical operations
import sqlite3  # For SQLite database operations

# the script immediately following is written for processing NCBI Refseq .fna

# Define constants for minimum and maximum lengths
MINLEN = 16
MAXLEN = 50

# Define file paths and directories
RESULTSDIRECTORY = '/filepath'
RNAFILENAME = '/filepath to transcriptome' 
# for human_refseq: GCF_000001405.40_GRCh38_rna.fna 
# for mouse_refseq: GCF_000001635.27_GRCm39_rna.fna

database = '/filepath to specific database created in the scripts above'
# for human_refseq: GRCh38_refseq.db
# for mouse_refseq: GRCm39_refseq.db

U2TYPE = '/filepath/u2_filename'
U5TYPE = '/filepath/u5_filename'
U1TYPE = '/filepath/u1_filename'
U7TYPE = '/filepath/u7_filename'
NONCANONICAL = '/filepath/noncanonical_filename'

# Change current working directory to the results directory
os.chdir(RESULTSDIRECTORY)

# Connect to the SQLite database and create a cursor
conn = sqlite3.connect(database)
cur = conn.cursor()

# Function to generate candidate sequences
def generate_candidates(search_seq, output_file_prefix):
    region_count = 0
    num_mismatch_count = 0

    search_seq = str(search_seq)
    output_file_essential = open(output_file_prefix + '_essential.txt', "w")
    output_file_nonessential = open(output_file_prefix + '_nonessential.txt', "w")
    print("Generating Candidates for {}...".format(output_file_prefix))
    num_found = 0

    # Iterate through the RNA file
    for seq_record in SeqIO.parse(RNAFILENAME, "fasta"):
        str_trans_seq = str(seq_record.seq)
        transcript_matches = SeqUtils.nt_search(str_trans_seq, search_seq)

        # If there are matches for the search sequence in the transcript
        if len(transcript_matches) > 1:
            trans_id = seq_record.id
            desc = seq_record.description.split(",")
            type_of_rna_by_rna_file = desc[-1].split(";")[0].strip()  # Get the type of RNA

            # Query the database for transcript information
            cur.execute("SELECT TRANSCRIPT_TYPE, GENE_NAME, PARENT_GENE_NAME, CDS_START_INDEX, CDS_LEN, EXON_LENS FROM TRANSCRIPTS WHERE TRANSCRIPT_ID = ?", [trans_id])
            trans_results = cur.fetchone()

            transcript_type = trans_results[0]
            gene_name = trans_results[1]
            exact_gene_name = trans_results[2]

            # Determine if the sequence is predicted
            predict = 'predicted' if 'PREDICTED' in seq_record.description else 'known'

            # Check if CDS is annotated and if transcript length matches exon length
            has_regions = False
            if trans_results[3] is not None and trans_results[4] is not None and len(seq_record.seq) == trans_results[5]:
                has_regions = True
                cds_start = trans_results[3]
                cds_end = trans_results[4] + cds_start
                region_count += 1

            if len(seq_record.seq) != trans_results[5]:
                num_mismatch_count += 1

            # Get gene type using gene name or exact gene name
            if trans_results is not None:
                gene_type_query = exact_gene_name if exact_gene_name else gene_name
                cur.execute("SELECT GENE_TYPE FROM GENES WHERE GENE_ID = ?", [gene_type_query])
                gene_type = cur.fetchone()

            gene_type = gene_type[0] if gene_type is not None else 'N/A'

            # Get all transcript matches
            transcript_matches = transcript_matches[1:]
            transcript_matches = [x + len(search_seq) for x in transcript_matches]
            str_trans_seq = str_trans_seq.replace('T', 'U').replace('t', 'u')  # Replace T with U
            trans_seq_len = len(str_trans_seq)

            for match_index in transcript_matches:
                num_found += 1

                # Get sequence after the match
                if (trans_seq_len - match_index) < MINLEN:
                    continue
                elif (trans_seq_len - match_index) < MAXLEN:
                    candidate_seq = str_trans_seq[match_index:]
                else:
                    candidate_seq = str_trans_seq[match_index:match_index + MAXLEN]

                # Get sequence before the match
                if match_index - len(search_seq) < 9:
                    continue
                elif match_index - len(search_seq) < 200:
                    prior_seq = str_trans_seq[:match_index - len(search_seq)]
                else:
                    prior_seq = str_trans_seq[match_index - len(search_seq) - 200:match_index - len(search_seq)]

                # Determine the region of the sequence
                region = 'N/A'
                distance = 'N/A'
                if has_regions:
                    if match_index <= cds_start and cds_start != 0:
                        region = '5\'UTR'
                        distance = str(cds_start - match_index)
                    elif match_index - len(search_seq) <= cds_end and cds_end != cds_start:
                        region = 'CDS'
                    elif cds_end < trans_seq_len:
                        region = '3\'UTR'
                        distance = str(match_index - len(search_seq) - cds_end)
                elif len(seq_record.seq) != trans_results[5]:
                    region = 'Exon_Mismatch'

                header = '{}|{}|{}|{}|{}|{}|{}|{}|{}|{}'.format(trans_id, match_index, region, transcript_type, gene_name, exact_gene_name, gene_type, type_of_rna_by_rna_file, predict, distance)

                # Write the results to the output files
                # essential is for 3' stem loop, nonessential is for 5' stem loop 
                write_to_output_file(header, candidate_seq, output_file_essential)
                write_to_output_file(header, prior_seq, output_file_nonessential)

    print("{} Candidates Found for {}...".format(num_found, output_file_prefix))
    print('regions ' + str(region_count))
    print('exon mismatch ' + str(num_mismatch_count))

    output_file_essential.close()
    output_file_nonessential.close()

# Function to write to the output file
def write_to_output_file(header, seq, output_file):
    output_file.write('>' + header + '\n')
    output_file.write(seq + '\n')

# Generate candidate sequences
generate_candidates('ATTTTTG', U2TYPE)
generate_candidates('ATTTTTTG', U5TYPE)
generate_candidates('ATTTGTG', U1TYPE)
generate_candidates('ATTTGTCTAG', U7TYPE)
generate_candidates('ATNTKTN', NONCANONICAL)

# Close the database connection
conn.close()

#######################################################################################################################################################################################

# this script following is written to for processing Gencode .fa

# Define constants for minimum and maximum lengths
MINLEN = 16
MAXLEN = 50

# Define file paths and directories
RESULTSDIRECTORY = '/filepath'
RNAFILENAME = '/filepath to transcriptome' 
# for human_gencode: gencode.v46.transcripts.fa 
# for mouse_gencode: gencode.vM35.transcripts.fa

database = '/filepath to specific database created in the scripts above'
# for human_gencode: GRCh38_gencode.db
# for mouse_gencode: GRCm39_gencode.db

U2TYPE = '/filepath/u2_filename'
U5TYPE = '/filepath/u5_filename'
U1TYPE = '/filepath/u1_filename'
U7TYPE = '/filepath/u7_filename'
NONCANONICAL = '/filepath/noncanonical_filename'

# Change current working directory to the results directory
os.chdir(RESULTSDIRECTORY)

# Connect to the SQLite database and create a cursor
conn = sqlite3.connect(database)
cur = conn.cursor()

# Function to generate candidate sequences
def generate_candidates(search_seq, output_file_prefix):
    region_count = 0
    num_mismatch_count = 0

    search_seq = str(search_seq)
    output_file_essential = open(output_file_prefix + '_essential.txt', "w")
    output_file_nonessential = open(output_file_prefix + '_nonessential.txt', "w")
    print("Generating Candidates for {}...".format(output_file_prefix))
    num_found = 0

    # Iterate through the RNA file
    for seq_record in SeqIO.parse(RNAFILENAME, "fasta"):
        str_trans_seq = str(seq_record.seq)
        transcript_matches = SeqUtils.nt_search(str_trans_seq, search_seq)

        # If there are matches for the search sequence in the transcript
        if len(transcript_matches) > 1:
            # Split the header and remove '>'
            header_parts = seq_record.id.lstrip('>').split('|')
            trans_id = header_parts[0].strip()

            # Query the database for transcript information
            cur.execute("SELECT TRANSCRIPT_TYPE, GENE_NAME, PARENT_GENE_NAME, CDS_START_INDEX, CDS_LEN, EXON_LENS FROM TRANSCRIPTS WHERE TRANSCRIPT_ID = ?", [trans_id])
            trans_results = cur.fetchone()

            # Check if trans_results is None
            if trans_results is None:
                print(f"No transcript information found for transcript ID: {trans_id}")
                continue  # Skip processing this transcript

            transcript_type = trans_results[0]
            gene_name = trans_results[1]
            exact_gene_name = trans_results[2]

            # Determine if the sequence is predicted
            predict = 'predicted' if 'PREDICTED' in seq_record.description else 'known'

            # Check if CDS is annotated and if transcript length matches exon length
            has_regions = False
            if trans_results[3] is not None and trans_results[4] is not None and len(seq_record.seq) == trans_results[5]:
                has_regions = True
                cds_start = trans_results[3]
                cds_end = trans_results[4] + cds_start
                region_count += 1

            if len(seq_record.seq) != trans_results[5]:
                num_mismatch_count += 1

            # Get gene type using gene name or exact gene name
            if trans_results is not None:
                gene_type_query = exact_gene_name if exact_gene_name else gene_name
                cur.execute("SELECT GENE_TYPE FROM GENES WHERE GENE_ID = ?", [gene_type_query])
                gene_type = cur.fetchone()

            gene_type = gene_type[0] if gene_type is not None else 'N/A'

            # Get all transcript matches
            transcript_matches = transcript_matches[1:]
            transcript_matches = [x + len(search_seq) for x in transcript_matches]
            str_trans_seq = str_trans_seq.replace('T', 'U').replace('t', 'u')  # Replace T with U
            trans_seq_len = len(str_trans_seq)

            for match_index in transcript_matches:
                num_found += 1

                # Get sequence after the match
                if (trans_seq_len - match_index) < MINLEN:
                    continue
                elif (trans_seq_len - match_index) < MAXLEN:
                    candidate_seq = str_trans_seq[match_index:]
                else:
                    candidate_seq = str_trans_seq[match_index:match_index + MAXLEN]

                # Get sequence before the match
                if match_index - len(search_seq) < 9:
                    continue
                elif match_index - len(search_seq) < 200:
                    prior_seq = str_trans_seq[:match_index - len(search_seq)]
                else:
                    prior_seq = str_trans_seq[match_index - len(search_seq) - 200:match_index - len(search_seq)]

                # Determine the region of the sequence
                region = 'N/A'
                distance = 'N/A'
                if has_regions:
                    if match_index <= cds_start and cds_start != 0:
                        region = '5\'UTR'
                        distance = str(cds_start - match_index)
                    elif match_index - len(search_seq) <= cds_end and cds_end != cds_start:
                        region = 'CDS'
                    elif cds_end < trans_seq_len:
                        region = '3\'UTR'
                        distance = str(match_index - len(search_seq) - cds_end)
                elif len(seq_record.seq) != trans_results[5]:
                    region = 'Exon_Mismatch'

                header = '{}|{}|{}|{}|{}|{}|{}|{}|{}'.format(trans_id, match_index, region, transcript_type, gene_name, exact_gene_name, gene_type, predict, distance)

                # Write the results to the output files
                # essential is for 3' stem loop, nonessential is for 5' stem loop 
                write_to_output_file(header, candidate_seq, output_file_essential)
                write_to_output_file(header, prior_seq, output_file_nonessential)

    print("{} Candidates Found for {}...".format(num_found, output_file_prefix))
    print('regions ' + str(region_count))
    print('exon mismatch ' + str(num_mismatch_count))

    output_file_essential.close()
    output_file_nonessential.close()

# Function to write to the output file
def write_to_output_file(header, seq, output_file):
    output_file.write('>' + header + '\n')
    output_file.write(seq + '\n')

# Generate candidate sequences
generate_candidates('ATTTTTG', U2TYPE)
generate_candidates('ATTTTTTG', U5TYPE)
generate_candidates('ATTTGTG', U1TYPE)
generate_candidates('ATTTGTCTAG', U7TYPE)
generate_candidates('ATNTKTN', NONCANONICAL)

# Close the database connection
conn.close()

#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################

# rnafold

# the following is the script to perform RNAfold secondary structure prediction on sm_sequence.py output files.

# Summary
# The script prints status messages to the terminal.
# It runs RNAfold to analyze RNA sequences from specific input files.
# The --noPS option is used to skip generating PostScript files.
# Output files are generated in the specified directories.
# Some sections are commented out and can be uncommented if needed.

# Need to install RNAfold from http://www.tbi.univie.ac.at/RNA/index.html and upload into supercomputer
# install tutorial is here : https://www.tbi.univie.ac.at/RNA/tutorial/ 
# unzip using the following command : tar -zxf ViennaRNA-2.6.4.tar.gz
# check that it unzips using: ls
# Change directory to : cd ViennaRNA-2.6.4
# make new directory : mkdir -p ~/Tutorial/Progs/VRP
# configure the buld using : ./configure --prefix=/users/PAS0631/frankenchicken/identify_sm_sites/rnafold/programs
# make install
# combine into a new directory : cd ~/Tutorial/Progs/
# cp VRP/share/ViennaRNA/bin/* .

# Print a message to the terminal
echo "u2type rnafold essential"

# Run RNAfold on the essential 3' stem loop RNA sequences (5' is nonessential)
# The --noPS option prevents the generation of PostScript plots
# Input: .txt output from sm_sequence.py
# Output: rnafold_output_sm_sequence.txt
RNAfold --noPS < /filepath/sm_sequence_output.txt > /filepath/rnafold_output_sm_sequence.txt

# this will yield a text file providing secondary structure prediction notation for each 'essential' (3') and 'nonessential' (5') sequence flanking the identified Sm-site.

#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################

# stem_loops.py

# From RNAs containing Sm-site sequence, determine if stem loops are also present

# Need to run in python environment created to make the database in create_database.py 
# module load miniconda3
# conda smsite_env

# Bio, pandas, numpy, sqlite3, and os are imported for handling sequences, dataframes, numeric operations, databases, and directory operations respectively. 
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqUtils
import os
import numpy as np
import sqlite3

# Define input and output filenames for RNAfold files
# The paths can be updated to use Ensembl annotations if required
# Paths for input and output files are defined.
# this script can take all output files generated from the RNAfold output, regardless of species or annotation.
# this script takes both essential (3') and nonessential (5') RNAfold out[ut files for each Sm-site type, 
# to determine if folding occurs in both, and yields a single output file retaining only the essential (3') sequence for easily finding Sm-sites in cDNA sequences.

U2INPUTFILENAME = '/filepath/U2 RNAfold output filename'
U2OUTPUTFILENAME = '/filepath/U2 output filename .txt'
U5INPUTFILENAME = '/filepath/U5 RNAfold output filename'
U5OUTPUTFILENAME = '/filepath/U5 output filename .txt'
U1INPUTFILENAME = '/filepath/U1 RNAfold output filename'
U1OUTPUTFILENAME = '/filepath/U1 output filename .txt'
U7INPUTFILENAME = '/filepath/U7 RNAfold output filename'
U7OUTPUTFILENAME = '/filepath/U7 output filename .txt'
NCINPUTFILENAME = '/filepath/Noncanonical RNAfold output filename'
NCOUTPUTFILENAME = '/filepath/Noncanonical output filename .txt'

# Index positions for parsing the output files
# Index positions for various fields in the output files are defined.
trans_id_idx = 0
match_idx = 1
region_idx = 2
trans_type_idx = 3
gene_name_idx = 4
exact_gene_name_idx = 5
gene_type_idx = 6
trans_type_by_rna_file_idx = 7  # only in NIH
predicted_transcript_idx = 8  # only in NIH

# Function to parse RNAfold output files and identify stem-loop structures
# Reads the RNAfold output files, checks for stem-loop structures in the sequences, and writes sequences with stem-loop structures to an output file.
# The function handles essential (3' sequence) and nonessential (5' sequence) files separately but processes them in parallel for each sequence.
def parse_rnafold(input_file_prefix, output_file_name):
    input_file_essential = input_file_prefix + '_rnafold_gencode_essential.txt'
    input_file_nonessential = input_file_prefix + '_rnafold_gencode_nonessential.txt'

    line_num = 0
    print('Parsing output for ' + input_file_essential)
    output_file = open(output_file_name, "w")

    with open(input_file_nonessential, 'r') as nonessential_file, open(input_file_essential, 'r') as essential_file:
        for line_nonessential, line_essential in zip(nonessential_file, essential_file):
            line_num += 1

            # This is the header line
            if line_num % 3 == 1:
                header = line_essential
                if line_nonessential != header:
                    print("error: headers do not match")

            # This is the sequence line
            elif line_num % 3 == 2:
                seq_essential = line_essential

            # This is the structure line
            else:
                structure = line_essential.split(' ')[0]
                stem_loop = False

                # Check for stem-loop structures in the first 10 nucleotides
                for index in range(0, len(structure)):
                    if index == 9:
                        break
                    elif structure[index] == '(':
                        if structure[index:].startswith('((((') or structure[index:].startswith('(((.') or \
                                structure[index:].startswith('((.(') or structure[index:].startswith('(.(('):
                            stem_loop = True
                        break

                structure = line_nonessential.split(' ')[0]
                nonessential_stem_loop = False

                # Check for stem-loop structures in the first 20 nucleotides
                for index in range(0, len(structure)):
                    if index == 19:
                        break
                    elif structure[index] == '(':
                        if structure[index:].startswith('((((') or structure[index:].startswith('(((.') or \
                                structure[index:].startswith('((.(') or structure[index:].startswith('(.(('):
                            nonessential_stem_loop = True
                        break

                # Write to output file if both nonessential and essential sequences have stem-loop structures
                if stem_loop and nonessential_stem_loop:
                    write_to_output_file(header, seq_essential, output_file)

    output_file.close()

# Function to write sequences with stem-loop structures to the output file
# Writes the header and sequence to the output file if a stem-loop structure is found.
def write_to_output_file(header, seq, output_file):
    output_file.write(header)
    output_file.write(seq)

# Function to get statistics from the parsed RNAfold output
# Extracts statistics from the RNAfold output files regarding transcript types and regions.
# Uses helper function add_to_dict to manage counts of transcript types and regions.
def get_stats_nih(output_file_name):
    transcript_types = {}
    transcript_type_genes = []
    regions = {}
    regions_genes = []

    with open(output_file_name, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line.split('|')
                trans_type = '{}-{}'.format(line[gene_name_idx], line[predicted_transcript_idx]).replace('\n', '')

                if trans_type not in transcript_type_genes:
                    transcript_type_genes.append(trans_type)
                    add_to_dict(transcript_types, '{}-{}'.format(line[trans_type_by_rna_file_idx], line[predicted_transcript_idx]).replace('\n', ''))

                trans_type_region = '{}-{}-{}'.format(line[region_idx], line[gene_name_idx], line[predicted_transcript_idx]).replace('\n', '')
                if trans_type_region not in regions_genes:
                    regions_genes.append(trans_type_region)
                    add_to_dict(regions, '{}:{}-{}'.format(line[region_idx], line[trans_type_by_rna_file_idx], line[predicted_transcript_idx]).replace('\n', ''))

    print()
    print(output_file_name)
    print('---------------------------------------------------------------------------------------------------------------------------------------')
    print('transcript_types:')
    print(transcript_types)
    print('---------------------------------------------------------------------------------------------------------------------------------------')
    print('regions:')
    print(regions)
    print()

# Helper function to add keys to a dictionary and increment their count
# Adds keys to a dictionary and increments their counts.
def add_to_dict(dict, key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1

# Parse the RNAfold output files for Sm-site types
# parse_rnafold is called forSm-site type RNAfold output files.
parse_rnafold(U2INPUTFILENAME, U2OUTPUTFILENAME)
parse_rnafold(U5INPUTFILENAME, U5OUTPUTFILENAME)
parse_rnafold(U1INPUTFILENAME, U1OUTPUTFILENAME)
parse_rnafold(U7INPUTFILENAME, U7OUTPUTFILENAME)
parse_rnafold(NCINPUTFILENAME, NCOUTPUTFILENAME)

# these output files were then cross-referenced and combined using custom R scripts which can be made available upon request.

#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
