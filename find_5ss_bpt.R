
# Python script to query sequences, controlling for mismatches or deletions, from from fasta files 

# Import necessary modules
from Bio import SeqIO  # For parsing sequence files
import os  # For operating system dependent functionality
import sqlite3  # For SQLite database operations
import regex as re  # For fuzzy matching

# Define file paths and directories
RESULTSDIRECTORY = 'PATH TO RESULTS DIRECTORY'
RNAFILENAME = 'PATH TO TRANSCRIPTOME FASTA'  
    # for human gencode     gencode.v46.transcripts.fa
    # for human refseq      GRCh38_latest_rna.fna
    # for mouse gencode     gencode.vM35.transcripts.fa
    # for mouse refseq      GCF_000001635.27_GRCm39_rna.fna
database = 'PATHT TO DATABASE'  
    # for human gencode     GRCh38_gencode.db
    # for human refseq      hg38_nih_rna.db
    # for mouse gencode     GRCm39_gencode.db
    # for mouse refseq      GRCm39_nih_rna.db

# Change current working directory to the results directory
os.chdir(RESULTSDIRECTORY)

# Connect to the SQLite database and create a cursor
conn = sqlite3.connect(database)
cur = conn.cursor()

# Function to generate candidate sequences with mismatches or deletions
def generate_candidates_with_specific_edits(search_seq, output_file_prefix, max_mismatches=1, max_deletions=1):
    """
    Search for sequences matching the search pattern with allowed mismatches and deletions.
    Writes a single output file containing the matching sequences with relevant transcript details.
    """
    search_seq = search_seq.upper()
    ambiguous_pattern = search_seq.replace('M', '[AC]').replace('R', '[AG]').replace('U', 'T').replace('K', '[GT]').replace('Y', '[CT]').replace('W', '[AT]').replace('N', '[ATCG]')

    # Allow mismatches
    mismatch_pattern = f"({ambiguous_pattern}){{e<={max_mismatches}}}"
    
    # Allow single deletions
    deletion_pattern = f"({ambiguous_pattern}){{s<={max_deletions}}}"

    # Combine patterns to allow mismatches OR deletions
    regex_pattern = f"{mismatch_pattern}|{deletion_pattern}"

    # Output file for storing matches
    with open(output_file_prefix + '_specific_edits_matches.txt', "w") as output_file:
        output_file.write("Transcript_ID\tGene_Name\tExact_Gene_Name\tTranscript_Type\tRegion\tPredicted\tMatched_Sequence\n")

        print(f"Generating Candidates for {output_file_prefix} with {max_mismatches} mismatches and {max_deletions} deletions allowed...")

        num_found = 0
        region_count = 0
        num_mismatch_count = 0

        # Iterate through the RNA file
        for seq_record in SeqIO.parse(RNAFILENAME, "fasta"):
            str_trans_seq = str(seq_record.seq).upper()

            try:
                matches = [(m.start(), m.end(), m.group()) for m in re.finditer(regex_pattern, str_trans_seq, overlapped=True)]
            except re.error as e:
                print(f"Regex error: {e}")
                print(f"Regex pattern: {regex_pattern}")
                return

            if matches:
                # Extract transcript information
                header_parts = seq_record.id.lstrip('>').split('|')
                trans_id = header_parts[0].strip()

                cur.execute("SELECT TRANSCRIPT_TYPE, GENE_NAME, PARENT_GENE_NAME, CDS_START_INDEX, CDS_LEN, EXON_LENS FROM TRANSCRIPTS WHERE TRANSCRIPT_ID = ?", [trans_id])
                trans_results = cur.fetchone()

                if trans_results is None:
                    print(f"No transcript information found for transcript ID: {trans_id}")
                    continue

                transcript_type = trans_results[0]
                gene_name = trans_results[1]
                exact_gene_name = trans_results[2]

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
                gene_type_query = exact_gene_name if exact_gene_name else gene_name
                cur.execute("SELECT GENE_TYPE FROM GENES WHERE GENE_ID = ?", [gene_type_query])
                gene_type = cur.fetchone()

                gene_type = gene_type[0] if gene_type is not None else 'N/A'

                # Replace T with U for RNA sequence
                str_trans_seq = str_trans_seq.replace('T', 'U').replace('t', 'u')  # Replace T with U for RNA
                trans_seq_len = len(str_trans_seq)

                # Process matches
                for match_start, match_end, match_seq in matches:
                    num_found += 1

                    # Determine region of the sequence
                    region = 'N/A'
                    if has_regions:
                        if match_start <= cds_start and cds_start != 0:
                            region = '5\'UTR'
                        elif match_start <= cds_end and cds_end != cds_start:
                            region = 'CDS'
                        elif cds_end < trans_seq_len:
                            region = '3\'UTR'
                    elif len(seq_record.seq) != trans_results[5]:
                        region = 'Exon_Mismatch'

                    # Extract flanking sequences
                    flank_start = max(0, match_start - 20)  # Ensure start index is not negative
                    flank_end = min(trans_seq_len, match_end + 20)  # Ensure end index doesn't exceed sequence length
                    flanking_sequence = str_trans_seq[flank_start:flank_end]

                    # Write to output file with relevant information
                    output_file.write(f"{trans_id}\t{gene_name}\t{exact_gene_name}\t{transcript_type}\t{region}\t{predict}\t{flanking_sequence}\n")

        print(f"{num_found} Candidates Found for {output_file_prefix} with {max_mismatches} mismatches and {max_deletions} deletions allowed...")
        print(f"Regions: {region_count}")
        print(f"Exon Mismatch: {num_mismatch_count}")

# Example: Generate candidates with mismatches for the search sequence
generate_candidates('MAGGURAGK', 'human_gencode_fivess', max_mismatches=0, max_deletions=0)
generate_candidates('MAGGURAGK', 'human_gencode_fivess_1mis_1del', max_mismatches=1, max_deletions=1)
generate_candidates('YMYUNACW', 'human_gencode_U2U12N', max_mismatches=0, max_deletions=0)
generate_candidates('YUNAYYYYY', 'human_gencode_brchpt_short', max_mismatches=0, max_deletions=0)

# Close the database connection
conn.close()
