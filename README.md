# sm_ring_assembly_mrna

These are scripts used to identify and characterize Sm-site containing RNAs within the human and mouse transcriptomes.

a) to identify Sm-sites in human and mouse transcriptomes, use find_smsites.py

Custom R script used to combine the outputs from find_smsites.py into a single file can be accessed in combine_smsite_output.R, tables in smsite_combined.zip

Custom R scripts used to make a single master file giving the frequency of Sm-sites per gene_id -> human_smsite_utype.R, mouse_smsite_utype.R

Custom R scripts used to describe the RNAs identfiied with Sm-sites -> insilico.R, master tables in smsite_master.zip



b) to identify 5'ss and branchpoint sequences in supplied transcriptomes use -> find_5ss_bpt.py, output files in human_5ss_bpt.zip and mouse_5ss_bpt.zip



c) to create a Simplified Annotation File containing intron coordinates for counting intronic reads only -> intron_saf.R



d) to process raw sequencing files from publication, use sm_rip_seq_processing.sh

To visualize processing of reads -> multiqc_report_human.html, multiqc_report_mouse.html

To combine featureCounts outputs into a single table -> make_master.py

To obtain Wald test comparisons of featureCounts outputs -> deseq2_analysis.R (for SmB and Y12 comparison, built into DESeq_Y12_SmB.R)



e) characterize Deseq2 generated Wald test comparisons:

For Sm-RIP vs polyA-RNA -> sm_rip_vs_polyarna.R, tables in equilibrium.zip

For ATP vs no ATP -> atp_vs_noatp.R, tables in atp_vs_no_atp.zip

For analysis of 5'ss, bpt, and intron -> 5ss_bpt_intron.R, used to annotate master tables in smsite_master.zip, equilibrium.zip, atp_vs_no_atp.zip, h_Y12.tsv, h_SmB.tsv, and m_atp_stringent.tsv

For comparisons of Y12 and SmB RIP-seq -> Y12_SmB_enriched.R, tables in h_SmB.tsv and h_Y12.tsv

For processing ATP vs no ATP in higher stringency washing conditions -> stringent_Y12_ATP.R, table in m_atp_stringent.tsv

For identifying candidate mRNAs to test Sm-ring assembly directly -> candidates.R

For processing Deseq2 outputs from SMA models for Sm-site and Sm-RIP characterizations -> sma_models.R

For all other supplemental data (Lu Y12-RIP, Briese SmB/B'-CLIP, Nijssen axon vs soma, Todd alpha-COP) -> supplemental_data.R

     
