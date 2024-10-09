# sm_ring_assembly_mrna

These are scripts used to identify and characterize Sm-site containing RNAs within the human and mouse transcriptomes.

a) to identify Sm-sites in human and mouse transcriptomes, use find_smsites.py

Custom R script used to combine the outputs from find_smsites.py into a single file can be accessed in combine_smsite_output.R

Custom R scripts used to make a single master file giving the frequency of Sm-sites per gene_id -> human_smsite_utype.R, mouse_smsite_utype.R

Custom R scripts used to describe the RNAs identfiied with Sm-sites -> insilico.R

b) to process raw sequencing files from publication, use sm_rip_seq_processing.sh

To visualize processing of reads -> multiqc_report_human.html, multiqc_report_mouse.html

To combine featureCounts outputs into a single table -> make_master.py

To obtain Wald test comparisons of featureCounts outputs -> deseq2_analysis.R

c) characterize Deseq2 generated Wald test comparisons:

For Sm-RIP vs polyA-RNA -> sm_rip_vs_polyarna.R

For ATP vs no ATP -> atp_vs_noatp.R

For identifying candidate mRNAs to test Sm-ring assembly directly -> candidates.R

For processing Deseq2 outputs from SMA models for Sm-site and Sm-RIP characterizations -> sma_models.R

For all other supplemental data (Lu Y12-RIP, Briese SmB/B'-CLIP, Nijssen axon vs soma, Todd alpha-COP) -> supplemental_data.R

     
