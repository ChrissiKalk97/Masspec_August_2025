#!/bin/bash

# ----- This script analyses the HUVEC and CM Masspec Data from June 2025 ----- #
# ----- The PD results are first analyzed assumint that they are correct  ----- #
# ----- but actually the marking did not work, the second script performs ----- #
# ----- the tryptic digest of 2 possible references and really finds the  ----- #
# ----- unique peptides  ----- #


outdir="/Users/christina/Documents/own_data/Masspec/analysis_results"

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi


python splitorfs_validated_peptide_groups_PD_based_analysis.py \
 --peptides_file "/Users/christina/Documents/own_data/Masspec/HUVEC_July_2025/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt" \
 --so_id_mapping_file "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_unique_SO_ID_with_assembly_info_so_id_mapping_with_assembly_info.tsv" \
 --cell_type "huvec" \
 --outdir "$outdir"


python splitorfs_validated_peptide_groups_PD_based_analysis.py \
 --peptides_file "/Users/christina/Documents/own_data/Masspec/IPSC_CM_July_2025/48F_IPSC/20250505_AS_LC4_MAA_20050_01_VLD_iPSC_F_PeptideGroups.txt" \
 --so_id_mapping_file "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_unique_SO_ID_with_assembly_info_so_id_mapping_with_assembly_info.tsv" \
 --cell_type "cm" \
 --outdir "$outdir"


python splitorfs_validated_peptide_groups_fasta_digest.py \
 --filtered_peptide_df "/Users/christina/Documents/own_data/Masspec/analysis_results/huvec_PD_filtered_peptides.csv" \
 --unfiltered_peptide_df "/Users/christina/Documents/own_data/Masspec/analysis_results/huvec_PD_before_quan_filtering_peptides.csv" \
 --not_human_iso_filtered_df "/Users/christina/Documents/own_data/Masspec/analysis_results/huvec_PD_before_human_iso_filtering_peptides.csv" \
 --proteome_fasta_file "/Users/christina/Documents/own_data/Masspec/20250327_AS_LC4_MAA_VLD_12420_01_iPS_TS_PeptideGroups/Homo sapiens (sp_incl_isoforms TaxID=9606) - [Release=408].fasta" \
 --so_pipe_proteome_fasta "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Input2023/TSL_filtered/protein_coding_peptide_sequences_110_tsl_refseq_filtered.fa" \
 --cell_type "huvec" \
 --digest False


python splitorfs_validated_peptide_groups_fasta_digest.py \
 --filtered_peptide_df "/Users/christina/Documents/own_data/Masspec/analysis_results/cm_PD_filtered_peptides.csv" \
 --not_human_iso_filtered_df "/Users/christina/Documents/own_data/Masspec/analysis_results/cm_PD_before_human_iso_filtering_peptides.csv" \
 --proteome_fasta_file "/Users/christina/Documents/own_data/Masspec/20250327_AS_LC4_MAA_VLD_12420_01_iPS_TS_PeptideGroups/Homo sapiens (sp_incl_isoforms TaxID=9606) - [Release=408].fasta" \
 --so_pipe_proteome_fasta "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Input2023/TSL_filtered/protein_coding_peptide_sequences_110_tsl_refseq_filtered.fa" \
 --cell_type "cm" \
 --digest False

python compare_abundance_across_samples.py \
    --cm_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/analysis_results/cm_unique_so_peptides_after_proteome_filtering.csv" \
    --huvec_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/analysis_results/huvec_unique_so_peptides_after_proteome_filtering.csv"
