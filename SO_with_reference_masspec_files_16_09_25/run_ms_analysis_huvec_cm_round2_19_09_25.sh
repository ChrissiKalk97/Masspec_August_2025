#!/bin/bash

# ----- This script analyses the HUVEC and CM Masspec Data from June 2025 ----- #
# ----- The PD results are first analyzed assumint that they are correct  ----- #
# ----- but actually the marking did not work, the second script performs ----- #
# ----- the tryptic digest of 2 possible references and really finds the  ----- #
# ----- unique peptides  ----- #


outdir="/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25"

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi


python peptide_groups_pd_analysis_with_reference_19_09_25.py \
 --peptides_file "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt" \
 --so_id_mapping_file "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/tama_NMD_RI_masspec_files/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25_so_id_mapping_with_assembly_info.tsv" \
 --cell_type "huvec" \
 --outdir "$outdir"


python peptide_groups_pd_analysis_with_reference_19_09_25.py \
 --peptides_file "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/20250505_AS_LC4_MAA_20050_01_VLD_iPSC_F_PeptideGroups.txt" \
 --so_id_mapping_file "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/tama_NMD_RI_masspec_files/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25_so_id_mapping_with_assembly_info.tsv" \
 --cell_type "cm" \
 --outdir "$outdir"



python /Users/christina/Documents/own_data/Masspec/compare_abundance_across_samples.py \
    --cm_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/cm_PD_split_orf_only_peptides_quan_filtered.csv" \
    --huvec_unique_peps_csv ""/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_PD_split_orf_only_peptides_quan_filtered.csv"" \
    --scaled "Normalized"


python /Users/christina/Documents/own_data/Masspec/compare_abundance_across_samples.py \
    --cm_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/cm_PD_split_orf_only_peptides_quan_filtered.csv" \
    --huvec_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_PD_split_orf_only_peptides_quan_filtered.csv" \
    --scaled "Scaled"