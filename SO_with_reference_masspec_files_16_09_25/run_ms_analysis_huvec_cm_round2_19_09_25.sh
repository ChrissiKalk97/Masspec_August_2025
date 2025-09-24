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


# Only need the genomic positions for HUVEC right now as for CM there is no Riboseq data
# Prepare per assembly the unique peptides (no redundancy) to calculate genomic positions 
# with the Split-ORF pipeline
python prepare_unique_peptides_for_genomic_positions.py \
 --unique_peptide_information_csv /Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_validated_SO_protein_original_Ids_with_assembly.csv \
 --favorite_assembly TAMA_HUVEC \
 --cell_type HUVEC



# Get all .sh files into an array
unique_pep_coord_files=(
    "New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/TAMA_HUVEC_unique_peptides_of_HUVEC.bed" 
    "New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/TAMA_CM_unique_peptides_of_HUVEC.bed" 
    "New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/NMD_unique_peptides_of_HUVEC.bed" 
    "New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/RI_unique_peptides_of_HUVEC.bed" 
    )

coord_files=(
    "/Users/christina/Documents/own_data/PacBio_HiFi/compare_mando_stringtie/SO_pipe/run_12.09.2025-17.51.04_HUVEC_tama_merged/HUVEC_merged_tama_gene_id_ExonCoordsOfTranscriptsForSO.txt_transcript_positions.bed"
    "/Users/christina/Documents/own_data/PacBio_HiFi/compare_mando_stringtie/SO_pipe/run_12.09.2025-14.10.14_CM_tama_merged/CM_merged_tama_gene_id_ExonCoordsOfTranscriptsForSO.txt_transcript_positions.bed"
    "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_18.06.2025-09.35.29_NMD_transcripts/ExonCoordsWIthChr110_transcript_positions.bed"
    "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_18.06.2025-09.35.29_NMD_transcripts/ExonCoordsWIthChr110_transcript_positions.bed"
    )


source $(conda info --base)/etc/profile.d/conda.sh
conda activate Splitorf
# Iterate through them
for i in "${!unique_pep_coord_files[@]}"; do

    unique_pep_file=${unique_pep_coord_files[$i]}
    coord_file=${coord_files[$i]}
    outdir=$(dirname $unique_pep_file)
    outname=$(basename $unique_pep_file .bed)
    python /Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py \
     $unique_pep_file \
     $coord_file\
    $outdir/${outname}_genomic_coordinates.bed

done


# python /Users/christina/Documents/own_data/Masspec/compare_abundance_across_samples.py \
#     --cm_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/cm_PD_split_orf_only_peptides_quan_filtered.csv" \
#     --huvec_unique_peps_csv ""/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_PD_split_orf_only_peptides_quan_filtered.csv"" \
#     --scaled "Normalized"


# python /Users/christina/Documents/own_data/Masspec/compare_abundance_across_samples.py \
#     --cm_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/cm_PD_split_orf_only_peptides_quan_filtered.csv" \
#     --huvec_unique_peps_csv "/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_PD_split_orf_only_peptides_quan_filtered.csv" \
#     --scaled "Scaled"