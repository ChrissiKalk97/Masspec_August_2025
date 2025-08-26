
import os
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--peptides_file",
        required=True,
        help="Path to the peptides file."
    )

    parser.add_argument(
        "--so_id_mapping_file",
        required=True,
        help="Path to the SO ID mapping file."
    )

    parser.add_argument(
        "--cell_type",
        required=True,
        help="Cell type name."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Outdirectory for results."
    )

    return parser.parse_args()


################################################################################
# PATH DEFINITIONS
################################################################################
# peptides_file = '/Users/christina/Documents/own_data/Masspec/HUVEC_July_2025/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt'
# so_id_mapping_file = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_unique_SO_ID_with_assembly_info_so_id_mapping_with_assembly_info.tsv'
# cell_type = 'huvec'
# outdir = os.path.dirname(peptides_file)
# outdir = os.path.join(outdir, 'analysis_results')


################################################################################
# LOAD MS DATA
################################################################################
def main(peptides_file, so_id_mapping_file, cell_type, outdir):
    # read in peptides df
    peptides_df = pd.read_csv(peptides_file, sep='\t')

    peptides_df = peptides_df[['Protein Accessions',
                               'Sequence',
                               'Contaminant',
                               'Marked as',
                               'Number of Protein Groups',
                               'Number of Proteins',
                               'Number of PSMs',
                               'Master Protein Accessions',
                               'Positions in Master Proteins',
                               'Sequence Length',
                               'Quan Info',
                               'PEP',
                               'q-Value',
                               'Confidence',
                               'PSM Ambiguity'] + [col for col in peptides_df.columns if 'Abundance' in col]].copy()

    print('total number of peptides identified (unfiltered):',
          len(peptides_df.index))

    ################################################################################
    # FILTER OUT NON-UNIQUE AND NON-ABUNDANCE PEPTIDES
    ################################################################################
    peptides_df = peptides_df[peptides_df['Contaminant'] == False].copy()
    print('total number of peptides identified (no contaminants):',
          len(peptides_df.index))

    # save the unfiltered df
    peptides_df.to_csv(
        os.path.join(outdir, f'{cell_type}_PD_before_human_iso_filtering_peptides.csv'))

    peptides_df = peptides_df[peptides_df['Marked as'] != 'human_iso'].copy()
    print('total number of peptides identified (no contaminants, not human_iso):',
          len(peptides_df.index))

    peptides_df_with_quan_info = peptides_df[~peptides_df['Quan Info'].isin(
        ['NotUnique', 'NoQuanValues'])].copy()
    print('total number of peptides identified (no contaminants, not human_iso, with quan info):',
          len(peptides_df_with_quan_info.index))

    print('number of peptides uniquely attributable to one Split-ORF protein',
          sum(peptides_df_with_quan_info['Number of Proteins'] == 1))

    ################################################################################
    # CHECK WHICH ASSEMBLY THE PEPTIDES STEM FROM
    ################################################################################
    # read in the df with information of original SO name for each SO ID
    so_id_mapping_df = pd.read_csv(so_id_mapping_file, sep='\t')

    # generate list of the SO IDs (Accessions)
    peptides_df_with_quan_info['Protein Accessions list'] = peptides_df_with_quan_info.loc[:, 'Protein Accessions'].apply(
        lambda x: x.split('; ')).copy()
    # remove the sp| from the SO IDs
    peptides_df_with_quan_info['Protein Accessions list'] = peptides_df_with_quan_info.loc[:, 'Protein Accessions list'].apply(
        lambda x: [so_id.split('|')[1] for so_id in x]).copy()
    # assign assembly information: for all SO proteins the respective assembly is noted down
    # if the same protein is found in different assemblies: all are noted down
    # if different proteins are in the list from the same assembly: assembly is noted down a number of times
    peptides_df_with_quan_info['assembly of protein'] = peptides_df_with_quan_info.loc[:, 'Protein Accessions list'].apply(
        lambda x: so_id_mapping_df.loc[so_id_mapping_df['SO_unique_ID'].isin(x), 'assembly'].tolist()).copy()
    # how often is the cell type of the MS samples also the one in which the protein was predicted
    peptides_df_with_quan_info[f'assembly contains {cell_type}'] = peptides_df_with_quan_info.loc[:, 'assembly of protein'].apply(
        lambda x: cell_type in x).copy()

    ################################################################################
    # FILTER FOR VALIDATED PROTEINS
    ################################################################################
    peptides_df_with_quan_info['shared SO peptide'] = peptides_df_with_quan_info.loc[:,
                                                                                     'Protein Accessions list'].apply(lambda x: len(x) > 1).copy()
    peptides_df_with_quan_info_exploded = peptides_df_with_quan_info.explode(
        'Protein Accessions list').copy()

    # Regard shared SO peptides as evidence
    peptides_grouped_by_proteins_df = peptides_df_with_quan_info_exploded.groupby('Protein Accessions list').agg({'Sequence': 'nunique',
                                                                                                                  f'assembly contains {cell_type}': 'first'}).copy()

    # if at least two peptides are found a protein is kept, other wise not considered as validated
    peptides_grouped_by_proteins_filtered_df = peptides_grouped_by_proteins_df.loc[
        peptides_grouped_by_proteins_df['Sequence'] > 1, :].copy()

    so_with_evidence = peptides_grouped_by_proteins_filtered_df.index.to_list()

    peptides_df_with_quan_info['Filtered SO IDs'] = peptides_df_with_quan_info.loc[:,
                                                                                   'Protein Accessions list'].apply(
        lambda x: [
            so_id for so_id in x if so_id in so_with_evidence]).copy()
    peptides_df_with_quan_info['Keep Peptide'] = peptides_df_with_quan_info.loc[:,
                                                                                'Filtered SO IDs'].apply(lambda x: len(x) > 0).copy()
    peptides_df_with_quan_info_filtered = peptides_df_with_quan_info[
        peptides_df_with_quan_info['Keep Peptide'] == True].copy(
    )

    peptides_df_with_quan_info_filtered = peptides_df_with_quan_info_filtered.reset_index()

    # save the filtered df
    peptides_df_with_quan_info_filtered.to_csv(
        os.path.join(outdir, f'{cell_type}_PD_filtered_peptides.csv'))

    # save the unfiltered df
    peptides_df.to_csv(
        os.path.join(outdir, f'{cell_type}_PD_before_quan_filtering_peptides.csv'))


if __name__ == "__main__":
    args = parse_arguments()

    peptides_file = args.peptides_file
    so_id_mapping_file = args.so_id_mapping_file
    cell_type = args.cell_type
    outdir = args.outdir

    main(peptides_file, so_id_mapping_file, cell_type, outdir)
