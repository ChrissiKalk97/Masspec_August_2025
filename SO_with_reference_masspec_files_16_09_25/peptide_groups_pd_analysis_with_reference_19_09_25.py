
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
peptides_file = '/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt'
so_id_mapping_file = '/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/tama_NMD_RI_masspec_files/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25_so_id_mapping_with_assembly_info.tsv'
cell_type = 'huvec'
outdir = os.path.dirname(peptides_file)
outdir = os.path.join(outdir, 'analysis_results_with_ref_19_09_25')


def main(peptides_file, so_id_mapping_file, cell_type, outdir):
    ################################################################################
    # HELPER FUNCTION DEFINITIONS
    ################################################################################
    def nr_proteins_with_at_least_two_unique_peptides(split_orf_peptides_with_quan_info_df):
        exploded_df = split_orf_peptides_with_quan_info_df.explode(
            'Protein Accessions List')
        grouped_df = exploded_df.groupby(
            'Protein Accessions List').agg({'Sequence': 'nunique'}).copy()
        grouped_df = grouped_df[grouped_df['Sequence'] > 1]
        print('Number of proteins with two or more unique peptides:',
              len(grouped_df.index))

    def filter_peptides_contamination_and_reference(peptides_df):
        peptides_df = peptides_df[peptides_df['Contaminant'] == False].copy()
        print('total number of peptides identified (no contaminants):',
              len(peptides_df.index))

        peptides_df['Protein Accessions List'] = peptides_df['Protein Accessions'].apply(
            lambda x: x.split(';'))

        peptides_df['Nr Reference Proteins'] = peptides_df['Protein Accessions List'].apply(
            lambda x: len([protein for protein in x if 'ReferenceProtein' in protein]))

        split_orf_peptides_df = peptides_df[peptides_df['Nr Reference Proteins'] == 0]

        print('total number of split-orf peptides identified (no contaminants):',
              len(split_orf_peptides_df.index))

        # save the unfiltered df
        split_orf_peptides_df.to_csv(
            os.path.join(outdir, f'{cell_type}_PD_split_orf_only_peptides_unfiltered.csv'))

        split_orf_peptides_with_quan_info_df = split_orf_peptides_df[~split_orf_peptides_df['Quan Info'].isin(
            ['NotUnique', 'NoQuanValues', 'ExcludedByMethod'])].copy()
        print('total number of split-orf peptides identified (no contaminants, with quan info):',
              len(split_orf_peptides_with_quan_info_df.index))

        split_orf_peptides_with_quan_info_df['Positions in Proteins'] = split_orf_peptides_with_quan_info_df['Positions in Proteins'].apply(
            lambda x: x.split("; ")).apply(lambda x: [x[i] if x[i].startswith('Split') else x[i-1].split('[')[0]+x[i] for i in range(len(x))])

        split_orf_peptides_with_quan_info_df['Positions'] = split_orf_peptides_with_quan_info_df['Positions in Proteins'].apply(
            lambda x: [prot_info.split('[')[1].strip(']') for prot_info in x])

        split_orf_peptides_with_quan_info_df['Protein Accessions List'] = split_orf_peptides_with_quan_info_df['Positions in Proteins'].apply(
            lambda x: [prot_info.split('[')[0].strip() for prot_info in x])

        split_orf_peptides_with_quan_info_df.to_csv(
            os.path.join(outdir, f'{cell_type}_PD_split_orf_only_peptides_quan_filtered.csv'))

        return peptides_df, split_orf_peptides_df, split_orf_peptides_with_quan_info_df

    def get_validated_splitorf_proteins(split_orf_peptides_with_quan_info_df):
        so_proteins_with_unique_peptides_list = [
            protein for prot_list in split_orf_peptides_with_quan_info_df['Protein Accessions List'].to_list() for protein in prot_list]
        so_proteins_with_unique_peptides_set = set(
            so_proteins_with_unique_peptides_list)

        peptides_filtered_df = peptides_df[peptides_df['Protein Accessions List'].apply(lambda x: len(
            [protein for protein in x if protein in so_proteins_with_unique_peptides_set]) > 0)].copy()

        peptides_filtered_exploded_df = peptides_filtered_df.explode(
            'Protein Accessions List').copy()

        peptides_filtered_exploded_df = peptides_filtered_exploded_df[peptides_filtered_exploded_df['Protein Accessions List'].apply(
            lambda x: x in so_proteins_with_unique_peptides_set)]

        peptides_grouped_by_proteins_df = peptides_filtered_exploded_df.groupby(
            'Protein Accessions List').agg({'Sequence': 'nunique'}).copy()

        proteins_validated_df = peptides_grouped_by_proteins_df[
            peptides_grouped_by_proteins_df['Sequence'] > 1]

        proteins_validated_list = proteins_validated_df.index.to_list()

        proteins_validated_list = [prot.strip()
                                   for prot in proteins_validated_list]

        print('Number of SO proteins with at least two peptides of which one is unique', len(
            proteins_validated_list))

        with open(os.path.join(outdir, f'{cell_type}_validated_SO_protein_Ids.csv'), 'w') as fp:
            for protein in proteins_validated_list:
                # write each item on a new line
                fp.write("%s\n" % protein)

        return split_orf_peptides_with_quan_info_df, proteins_validated_df, proteins_validated_list

    ################################################################################
    # LOAD MS DATA
    ################################################################################
    # read in peptides df
    peptides_df = pd.read_csv(peptides_file, sep='\t')

    peptides_df = peptides_df[['Protein Accessions',
                               'Sequence',
                               'Contaminant',
                               'Number of Protein Groups',
                               'Number of Proteins',
                               'Number of PSMs',
                               'Master Protein Accessions',
                               'Positions in Proteins',
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

    peptides_df, split_orf_peptides_df, split_orf_peptides_with_quan_info_df = filter_peptides_contamination_and_reference(
        peptides_df)

    print('number of peptides uniquely attributable to one Split-ORF protein',
          sum(split_orf_peptides_with_quan_info_df['Number of Proteins'] == 1))

    nr_proteins_with_at_least_two_unique_peptides(
        split_orf_peptides_with_quan_info_df)

    ################################################################################
    # CHECK WHICH ASSEMBLY THE PEPTIDES STEM FROM
    ################################################################################
    # read in the df with information of original SO name for each SO ID
    so_id_mapping_df = pd.read_csv(so_id_mapping_file, sep='\t')

    # assign assembly information: for all SO proteins the respective assembly is noted down
    # if the same protein is found in different assemblies: all are noted down
    # if different proteins are in the list from the same assembly: assembly is noted down a number of times
    split_orf_peptides_with_quan_info_df['assembly of protein'] = split_orf_peptides_with_quan_info_df.loc[:, 'Protein Accessions List'].apply(
        lambda x: so_id_mapping_df.loc[so_id_mapping_df['SO_unique_ID'].isin(x), 'assembly'].tolist()).copy()
    # how often is the cell type of the MS samples also the one in which the protein was predicted
    split_orf_peptides_with_quan_info_df[f'assembly contains {cell_type}'] = split_orf_peptides_with_quan_info_df.loc[:, 'assembly of protein'].apply(
        lambda x: f'TAMA_{cell_type.upper()}' in x).copy()

    print(f'Number of peptides belonging to {cell_type} assembly', sum(
        split_orf_peptides_with_quan_info_df[f'assembly contains {cell_type}']))

    ################################################################################
    # FILTER FOR VALIDATED PROTEINS
    ################################################################################

    # need to filter peptides_df for the proteins that belong to the unique SO peptides
    # then can see valid proteins at least 2 peptides (1 unique is for sure)
    # can also check how many proteins with 2 or more unique peptides

    split_orf_peptides_with_quan_info_df, proteins_validated_df, proteins_validated_list = get_validated_splitorf_proteins(
        split_orf_peptides_with_quan_info_df)

    so_id_mapping_val_splitorfs_df = so_id_mapping_df[so_id_mapping_df['SO_unique_ID'].isin(
        proteins_validated_list)]

    # get start and end positions of the unique peptides in the respective protein
    so_valid_peptides_exploded_df = split_orf_peptides_with_quan_info_df.explode(
        ['Protein Accessions List', 'Positions'])
    # substract one from positions for 0-based coordinates
    so_valid_peptides_exploded_df['Prot_start_position'] = so_valid_peptides_exploded_df['Positions'].apply(
        lambda x: int(x.split('-')[0]) - 1)
    so_valid_peptides_exploded_df['Prot_end_position'] = so_valid_peptides_exploded_df['Positions'].apply(
        lambda x: int(x.split('-')[1]) - 1)

    # map start and end positions to the respective dataframe with the original ID
    so_id_mapping_val_splitorfs_df
    prot_start_dict = so_valid_peptides_exploded_df.groupby(
        'Protein Accessions List')['Prot_start_position'].apply(list).to_dict()
    prot_end_dict = so_valid_peptides_exploded_df.groupby(
        'Protein Accessions List')['Prot_end_position'].apply(list).to_dict()

    so_id_mapping_val_splitorfs_df['Prot_start_position'] = so_id_mapping_val_splitorfs_df['SO_unique_ID'].map(
        prot_start_dict)
    so_id_mapping_val_splitorfs_df['Prot_end_position'] = so_id_mapping_val_splitorfs_df['SO_unique_ID'].map(
        prot_end_dict)

    so_id_mapping_val_splitorfs_df.to_csv(
        os.path.join(outdir, f'{cell_type}_validated_SO_protein_original_Ids_with_assembly.csv'))


if __name__ == "__main__":
    args = parse_arguments()

    peptides_file = args.peptides_file
    so_id_mapping_file = args.so_id_mapping_file
    cell_type = args.cell_type
    outdir = args.outdir

    main(peptides_file, so_id_mapping_file, cell_type, outdir)
