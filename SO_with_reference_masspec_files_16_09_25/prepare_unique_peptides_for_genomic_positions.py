
import os
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--unique_peptide_information_csv",
        required=True,
        help="Path to the unique_peptide_information_csv."
    )

    parser.add_argument(
        "--favorite_assembly",
        required=True,
        help="Favorite assembly in which to search for the Split-ORF proteins of interest."
    )

    parser.add_argument(
        "--cell_type",
        required=True,
        help="cell_type of which the unique SO proteins were identified."
    )

    return parser.parse_args()


################################################################################
# PATH DEFINITIONS
################################################################################
unique_peptide_information_csv = '/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_validated_SO_protein_original_Ids_with_assembly.csv'
favorite_assembly = 'TAMA_HUVEC'
cell_type = 'HUVEC'


def main(unique_peptide_information_csv, favorite_assembly, cell_type):
    ################################################################################
    # READ IN DATA
    ################################################################################
    unique_peptide_information_df = pd.read_csv(unique_peptide_information_csv)
    outdir = os.path.dirname(unique_peptide_information_csv)

    ################################################################################
    # Get peptides perferably for fav assembly, the others assign to first assembl
    ################################################################################
    unique_peptide_information_fav_assembly_df = unique_peptide_information_df[
        unique_peptide_information_df['assembly'] == f'{favorite_assembly}']
    split_orf_protein_list = unique_peptide_information_fav_assembly_df['SO_unique_ID'].to_list(
    )

    other_unique_peptides_df = unique_peptide_information_df[~unique_peptide_information_df['SO_unique_ID'].isin(
        split_orf_protein_list)]
    other_unique_peptides_df = other_unique_peptides_df.groupby('SO_unique_ID').agg(
        {'Unnamed: 0': 'first', 'SO_unique_ID': 'first', 'assembly': 'first', 'Prot_start_position': 'first', 'Prot_end_position': 'first'})

    val_protein_no_redundancy_df = pd.concat(
        [unique_peptide_information_fav_assembly_df, other_unique_peptides_df], axis=0)
    val_protein_no_redundancy_df = val_protein_no_redundancy_df.reset_index(
        drop=True)
    val_protein_no_redundancy_df = val_protein_no_redundancy_df.drop(
        'Unnamed: 0.1', axis=1)

    ################################################################################
    # Convert to SO pipe format and write one bed file per assembly
    ################################################################################
    # convert string lists into real lists
    val_protein_no_redundancy_df['Prot_start_position'] = val_protein_no_redundancy_df['Prot_start_position'].apply(
        lambda x: x.strip('[]').split(', '))
    val_protein_no_redundancy_df['Prot_end_position'] = val_protein_no_redundancy_df['Prot_end_position'].apply(
        lambda x: x.strip('[]').split(', '))

    val_protein_no_redundancy_exploded_df = val_protein_no_redundancy_df.explode(
        ['Prot_start_position', 'Prot_end_position'])

    # remove duplicates in case there were several same peptides identified
    val_protein_no_redundancy_exploded_df = val_protein_no_redundancy_exploded_df.drop_duplicates()

    assembly_peptide_df_dict = {
        assembly: subdf for assembly, subdf in val_protein_no_redundancy_exploded_df.groupby('assembly')}

    for assembly, df in assembly_peptide_df_dict.items():
        df['gID|tID'] = df['Unnamed: 0'].apply(lambda x: x.split(':')[0])
        df['OrfID'] = df['Unnamed: 0'].apply(
            lambda x: ':'.join(x.split(':')[1:]))
        df['OrfStart'] = df['OrfID'].apply(lambda x: int(x.split(':')[1]))
        df['Prot_start_position'] = df.apply(lambda row: int(
            row['Prot_start_position']) * 3 + row['OrfStart'], axis=1)
        df['Prot_end_position'] = df.apply(lambda row: int(
            row['Prot_end_position']) * 3 + row['OrfStart'], axis=1)

        df[['gID|tID', 'Prot_start_position', 'Prot_end_position', 'OrfID']].to_csv(os.path.join(
            outdir, f'{assembly}_unique_peptides_of_{cell_type}.bed'), sep='\t', index=False, header=False)

    ################################################################################
    # LOAD MS DATA
    ################################################################################


if __name__ == "__main__":
    args = parse_arguments()

    unique_peptide_information_csv = args.unique_peptide_information_csv
    favorite_assembly = args.favorite_assembly
    cell_type = args.cell_type

    main(unique_peptide_information_csv, favorite_assembly, cell_type)
