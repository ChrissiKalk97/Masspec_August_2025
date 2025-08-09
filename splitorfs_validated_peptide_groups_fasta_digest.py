import os
import pandas as pd
from Bio import SeqIO
from pyteomics import parser
import pickle
import argparse
import ast


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--filtered_peptide_df",
        required=True,
        help="."
    )

    parser.add_argument(
        "--unfiltered_peptide_df",
        help="."
    )

    parser.add_argument(
        "--not_human_iso_filtered_df",
        required=True,
        help="."
    )

    parser.add_argument(
        "--proteome_fasta_file",
        required=True,
        help="."
    )

    parser.add_argument(
        "--so_pipe_proteome_fasta",
        required=True,
        help="."
    )

    parser.add_argument(
        "--digest",
        required=True,
        help="."
    )

    parser.add_argument(
        "--cell_type",
        required=True,
        help="."
    )
    return parser.parse_args()


################################################################################
# PATH DEFINITIONS
################################################################################
# filtered_peptide_df = '/Users/christina/Documents/own_data/Masspec/HUVEC_July_2025/analysis_results/huvec_PD_filtered_peptides.csv'
# unfiltered_peptide_df = '/Users/christina/Documents/own_data/Masspec/HUVEC_July_2025/analysis_results/huvec_PD_before_quan_filtering_peptides.csv'
# not_human_iso_filtered_df = '/Users/christina/Documents/own_data/Masspec/HUVEC_July_2025/analysis_results/huvec_PD_before_human_iso_filtering_peptides.csv'
# proteome_fasta_file = '/Users/christina/Documents/own_data/Masspec/20250327_AS_LC4_MAA_VLD_12420_01_iPS_TS_PeptideGroups/Homo sapiens (sp_incl_isoforms TaxID=9606) - [Release=408].fasta'
# so_pipe_proteome_fasta = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Input2023/TSL_filtered/protein_coding_peptide_sequences_110_tsl_refseq_filtered.fa'
# cell_type = 'huvec'
# outpath = os.path.dirname(filtered_peptide_df)
# digest = False


def main(filtered_peptide_df, not_human_iso_filtered_df,
         proteome_fasta_file, so_pipe_proteome_fasta, digest, outpath, cell_type, unfiltered_peptide_df=''):

    ################################################################################
    # LOAD FILTERED MS DATA
    ################################################################################
    peptides_df = pd.read_csv(filtered_peptide_df, index_col=0)
    peptides_df = peptides_df.drop('index', axis=1)
    min_peptide_len = peptides_df['Sequence Length'].min()
    max_peptide_length = peptides_df['Sequence Length'].max()

    ################################################################################
    # Unique peptides compared to proteome FASTA with isoforms
    ################################################################################
    # get all possible peptides from Uniprot to check the uniqueness of the peptides in SO proteins observed
    if digest:
        proteome_peptide_list = []
        for record in SeqIO.parse(proteome_fasta_file, 'fasta'):
            proteome_peptide_list = proteome_peptide_list + \
                list(parser.cleave(str(record.seq),
                                   parser.expasy_rules['trypsin'], min_length=min_peptide_len, missed_cleavages=2, max_length=max_peptide_length))

        proteome_peptide_set = set(proteome_peptide_list)

        with open(os.path.join(outpath, 'Homo_sp_incl_isoforms_9606_realease_408.pkl'), 'wb') as f:
            pickle.dump(proteome_peptide_set, f)

    else:
        with open(os.path.join(outpath, 'Homo_sp_incl_isoforms_9606_realease_408.pkl'), 'rb') as f:
            proteome_peptide_set = pickle.load(f)

    peptides_df['Not unique'] = peptides_df['Sequence'].apply(
        lambda x: x in proteome_peptide_set)

    ################################################################################
    # Unique peptides compared to SO protein FASTA
    ################################################################################
    # get all possible peptides from Uniprot to check the uniqueness of the peptides in SO proteins observed
    if digest:
        so_proteome_peptide_list = []
        for record in SeqIO.parse(so_pipe_proteome_fasta, 'fasta'):
            so_proteome_peptide_list = so_proteome_peptide_list + \
                list(parser.cleave(str(record.seq),
                                   parser.expasy_rules['trypsin'], min_length=min_peptide_len, missed_cleavages=2, max_length=max_peptide_length))

        so_proteome_peptide_set = set(so_proteome_peptide_list)

        with open(os.path.join(outpath, 'SO_proteome_Ens_110_TSL_and_Refseq_filtered.pkl'), 'wb') as f:
            pickle.dump(so_proteome_peptide_set, f)

    else:
        with open(os.path.join(outpath, 'SO_proteome_Ens_110_TSL_and_Refseq_filtered.pkl'), 'rb') as f:
            so_proteome_peptide_set = pickle.load(f)

    peptides_df['Not unique SO input'] = peptides_df['Sequence'].apply(
        lambda x: x in so_proteome_peptide_set)

    ################################################################################
    # Get the SO proteins that have a peptide found
    ################################################################################
    peptides_df_with_unique_evidence = peptides_df[peptides_df['Not unique'] == False].reset_index(
        drop=True).copy()

    peptides_df_with_unique_evidence['Protein Accessions list'] = peptides_df_with_unique_evidence.loc[:, 'Protein Accessions'].apply(
        lambda x: x.split('; ')).copy()

    peptides_with_unique_evidence_exploded_df = peptides_df_with_unique_evidence.explode(
        'Protein Accessions list')

    # these peptides have some kind of unique evidence, can be shared with other SO proteins
    peptides_with_unique_evidence_list = peptides_with_unique_evidence_exploded_df['Protein Accessions list'].to_list(
    )

    # read in df with all but contaminant peptides
    peptides_incl_human_iso_df = pd.read_csv(
        not_human_iso_filtered_df, index_col=0)

    # filter for proteins that have some kind of unique evidence (also shared amongst SOs)
    peptides_incl_human_iso_df['Protein Accessions list'] = peptides_incl_human_iso_df.loc[:, 'Protein Accessions'].apply(
        lambda x: x.split('; ')).copy()

    peptides_incl_human_iso_df['Proteins with unique evidence'] = peptides_incl_human_iso_df['Protein Accessions list'].apply(
        lambda x: [protein for protein in x if protein in peptides_with_unique_evidence_list])

    peptides_of_proteins_with_unique_df = peptides_incl_human_iso_df[peptides_incl_human_iso_df['Proteins with unique evidence'].apply(
        lambda x: len(x) > 0)]

    peptides_of_proteins_with_unique_exploded_df = peptides_of_proteins_with_unique_df.explode(
        'Protein Accessions list')
    peptides_of_proteins_with_unique_exploded_df = peptides_of_proteins_with_unique_exploded_df[
        peptides_of_proteins_with_unique_exploded_df['Protein Accessions list'].isin(peptides_with_unique_evidence_list)]

    proteins_with_unique_grouped = peptides_of_proteins_with_unique_exploded_df.groupby(
        'Protein Accessions list').agg({'Protein Accessions': list,
                                        'Sequence': list,
                                        'Number of PSMs': list,
                                        'Positions in Master Proteins': list,
                                        'Sequence Length': list,
                                        'Proteins with unique evidence': list,
                                        'PEP': 'nunique',
                                        'q-Value': 'nunique'
                                        })

    proteins_with_unique_grouped['nr all peptides'] = proteins_with_unique_grouped['Sequence'].apply(
        lambda x: len(x))

    # no proteins filtered out, all have at least 2 peptides
    proteins_with_unique_grouped = proteins_with_unique_grouped[
        proteins_with_unique_grouped['nr all peptides'] > 1]

    proteins_with_unique_grouped['Unique Peptide Sequences'] = [
        [] for protein in proteins_with_unique_grouped.index]
    proteins_with_unique_grouped['Nr Unique Peptides'] = 0
    for protein in proteins_with_unique_grouped.index:
        unique_sequences = peptides_with_unique_evidence_exploded_df[
            peptides_with_unique_evidence_exploded_df['Protein Accessions list'] == protein]['Sequence'].to_list()
        proteins_with_unique_grouped.loc[protein, 'Unique Peptide Sequences'].append(
            unique_sequences)

        proteins_with_unique_grouped.loc[protein,
                                         'Nr Unique Peptides'] = len(unique_sequences)

        print('Number of SO proteins with 2 or more unique peptides:',
              sum(proteins_with_unique_grouped['Nr Unique Peptides'] > 1))

        proteins_with_unique_grouped.to_csv(os.path.join(
            outpath, f'{cell_type}_proteome_fasta_uniqueness_df.csv'))


if __name__ == "__main__":
    args = parse_arguments()

    filtered_peptide_df = args.filtered_peptide_df
    unfiltered_peptide_df = args.unfiltered_peptide_df
    not_human_iso_filtered_df = args.not_human_iso_filtered_df
    proteome_fasta_file = args.proteome_fasta_file
    so_pipe_proteome_fasta = args.so_pipe_proteome_fasta
    digest = ast.literal_eval(args.digest)
    cell_type = args.cell_type
    outpath = os.path.dirname(filtered_peptide_df)

    main(filtered_peptide_df, not_human_iso_filtered_df,
         proteome_fasta_file, so_pipe_proteome_fasta, digest, outpath, cell_type, unfiltered_peptide_df)
