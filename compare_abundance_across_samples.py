import os
import pandas as pd
import argparse
from scipy import stats
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--cm_unique_peps_csv",
        required=True,
        help="."
    )

    parser.add_argument(
        "--huvec_unique_peps_csv",
        help="."
    )
    return parser.parse_args()


def main(cm_unique_peps_csv, cm_dict, huvec_unique_peps_csv, huvec_dict):
    cm_unique_peptdies_df = pd.read_csv(cm_unique_peps_csv, index_col=0)
    huvec_unique_peptdies_df = pd.read_csv(huvec_unique_peps_csv, index_col=0)

    ################################################################################
    # Compare CM NMD vs DMSO by wilcoxon signed rank test
    ################################################################################
    print('Test differences in expression of unique peptides for CM NMD vs DMSO')
    wilcoxon_signed_test_ps = []
    for unique_peptide in cm_unique_peptdies_df.index:

        wilcox_signed_result = stats.wilcoxon(
            x=np.array(
                cm_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 130C Sample',
                                                           f'Abundances Normalized F1 131N Sample',
                                                           f'Abundances Normalized F1 131C Sample']].copy(),
                dtype=float),
            y=np.array(
                cm_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 {cm_dict['130C']} Sample',
                                                           f'Abundances Normalized F1 {cm_dict['131N']} Sample',
                                                           f'Abundances Normalized F1 {cm_dict['131C']} Sample']].copy(),
                dtype=float),
            alternative='greater')
        # print(wilcox_signed_result)
        wilcoxon_signed_test_ps.append(wilcox_signed_result.pvalue)

    # the lowest possible pvalue is 0.125
    avg_p = sum(wilcoxon_signed_test_ps)/len(cm_unique_peptdies_df.index)
    print('The average p-value amongst all tests is:', avg_p)

    print('The number of tests with a p-value of 0.25 or smaller:',
          len([p for p in wilcoxon_signed_test_ps if p <= 0.25]))

    print('The number of tests with a p-value smaller than 0.5:',
          len([p for p in wilcoxon_signed_test_ps if p < 0.5]))

    ################################################################################
    # Compare CM Hypo vs DMSO by wilcoxon signed rank test
    ################################################################################
    print('Test differences in expression of unique peptides for CM Hypo vs DMSO')
    wilcoxon_signed_test_ps = []
    for unique_peptide in cm_unique_peptdies_df.index:

        wilcox_signed_result = stats.wilcoxon(
            x=np.array(
                cm_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 127N Sample',
                                                           f'Abundances Normalized F1 128N Sample',
                                                           f'Abundances Normalized F1 129N Sample']].copy(),
                dtype=float),
            y=np.array(
                cm_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 {cm_dict['127N']} Sample',
                                                           f'Abundances Normalized F1 {cm_dict['128N']} Sample',
                                                           f'Abundances Normalized F1 {cm_dict['129N']} Sample']].copy(),
                dtype=float),
            alternative='greater')
        # print(wilcox_signed_result)
        wilcoxon_signed_test_ps.append(wilcox_signed_result.pvalue)

    # the lowest possible pvalue is 0.125
    avg_p = sum(wilcoxon_signed_test_ps)/len(cm_unique_peptdies_df.index)
    print('The average p-value amongst all tests is:', avg_p)

    print('The number of tests with a p-value of 0.25 or smaller:',
          len([p for p in wilcoxon_signed_test_ps if p <= 0.25]))

    print('The number of tests with a p-value smaller than 0.5:',
          len([p for p in wilcoxon_signed_test_ps if p < 0.5]))

    ################################################################################
    # Compare HUVEC hypo vs DMSO by wilcoxon signed rank test
    ################################################################################
    print('Test differences in expression of unique peptides for HUVEC Hypo vs DMSO')
    wilcoxon_signed_test_ps = []
    for unique_peptide in huvec_unique_peptdies_df.index:

        wilcox_signed_result = stats.wilcoxon(
            x=np.array(
                huvec_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 127 Sample',
                                                              f'Abundances Normalized F1 129 Sample',
                                                              f'Abundances Normalized F1 131 Sample']].copy(),
                dtype=float),
            y=np.array(
                huvec_unique_peptdies_df.loc[unique_peptide, [f'Abundances Normalized F1 {huvec_dict['127']} Sample',
                                                              f'Abundances Normalized F1 {huvec_dict['129']} Sample',
                                                              f'Abundances Normalized F1 {huvec_dict['131']} Sample']].copy(),
                dtype=float),
            alternative='greater')
        # print(wilcox_signed_result)
        wilcoxon_signed_test_ps.append(wilcox_signed_result.pvalue)

    # the lowest possible pvalue is 0.125
    avg_p = sum(wilcoxon_signed_test_ps)/len(huvec_unique_peptdies_df.index)
    print('The average p-value amongst all tests is:', avg_p)

    print('The number of tests with a p-value of 0.25 or smaller:',
          len([p for p in wilcoxon_signed_test_ps if p <= 0.25]))

    print('The number of tests with a p-value smaller than 0.5:',
          len([p for p in wilcoxon_signed_test_ps if p < 0.5]))


if __name__ == "__main__":
    args = parse_arguments()
    cm_unique_peps_csv = args.cm_unique_peps_csv
    huvec_unique_peps_csv = args.huvec_unique_peps_csv

    # cm_unique_peps_csv = '/Users/christina/Documents/own_data/Masspec/analysis_results/cm_unique_so_peptides_after_proteome_filtering.csv'
    # huvec_unique_peps_csv = '/Users/christina/Documents/own_data/Masspec/analysis_results/huvec_unique_so_peptides_after_proteome_filtering.csv'

    cm_dict = {
        # hypoxia vs dmso by clone
        '127N': '126', '128N': '127C', '129N': '128C',
        # nmd vs dmso by clone
        '130C': '126', '131N': '127C', '131C': '128C',
    }

    huvec_dict = {
        # hypoxia vs dmso by clone
        '127': '126', '129': '128', '131': '130',
    }

    main(cm_unique_peps_csv, cm_dict, huvec_unique_peps_csv, huvec_dict)
