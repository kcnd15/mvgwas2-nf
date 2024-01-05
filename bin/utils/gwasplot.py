# gwasplot.py
# create GWAS plots

import argparse
import pandas as pd
from pandas import DataFrame

import os
import glob

import numpy as np
import matplotlib.pyplot as plt


def manhattan_plot_sample():

    # import libraries
    # from pandas import DataFrame
    from scipy.stats import uniform
    from scipy.stats import randint

    # sample data
    df = DataFrame({'gene' : ['gene-%i' % i for i in np.arange(10000)],
                    'pvalue' : uniform.rvs(size=10000),
                    'chromosome' : ['ch-%i' % i for i in randint.rvs(0,12,size=10000)]})

    # -log_10(pvalue)
    df['minuslog10pvalue'] = -np.log10(df.pvalue)
    df.chromosome = df.chromosome.astype('category')
    df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(12)], ordered=True)
    df = df.sort_values('chromosome')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('chromosome'))

    # manhattan plot
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, 3.5])

    # x axis label
    ax.set_xlabel('Chromosome')

    # show the graph
    plt.show()


def manhattan_plot(gwas_data: DataFrame, title: str = None):

    df = gwas_data.copy()

    # rename columns
    df = df.set_axis(["chromosome", "position", "pvalue"], axis=1)

    # -log_10(pvalue)
    df['minuslog10pvalue'] = -np.log10(df.pvalue)
    print(df.describe())

    # create category for chromosome field
    df.chromosome = df.chromosome.map(lambda x: "chr" + str(x))
    chromosomes = df.chromosome.unique().tolist()
    df.chromosome = df.chromosome.astype('category')
    # range must represent actual number of chromosomes; here we have only chromosome 1
    df.chromosome = df.chromosome.cat.set_categories(chromosomes, ordered=True)

    df = df.sort_values(['chromosome', 'position'])

    # create an index for the x-axis and group by chromosome
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('chromosome'))

    # manhattan plot
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    ax = fig.add_subplot(111)

    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    ax.set_xlim([0, len(df)])

    max_y = df["minuslog10pvalue"].max()
    ax.set_ylim([0, max_y])

    # x axis label
    ax.set_xlabel('Chromosome')

    # set title
    if title is None:
        title_string = "GWAS plot"
    else:
        title_string = f"GWAS plot for {title}"

    plt.title(title_string)

    # show the graph
    plt.show()


def process_result_file(glob_files, col_chrom: int, col_pos: int, col_p: int) -> DataFrame:

    df = DataFrame()

    show_once = True

    for result_file in glob_files:
        print(f"processing result file: {result_file}")

        try:
            result_data = pd.read_csv(result_file, sep=args.sep)
        except Exception as e:
            print(e)
            exit(1)

        # process result file
        print(f"rows: {len(result_data)}")

        relevant_columns = result_data.columns[[col_chrom, col_pos, col_p]]

        if show_once:
            print(f"columns to be processed: {relevant_columns.to_list()}")
            show_once = False

        result_data_relevant = result_data[relevant_columns].copy()

        # rename columns to standard names
        result_data_relevant = result_data_relevant.set_axis(["CHROM", "POS", "P"], axis=1)

        df = pd.concat([df, result_data_relevant])

    # return dataframe containing all results
    return df


# ------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='create plots for GWAS results.')

parser.add_argument('--resultdir', action='store',
                    help='directory containing GWAS results', required=True)
parser.add_argument('--resultfile', action='store',
                    help='regular expression for result files, e.g. manta*.txt', required=True)
parser.add_argument('--sep', action='store',
                    help='separator', default="\t")
parser.add_argument('--col_chrom', action='store',
                    help='column number of CHROM', required=True, type=int)
parser.add_argument('--col_pos', action='store',
                    help='column number of POS', required=True, type=int)
parser.add_argument('--col_p', action='store',
                    help='column number of p-value', required=True, type=int)
parser.add_argument('--title', action='store',
                    help='title of plot')
parser.add_argument('--sample', action='store_true',
                    help='sample Manhattan plot')

args = parser.parse_args()

print("GWAS plots v0.1")

# show sample
if args.sample:
    manhattan_plot_sample()
    exit(0)

# print arguments
print(f"result directory: {args.resultdir}")
print(f"result file: {args.resultfile}")

# get relevant columns
column_chrom = args.col_chrom - 1
column_pos = args.col_pos - 1
column_p = args.col_p - 1

# read result file
try:
    result_file_glob = os.path.join(args.resultdir, args.resultfile)

    g = glob.glob(result_file_glob)
    len_g = len(g)
    if len_g == 0:
        print("no files found")
        exit(1)
    else:
        print(f"{len_g} files found")

    result_df = process_result_file(glob_files=g,
                                    col_chrom=column_chrom, col_pos=column_pos, col_p=column_p)

except Exception as e:
    print(e)
    exit(1)


# create manhattan plot
manhattan_plot(result_df, title=args.title)
