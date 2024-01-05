# gwasplot.py
# create GWAS plots
# sample arguments:
# GEMMA: --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI --resultfile "gemma.*.assoc.txt" --col_chrom 1 --col_pos 3 --col_p 73 --title GEMMA --showplot --saveplot gemma_adni
# MANTA: --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI --resultfile "manta*.txt" --col_chrom 1 --col_pos 2 --col_p 8 --title MANTA --saveplot manta_adni
# MOSTest: --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI --resultfile "mostest_results_genotypes.*.most_perm.sumstats" --col_chrom 1 --col_pos 3 --col_p 6 --title MOSTest --saveplot mostest_adni

import argparse
import pandas as pd
from pandas import DataFrame

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import math

import statsmodels.api as sm
import seaborn as sns


def manhattan_plot_sample():

    # import libraries
    # from pandas import DataFrame
    from scipy.stats import uniform
    from scipy.stats import randint

    # sample data
    df = DataFrame({'gene': ['gene-%i' % i for i in np.arange(10000)],
                    'pvalue': uniform.rvs(size=10000),
                    'chromosome': ['ch-%i' % i for i in randint.rvs(0, 12, size=10000)]})

    # -log_10(pvalue)
    df['minuslog10pvalue'] = -np.log10(df.pvalue)
    df.chromosome = df.chromosome.astype('category')
    df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(12)], ordered=True)
    df = df.sort_values('chromosome')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    df['ind'] = range(len(df))
    df_grouped = df.groupby('chromosome')

    # manhattan plot
    fig = plt.figure(figsize=(14, 8))  # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['darkred', 'darkgreen', 'darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue', color=colors[num % len(colors)], ax=ax)
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


def manhattan_plot(gwas_data: DataFrame, title: str = None,
                   show_plot: bool = True, save_path: str = None):

    df = gwas_data.copy()

    # rename columns
    df = df.set_axis(["chromosome", "position", "pvalue"], axis=1)

    # -log_10(pvalue)
    df['minuslog10pvalue'] = -np.log10(df.pvalue)

    # create category for chromosome field
    df.chromosome = df.chromosome.map(lambda x: "chr" + str(x))
    chromosomes = df.chromosome.unique().tolist()
    df.chromosome = df.chromosome.astype('category')
    # range must represent actual number of chromosomes; here we have only chromosome 1
    df.chromosome = df.chromosome.cat.set_categories(chromosomes, ordered=True)

    df = df.sort_values(['chromosome', 'position'])

    # create an index for the x-axis and group by chromosome
    df['ind'] = range(len(df))
    df_grouped = df.groupby('chromosome')

    # manhattan plot
    fig = plt.figure(figsize=(14, 8))  # Set the figure size
    ax = fig.add_subplot(111)

    colors = ['darkred', 'darkgreen', 'darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue', color=colors[num % len(colors)], ax=ax)
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
    if save_path:
        plt.savefig(save_path)
        print(f"{save_path} saved.")

    if show_plot:
        plt.show()


def qqplot(gwas_data: DataFrame, save_path: str = None):
    # pvals <- read.table("DGI_chr3_pvals.txt", header=T)

    df = gwas_data.copy()

    # rename columns
    df = df.set_axis(["chromosome", "position", "pvalue"], axis=1)

    # observed <- sort(pvals$PVAL)
    df_sorted = df.sort_values("pvalue")

    # -log_10(pvalue)
    # lobs <- -(log10(observed))
    df_sorted['minuslog10pvalue'] = -np.log10(df_sorted.pvalue)

    # expected <- c(1:length(observed))
    number_of_rows = len(df_sorted)
    df_sorted['expected'] = range(1, number_of_rows + 1)

    # lexp <- -(log10(expected / (length(expected)+1)))
    df_sorted['minuslog10expected'] = -np.log10(df_sorted.expected / (number_of_rows + 1))

    # pdf("qqplot.pdf", width=6, height=6)
    # plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)",
    # ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
    # points(lexp, lobs, pch=23, cex=.4, bg="black")

    # create Q-Q plot with 45-degree line added to plot
    # fig = sm.qqplot(data, line='45')

    # sns.lineplot(x="expected", y="minuslog10expected", data=df_sorted)
    # sm.qqplot(df_sorted.minuslog10pvalue)

    # now plot these two arrays against each other using matplotlib
    # retrive pmin, the most significant (i.e. min) p value (for defining the axes)
    # pmin = df_sorted['minuslog10pvalue'][0]

    first = df_sorted['minuslog10pvalue'].iloc[0]
    axisMax = math.ceil(first)

    # fig = plt.figure()

    # plt.xlim([0, axisMax])
    # plt.xlabel("Expected")

    # plt.ylim([0, axisMax])
    # plt.ylabel("Observed")

    # plt.title("QQ Plot: Observed vs. Expected distribution of p values (-log10)")

    # the observed vs. expected data
    # dataAx = fig.add_subplot(111)

    # dataAx.set_xlim([0, axisMax])
    # dataAx.set_ylim([0, axisMax])
    # dataAx.set_xlabel('Expected')
    # dataAx.set_ylabel('Observed')

    # dataAx.plot(df_sorted.minuslog10pvalue, df_sorted.minuslog10expected, 'r.')  # red dots

    # a diagonal line for comparison
    # lineAx = fig.add_subplot(1, 1, 1)
    # dataAx.plot([0, axisMax], [0, axisMax], 'b-')  # blue line

    # qq = sns.lineplot(x="minuslog10expected", y="minuslog10expected", data=df_sorted)
    qq = sns.scatterplot(x="minuslog10pvalue", y="minuslog10expected", data=df_sorted, color='red')
    sns.lineplot(x="minuslog10expected", y="minuslog10expected", data=df_sorted, color='blue')
    qq.set(xlabel="Expected", ylabel="Observed",
           title='Observed vs. Expected distribution of p-values (-log10)')

    plt.legend(labels=['MANTA'])

    if save_path:
        qq.figure.savefig(save_path)

    plt.show()

    pass


def process_result_file(glob_files, sep: str, col_chrom: int, col_pos: int, col_p: int,
                        verbose: bool = False) -> DataFrame:

    df = DataFrame()

    show_once = True

    for result_file in glob_files:
        if verbose:
            print(f"processing result file: {result_file}")

        try:
            result_data = pd.read_csv(result_file, sep=sep)
        except Exception as e:
            print(e)
            exit(1)

        # process result file
        if verbose:
            print(f"rows: {len(result_data)}")

        relevant_columns = result_data.columns[[col_chrom, col_pos, col_p]]

        if show_once:
            if verbose:
                print(f"columns to be processed: {relevant_columns.to_list()}")
            show_once = False

        result_data_relevant = result_data[relevant_columns].copy()

        # rename columns to standard names
        result_data_relevant = result_data_relevant.set_axis(["CHROM", "POS", "P"], axis=1)

        df = pd.concat([df, result_data_relevant])

    # return dataframe containing all results
    return df


def read_result_file(result_dir, result_row):

    # get relevant columns
    row_method = result_row[['method']].item()
    row_col_chrom = result_row[['col_chrom']].item()
    row_col_pos = result_row[['col_pos']].item()
    row_col_p = result_row[['col_p']].item()
    row_resultfile = result_row[['resultfile']].item()

    column_chrom = row_col_chrom - 1
    column_pos = row_col_pos - 1
    column_p = row_col_p - 1

    # read result file
    try:
        result_file_glob = os.path.join(result_dir, row_resultfile)

        g = glob.glob(result_file_glob)
        len_g = len(g)
        if len_g == 0:
            print("no files found")
            exit(1)
        else:
            print(f"{len_g} files found")

        result_df = process_result_file(glob_files=g, sep=args.sep,
                                        col_chrom=column_chrom, col_pos=column_pos, col_p=column_p)
        result_df["method"] = row_method

    except Exception as e:
        print(e)
        exit(1)

    return result_df


# ------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='create plots for GWAS results.',
                                 prog="gwasplot",
                                 epilog="Columns (CHROM,POS,P):" 
                                        "MANTA: 1,2,8 GEMMA: 1,3,73 MOSTest: 1,3,6",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--resultdir', action='store', default='.',
                    help='directory containing GWAS results', required=True)
parser.add_argument('--resultfile', action='store',
                    help='regular expression for result files, e.g. "manta*.txt"')
parser.add_argument('--sep', action='store',
                    help='separator', default="\t")
parser.add_argument('--col_chrom', action='store',
                    help='column number of CHROM', type=int)
parser.add_argument('--col_pos', action='store',
                    help='column number of POS', type=int)
parser.add_argument('--col_p', action='store',
                    help='column number of p-value', type=int)
parser.add_argument('--title', action='store',
                    help='title of plot')
parser.add_argument('--sample', action='store_true',
                    help='sample Manhattan plot')
parser.add_argument('--showplot', action='store_true', default=False,
                    help='show plot')
parser.add_argument('--saveplot', action='store',
                    help='save plot with given prefix in the result directory')
parser.add_argument('--verbose', action='store_true', default=False,
                    help='verbose output')
parser.add_argument('--input', action='store',
                    help='input csv for all method results')

args = parser.parse_args()

print("GWAS plots v0.1")

if args.input:
    print(f"input csv: {args.input}")

    try:
        gwas_input = pd.read_csv(args.input, sep=",")

        for index, row in gwas_input.iterrows():
            result_df = read_result_file(result_dir=args.resultdir, result_row=row)

    except Exception as e:
        print(e)

    exit(0)

# show sample
if args.sample:
    manhattan_plot_sample()
    exit(0)

# print arguments
if args.verbose:
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

    result_df = process_result_file(glob_files=g, sep=args.sep,
                                    col_chrom=column_chrom, col_pos=column_pos, col_p=column_p)

except Exception as e:
    print(e)
    exit(1)


# create manhattan plot
if args.saveplot:
    save_plot_path = str(os.path.join(args.resultdir, args.saveplot)) + ".png"
    save_qqplot_path = str(os.path.join(args.resultdir, args.saveplot)) + "_qq.png"
else:
    save_plot_path = None
    save_qqplot_path = None

# manhattan_plot(result_df, title=args.title, show_plot=args.showplot, save_path=save_plot_path)
qqplot(result_df, save_path=save_qqplot_path)
