# gwasplot.py
# create GWAS plots
# sample arguments:
# Miami plot for all GWAS methods
# --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI
# --input /home/kcan/UOC/tfm/mvgwas2-nf/bin/gwasplot/gwasplot_input.csv --plot miami --saveplot adni --showplot

import argparse
import pandas as pd
from pandas import DataFrame

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns


def preprocess_result(single_gwas_data: dict):

    # ----------------------------------
    # preprocessing of the result data
    # ----------------------------------
    df = single_gwas_data["result_df"].copy()
    gwas_method = single_gwas_data["method"]

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

    return df, df_grouped


def manhattan_plot(single_gwas_data: dict, method_colors: dict, figsize: tuple,
                   show_plot: bool = True, save_path: str = None):

    # ----------------------------------
    # preprocessing of the result data
    # ----------------------------------
    df = single_gwas_data["result_df"].copy()
    gwas_method = single_gwas_data["method"]

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

    # ----------------------------
    # display the manhattan plot
    # ----------------------------

    fig = plt.figure(figsize=(figsize[0], figsize[1]))  # Set the figure size for width, height in inches
    plt.clf()  # clear any previous figures

    ax = fig.add_subplot(111)  # nrows, ncols, index

    # create plot
    plot_single_manhattan(gwas_method=gwas_method, df=df, df_grouped=df_grouped,
                          colors=method_colors[gwas_method], ax=ax, y_orientation_up=True)
    # set title
    title_string = f"GWAS plot for {gwas_method}"

    plt.title(title_string, fontweight="bold")

    # save the graph
    if save_path:
        png_file = save_path + "_" + gwas_method + "_manh.png"
        plt.savefig(png_file)
        print(f"{png_file} saved.")

    # show the graph
    if show_plot:
        plt.show()


def miami_plot_all(gwas_data: list, method_colors: dict, figsize: tuple,
                   show_plot: bool = False, save_path: str = None):

    number_of_methods = len(gwas_data)

    # get all comparison pairs
    comparisons = set()

    for method_number1 in range(number_of_methods):
        for method_number2 in range(number_of_methods):
            new_set = None
            if method_number1 < method_number2:
                new_set = frozenset([method_number1, method_number2])
            elif method_number1 > method_number2:
                new_set = frozenset([method_number2, method_number1])

            if new_set is not None:
                comparisons.add(new_set)

    # get a list of all pairs
    pair_list = list()
    for pair in comparisons:
        el_list = []
        for el in pair:
            el_list.append(el)
        pair_list.append(el_list)

    # compare all pairs with a Miami plot

    for pair in pair_list:
        method1 = gwas_data[pair[0]]
        method2 = gwas_data[pair[1]]

        method1["colors"] = method_colors[method1["method"]]
        method2["colors"] = method_colors[method2["method"]]

        miami_plot(gwas_data1=method1, gwas_data2=method2, figsize=figsize, show_plot=show_plot, save_path=save_path)

    pass


def plot_single_manhattan(gwas_method, df, df_grouped, colors, ax, y_orientation_up: bool = True):

    # list of labels and their positions
    x_labels = []
    x_labels_pos = []

    # process all groups, here chromosomes
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='position', y='minuslog10pvalue', color=colors[num % len(colors)], ax=ax)

        group_min = group["position"].min()
        group_max = group["position"].max()

        x_labels.append(str(group_min))
        x_labels_pos.append(group_min)

        # label and position of chromosome, e.g. CHR22
        x_labels.append(name)
        group_x = group_max - (group_max - group_min) / 2
        x_labels_pos.append(group_x)

        x_labels.append(str(group_max))
        x_labels_pos.append(group_max)

    # -----------------------
    # final layout settings
    # -----------------------
    # set x-labels and ticks
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    # min, max values of positions
    position_min = df["position"].min()
    position_max = df["position"].max()
    ax.set_xlim([position_min, position_max])

    max_y = df["minuslog10pvalue"].max()
    if y_orientation_up:
        ax.set_ylim([0, max_y])
    else:
        ax.set_ylim([max_y, 0])

    # axis labels
    number_of_variants = len(df)
    ax.set_xlabel(f"{gwas_method} ({number_of_variants} variants)", fontweight="bold")
    ax.set_ylabel("-log10(pvalue)")

def miami_plot(gwas_data1: dict, gwas_data2: dict,  figsize: tuple,
               show_plot: bool = False, save_path: str = None):

    method1 = gwas_data1['method']
    method2 = gwas_data2['method']

    # preprocess data
    df1, df1_grouped = preprocess_result(gwas_data1)
    df1_position_min = df1["position"].min()
    df1_position_max = df1["position"].max()

    df2, df2_grouped = preprocess_result(gwas_data2)
    df2_position_min = df2["position"].min()
    df2_position_max = df2["position"].max()

    position_min = min(df1_position_min, df2_position_min)
    position_max = max(df1_position_max, df2_position_max)

    # create graph with 2 rows and 1 column
    fig, axs = plt.subplots(2, 1, layout='constrained')
    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])

    # create first manhattan plot, upwards
    plot_single_manhattan(gwas_method=method1, df=df1, df_grouped=df1_grouped,
                          colors=gwas_data1["colors"], ax=axs[0], y_orientation_up=True)

    # create second manhattan plot, downwards
    plot_single_manhattan(gwas_method=method2, df=df2, df_grouped=df2_grouped,
                          colors=gwas_data2["colors"], ax=axs[1], y_orientation_up=False)

    # super title
    fig.suptitle("Miami plot", fontweight="bold")

    # save the graph
    if save_path:
        png_file = save_path + "_" + method1 + "_" + method2 + "_miami.png"
        plt.savefig(png_file)
        print(f"{png_file} saved.")

    # show the plot
    if show_plot:
        plt.show()

    pass


def qqplot(gwas_data: list, method_colors: dict, figsize: tuple, save_path: str = None):

    legend_labels = list()

    plt.clf()  # clear any previous figures

    for single_gwas_data in gwas_data:

        df = single_gwas_data["result_df"].copy()
        gwas_method = single_gwas_data["method"]
        legend_labels.append(gwas_method)

        # preprocess input data
        df = df.set_axis(["chromosome", "position", "pvalue"], axis=1)

        df_sorted = df.sort_values("pvalue")

        df_sorted['minuslog10pvalue'] = -np.log10(df_sorted.pvalue)

        # expected values
        number_of_rows = len(df_sorted)
        df_sorted['expected'] = range(1, number_of_rows + 1)

        df_sorted['minuslog10expected'] = -np.log10(df_sorted.expected / (number_of_rows + 1))

        first = df_sorted['minuslog10pvalue'].iloc[0]
        axisMax = math.ceil(first)

        single_plot = sns.scatterplot(x="minuslog10pvalue", y="minuslog10expected",
                                      data=df_sorted, color=method_colors[gwas_method][0])

    # show legend
    plt.legend(labels=legend_labels)

    # get current figure / Axis
    ax = plt.gca() # get Axes reference created by seaborn
    ax.set(xlabel="Expected", ylabel="Observed",
           title='Observed vs. Expected distribution of p-values (-log10)')
    ax.grid(visible=True)

    ax.figure.set_figwidth(figsize[0])
    ax.figure.set_figheight(figsize[1])

    # plot reference line
    left, right = plt.xlim()
    bottom, top = plt.ylim()
    top_floor = math.floor(top)
    point1 = [0, 0]
    point2 = [top_floor, top_floor]
    x_values = [point1[0], point2[0]]
    y_values = [point1[1], point2[1]]
    plt.plot(x_values, y_values, linestyle="-", color="black")

    if save_path:
        png_file = save_path + "_qq.png"
        ax.figure.savefig(png_file)
        print(f"{png_file} saved.")

    plt.show()


def process_result_file(glob_files: object, sep: str, col_chrom: int, col_pos: int, col_p: int,
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


def read_result_file(result_dir, result_row) -> (str, DataFrame):

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

    except Exception as e:
        print(e)
        exit(1)

    return row_method, result_df


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
parser.add_argument('--sep', action='store',
                    help='separator', default="\t")
parser.add_argument('--showplot', action='store_true', default=False,
                    help='show plot')
parser.add_argument('--saveplot', action='store',
                    help='save plot with given prefix in the result directory')
parser.add_argument('--verbose', action='store_true', default=False,
                    help='verbose output')
parser.add_argument('--input', action='store',
                    help='input csv for all method results')
parser.add_argument('--plot', action='store', default=[], nargs='+',
                    help='select plots of manh, qq, miami')

args = parser.parse_args()

print("GWAS plots v0.2")

result_df = None
all_results = list()

# show sample
if args.input:
    print(f"input csv: {args.input}")

    try:
        gwas_input = pd.read_csv(args.input, sep=",")

        # collect all results
        for index, row in gwas_input.iterrows():
            method, result_df = read_result_file(result_dir=args.resultdir, result_row=row)

            single_result = dict()
            single_result["method"] = method
            single_result["result_df"] = result_df
            all_results.append(single_result)

    except Exception as e:
        print(e)

# create manhattan plot
if args.saveplot:
    save_plot_path = str(os.path.join(args.resultdir, args.saveplot))
else:
    save_plot_path = None

# general plotting settings
gwas_method_colors = {
    "GEMMA": ['peru', 'darkgreen', 'darkblue', 'gold'],
    "MANTA": ['green', 'blue', 'grey', 'red'],
    "MOSTest": ['royalblue', 'pink', 'lightgreen', 'yellow'],
}

figure_size = (18, 10)  # width, height

# default settings for the plots
if args.plot:
    plt.rc('font', size=16)          # controls default text sizes

if args.input:
    if "manh" in args.plot:
        for single_result in all_results:
            manhattan_plot(single_result, method_colors=gwas_method_colors, figsize=figure_size,
                           show_plot=args.showplot, save_path=save_plot_path)

    if "miami" in args.plot:
        miami_plot_all(all_results, method_colors=gwas_method_colors, figsize=figure_size,
                       show_plot=args.showplot, save_path=save_plot_path)

    if "qq" in args.plot:
        # qqplot for multiple result data
        qqplot(all_results, method_colors=gwas_method_colors, figsize=figure_size, save_path=save_plot_path)
