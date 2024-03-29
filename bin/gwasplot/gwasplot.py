# gwasplot.py
# create GWAS plots
# sample arguments:
# Miami plot for all GWAS methods
# --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI
# --input /home/kcan/UOC/tfm/mvgwas2-nf/bin/gwasplot/gwasplot_input.csv --plot miami --saveplot adni --showplot
#
# --resultdir /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI --input /home/kcan/UOC/tfm/mvgwas2-nf/bin/gwasplot/gwasplot_input.csv
# --plot manh miami qq box --saveplot adni --showplot --remove_outliers
# --only_common --ntop 10 --genotypes /home/kcan/UOC/tfm/mvgwas2-nf/data/ADNI/genotypes.chr22.vcf.gz
# --phenotypes /home/kcan/UOC/tfm/mvgwas2-nf/data/ADNI/phenotypes.tsv

import argparse
import pandas as pd
from pandas import DataFrame

import os
from datetime import datetime
import glob
import gzip

import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns

import matplotlib_venn


def preprocess_result(single_gwas_data: dict):

    # ----------------------------------
    # preprocessing of the result data
    # ----------------------------------
    df = single_gwas_data["result_df"].copy()
    gwas_method = single_gwas_data["method"]

    # rename columns
    df = df[["CHROM", "POS", "P"]]  # consider only these columns
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
    df_grouped = df.groupby('chromosome', observed=False)

    return df, df_grouped


def manhattan_plot(single_gwas_data: dict, method_colors: dict, figsize: tuple,
                   show_plot: bool = True, save_path: str = None, significance_level_log: float = None,
                   adjusted_threshold_log: float = None):

    # ----------------------------------
    # preprocessing of the result data
    # ----------------------------------
    df = single_gwas_data["result_df"].copy()
    gwas_method = single_gwas_data["method"]

    # rename columns
    df = df[["CHROM", "POS", "P"]]  # consider only these columns
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
    df_grouped = df.groupby('chromosome', observed=False)

    # ----------------------------
    # display the manhattan plot
    # ----------------------------

    fig = plt.figure(figsize=(figsize[0], figsize[1]))  # Set the figure size for width, height in inches
    plt.clf()  # clear any previous figures

    ax = fig.add_subplot(111)  # nrows, ncols, index

    # create plot
    plot_single_manhattan(gwas_method=gwas_method, df=df, df_grouped=df_grouped,
                          colors=method_colors[gwas_method], ax=ax, y_orientation_up=True,
                          significance_level_log=significance_level_log,
                          adjusted_threshold_log=adjusted_threshold_log)
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
                   show_plot: bool = False, save_path: str = None,
                   significance_level_log: float = None,
                   adjusted_threshold_log: float = None):

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

        miami_plot(gwas_data1=method1, gwas_data2=method2, figsize=figsize,
                   show_plot=show_plot, save_path=save_path,
                   significance_level_log=significance_level_log,
                   adjusted_threshold_log=adjusted_threshold_log)

    pass


def plot_single_manhattan(gwas_method, df, df_grouped, colors, ax, y_orientation_up: bool = True,
                          position_min: int = 0, position_max: int = 0,
                          significance_level_log: float = None,
                          adjusted_threshold_log: float = None):

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


    # set axis limits
    # min, max values of positions
    df_position_min = df["position"].min()
    df_position_max = df["position"].max()
    if (position_min > 0) and (position_max > 0):
        # common min/max position of both GWAS method results is given
        ax.set_xlim([position_min, position_max])

        # correct first and last xlabel and its position
        x_labels_pos[0] = position_min
        x_labels[0] = str(position_min)

        x_labels_pos[-1] = position_max
        x_labels[-1] = str(position_max)

    else:
        # use min/max of this GWAS method
        ax.set_xlim([df_position_min, df_position_max])

    # set x-labels and ticks
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    max_y_data = df["minuslog10pvalue"].max()
    if significance_level_log is None:
        max_y = max_y_data
    else:
        max_y = max(max_y_data, significance_level_log)
        max_y = math.ceil(max_y)

        # plot line of significance level and adjusted significance level
        if significance_level_log is not None:
            ax.axhline(y=significance_level_log, color="red")
        if adjusted_threshold_log is not None:
            ax.axhline(y=adjusted_threshold_log, color="black")

    if y_orientation_up:
        ax.set_ylim([0, max_y])
    else:
        ax.set_ylim([max_y, 0])

    # axis labels
    number_of_variants = len(df)
    ax.set_xlabel(f"{gwas_method} ({number_of_variants} variants)", fontweight="bold")
    ax.set_ylabel("-log10(pvalue)")

def miami_plot(gwas_data1: dict, gwas_data2: dict,  figsize: tuple,
               show_plot: bool = False, save_path: str = None,
               significance_level_log: float = None, adjusted_threshold_log: float = None):

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
                          colors=gwas_data1["colors"], ax=axs[0], y_orientation_up=True,
                          position_min=position_min, position_max=position_max,
                          significance_level_log=significance_level_log,
                          adjusted_threshold_log=adjusted_threshold_log)

    # create second manhattan plot, downwards
    plot_single_manhattan(gwas_method=method2, df=df2, df_grouped=df2_grouped,
                          colors=gwas_data2["colors"], ax=axs[1], y_orientation_up=False,
                          position_min=position_min, position_max=position_max,
                          significance_level_log=significance_level_log,
                          adjusted_threshold_log=adjusted_threshold_log)

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


def plot_venn_diagram(a: set, b: set, c: set, labels: list, title: str, save_path: str = None,):

    # set operations: difference and intersections
    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)

    only_a_b = len(a & b - c)
    only_a_c = len(a & c - b)
    only_b_c = len(b & c - a)

    a_b_c = len(a & b & c)

    matplotlib_venn.venn3(subsets=(only_a, only_b, only_a_b, only_c, only_a_c, only_b_c, a_b_c),
                          set_labels=labels)

    # set title
    plt.title(title, fontweight="bold")

    # save the graph
    if save_path:
        png_file = save_path + "_venn.png"
        plt.savefig(png_file)
        print(f"{png_file} saved.")

    # show plot
    plt.show()

def venn_diagram(gwas_data: list, significance_threshold: float, method_colors: dict, figsize: tuple, save_path: str = None):

    significant_results = dict()

    for method_result in gwas_data:
        df = method_result["result_df"].copy()
        df['minuslog10pvalue'] = -np.log10(df["P"])
        if significance_threshold is not None:
            df_filtered = df[df['minuslog10pvalue'] >= significance_threshold]
        else:
            df_filtered = df
        df_filtered = df_filtered[["CHROM", "POS", "minuslog10pvalue"]]  # consider only these columns

        significant_results[method_result["method"]] = df_filtered.copy()

    # intersect the results
    all_chrom_pos_set = set()
    method_locus = dict()
    for method in significant_results: # eg GEMMA
        df = significant_results[method]
        method_locus[method] = set()
        for index, row in df.iterrows():
            row_chrom = int(row["CHROM"].item())
            row_pos = int(row["POS"].item())
            locus = (row_chrom, row_pos)
            all_chrom_pos_set.add(locus)
            method_locus[method].add(locus)

    all_chrom_pos_df = pd.DataFrame(all_chrom_pos_set)
    all_chrom_pos_df = all_chrom_pos_df.set_axis(["CHROM", "POS"], axis=1)
    all_chrom_pos_sorted_df = all_chrom_pos_df.sort_values(["CHROM", "POS"])

    # join all method results
    all_chrom_pos_join_df = all_chrom_pos_sorted_df.copy()

    for method_result in gwas_data:
        df = method_result["result_df"]
        df = df[["CHROM", "POS", "P"]].copy()  # consider only these columns
        df["method"] = method_result["method"]
        method_p = "P_" + method_result["method"]
        df = df.set_axis(["CHROM", "POS", method_p, "method"], axis=1)
        # all_chrom_pos_join_df = all_chrom_pos_join_df.join(df, on=["CHROM","POS"], how="left")
        all_chrom_pos_join_df = pd.merge(all_chrom_pos_join_df, df, on=["CHROM","POS"], how="left")
        pass

    if significance_threshold is None:
        title_string = "Venn diagram for GWAS methods"
    else:
        title_string = f"Venn diagram for threshold {significance_threshold}"

    plot_venn_diagram(method_locus["GEMMA"], method_locus["MANTA"], method_locus["MOSTest"],
                      labels=["GEMMA", "MANTA", "MOSTest"], title=title_string, save_path=save_path)
    pass


def qqplot(gwas_data: list, method_colors: dict, figsize: tuple, save_path: str = None):

    legend_labels = list()
    number_of_gwas_data = len(gwas_data)

    plt.clf()  # clear any previous figures

    for single_gwas_data in gwas_data:

        df = single_gwas_data["result_df"].copy()
        gwas_method = single_gwas_data["method"]
        legend_labels.append(gwas_method)

        # preprocess input data
        df = df[["CHROM", "POS", "P"]]  # consider only these columns
        df = df.set_axis(["chromosome", "position", "pvalue"], axis=1)

        df_sorted = df.sort_values("pvalue")

        df_sorted['minuslog10pvalue'] = -np.log10(df_sorted.pvalue)

        # expected values
        number_of_rows = len(df_sorted)
        df_sorted['expected'] = range(1, number_of_rows + 1)

        df_sorted['minuslog10expected'] = -np.log10(df_sorted.expected / (number_of_rows + 1))

        first = df_sorted['minuslog10pvalue'].iloc[0]
        axisMax = math.ceil(first)

        single_plot = sns.scatterplot(x="minuslog10expected", y="minuslog10pvalue", s=150, alpha=1.0,
                                      data=df_sorted, color=method_colors[gwas_method][0])

    # show legend
    default_fontsize = 40
    legend = plt.legend(labels=legend_labels, fontsize=default_fontsize)

    # show larger legend markers
    for handle in legend.legend_handles:
        handle.set_sizes([500.0])

    # get current figure / Axis
    ax = plt.gca() # get Axes reference created by seaborn
    ax.set_xlabel("Expected $-log_{10}(p)$", fontweight="bold", fontsize=default_fontsize)
    ax.set_ylabel("Observed $-log_{10}(p)$", fontweight="bold", fontsize=default_fontsize)
    ax.set_title("Observed vs. Expected distribution of p-values",
                 fontweight="bold", fontsize=default_fontsize)
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
        if number_of_gwas_data == 1:
            png_file = save_path + "_" + gwas_data[0]["method"] + "_qq.png"
        else:
            png_file = save_path + "_qq.png"

        ax.figure.savefig(png_file)
        print(f"{png_file} saved.")

    plt.show()


def top_snp_boxplots(snp_genotype_volumes: dict, figsize: tuple, number_of_top_snps: int = 3, save_path: str = None):

    i = 0
    for snp in snp_genotype_volumes:
        boxplot(snp_genotype_volumes=snp_genotype_volumes, snp=snp, figsize=figsize,  save_path=save_path)
        i += 1
        if i == number_of_top_snps:
            break


def compare_p(all_results: dict, figsize: tuple, save_path: str = None):
    # compare p-values

    all_pvalues_df = None

    for single_result in all_results:

        method_p_values = single_result["result_df"].copy()
        method_p_values = method_p_values[["CHROM","POS","P_MINUS_LOG10"]]
        method_p_values["method"] = single_result["method"]

        if all_pvalues_df is None:
            all_pvalues_df = method_p_values.copy()
        else:
            all_pvalues_df = pd.concat([all_pvalues_df, method_p_values])
        pass

    default_fontsize = 30
    ax = plt.gca() # get Axes reference created by seaborn
    ax.set_title("Comparison of -log10(p)", fontweight="bold", fontsize=default_fontsize)
    ax.set_xlabel("GWAS method", fontsize=default_fontsize)
    ax.set_ylabel("$-log_{10}(p)$", fontweight="bold", fontsize=default_fontsize)
    ax.figure.set_figwidth(figsize[0])
    ax.figure.set_figheight(figsize[1])

    ax_boxplot = sns.boxplot(x="method", y="P_MINUS_LOG10", hue="CHROM", data=all_pvalues_df)

    # save plot
    if save_path:
        png_file = save_path + "_pcomp_box.png"
        ax.figure.savefig(png_file)
        print(f"{png_file} saved.")

    # show plot
    plt.show()

    # show distribution of p-values
    all_pvalues_df["P_MINUS_LOG10_ROUNDED"] = np.round(all_pvalues_df["P_MINUS_LOG10"])
    all_pvalues_df = all_pvalues_df.astype({'P_MINUS_LOG10_ROUNDED': 'int'})
    all_pvalues_grouped_df = all_pvalues_df.groupby(["method","P_MINUS_LOG10_ROUNDED"]).agg('count').reset_index()
    all_pvalues_grouped_df = all_pvalues_grouped_df.drop(["CHROM", "POS"], axis=1)
    all_pvalues_grouped_df.rename(columns={"P_MINUS_LOG10": "Count"}, inplace=True)

    plt.clf()  # clear any previous figures
    default_fontsize = 30
    ax = plt.gca() # get Axes reference created by seaborn
    ax.set_title("Comparison of -log10(p) counts", fontweight="bold", fontsize=default_fontsize)
    ax.set_xlabel("$-log_{10}(p)$", fontsize=default_fontsize)
    ax.set_ylabel("count", fontsize=default_fontsize)
    ax.figure.set_figwidth(figsize[0])
    ax.figure.set_figheight(figsize[1])

    sns.barplot(data=all_pvalues_grouped_df, x="P_MINUS_LOG10_ROUNDED", y="Count", hue="method")

    # save plot
    if save_path:
        png_file = save_path + "_pcomp_bar.png"
        ax.figure.savefig(png_file)
        print(f"{png_file} saved.")

    plt.show()

    pass

def boxplot(snp_genotype_volumes: dict, figsize: tuple, snp: str, save_path: str = None):

    # SNP, genotype, hippocampus subfield, sample volumes
    df = snp_genotype_volumes[snp]

    all_subfields_df = None

    genotype_mapping =  {
        "genotype_0_0": "0/0",
        "genotype_0_1": "0/1",
        "genotype_1_1": "1/1",
        "genotype_oth": "oth"
    }

    for genotype in df:
        for subfield in df[genotype]:
            sample_volumes = snp_genotype_volumes[snp][genotype][subfield]["sample_volumes"]
            df_subfield_volumes = pd.DataFrame({"volume": sample_volumes})
            df_subfield_volumes["subfield"] = subfield
            df_subfield_volumes["genotype"] = genotype_mapping[genotype]

            if all_subfields_df is None:
                all_subfields_df = df_subfield_volumes.copy()
            else:
                all_subfields_df = pd.concat([all_subfields_df, df_subfield_volumes])
            pass

    default_fontsize = 30
    ax = plt.gca() # get Axes reference created by seaborn
    ax.set_title("Volumes of hippocampus subfields of " + snp, fontweight="bold", fontsize=default_fontsize)
    ax.figure.set_figwidth(figsize[0])
    ax.figure.set_figheight(figsize[1])

    ax_boxplot = sns.boxplot(x="subfield", y="volume", hue="genotype", data=all_subfields_df)

    # rotate the x-axis labels
    ax_boxplot.set_xticks(ax_boxplot.get_xticks())
    ax_boxplot.set_xticklabels(ax.get_xticklabels(), rotation=30)

    # save plot
    if save_path:
        png_file = save_path + "_" + snp + "_box.png"
        ax.figure.savefig(png_file)
        print(f"{png_file} saved.")

    # show plot
    plt.show()

    pass


def process_result_file(glob_files, sep: str, result_dir: str,
                        outlier_threshold: float, remove_outliers: bool,
                        col_chrom: int, col_pos: int, col_p: int,
                        col_snp: int, col_ref: int, col_alt: int,
                        verbose: bool = False) -> tuple:

    df = DataFrame()

    show_once = True
    first_outlier_found = False
    all_outliers_df = DataFrame()
    outliers_removed = 0

    number_of_na_found = 0

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

        # gemma: rs:2, allele1:5, allele0:6, af:7; manta: ID:3, REF:4, ALT:5, F:6; MOSTest: SNP:2, A1:4, A2:5
        relevant_columns = result_data.columns[[col_chrom, col_pos, col_p, col_snp, col_ref, col_alt]]

        if show_once:
            if verbose:
                print(f"columns to be processed: {relevant_columns.to_list()}")
            show_once = False

        result_data_relevant = result_data[relevant_columns].copy()

        # rename columns to standard names
        result_data_relevant = result_data_relevant.set_axis(["CHROM", "POS", "P", "SNP", "REF", "ALT"], axis=1)

        # check for NANs
        count_na = result_data_relevant["P"].isna().sum()
        if count_na > 0:
            if verbose:
                print(f"{count_na} NAs in {result_file}")
            number_of_na_found += count_na
            result_data_relevant = result_data_relevant.dropna()

        # check for outliers
        outliers_df = result_data_relevant[result_data_relevant["P"] < outlier_threshold]
        outliers_len = len(outliers_df.index)
        if outliers_len > 0:
            outliers_copy_df = outliers_df.copy()
            outliers_copy_df["file"] = result_file
            if not first_outlier_found:
                first_outlier_found = True
                all_outliers_df = outliers_copy_df.copy()
            else:
                all_outliers_df = pd.concat([all_outliers_df, outliers_copy_df])

            for df_row in range(outliers_len):
                outlier_chrom = outliers_df["CHROM"].values[df_row]
                outlier_pos = outliers_df["POS"].values[df_row]
                outlier_p = outliers_df["P"].values[df_row]
                if verbose:
                    print(f"outlier found for {result_file}, CHROM: {outlier_chrom}, POS: {outlier_pos}, P: {outlier_p}")

            if remove_outliers:
                len_result_data_relevant_before = len(result_data_relevant.index)
                result_data_relevant = result_data_relevant[result_data_relevant["P"] >= outlier_threshold]
                # result_data_relevant_ok = result_data_relevant[result_data_relevant["P"] >= outlier_threshold]
                # result_data_relevant_notok = result_data_relevant[result_data_relevant["P"] < outlier_threshold]

                # result_data_relevant = result_data_relevant_ok.copy()
                len_result_data_relevant_after = len(result_data_relevant.index)

                number_of_outliers_removed = len_result_data_relevant_before - len_result_data_relevant_after

                outliers_removed += number_of_outliers_removed
                if verbose and number_of_outliers_removed > 0:
                    print(f"{number_of_outliers_removed} outliers removed")

        # append result data of current file
        df = pd.concat([df, result_data_relevant])

    # save outlier report
    all_outliers_df_len = len(all_outliers_df)
    if all_outliers_df_len > 0:
        now = datetime.now()
        now_formatted =  now.strftime("%Y%m%d_%H%M%S")
        outlier_file = f"outlier_{now_formatted}.csv"

        outlier_path = os.path.join(result_dir, outlier_file)
        all_outliers_df.to_csv(outlier_path, index=False)

        # print(f"{all_outliers_df_len} outliers found with p-values < {outlier_threshold} and saved to {outlier_path}")
        # print(f"{outliers_removed} outliers removed")

    # print(f"{number_of_na_found} NAs found")

    # return dataframe containing all results
    return df, number_of_na_found, all_outliers_df_len, outliers_removed


def read_result_file(result_dir, result_row,
                     outlier_threshold, remove_outliers: bool = False) -> tuple:

    # get relevant columns
    row_method = result_row[['method']].item()
    row_col_chrom = result_row[['col_chrom']].item()
    row_col_pos = result_row[['col_pos']].item()
    row_col_p = result_row[['col_p']].item()
    row_resultfile = result_row[['resultfile']].item()

    # additional columns for variant/SNP, Reference allele and Alternative allele
    row_col_snp = result_row[['col_snp']].item()
    row_col_ref = result_row[['col_ref']].item()
    row_col_alt = result_row[['col_alt']].item()

    column_chrom = row_col_chrom - 1
    column_pos = row_col_pos - 1
    column_p = row_col_p - 1

    column_snp = row_col_snp - 1
    column_ref = row_col_ref - 1
    column_alt = row_col_alt - 1

    # read result file
    try:
        result_file_glob = os.path.join(result_dir, row_resultfile)

        g = glob.glob(result_file_glob)
        len_g = len(g)

        result_df, number_of_na_found, all_outliers, outliers_removed =\
            process_result_file(glob_files=g, sep=args.sep, result_dir=result_dir,
                                outlier_threshold=outlier_threshold, remove_outliers=remove_outliers,
                                col_chrom=column_chrom, col_pos=column_pos, col_p=column_p,
                                col_snp=column_snp, col_ref=column_ref, col_alt=column_alt
                                )

    except Exception as e:
        print(e)
        exit(1)

    return row_method, result_df, number_of_na_found, all_outliers, outliers_removed, len_g


# keep only variants which are common to all results
def intersect_all_results(all_results_list: list, top_n: int = 10) -> tuple:

    all_positions = []

    # create list of chromosome-position tuples of all results / methods
    for single_result in all_results_list:
        result_df = single_result["result_df"]

        chrom_pos_set = set()

        for index, row in result_df.iterrows():
            row_chrom = int(row["CHROM"])
            row_pos = int(row["POS"])
            chromosome_position = (row_chrom, row_pos)
            chrom_pos_set.add(chromosome_position)

        all_positions.append(chrom_pos_set)

    # intersect all chromosome-position sets
    common_positions_set = None
    for single_position_set in all_positions:
        if common_positions_set is None:
            common_positions_set = single_position_set
        else:
            common_positions_set = common_positions_set.intersection(single_position_set)

    # convert set to dataframe
    common_positions_df = pd.DataFrame(common_positions_set, columns =['CHROM', 'POS'])


    # keep only results which have common chrom/pos rows
    common_results_list = list()
    all_snps_df = None
    p_columns = list()

    for single_result in all_results_list:
        common_single_result = dict()
        gwas_method = single_result["method"]
        common_single_result["method"] = single_result["method"]

        int_df = pd.merge(single_result["result_df"], common_positions_df, how='inner', on=['CHROM', 'POS'])
        common_single_result["result_df"] = int_df
        common_results_list.append(common_single_result)

        # primary key is CHROM, POS; rename the remaining columns which are specific to the GWAS method
        # create a dictionary: key = old name, value = new name
        rename_columns = {
            'P': 'P_' + gwas_method,
            'SNP': 'SNP_' + gwas_method,
            'REF': 'REF_' + gwas_method,
            'ALT': 'ALT_' + gwas_method
        }
        int_method_df = int_df.rename(columns=rename_columns)

        # store the column names of all p-value columns for later minimum processing
        p_columns.append('P_' + gwas_method)

        if all_snps_df is None:
            all_snps_df = int_method_df.copy()
        else:
            # merge
            all_snps_df = pd.merge(all_snps_df, int_method_df, how='inner', on=['CHROM', 'POS'])

    # top_snps = common_results_list.nlargest(5, "c")

    # all_snps_df["P_MIN"] = np.min(all_snps_df[['flow_h','flow_c']],axis=1)
    all_snps_df["P_MIN"] = all_snps_df[p_columns].min(axis=1)  # per row
    all_snps_df['P_MIN_MINUS_LOG10'] = -np.log10(all_snps_df["P_MIN"])
    top_snps_df = all_snps_df.nsmallest(top_n, "P_MIN")

    return common_results_list, top_snps_df


def process_top_snps(top_snps_df: DataFrame, genotypes: str, phenotypes: str) -> dict:

    # create list of (CHROM,POS) of the top-SNPs
    top_snps_set = set()
    for index, row in top_snps_df.iterrows():
        snp_tuple = (row["CHROM"], row["POS"])
        top_snps_set.add(snp_tuple)
    number_of_top_snps = len(top_snps_df)

    # read phenotypes
    phenotypes_df = pd.read_csv(phenotypes, sep="\t")
    phenotypes_df.set_index("ID", inplace=True)

    # sum all volumes per region
    # pheno_legend.txt:
    # 1       ID      Individual_ID
    # 2       ST131HS LeftCA1
    # 3       ST132HS LeftCA2_3
    # 4       ST133HS LeftCA4_DG
    # 5       ST134HS LeftFimbria
    # 6       ST135HS LeftHippocampalFissure
    # 7       ST136HS LeftPresubiculum
    # 8       ST137HS LeftSubiculum
    # 9       ST138HS LeftTail
    # 10      ST139HS RightCA1
    # 11      ST140HS RightCA2_3
    # 12      ST141HS RightCA4_DG
    # 13      ST142HS RightFimbria
    # 14      ST143HS RightHippocampalFissure
    # 15      ST144HS RightPresubiculum
    # 16      ST145HS RightSubiculum
    # 17      ST146HS RightTail
    phenotypes_df['volumes'] = phenotypes_df.sum(axis=1, numeric_only=True)

    # sum the 8 regions of the hippocampus (left / right)
    phenotypes_df['CA1'] = phenotypes_df['ST131HS'] + phenotypes_df['ST139HS']
    phenotypes_df['CA2_3'] = phenotypes_df['ST132HS'] + phenotypes_df['ST140HS']
    phenotypes_df['CA4_DG'] = phenotypes_df['ST133HS'] + phenotypes_df['ST141HS']
    phenotypes_df['Fimbria'] = phenotypes_df['ST134HS'] + phenotypes_df['ST142HS']
    phenotypes_df['HippocampalFissure'] = phenotypes_df['ST135HS'] + phenotypes_df['ST143HS']
    phenotypes_df['Presubiculum'] = phenotypes_df['ST136HS'] + phenotypes_df['ST144HS']
    phenotypes_df['Subiculum'] = phenotypes_df['ST137HS'] + phenotypes_df['ST145HS']
    phenotypes_df['Tail'] = phenotypes_df['ST138HS'] + phenotypes_df['ST146HS']

    hippocampus_subfields = ['CA1', 'CA2_3', 'CA4_DG', 'Fimbria', 'HippocampalFissure',
                             'Presubiculum', 'Subiculum', 'Tail']

    # read relevant genotypes entries from gz-file
    chrom_found = False
    line_count = 0
    sample_columns_start = None
    snps_found = 0

    count_0_0 = 0
    count_0_1 = 0
    count_1_1 = 0
    count_oth = 0

    snp_genotype_samples = dict()


    with gzip.open(genotypes,'r') as f:
        for line in f:
            line_count += 1
            str_line = line.decode()

            # process #CHROM-line
            if not chrom_found:
                if str_line.startswith("#CHROM"):
                    chrom_found = True
                    chrom_line = str_line
                    chrom_line_number = line_count
                    print(f"#CHROM found in line {chrom_line_number}")
                    chrom_line_split = chrom_line.split("\t")

                    # get sample IDs (start after the FORMAT columns)
                    index_format = chrom_line_split.index("FORMAT")
                    sample_columns_start = index_format+1
                    sample_ids = chrom_line_split[sample_columns_start:]

                    continue

            # compare loci after #CHROM line has been found
            if chrom_found:
                locus_line = str_line.split("\t")
                line_chrom = int(locus_line[0])
                line_pos = int(locus_line[1])
                line_locus = (line_chrom, line_pos)

                # compare CHROM and POS values with those found in the top-SNPs
                if line_locus in top_snps_set:
                    snps_found +=  1
                    print(f"top SNP {snps_found} found: locus {line_locus} in line {line_count}")
                    locus_samples_list = locus_line[sample_columns_start:]

                    locus_snp = locus_line[2]
                    if locus_snp == ".":  # SNP not in dbSNP-catalog
                        snp_reference = f"chr{line_chrom}_{line_pos}"
                    else:
                        snp_reference = locus_snp  # e.g. rs12345

                    # group locus samples by genotype 0/0, 0/1, 1/1
                    snp_genotype_samples[snp_reference] = {
                        "genotype_0_0": {"sample_id": list()},
                        "genotype_0_1": {"sample_id": list()},
                        "genotype_1_1": {"sample_id": list()},
                        "genotype_oth": {"sample_id": list()},
                    }

                    locus_sample_number = 0
                    for locus_sample in locus_samples_list:
                        sample_split = locus_sample.split(":")
                        geno_type = sample_split[0]  # 0/0, 0/1, 1/1, or unknown ./.
                        if geno_type == "0/0":
                            count_0_0 += 1
                            snp_genotype_samples[snp_reference]["genotype_0_0"]["sample_id"]\
                                .append(sample_ids[locus_sample_number])
                        elif geno_type == "0/1":
                            count_0_1 += 1
                            snp_genotype_samples[snp_reference]["genotype_0_1"]["sample_id"]\
                                .append(sample_ids[locus_sample_number])
                        elif  geno_type == "1/1":
                            count_1_1 += 1
                            snp_genotype_samples[snp_reference]["genotype_1_1"]["sample_id"]\
                                .append(sample_ids[locus_sample_number])
                        else:
                            snp_genotype_samples[snp_reference]["genotype_oth"]["sample_id"]\
                                .append(sample_ids[locus_sample_number])
                            count_oth += 1

                        locus_sample_number += 1

                    if snps_found == number_of_top_snps:
                        # no need to look for further SNPs since all top-SNPs already found
                        break
                pass

    # combine locus / genotype with phenotype values (volumes of hippocampus subfields)
    # for all SNPs
    snp_genotype_volumes = dict()
    for snp in snp_genotype_samples:  # rs1018834
        snp_genotype_volumes[snp] = dict()

        # for all genotypes
        for genotype in snp_genotype_samples[snp]:  # genotype_0_0
            snp_genotype_volumes[snp][genotype] = dict()

            # get values for all hippocampus subfields
            for subfield in hippocampus_subfields:
                snp_genotype_volumes[snp][genotype][subfield] = dict()

                # get values
                sample_count = 0
                sample_sum = 0.0
                sample_volume_list = list()

                for sample in snp_genotype_samples[snp][genotype]['sample_id']:  # 136_S_0873
                    try:
                        subfield_samples = phenotypes_df[subfield]
                        subfield_vol_numpy = subfield_samples[sample]
                        sample_count += 1

                        try:
                            subfield_vol = subfield_vol_numpy.item()
                        except ValueError as e:
                            # entry has duplicates, take the first value
                            subfield_vol = subfield_vol_numpy.iloc[0].item()

                        sample_sum += subfield_vol
                        sample_volume_list.append(subfield_vol)
                    except KeyError:
                        pass

                snp_genotype_volumes[snp][genotype][subfield]["sample_count"] = sample_count
                snp_genotype_volumes[snp][genotype][subfield]["sample_sum"] = sample_sum
                snp_genotype_volumes[snp][genotype][subfield]["sample_volumes"] = sample_volume_list
                if sample_count > 0:
                    snp_genotype_volumes[snp][genotype][subfield]["sample_avg"] = sample_sum / sample_count
                else:
                    snp_genotype_volumes[snp][genotype][subfield]["sample_avg"] = 0.0

    # create box plots
    return snp_genotype_volumes


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
parser.add_argument('--outlier_threshold', action='store',
                    help='threshold for outlier check of p-values', default="1e-10")
parser.add_argument('--remove_outliers', action='store_true', default=False,
                    help='remove outliers')
parser.add_argument('--plot', action='store', default=[], nargs='+',
                    help='select plots of qq manh miami box venn pcomp')
parser.add_argument('--only_common', action='store_true', default=False,
                    help='keep only common variants of all methods')
parser.add_argument('--ntop', action='store', type=int, default=10,
                    help='select top n SNPs; used with --only_common')
parser.add_argument('--threshold', action='store', type=float,
                    help='adjusted GWAS-threshold')
parser.add_argument('--vennthreshold', action='store', type=float,
                    help='adjusted GWAS-threshold for Venn diagram')
parser.add_argument('--genotypes', action='store',
                    help='genotypes gz-file')
parser.add_argument('--phenotypes', action='store',
                    help='phenotypes file')

args = parser.parse_args()

print("GWAS plots v1.1")

print("\nparameters used:\n")

sep = "{:02x}x".format(ord(args.sep))
if sep == "09x":
    sep = "\\t"

outlier_threshold = float(args.outlier_threshold)

print(f"resultdir        : {args.resultdir}")
print(f"input            : {args.input}")
print(f"sep              : {sep}")
print(f"outlier_threshold: {outlier_threshold}")
print(f"remove_outliers  : {args.remove_outliers}")
print(f"only_common      : {args.only_common}")
print(f"ntop             : {args.ntop}")
print(f"showplot         : {args.showplot}")
print(f"saveplot         : {args.saveplot}")
print(f"plot             : {args.plot}")
print(f"verbose          : {args.verbose}")
print(f"genotypes        : {args.genotypes}")
print(f"phenotypes       : {args.phenotypes}")
print(f"threshold        : {args.threshold}")
print(f"vennthreshold    : {args.vennthreshold}")
print()


result_df = None
all_results = list()
snp_genotype_volumes = None

significance_level = 5e-8
significance_level_log = - math.log10(significance_level)

if args.threshold:
    adjusted_threshold = args.threshold
    adjusted_threshold_log = - math.log10(adjusted_threshold)
else:
    adjusted_threshold = None
    adjusted_threshold_log = None

if args.vennthreshold:
    venn_threshold = args.vennthreshold
    venn_threshold_log = - math.log10(venn_threshold)
else:
    venn_threshold = None
    venn_threshold_log = None

print(f"GWAS significance level  : {significance_level}, -log10: {significance_level_log}")
print(f"GWAS adjusted sign. level: {adjusted_threshold}, -log10: {adjusted_threshold_log}\n")
print(f"Venn threshold           : {venn_threshold}, -log10: {venn_threshold_log}\n")

# read GWAS result data
if args.input:
    print(f"input csv: {args.input}")

    try:
        gwas_input = pd.read_csv(args.input, sep=",")

        # collect all results
        for index, row in gwas_input.iterrows():

            method, result_df, number_na_found, number_outliers, number_outliers_removed, files_found =\
                read_result_file(result_dir=args.resultdir, result_row=row,
                                                 outlier_threshold=outlier_threshold,
                                                 remove_outliers=args.remove_outliers)

            print(f"{method:10}: NAs: {number_na_found}, files: {files_found} " +
                  f"outlier: {number_outliers}, outlier removed: {number_outliers_removed}")

            # add -log10(p) to the results
            result_df['P_MINUS_LOG10'] = -np.log10(result_df["P"])

            single_result = dict()
            single_result["method"] = method
            single_result["result_df"] = result_df
            all_results.append(single_result)

        if args.only_common:
            all_results, top_snps = intersect_all_results(all_results, top_n=args.ntop)

            # save top SNPs
            now = datetime.now()
            now_formatted = now.strftime("%Y%m%d_%H%M%S")
            top_snp_file = f"top_snps_{now_formatted}.csv"

            top_snp_path = os.path.join(args.resultdir, top_snp_file)
            top_snps.to_csv(top_snp_path, index=False)

            snp_genotype_volumes = process_top_snps(top_snps, args.genotypes, args.phenotypes)

        print()

    except Exception as e:
        print(e)

# settings for saving the plots
if args.saveplot:
    if args.only_common:
        saveplot_prefix = args.saveplot + "_common"
    else:
        saveplot_prefix = args.saveplot
    save_plot_path = str(os.path.join(args.resultdir, saveplot_prefix))
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

# create the plots: Manhattan, Miami, QQ, Venn, boxplots
if args.input:
    if "pcomp" in args.plot:
        # compare p-values
        compare_p(all_results, figsize=figure_size, save_path=save_plot_path)

    if "manh" in args.plot:
        for single_result in all_results:
            manhattan_plot(single_result, method_colors=gwas_method_colors, figsize=figure_size,
                           show_plot=args.showplot, save_path=save_plot_path,
                           significance_level_log=significance_level_log,
                           adjusted_threshold_log=adjusted_threshold_log)

    if "miami" in args.plot:
        miami_plot_all(all_results, method_colors=gwas_method_colors, figsize=figure_size,
                       show_plot=args.showplot, save_path=save_plot_path,
                       significance_level_log=significance_level_log,
                       adjusted_threshold_log=adjusted_threshold_log)

    if "qq" in args.plot:
        # qqplot for multiple result data
        qqplot(all_results, method_colors=gwas_method_colors, figsize=figure_size, save_path=save_plot_path)

    if "box" in args.plot and args.only_common:
        # qqplot for multiple result data
        top_snp_boxplots(snp_genotype_volumes, figsize=figure_size, save_path=save_plot_path, number_of_top_snps=3)

    if "venn" in args.plot:
        # qqplot for multiple result data
        venn_diagram(all_results, significance_threshold=venn_threshold_log,
                     method_colors=gwas_method_colors, figsize=figure_size, save_path=save_plot_path)
