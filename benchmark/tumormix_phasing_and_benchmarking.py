"""
# —————————————————————————————————
# encoding: utf-8
# tumormix_phasing_and_benchmarking.py
# Mingyun Bae, Seunghyun Wang
# Last Modified: 2026.02
# —————————————————————————————————
"""

from all_common_functions import *
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")

#====================================================================
def tumormix_add_true_positive(df_callset):

    # Please note that this is the simplified version of SMaHT somatic MEI benchmarking set (LINE1, tier2) for tutorial.
    # The more detailed benchmarking sets are available in the supplementary table of original publication.

    df_benchmark_soma_tier2 = pd.read_csv("./data/benchmarking-set/benchmarking-set-SMaHT-sMEI-LINE1-tier2-TumorMix.tsv", sep="\t")
    df_benchmark_soma_tier3 = pd.read_csv("./data/benchmarking-set/benchmarking-set-SMaHT-sMEI-LINE1-tier3-TumorMix.tsv", sep="\t")
    num_all_soma_tier2 = df_benchmark_soma_tier2.shape[0]

    df_callset["SOMA(tier2)"] = -1
    df_callset["SOMA(tier3)"] = -1
    for idx, row in df_callset.iterrows():
        chrom = row["#CHROM"]
        pos = row["POS"]
        somatic_position_benchmarking_set_tier2 = set(df_benchmark_soma_tier2[df_benchmark_soma_tier2["#CHROM"] == chrom]["POS"])
        somatic_position_benchmarking_set_tier3 = set(df_benchmark_soma_tier3[df_benchmark_soma_tier3["#CHROM"] == chrom]["POS"])
        call_position_range = set([pos for pos in range(pos - 100, pos + 101)])
        intersection_tier2 = list(somatic_position_benchmarking_set_tier2.intersection(call_position_range))
        intersection_tier3 = list(somatic_position_benchmarking_set_tier3.intersection(call_position_range))

        if len(intersection_tier2) == 1:
            tp_position_tier2 = intersection_tier2[0]
            soma_id_tier2 = df_benchmark_soma_tier2[(df_benchmark_soma_tier2["#CHROM"] == chrom) & (df_benchmark_soma_tier2["POS"] == tp_position_tier2)]["ID"].values[0]
            df_callset.loc[idx, "SOMA(tier2)"] = soma_id_tier2
        elif len(intersection_tier2) > 1:
            diff_min = 300
            for position in intersection_tier2:
                diff = abs(position - pos)
                if diff < diff_min:
                    diff_min = diff
                    tp_position_tier2 = position
            soma_id_tier2 = df_benchmark_soma_tier2[(df_benchmark_soma_tier2["#CHROM"] == chrom) & (df_benchmark_soma_tier2["POS"] == tp_position_tier2)]["ID"].values[0]
            df_callset.loc[idx, "SOMA(tier2)"] = soma_id_tier2

        if len(intersection_tier3) == 1:
            tp_position_tier3 = intersection_tier3[0]
            soma_id_tier3 = df_benchmark_soma_tier3[(df_benchmark_soma_tier3["#CHROM"] == chrom) & (df_benchmark_soma_tier3["POS"] == tp_position_tier3)]["ID"].values[0]
            df_callset.loc[idx, "SOMA(tier3)"] = soma_id_tier3
        elif len(intersection_tier3) > 1:
            diff_min = 300
            for position in intersection_tier3:
                diff = abs(position - pos)
                if diff < diff_min:
                    diff_min = diff
                    tp_position_tier3 = position
            soma_id_tier3 = df_benchmark_soma_tier3[(df_benchmark_soma_tier3["#CHROM"] == chrom) & (df_benchmark_soma_tier3["POS"] == tp_position_tier3)]["ID"].values[0]
            df_callset.loc[idx, "SOMA(tier3)"] = soma_id_tier3

    return df_callset, num_all_soma_tier2

#====================================================================
def tumormix_evaluate_detection_performance(dataframe_name):

    # Function purpose: To evaluate F1, recall, precision of somatic MEI detection method in HapMap mixture
    # Input: 1) HapMap mixture benchmarking set / 2) Call set (Dataframe converted from VCF)
    # Output: F1, recall, precision score

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset = df_callset[~df_callset["PHASING"].isin(["germline", "insertion_in_multiple_hap"])]
    df_callset = df_callset[df_callset["POLYMORPHIC"]=="No"]

    df_callset, num_all_soma_tier2 = tumormix_add_true_positive(df_callset)

    tp = df_callset[df_callset["SOMA(tier2)"] != -1].shape[0]
    ## Except the potential true positive calls that existing lower-level benchmarking set in precision calculation
    undecisional = df_callset[(df_callset["SOMA(tier3)"]!= -1)&(df_callset["SOMA(tier2)"]== -1)].shape[0]

    fp = (df_callset.shape[0] - tp) - undecisional

    recall = tp / num_all_soma_tier2
    precision = tp / (tp + fp)
    f1 = (2 * precision * recall) / (precision + recall)

    print("F1: ", round(f1,2))
    print("Recall: ", round(recall,2))
    print("Precision: ", round(precision,2))

#====================================================================
if __name__ == '__main__':

    # Preprocess call set
    long = "PALMER-SMaHT-sMEI-LINE1-TumorMix-PacBio-60x"
    long_bam = "TumorMix-PacBio-60x"

    convert_vcf_to_dataframe(long)
    remove_duplicated_calls(long)

    # Identify germline calls
    run_haplotypcaller(long, long_bam)
    do_phasing(long, long_bam, 5)
    find_polymorphic(long)

    #Evaluate detection performance
    tumormix_evaluate_detection_performance(long)
