"""
# —————————————————————————————————
# encoding: utf-8
# hapmapmix_benchmarking_and_integration.py
# Seunghyun Wang
# Last Modified: 2026.02
# —————————————————————————————————
"""

from all_common_functions import *
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

#====================================================================
def hapmapmix_add_germline(dataframe_name):

    # Function purpose: To exclude germline call in HapMap mixture call set before somatic detection performance evaluation
    # Input: Call set dataframe
    # Output: Call set dataframe with "GERMLINE" column

    df_benchmark_germ_tier3 = pd.read_csv("./data/benchmarking-set/benchmarking-set-SMaHT-gMEI-LINE1-tier3-HapMapMix.tsv", sep="\t")
    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")

    germlines = []
    for idx, row in df_callset.iterrows():
        chrom = row["#CHROM"]
        pos = row["POS"]
        germline_position_benchmarking_set = set(df_benchmark_germ_tier3[df_benchmark_germ_tier3["#CHROM"] == chrom]["POS"])
        call_position_range = set([pos for pos in range(pos-100, pos+101)])
        intersection = germline_position_benchmarking_set.intersection(call_position_range)
        if len(intersection) == 0:
            germlines.append("No")
        else:
            germlines.append("Yes")

    df_callset["GERMLINE"] = germlines
    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)

#====================================================================
def hapmapmix_add_true_positive(df_callset):

    # Function purpose: To label true positives 
    # Input: Call set dataframe
    # Output: Call set dataframe with "SOMA(tier)" column

    df_benchmark_soma_tier2 = pd.read_csv("./data/benchmarking-set/benchmarking-set-SMaHT-sMEI-LINE1-tier2-HapMapMix.tsv", sep="\t")
    df_benchmark_soma_tier3 = pd.read_csv("./data/benchmarking-set/benchmarking-set-SMaHT-sMEI-LINE1-tier3-HapMapMix.tsv", sep="\t")
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
def hapmapmix_evaluate_detection_performance(dataframe_name):

    # Function purpose: To evaluate F1, recall, precision of somatic MEI detection method in HapMap mixture
    # Input: Call set dataframe 
    # Output: F1, recall, precision score

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset = df_callset[df_callset["GERMLINE"]=="No"]
    df_callset, num_all_soma_tier2 = hapmapmix_add_true_positive(df_callset)

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
def hapmapmix_multi_platform_integration(dataframe_name1, dataframe_name2):

    # Function purpose: To integrate short-read and long-read based call set 
    # Input: Short-read nad long-read call set dataframe 
    # Output: Integrated call set 

    df_callset1 = pd.read_csv("./data/call-set-processed/" + dataframe_name1 + ".tsv", sep="\t")
    df_callset2 = pd.read_csv("./data/call-set-processed/" + dataframe_name2 + ".tsv", sep="\t")

    df_callset1 = df_callset1[df_callset1["GERMLINE"] == "No"]
    df_callset2 = df_callset2[df_callset2["GERMLINE"] == "No"]

    # S1
    df_callset1_s1 = df_callset1[df_callset1["CROSS-CALL"]=="Yes"]
    df_callset1_s1, _ = hapmapmix_add_true_positive(df_callset1_s1)
    df_callset1_s1_tp = df_callset1_s1[df_callset1_s1["SOMA(tier2)"]!=-1]
    df_callset1_s1_fp = df_callset1_s1[df_callset1_s1["SOMA(tier2)"]==-1]
    tp_s1 = df_callset1_s1_tp.shape[0]
    fp_s1 = df_callset1_s1_fp.shape[0]
    print("S1: ", "TP ", tp_s1, "FP ", fp_s1)
    df_callset1_s1["SUPPORT"] = "S1"
    df_callset1_s1 = df_callset1_s1[["#CHROM", "POS", "SUPPORT", "SOMA(tier2)"]]

    #S2
    df_callset1_s2 = df_callset1[df_callset1["CROSS-CALL"] == "No"]
    df_callset1_s2 = df_callset1_s2[df_callset1_s2["CROSS-RAW"] == "Yes"]
    df_callset1_s2, _ = hapmapmix_add_true_positive(df_callset1_s2)
    df_callset1_s2_tp = df_callset1_s2[df_callset1_s2["SOMA(tier2)"] != -1]
    df_callset1_s2_fp = df_callset1_s2[df_callset1_s2["SOMA(tier2)"] == -1]
    tp_s2 = df_callset1_s2_tp.shape[0]
    fp_s2 = df_callset1_s2_fp.shape[0]
    print("S2: ", "TP ", tp_s2, "FP ", fp_s2)
    df_callset1_s2["SUPPORT"] = "S2"
    df_callset1_s2 = df_callset1_s2[["#CHROM", "POS", "SUPPORT", "SOMA(tier2)"]]

    # S3'
    df_callset2_s3 = df_callset2[df_callset2["CROSS-CALL"] == "No"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["CROSS-RAW"] == "Yes"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["SUBFAMILY"]=="Young"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["PRIMARY"]!="NoPrimary"]
    df_callset2_s3, _ = hapmapmix_add_true_positive(df_callset2_s3)
    df_callset2_s3_tp = df_callset2_s3[df_callset2_s3["SOMA(tier2)"] != -1]
    df_callset2_s3_fp = df_callset2_s3[df_callset2_s3["SOMA(tier2)"] == -1]
    tp_s3 = df_callset2_s3_tp.shape[0]
    fp_s3 = df_callset2_s3_fp.shape[0]
    print("S3: ", "TP ", tp_s3, "FP ", fp_s3)
    df_callset2_s3["SUPPORT"] = "S3'"
    df_callset2_s3 = df_callset2_s3[["#CHROM", "POS", "SUPPORT", "SOMA(tier2)"]]

    # S4'
    df_callset1_s4 = df_callset1[df_callset1["CROSS-CALL"] == "No"]
    df_callset1_s4 = df_callset1_s4[df_callset1_s4["CROSS-RAW"] == "No"]
    df_callset1_s4 = df_callset1_s4[df_callset1_s4["NON-UNIQUE"]=="No"]
    df_callset1_s4, _ = hapmapmix_add_true_positive(df_callset1_s4)
    df_callset1_s4_tp = df_callset1_s4[df_callset1_s4["SOMA(tier2)"] != -1]
    df_callset1_s4_fp = df_callset1_s4[df_callset1_s4["SOMA(tier2)"] == -1]
    tp_s4 = df_callset1_s4_tp.shape[0]
    fp_s4 = df_callset1_s4_fp.shape[0]
    print("S4: ", "TP ", tp_s4, "FP ", fp_s4)
    df_callset1_s4["SUPPORT"] = "S4'"
    df_callset1_s4 = df_callset1_s4[["#CHROM", "POS", "SUPPORT", "SOMA(tier2)"]]

    # S5'
    df_callset2_s5 = df_callset2[df_callset2["CROSS-CALL"] == "No"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["CROSS-RAW"] == "No"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["SUBFAMILY"] == "Young"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["PRIMARY"] != "NoPrimary"]
    df_callset2_s5, num_all = hapmapmix_add_true_positive(df_callset2_s5)
    df_callset2_s5_tp = df_callset2_s5[df_callset2_s5["SOMA(tier2)"] != -1]
    df_callset2_s5_fp = df_callset2_s5[df_callset2_s5["SOMA(tier2)"] == -1]
    tp_s5 = df_callset2_s5_tp.shape[0]
    fp_s5 = df_callset2_s5_fp.shape[0]
    print("S5: ", "TP ", tp_s5, "FP ", fp_s5)
    df_callset2_s5["SUPPORT"] = "S5'"
    df_callset2_s5 = df_callset2_s5[["#CHROM", "POS", "SUPPORT", "SOMA(tier2)"]]

    tp = tp_s1 + tp_s2 + tp_s3 + tp_s4 + tp_s5
    fp = fp_s1 + fp_s2 + fp_s3 + fp_s4 + fp_s5

    recall = tp / num_all
    precision = tp / (tp+fp)
    f1 = 2 * (precision * recall) / (precision + recall)

    print("F1: ", round(f1, 2))
    print("Recall: ", round(recall, 2))
    print("Precision: ", round(precision, 2))

    df_whole = pd.concat([df_callset1_s1, df_callset1_s2, df_callset2_s3, df_callset1_s4, df_callset2_s5],
                         ignore_index=True)

    df_whole.to_csv("./data/call-set-processed/hapmapmix_integrated_after.tsv", sep="\t", index=False)

#====================================================================
if __name__ == '__main__':

    # Preprocess call set
    short = "xTea_mosaic-SMaHT-sMEI-LINE1-HapMapMix-Illumina-200x"
    long = "PALMER-SMaHT-sMEI-LINE1-HapMapMix-PacBio-60x"

    convert_vcf_to_dataframe(short)
    remove_duplicated_calls(short)
    hapmapmix_add_germline(short)
    add_genomic_region(short)

    convert_vcf_to_dataframe(long)
    remove_duplicated_calls(long)
    hapmapmix_add_germline(long)
    add_insertion_sequence_tmp(long)

    # Evaluate the detection performance of individual call set
    hapmapmix_evaluate_detection_performance(short)
    hapmapmix_evaluate_detection_performance(long)

    #Add information for platform-specific FP filtering
    add_genomic_region(short)
    add_insertion_sequence(long)
    run_blast(long)
    find_primary(long)

    # Multi-platform integration
    short_bam = "HapMapMix-Illumina-200x"
    long_bam = "HapMapMix-PacBio-60x"
    find_common_calls(short, long)
    find_short_read_signal(long, short_bam)
    find_long_read_signal(short, short_bam, long_bam)
    hapmapmix_multi_platform_integration(short, long)
