"""
# —————————————————————————————————
# encoding: utf-8
# tissue_integration.py
# Seunghyun Wang
# Last Modified: 2026.02
# —————————————————————————————————
"""

from all_common_functions import *
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

#====================================================================
def multi_platform_integration(dataframe_name1, dataframe_name2):

    df_callset1 = pd.read_csv("./data/call-set-processed/" + dataframe_name1 + ".tsv", sep="\t")
    df_callset2 = pd.read_csv("./data/call-set-processed/" + dataframe_name2 + ".tsv", sep="\t")

    # Filter germline calls
    df_callset1 = df_callset1[df_callset1["POLYMORPHIC"] == "No"]

    df_callset2 = df_callset2[~df_callset2["PHASING"].isin(["germline", "insertion_in_multiple_hap"])]
    df_callset2 = df_callset2[df_callset2["POLYMORPHIC"] == "No"]

    # S1
    df_callset1_s1 = df_callset1[df_callset1["CROSS-CALL"]=="Yes"]
    df_callset1_s1["SUPPORT"] = "S1"
    df_callset1_s1 = df_callset1_s1[["#CHROM", "POS", "SUPPORT"]]
    num_s1 = df_callset1_s1.shape[0]
    print("S1: ", num_s1)
    print(df_callset1_s1[["#CHROM", "POS"]])

    #S2
    df_callset1_s2 = df_callset1[df_callset1["CROSS-CALL"] == "No"]
    df_callset1_s2 = df_callset1_s2[df_callset1_s2["CROSS-RAW"] == "Yes"]
    df_callset1_s2["SUPPORT"] = "S2"
    df_callset1_s2 = df_callset1_s2[["#CHROM", "POS", "SUPPORT"]]
    num_s2 = df_callset1_s2.shape[0]
    print("S2: ", num_s2)
    print(df_callset1_s2[["#CHROM", "POS"]])

    # S3'
    df_callset2_s3 = df_callset2[df_callset2["CROSS-CALL"] == "No"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["CROSS-RAW"] == "Yes"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["SUBFAMILY"]=="Young"]
    df_callset2_s3 = df_callset2_s3[df_callset2_s3["PRIMARY"]!="NoPrimary"]
    df_callset2_s3["SUPPORT"] = "S3'"
    df_callset2_s3 = df_callset2_s3[["#CHROM", "POS", "SUPPORT"]]
    num_s3 = df_callset2_s3.shape[0]
    print("S3': ", num_s3)
    print(df_callset2_s3[["#CHROM", "POS"]])

    # S4'
    df_callset1_s4 = df_callset1[df_callset1["CROSS-CALL"] == "No"]
    df_callset1_s4 = df_callset1_s4[df_callset1_s4["CROSS-RAW"] == "No"]
    df_callset1_s4 = df_callset1_s4[df_callset1_s4["NON-UNIQUE"]=="No"]
    df_callset1_s4["SUPPORT"] = "S4'"
    df_callset1_s4 = df_callset1_s4[["#CHROM", "POS", "SUPPORT"]]
    num_s4 = df_callset1_s4.shape[0]
    print("S4': ", num_s4)
    print(df_callset1_s4[["#CHROM", "POS"]])

    # S4'
    df_callset2_s5 = df_callset2[df_callset2["CROSS-CALL"] == "No"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["CROSS-RAW"] == "No"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["SUBFAMILY"] == "Young"]
    df_callset2_s5 = df_callset2_s5[df_callset2_s5["PRIMARY"] != "NoPrimary"]
    df_callset2_s5["SUPPORT"] = "S5'"
    df_callset2_s5 = df_callset2_s5[["#CHROM", "POS", "SUPPORT"]]
    num_s5 = df_callset2_s5.shape[0]
    print("S5': ", num_s5)
    print(df_callset2_s5[["#CHROM", "POS"]])

    df_whole = pd.concat([df_callset1_s1, df_callset1_s2, df_callset2_s3, df_callset1_s4, df_callset2_s5],
                         ignore_index=True)

    df_whole.to_csv("./data/call-set-processed/tissue_integrated.tsv", sep="\t", index=False)

#====================================================================
if __name__ == '__main__':

    # Running example 
    # Please change dataframe and bam name 
    
    # Preprocessing short-read call set 
    short = "xTea_mosaic-SMaHT-sMEI-LINE1-ST003-1Q-NYGC-Illumina-200x"
    short_bam = "ST003-1Q-NYGC-Illumina-200x"
    convert_vcf_to_dataframe(short)
    remove_duplicated_calls(short)
    add_genomic_region(short)
    find_polymorphic(short)

    # Preprocessing long-read call set 
    long = "PALMER-SMaHT-sMEI-LINE1-ST003-1Q-UW-PacBio-49x"
    long_bam = "ST003-1Q-UW-PacBio-49x"
    convert_vcf_to_dataframe(long)
    remove_duplicated_calls(long)
    add_insertion_sequence(long)
    run_blast(long)
    find_primary(long)
    run_haplotypcaller(long, long_bam)
    do_phasing(long, long_bam, 5)
    find_polymorphic(long)

    # Multi-platform integration
    find_common_calls(short, long)
    find_short_read_signal(long, short_bam)
    find_long_read_signal(short, short_bam, long_bam)
    multi_platform_integration(short, long)
