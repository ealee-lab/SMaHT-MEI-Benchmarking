"""
# —————————————————————————————————
# encoding: utf-8
# all_common_functions.py
# Seunghyun Wang, Mingyun Bae
# Last Modified: 2026.02
# —————————————————————————————————
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
import vcfpy
import pyranges as pr
import pysam
from collections import Counter
import subprocess
import regex
import ast
import os

#====================================================================
def convert_vcf_to_dataframe(vcf_name):

    # Function purpose: To convert vcf format to dataframe (.tsv)
    # Input: Caller VCF
    # Output: Call set dataframe (column: #CHROM, POS)

    reader = vcfpy.Reader.from_path("./data/call-set/" + vcf_name + ".vcf")
    ## Convert VCF to DataFrame
    rows = []
    for record in reader:
        rows.append({"#CHROM": record.CHROM, "POS": record.POS, "INFO": record.INFO})
    df_callset = pd.DataFrame(rows)

    df_callset.to_csv("./data/call-set-processed/" + vcf_name+ ".tsv", sep="\t", index=False)


#====================================================================
def remove_duplicated_calls(dataframe_name):

    # Function purpose: To remove duplicated calls
    # Input: Call set dataframe
    # Output: Call set dataframe after removing duplicated calls

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    print("The number of whole calls in call set before removing duplicated calls: ", df_callset.shape[0])

    keep_index = []
    for idx, row in df_callset.iterrows():
        chrom = row["#CHROM"]
        pos = row["POS"]

        call_position_call_set = set(list(df_callset[df_callset["#CHROM"] == chrom]["POS"]))
        call_position_range = set([pos for pos in range(pos - 100, pos + 101)])
        intersection = call_position_call_set.intersection(call_position_range)

        df_tmp = df_callset[(df_callset["#CHROM"] == chrom) & (df_callset["POS"].isin(list(intersection)))]
        all_idx = list(df_tmp.index)

        if len(all_idx) == 1:
            keep_index.append(all_idx[0])
        elif len(all_idx) > 1 and all_idx[0] not in keep_index:
            keep_index.append(all_idx[0])

    df_callset = df_callset.iloc[keep_index]
    df_callset = df_callset.sort_values(by=["#CHROM", "POS"]).reset_index(drop=True)

    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)
    print("The number of whole calls in call set after removing duplicated calls: ", df_callset.shape[0])

# ====================================================================
def add_genomic_region(dataframe_name):

    # Function purpose: To add "NON-UNIQUE" (values: 'True' or 'False') column in call set dataframe
    # Input: Call set dataframe
    # Output: Call set dataframe with "NON-UNIQUE" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset_copy = df_callset.copy()

    df_callset_copy["Start"] = df_callset_copy["POS"].astype(int)
    df_callset_copy["End"] = df_callset_copy["POS"]

    df_callset_copy = df_callset_copy.rename(columns={"#CHROM": "Chromosome"})
    candidates_ranges = pr.PyRanges(df_callset_copy[["Chromosome", "Start", "End"]].assign(_idx=df_callset_copy.index))

    # Please download genomic region SMaHT_extreme_v2.bed from SMaHT data portal
    # Extreme --> non-unique
    reference_ranges = pr.read_bed("./data/benchmarking-set/SMaHT_extreme_v2.bed")
    hit = candidates_ranges.join(reference_ranges)

    hit_idx = set(hit.df["_idx"])
    df_callset["NON-UNIQUE"] = df_callset.index.map(lambda i: i in hit_idx)

    converted = []
    for value in list(df_callset["NON-UNIQUE"]):
        if value==False:
            converted.append("No")
        else:
            converted.append("Yes")

    df_callset["NON-UNIQUE"] = converted
    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)

# ====================================================================
def add_insertion_sequence(dataframe_name):

    # Function purpose: To add insertion sequence for PALMER call set
    # Input: Call set dataframe (PALMER)
    # Output: Call set dataframe with "SEQ" column (PALMER)

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset["READ"] = ""
    df_callset["SEQ"] = ""

    raw_calls = "./data/call-set/" + dataframe_name + ".TSD_reads.txt"
    with open(raw_calls) as indata:
        for i, l in enumerate(tqdm(indata)):
            ln = l.rstrip().split("\t")
            chrom = ln[0].split("_")[1]
            pos = int(ln[0].split("_")[2])
            if str(ln[2]) == "N":
                tsd5 = ""
            else:
                tsd5 = str(ln[2])
            seq = str(ln[-1])

            read_name_info = ln[1].split(".")
            read_name_info = read_name_info[0].split("_")
            read_name_info = read_name_info[:-2]
            read_name = "_".join(read_name_info)

            pos_candidates = set([i for i in range(pos-100, pos+100)])
            callset_positions  = set(df_callset[df_callset["#CHROM"] == chrom]["POS"])

            intersection = pos_candidates.intersection(callset_positions)
            intersection = list(intersection)

            if len(intersection) == 1:

                call_idx = df_callset[(df_callset["#CHROM"] == chrom) & (df_callset["POS"] == intersection[0])].index[0]
                if df_callset.loc[call_idx, "READ"] == "":
                    df_callset.loc[call_idx, "READ"] = read_name
                    df_callset.loc[call_idx, "SEQ"] = tsd5+seq
                else:
                    df_callset.loc[call_idx, "READ"] = df_callset.loc[call_idx, "READ"] + ";" +read_name
                    df_callset.loc[call_idx, "SEQ"] = df_callset.loc[call_idx, "SEQ"] + ";" +tsd5+seq

    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)


# ====================================================================
def run_blast(dataframe_name):

    # Function purpose: To annotate repeat element subfamily information
    # Input: Call set dataframe
    # Output: Call set dataframe with "SUBFAMILY" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    blast_result = []

    for idx, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):

        try:
            seqs = row["SEQ"].split(";")
            subtypes = []
            for seq in seqs:
                out = open("./data/blastndb/seq.fa", 'w')
                out.write(">SEQ"+"\n")
                out.write(seq + "\n")
                out.close()
                command = "blastn -db ./data/blastndb/hg38_repeat.fa -query ./data/blastndb/seq.fa -outfmt 6"
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
                if len(result.stdout) == 0:
                    subtypes.append(".")

                else:
                    df_tmp = pd.DataFrame([row.split("\t") for row in result.stdout.split("\n")])
                    subtypes.append(df_tmp.loc[0][1])

            most_common = Counter(subtypes).most_common(1)[0][0]
            blast_result.append(most_common)
        except:
            blast_result.append(".")

    df_callset["SUBFAMILY"] = blast_result
    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)


# ====================================================================
def find_primary(dataframe_name):

    # Function purpose: To annotate read alignment information
    # Input: Call set dataframe
    # Output: Call set dataframe with "PRIMARY" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    bamfile = pysam.AlignmentFile("./data/bam/HapMapMix-PacBio-60x.bam", "r")
    df_callset["PRIMARY"] = "NoPrimary"

    for idx, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):
        try:
            chrom = row["#CHROM"]
            pos = row["POS"]
            reads = row["READ"].split(";")

            primaries = []
            for read in bamfile.fetch(chrom, pos, pos+1):
                if read.query_name in reads:
                    if read.is_supplementary == False:
                        primaries.append(read.query_name)

            if len(primaries) > 0:
                df_callset.loc[idx, "PRIMARY"] = ";".join(primaries)

        except:
            pass

    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)


# ====================================================================
def find_common_calls(dataframe_name1, dataframe_name2):

    # Function purpose: To find intersected calls between two call set (short-read and long-read call set)
    # Input: Each of call set dataframe
    # Output: Each of call set dataframe with "CROSS-CALL" column

    df_callset1 = pd.read_csv("./data/call-set-processed/" + dataframe_name1 + ".tsv", sep="\t")
    df_callset2 = pd.read_csv("./data/call-set-processed/" + dataframe_name2 + ".tsv", sep="\t")

    df_callset1["CROSS-CALL"] = "No"
    df_callset2["CROSS-CALL"] = "No"

    for idx, row in tqdm(df_callset1.iterrows(), total=df_callset1.shape[0]):

        chrom = row["#CHROM"]
        pos = row["POS"]

        pos_candidates = set([i for i in range(pos-100, pos+101)])
        callset2_positions = set(list(df_callset2[df_callset2["#CHROM"] == chrom]["POS"]))

        intersection = pos_candidates.intersection(callset2_positions)
        if len(intersection) >0:
            df_callset1.loc[idx, "CROSS-CALL"] = "Yes"

    for idx, row in tqdm(df_callset2.iterrows(), total=df_callset2.shape[0]):

        chrom = row["#CHROM"]
        pos = row["POS"]

        pos_candidates = set([i for i in range(pos - 100, pos + 101)])
        callset1_positions = set(list(df_callset1[df_callset1["#CHROM"] == chrom]["POS"]))

        intersection = pos_candidates.intersection(callset1_positions)
        if len(intersection) > 0:
            df_callset2.loc[idx, "CROSS-CALL"] = "Yes"

    df_callset1.to_csv("./data/call-set-processed/" + dataframe_name1 + ".tsv", sep="\t", index=False)
    df_callset2.to_csv("./data/call-set-processed/" + dataframe_name2 + ".tsv", sep="\t", index=False)


# ====================================================================
def find_short_read_signal(dataframe_name, bam_name):

    # Function purpose: To find short-read raw-alignment signal at each long-read call set
    # Input: Each of call set dataframe
    # Output: Each of call set dataframe with "CROSS-RAW" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset["CROSS-RAW"] = "No"

    bamfile = pysam.AlignmentFile("./data/bam/"+bam_name+".bam", "r")

    for idx, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):
        cross_call = row["CROSS-CALL"]

        if cross_call == "Yes":
            df_callset.loc[idx, "CROSS-RAW"] = "Yes"

        else:
            try:
                chrom = row["#CHROM"]
                pos = row["POS"]
                insertions = row["SEQ"].split(";")
                info = ast.literal_eval(row["INFO"])
                tsd = int(info["TSD_LEN"][0])

                matches = []
                for insertion in insertions:
                    for read in bamfile.fetch(chrom, pos - 25, pos + 25):
                        if read.mapping_quality > 0:

                            ref_start = read.reference_start
                            ref_end = read.reference_end
                            cigar_dict = dict(read.cigartuples)
                            cigar_dict_keys = list(cigar_dict.keys())

                            if cigar_dict_keys[0] == 4:
                                if abs(ref_start - pos) < 25:
                                    length = cigar_dict[4]
                                    if length >= 7:
                                        insertion_start = insertion[tsd:tsd + length]
                                        insertion_end = insertion[-length:]
                                        seq = read.query_sequence[0:length]

                                        m = regex.findall("(" + seq + "){e<=2}", insertion_start)
                                        matches = matches + m
                                        m = regex.findall("(" + seq + "){e<=2}", insertion_end)
                                        matches = matches + m

                            elif cigar_dict_keys[-1] == 4:
                                if abs(ref_end - pos) < 25:
                                    length = cigar_dict[4]
                                    if length >= 7:
                                        insertion_start = insertion[tsd:tsd + length]
                                        insertion_end = insertion[-length:]
                                        seq = read.query_sequence[-length:]
                                        m = regex.findall("(" + seq + "){e<=2}", insertion_start)
                                        matches = matches + m
                                        m = regex.findall("(" + seq + "){e<=2}", insertion_end)
                                        matches = matches + m

                matches_num = [len(match) for match in matches]
                matches_num = np.array(matches_num)
                matches_num = matches_num >= 7
                matches_num = matches_num * 1
                matches_num = np.sum(matches_num)

                if len(matches) > 0 and matches_num > 0:
                    df_callset.loc[idx, "CROSS-RAW"] = "Yes"

            except:
                pass

    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)

# ====================================================================
def get_long_insertions(bamfile, chrom, pos, min_len=50):

    # Function purpose: To get insertion seqeunce (raw read-alignment signal) at short-read based WGS call set
    # Input: long-read WGS BAM / chromosome and position from short-read based WGS call set
    # Output: List including insertion sequences

    insertion_seqs = []

    for read in bamfile.fetch(chrom, pos-25, pos + 25):
        if read.is_unmapped or read.cigartuples is None:
            continue

        read_pos = 0
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                if ref_pos <= pos < ref_pos + length:
                    break
                ref_pos += length
                read_pos += length
            elif op == 1:
                if abs(ref_pos-pos)<5 and length >=min_len:
                    seq = read.query_sequence[read_pos:read_pos + length]
                    insertion_seqs.append(seq)
                read_pos += length
            elif op in (2, 3):
                ref_pos += length
            elif op in (4, 5):
                read_pos += length
    return insertion_seqs


# ====================================================================
def get_clipped_sequences(samfile, chrom, pos, min_clip_len=50):

    # Function purpose: To get clipped sequence (raw read-alignment signal) at short-read based WGS call set
    # Input: long-read WGS BAM / chromosome and position from short-read based WGS call set
    # Output: List including clipped sequences

    clipped_seqs = []

    for read in samfile.fetch(chrom, pos-25, pos + 25):
        if read.is_unmapped or read.cigartuples is None:
            continue

        cigars = read.cigartuples
        seq = read.query_sequence
        ref_start = read.reference_start
        ref_end = read.reference_end
        # Left (5') soft clip
        if cigars[0][0] == 4 and cigars[0][1] >= min_clip_len and abs(ref_start-pos)<25:
            clip_len = cigars[0][1]
            clipped_seq = seq[:clip_len]
            clipped_seqs.append(clipped_seq)

        # Right (3') soft clip
        if cigars[-1][0] == 4 and cigars[-1][1] >= min_clip_len and abs(ref_end-pos)<25:

            clip_len = cigars[-1][1]
            clipped_seq = seq[-clip_len:]
            clipped_seqs.append(clipped_seq)

    return clipped_seqs



# ====================================================================
def find_raw_signal_per_insertion_sequence(samfile, chrom, pos, insertion):

    # Function purpose: To find matched long-read raw-aligmented signal with short-read call
    # Input: Short-read BAM and insertion or clipped sequence from long-read
    # Output: Matched sequence candidates

    matches = []
    for read in samfile.fetch(chrom, pos - 25, pos + 25):
        if read.mapping_quality > 0:
            try:
                cigar_dict = dict(read.cigartuples)
                cigar_dict_keys = list(cigar_dict.keys())

                if cigar_dict_keys[0] == 4 or cigar_dict_keys[-1] == 4:
                    length = cigar_dict[4]
                if cigar_dict_keys[0] == 4:

                    if length > 7:
                        seq = read.query_sequence[0:length]
                        m = regex.findall("(" + seq + "){e<=2}", insertion)
                        matches = matches + m

                elif cigar_dict_keys[-1] == 4:

                    if length > 7:
                        seq = read.query_sequence[-length:]
                        m = regex.findall("(" + seq + "){e<=2}", insertion)
                        matches = matches + m
            except:
                pass
    return matches

# ====================================================================
def find_long_read_signal(dataframe_name, bam_name1, bam_name2):

    # Function purpose: To find long-read raw-alignment signal at each short-read call set
    # Input: Each of call set dataframe
    # Output: Each of call set dataframe with "CROSS-RAW" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    df_callset["CROSS-RAW"] = "No"

    bamfile1 = pysam.AlignmentFile("./data/bam/"+bam_name1+".bam", "r")
    bamfile2 = pysam.AlignmentFile("./data/bam/"+bam_name2+".bam", "r")

    for idx, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):
        cross_call = row["CROSS-CALL"]

        if cross_call == "Yes":
            df_callset.loc[idx, "CROSS-RAW"] = "Yes"

        else:
            try:
                chrom = row["#CHROM"]
                pos = row["POS"]
                matches_final = []

                insertions = get_long_insertions(bamfile2, chrom, pos, min_len=50)
                insertions = list(set(insertions))

                if len(insertions) >0:
                    matches_final = []
                    for insertion in insertions:
                        matches = find_raw_signal_per_insertion_sequence(bamfile1, chrom, pos, insertion)
                        matches_final = matches_final + matches

                matches_num = [len(match) for match in matches_final]
                matches_num = np.array(matches_num)
                matches_num = matches_num >= 7
                matches_num = matches_num * 1
                matches_num = np.sum(matches_num)
                if len(matches) > 0 and matches_num > 0:
                    df_callset.loc[idx, "CROSS-RAW"] = "Yes"

                else:
                    clipped_seqs = get_clipped_sequences(bamfile2, chrom, pos, min_clip_len=50)
                    clipped_seqs = list(set(clipped_seqs))

                    if len(clipped_seqs) >0:

                        matches_final = []
                        for insertion in clipped_seqs:
                            matches = find_raw_signal_per_insertion_sequence(bamfile1, chrom, pos, insertion)
                            matches_final = matches_final + matches
                        matches_num = [len(match) for match in matches_final]
                        matches_num = np.array(matches_num)
                        matches_num = matches_num >= 7
                        matches_num = matches_num * 1
                        matches_num = np.sum(matches_num)
                        if len(matches) > 0 and matches_num > 0:
                            df_callset.loc[idx, "CROSS-RAW"] = "Yes"

            except:
                pass
    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)


# ====================================================================
def run_haplotypcaller(dataframe_name, bam_name):

    # Function purpose: To find long-read raw-alignment signal at each short-read call set
    # Input: Each of call set dataframe
    # Output: Each of call set dataframe with "CROSS-RAW" column

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    bamfile = "./data/bam/" + bam_name + ".bam"
    ref = "./data/ref/human_GRCh38_no_alt_analysis_set.fasta"
    snp = "./data/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

    df_callset["start"] = df_callset["POS"]-20000
    df_callset["end"] = df_callset["POS"]+20000
    df_callset = df_callset[["#CHROM","start", "end"]]
    df_callset.to_csv("./data/call-set-processed/candidates_bed/"+dataframe_name+"_candidates.bed", sep="\t", index=False, header=False)

    os.system("gatk HaplotypeCaller -R %s -I %s -O %s -L %s --dbsnp %s"%(ref, bamfile, "./data/call-set-processed/candidates_vcf/"+dataframe_name+"_candidates.vcf", "./data/call-set-processed/candidates_bed/"+dataframe_name+"_candidates.bed", snp))

    os.system("bgzip -c "+"./data/call-set-processed/candidates_vcf/"+dataframe_name+"_candidates.vcf"+" >"+"./data/call-set-processed/candidates_vcf/"+dataframe_name+"_candidates.vcf.gz")
    os.system("tabix -p vcf "+"./data/call-set-processed/candidates_vcf/"+dataframe_name+"_candidates.vcf.gz")


# ====================================================================
def get_closest_variants(vcf, chrom, position, n):
    variants = []
    variants_pos = []
    variants_neg = []
    for i in vcf.fetch(chrom):
        if i.qual > 1000:
            dist_abs = abs(i.pos - position)
            dist_real = i.pos - position
            if len(i.ref) == 1 and len(i.alts[0]) == 1 and i.info["AF"][0]==0.5:
                variants.append((dist_abs, i))
                if dist_real > 0:
                    variants_pos.append((dist_abs, i))
                else:
                    variants_neg.append((dist_abs, i))
    closest = sorted(variants, key=lambda x: x[0])[:n]
    return closest

# ====================================================================
def get_snp_read(bam, chrom, position):
    info_base = {}
    info = {}
    for read in bam.fetch(chrom, position-1, position):
        if read.mapq >= 30:
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == position-1:
                    base = read.query_sequence[qpos]
                    info_base[read.query_name] = read
                    info[read.query_name] = base
    return info_base, info

# ====================================================================
def get_insertion_pos_read(bam, chrom, position):
    info = {}
    for read in bam.fetch(chrom, position-1, position):
        if read.mapq >= 30:
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == position:
                    base = read.query_sequence[qpos]
                    info[read.query_name] = base

    return info

# ====================================================================
def get_hap(bam, vcf, chrom, pos, n):
    snps = get_closest_variants(vcf, chrom, pos, n)
    if len(snps) == 0:
        return ["no_hSNPs", "no_hSNPs"]
    dic_snp = {}
    dic_snp_base = {}
    for i in range(len(snps)):
        dic_snp_base[i], dic_snp[i] = get_snp_read(bam, chrom, snps[i][1].pos)
    insertion_pos_reads = get_insertion_pos_read(bam, chrom, pos)
    insertion = get_long_insertions(bam, chrom, pos)
    phasing = []
    result = []
    for i in range(len(snps)):
        common_reads = list(set(dic_snp[i].keys()) & set(insertion_pos_reads))
        df_tmp = pd.DataFrame.from_dict(dic_snp[i], orient="index")
        df_tmp = df_tmp.loc[common_reads]
        if df_tmp.shape[0] == 0:
            result.append("cannot_phasing")
            phasing.append("cannot_phasing")
            continue
        df_tmp["insertion"] = df_tmp.index.isin(insertion).astype(int)
        if 1 not in df_tmp["insertion"].values:
            result.append("insertion_read_not_cover")
            phasing.append("insertion_read_not_cover")
            continue
        insertion_snp = df_tmp.loc[df_tmp["insertion"]==1][0].value_counts().idxmax()
        num_no_insertion_snp = df_tmp.loc[(df_tmp[0]==insertion_snp) & (df_tmp["insertion"]==0)].shape[0]
        if len(set(df_tmp.loc[df_tmp["insertion"]==1][0].values)) > 1:
            result.append("insertion_in_multiple_hap")
            phasing.append("insertion_in_multiple_hap")
        else:
            if num_no_insertion_snp > 0:
                result.append("sMEI")
                phasing.append("sMEI")
            else:
                phasing.append("germline")
                result.append("germline")

    result = "/".join(result)
    phasing = "/".join(phasing)
    return [result, phasing]

# ====================================================================
def do_phasing(dataframe_name, bam_name, n):

    vcf = pysam.VariantFile("./data/call-set-processed/candidates_vcf/" + dataframe_name + "_candidates.vcf.gz")
    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    bamfile = pysam.AlignmentFile("./data/bam/"+bam_name+".bam", "r")

    results = []
    for i, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):
        chrom = row["#CHROM"]
        pos = row["POS"]
        r, combi = get_hap(bamfile, vcf, chrom, pos, n)
        rr = r.split("/")
        rr = [val for val in rr if val != "cannot_phasing"]
        if len(rr) == 0:
            results.append(["cannot_phasing"])
        else:
            counter = Counter(rr)
            result = counter.most_common(1)[0][0]
            results.append([result])

    df_tmp = pd.DataFrame(results, columns=["PHASING"])
    df_combined = pd.concat([df_callset, df_tmp], axis=1)

    df_combined.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)

# ====================================================================
def find_polymorphic(dataframe_name):

    df_callset = pd.read_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t")
    
    #Please use your polymorphic MEI collections 
    df_polymorphic = pd.read_csv("./data/polymorphic/LINE1-polymorphic.tsv", sep="\t")
    df_callset["POLYMORPHIC"] = "No"

    for idx, row in tqdm(df_callset.iterrows(), total=df_callset.shape[0]):
        chrom = row["#CHROM"]
        pos = row["POS"]
      
        polymorphic_sites = list(df_polymorphic[df_polymorphic["#CHROM"] == chrom]["POS"])
        polymorphic_sites = [int(pos) for pos in polymorphic_sites]

        pos_candidates = [i for i in range(pos-100, pos+100)]
        intersection = set(polymorphic_sites).intersection(pos_candidates)

        if len(intersection) > 0:
            df_callset.loc[idx, "POLYMORPHIC"] = "Yes"

    df_callset.to_csv("./data/call-set-processed/" + dataframe_name + ".tsv", sep="\t", index=False)
