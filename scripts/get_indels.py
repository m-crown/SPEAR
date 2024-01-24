#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import re
from summarise_snpeff import parse_vcf, write_vcf
import pandas as pd
import os
import multiprocessing

#need to have a flag to exclude ambiguous indels in get_indels to go alongside excluding ambiguous snps from fatovcf
def mask_trimmed_sequence(sample):
    sample_seq = str(sample.seq)
    def repl(m):
        return 'N' * len(m.group())
    sample_seq = re.sub(r'^-+|-+$', repl, sample_seq)
    sample.seq = Seq(sample_seq)
    return sample

def get_indels(reference, sample, window, allow_ambiguous):
    vcf = [f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample.id}'] #set indels vcf header.
    offset = 0 #offset used when insertions detected, subsequent indels should be shifted by this offset
    indels = []
    #find insertions from reference sequence
    for match in re.finditer("-+", str(reference.seq)):
        ref = match.start() - 1 #the 0 based index position of the last ref base before indel
        length = len(match.group()) #length of the indel itself
        ins_end = ref + length + 1
        ref_base = reference.seq[ref]
        alt_base = sample.seq[ref:ins_end]
        indel = {"type": "insertion","pos": ref, "ref_base": ref_base, "alt_base": alt_base, "length": length}
        indels.append(indel)
    #find deletions from sample
    for match in re.finditer("-+", str(sample.seq)):
        ref = match.start() - 1
        length = len(match.group())
        ins_end = ref + length + 1
        ref_base = reference.seq[ref:ins_end]
        alt_base = sample.seq[ref]
        indel = {"type": "deletion","pos": ref, "ref_base": ref_base, "alt_base": alt_base, "length": length}
        indels.append(indel)
    indels = sorted(indels, key=lambda d: d['pos'])
    if allow_ambiguous:
        iupac_chars = ["N"]
    else:
        iupac_chars = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    for indel in indels:
        current_offset = offset
        if indel["type"] == "insertion": #filtering insertions that consist of only N characters , as these are interpreted as any nucletide downstream and are not reliable for annotation anyway. 
            offset += indel["length"] #still have to iterate the offset but dont mark the variant
            if any(char in indel["alt_base"][1:] for char in iupac_chars): 
                continue
            if indel["alt_base"][0] in iupac_chars: #if the reference base is N in the sample convert to ref in the alt description (ignore N snp)
                indel["alt_base"] = indel["ref_base"][0] + indel["alt_base"][1:]
        if indel["type"] == "deletion":
            if indel["alt_base"][0] in iupac_chars: #if the reference base is IUPAC ambiguous in the sample convert to ref in the alt description (ignore N snp)
                indel["alt_base"] = indel["ref_base"][0]

        if (window != 0 and ((sample.seq[indel["pos"] +1 - window: indel["pos"] + 1 ] == "N" * window) or (sample.seq[indel["pos"] +1 + indel["length"]: indel["pos"] +1 + indel["length"] + window] == "N" * window))):
            continue
        else:
            count = 0
            p = 1
            while count <= window and sample.seq[min((max(0, indel["pos"]) + indel["length"] + p), len(sample.seq)-1)] == "N":
                count += 1
                p+=1
            p = 0
            while count <= window and sample.seq[max(0, indel["pos"]) - p] == "N":
                count +=1
                p+=1
            if (count < window) or (window == 0):
                variant = f'{reference.id}\t{indel["pos"] + 1 - current_offset}\t.\t{indel["ref_base"]}\t{indel["alt_base"]}\t.\t.\tAC=1;AN=1\tGT\t1'
                vcf.append(variant)

    vcf = pd.DataFrame.from_records([sub.split("\t") for sub in vcf[1:]], columns = vcf[0].split(sep="\t"))
    return vcf

def calculate_n_coverage(ref, sample, sample_name):
    #first need to work out where S gene begins and ends in the sample (may be different due to insertions compared to ref length)
    ref_seq = str(ref.seq)
    sample_seq = str(sample.seq)
    pre_s = ref_seq[:21562]
    pre_s_ref_indels = pre_s.count("-")
    ref_s = ref_seq[21563 + pre_s_ref_indels : 25385] #position of s gene in 0-index
    ref_indels = ref_s.count("-")
    refdiff = ref_indels
    s_start = 21563 + pre_s_ref_indels
    s_end = 25385 + ref_indels
    while refdiff != 0:
        ref_s = ref_seq[s_start : s_end]
        ref_indels = ref_s.count("-")
        refdiff -= ref_indels
        s_end = 25385 + refdiff
    sample_s = sample_seq[s_start: s_end]
    global_n_perc = sample_seq.count("N") / len(sample_seq)
    s_n_perc = sample_s.count("N") / len(sample_s)
    contig_S_n_len = len(max(re.compile("N+").findall(str(sample_s)),default=""))
    sample_rbd = sample_seq[s_start + (319 * 3): s_start + (541*3)]
    rbd_n_nts = sample_rbd.count("N")
    ncov = {"sample_id" : sample_name , "global_n" : global_n_perc, "s_n" : s_n_perc, "s_n_contig": contig_S_n_len, "rbd_n" : rbd_n_nts}
    n_cov_info = pd.DataFrame([ncov])
    return n_cov_info
    #22520 - 23186

def process_sample(alignment_file, ref, outdir, vcf_dir, window, allow_ambiguous, out_suffix):
    ref_aliases = ["NC_045512.2", "MN908947.3"]
    if ref in ref_aliases:
        ref = "NC_045512.2"
    else:
        ref = ref

    for record in SeqIO.parse(alignment_file, "fasta"):
        if record.id == ref:
            reference = record
        else:
            sample = record
    outfile = os.path.join(outdir, f"{sample.id}{out_suffix}")
    sample = mask_trimmed_sequence(sample)
    indels = get_indels(reference, sample, window, allow_ambiguous)
    sample_n_cov_info = calculate_n_coverage(reference, sample, sample.id)

    if vcf_dir is not None:
        vcf_file = os.path.join(vcf_dir, f"{sample.id}.vcf")
        snps_header, snps = parse_vcf(vcf_file, split_info_cols=False)
        snps["POS"] = snps["POS"].astype(int)
        snps["group"] = (snps["POS"].diff() != 1).cumsum()
        agg_cols = {"POS": "first", "REF": lambda x: ''.join(x), "ALT": lambda x: ''.join(x),
                    "QUAL": "first", "INFO": "first", sample.id: "first", "FILTER": "first", "FORMAT": "first"}
        snps = snps.groupby("group").agg(agg_cols)
        snps["ID"] = snps["REF"] + snps["POS"].astype("str") + snps["ALT"]

        snps = snps.loc[snps["POS"].astype("int").isin(indels["POS"].astype("int").values.tolist()) == False]
        fatovcf = pd.concat([indels, snps])
        fatovcf["POS"] = fatovcf["POS"].astype(int)
        fatovcf = fatovcf.sort_values(by=["POS"], ascending=True)
        fatovcf["#CHROM"] = ref
        write_vcf(snps_header, fatovcf, outfile)
    else:
        indels["POS"] = indels["POS"].astype(int)
        indels = indels.sort_values(by=["POS"], ascending=True)
        csv_vcf = indels
        csv_vcf["#CHROM"] = sample.id
        write_vcf([], indels, outfile)
    return sample_n_cov_info

def process_alignment_file(alignment_file, ref, out_dir, vcf_dir, window, allow_ambiguous, out_suffix):
    return process_sample(alignment_file, ref, out_dir, vcf_dir,
                          window, allow_ambiguous, out_suffix)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--alignments-dir', metavar='path/to/alignments', type=str,
                        help='Directory containing alignment files in fasta format')
    parser.add_argument('--vcf-dir', metavar='path/to/vcfs', type=str,
                        help='Directory containing VCF files for merging SNPs and indels')
    parser.add_argument('--ref', metavar='MN908947.3', type=str,
                        help='Reference sequence')
    parser.add_argument('--out-dir', metavar='output_dir', type=str,
                        help="Output directory for processed files")
    parser.add_argument('--window', metavar="", type=int, default=2,
                        help="Flanking N filter for indels, set to 0 for off")
    parser.add_argument('--nperc', metavar="sample.nperc.csv", type=str,
                        help="Output file path for n-percentage csv")
    parser.add_argument('--allowAmbiguous', default=False, action='store_true',
                        help="Toggle whether to exclude ambiguous bases in SNPs and insertions")
    parser.add_argument('--out-suffix', metavar="", type=str, default = ".indels.vcf")
    parser.add_argument('--threads', metavar="", type=int, default = 1)
    args = parser.parse_args()

    alignments_dir = args.alignments_dir
    vcf_dir = args.vcf_dir
    ref = args.ref
    out_dir = args.out_dir
    all_sample_n_cov_info = []

    pool = multiprocessing.Pool(processes=args.threads)

    results = pool.starmap(process_alignment_file, [(os.path.join(alignments_dir, alignment_file), ref, out_dir, vcf_dir,
                                                     args.window, args.allowAmbiguous, args.out_suffix)
                                                    for alignment_file in os.listdir(alignments_dir)
                                                    if alignment_file.endswith(".fasta") or alignment_file.endswith(".fa")])

    all_sample_n_cov_info.extend(results)

    pool.close()
    pool.join()

    n_cov_info_df = pd.concat(all_sample_n_cov_info)

    n_cov_info_df.to_csv(args.nperc, index=False)
            
if __name__ == "__main__":
    main()