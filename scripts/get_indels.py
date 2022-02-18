#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from pathlib import Path
from re import finditer
import re
from summarise_snpeff import parse_vcf, write_vcf
import pandas as pd

#need to have a flag to exclude ambiguous indels in get_indels to go alongside excluding ambiguous snps from fatovcf
def mask_trimmed_sequence(sample):
    def repl(m):
        return 'N' * len(m.group())
    sample.seq = re.sub(r'^-+|-+$', repl, str(sample.seq))
    return sample

def get_indels(reference, sample, window, allow_ambiguous):
    vcf = [f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample.id}'] #set indels vcf header.
    offset = 0 #offset used when insertions detected, subsequent indels should be shifted by this offset
    indels = []
    #find insertions from reference sequence
    for match in finditer("-+", str(reference.seq)):
        ref = match.start() - 1 #the 0 based index position of the last ref base before indel
        length = len(match.group()) #length of the indel itself
        ins_end = ref + length + 1
        ref_base = reference.seq[ref]
        alt_base = sample.seq[ref:ins_end]
        indel = {"type": "insertion","pos": ref, "ref_base": ref_base, "alt_base": alt_base, "length": length}
        indels.append(indel)
    #find deletions from sample
    for match in finditer("-+", str(sample.seq)):
        ref = match.start() - 1
        length = len(match.group())
        ins_end = ref + length + 1
        ref_base = reference.seq[ref:ins_end]
        alt_base = sample.seq[ref]
        indel = {"type": "deletion","pos": ref, "ref_base": ref_base, "alt_base": alt_base, "length": length}
        indels.append(indel)
    indels = sorted(indels, key=lambda d: d['pos'])
    for indel in indels:
        current_offset = offset
        if allow_ambiguous:
            iupac_chars = ["N"]
            if indel["type"] == "insertion": #filtering insertions that consist of only N characters , as these are interpreted as any nucletide downstream and are not reliable for annotation anyway. 
                offset += indel["length"] #still have to iterate the offset but dont mark the variant
                if any(char in indel["alt_base"][1:] for char in iupac_chars): 
                    continue
                if indel["alt_base"][0] in iupac_chars: #if the reference base is N in the sample convert to ref in the alt description (ignore N snp)
                    indel["alt_base"] = indel["ref_base"][0] + indel["alt_base"][1:]
            if indel["type"] == "deletion":
                if indel["alt_base"][0] in iupac_chars: #if the reference base is IUPAC ambiguous in the sample convert to ref in the alt description (ignore N snp)
                    indel["alt_base"] = indel["ref_base"][0]
        else:
            iupac_chars = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
            if indel["type"] == "insertion": #filtering insertions that consist of only N characters , as these are interpreted as any nucletide downstream and are not reliable for annotation anyway. 
                offset += indel["length"] #still have to iterate the offset but dont mark the variant
                if any(char in indel["alt_base"][1:] for char in iupac_chars): 
                    continue
                if indel["alt_base"][0] in iupac_chars: #if the reference base is IUPAC ambiguous in the sample convert to ref in the alt description (ignore N snp)
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

def calculate_n_coverage(ref, sample, sample_name, outpath):
    #first need to work out where S gene begins and ends in the sample (may be different due to insertions compared to ref length)
    pre_s = str(ref.seq[:21562])
    pre_s_ref_indels = pre_s.count("-")
    ref_s = str(ref.seq[21563 + pre_s_ref_indels : 25385]) #position of s gene in 0-index
    ref_indels = ref_s.count("-")
    refdiff = ref_indels
    s_start = 21563 + pre_s_ref_indels
    s_end = 25385 + ref_indels
    while refdiff != 0:
        ref_s = str(ref.seq[s_start : s_end])
        ref_indels = ref_s.count("-")
        refdiff -= ref_indels
        s_end = 25385 + refdiff
    sample_s = sample.seq[s_start: s_end]
    global_n_perc = sample.seq.count("N") / len(sample.seq)
    s_n_perc = sample_s.count("N") / len(sample_s)
    contig_S_n_len = len(max(re.compile("N+").findall(str(sample_s)),default=""))
    sample_rbd = sample.seq[s_start + (319 * 3): s_start + (541*3)]
    rbd_n_nts = sample_rbd.count("N")
    ncov = {"sample" : sample_name , "global_n" : global_n_perc, "s_n" : s_n_perc, "longest_continuous_s_n": contig_S_n_len, "rbd_n_nts" : rbd_n_nts}
    n_cov_info = pd.DataFrame([ncov])
    n_cov_info.to_csv(outpath, index = False, header = False)
    #22520 - 23186
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('alignment', metavar='sample.muscle.aln', type=str,
        help='input alignment in fasta format')
    parser.add_argument('ref', metavar='MN908947.3', type=str,
        help='reference sequence')
    parser.add_argument('outfile', metavar="", type = str,
        help="output filename") #HAVE AN INITIAL CHECK AFTER READING IN SAMPLE, default is sample name as read from muscle alignment
    parser.add_argument('--vcf', metavar="", type = str,
        help="input vcf file for merging snps and indels")
    parser.add_argument('--window', metavar="", type = int, default = 2,
        help="flanking N filter for indels, set to 0 for off")
    parser.add_argument('--tsv', metavar="sample.indels.tsv", type = str,
        help="output file path for indels in tsv format")
    parser.add_argument('--nperc', metavar="sample.nperc.csv", type = str,
        help="output file path for n-percentage csv")   
    parser.add_argument('--allowAmbiguous', default=False, action='store_true',
        help = "Toggle whether to exclude ambiguous bases in SNPs and insertions")
    args = parser.parse_args()

    ref_aliases = ["NC_045512.2", "MN908947.3"]
    if args.ref in ref_aliases:
        ref = "NC_045512.2"
    else:
        ref = args.ref

    for record in SeqIO.parse(args.alignment, "fasta"):
        if record.id == args.ref:
            reference = record
        else:
            sample = record

    sample = mask_trimmed_sequence(sample)
    indels = get_indels(reference, sample, args.window, args.allowAmbiguous)
    calculate_n_coverage(reference,sample, sample.id, args.nperc)
    if args.tsv:
        indels.to_csv(args.tsv, mode='w', index = False, sep = "\t")
    if args.vcf:
        snps_header, snps = parse_vcf(args.vcf, split_info_cols=False)
        snps["POS"] = snps["POS"].astype(int)
        mnps_index = []
        mnps_pos = []
        for i in range(1, len(snps)):
            x = 0
            mnp_index = []
            mnp_pos = []
            if any((i + x) in sl for sl in mnps_index): #check if 
                continue
            else:
                while ((i + x) < len(snps)) and (snps.loc[i + x, "POS"] - snps.loc[i + x -1, "POS"] == 1): #while not at end of snps and snps are consecutive positions
                    if x == 0: #if this is the first snp in an mnp
                        mnp_index.append(i + x -1)
                        mnp_pos.append(snps.loc[i + x -1, "POS"])
                        mnp_index.append(i + x)
                        mnp_pos.append(snps.loc[i + x, "POS"])
                    else: #if not the first snp in the mnp
                        mnp_index.append(i + x)
                        mnp_pos.append(snps.loc[i + x, "POS"])
                    x += 1
            if len(mnp_index) > 0:
                mnps_index.append(mnp_index)
                mnps_pos.append(mnp_pos)
        for index, pos in zip(mnps_index, mnps_pos):
            #set the first snp to be the MNP and then drop all the others
            snps.loc[index[0],["ID"]] = "."
            #join together the ref and alt alleles from multiple rows in VCF using list comprehension after converting the pandas subset to list of lists. Then turn into string. 
            snps.loc[index[0],["REF"]] = ''.join([item for sublist in snps.loc[[ind for ind in index],["REF"]].values.tolist() for item in sublist])
            snps.loc[index[0],["ALT"]] = ''.join([item for sublist in snps.loc[[ind for ind in index],["ALT"]].values.tolist() for item in sublist])
            snps = snps.drop(index[1:])

        snps = snps.loc[snps["POS"].astype("int").isin(indels["POS"].astype("int").values.tolist()) == False]              
        fatovcf = pd.concat([indels,snps])
        fatovcf["POS"] = fatovcf["POS"].astype(int)
        fatovcf = fatovcf.sort_values(by = ["POS"], ascending = True)
        fatovcf["#CHROM"] = ref
        write_vcf(snps_header, fatovcf,args.outfile) #output the vcf header and body to file 
    else:
        indels["POS"] = indels["POS"].astype(int)
        indels = indels.sort_values(by = ["POS"], ascending = True)
        csv_vcf = indels
        csv_vcf["#CHROM"] = sample.id
        csv_vcf.to_csv(args.outfile, mode='w', index = False, sep = "\t") #this file is only really necessary for comparison exercise
        write_vcf([], indels, args.outfile) #output the vcf header and body to file 

if __name__ == "__main__":
    main()