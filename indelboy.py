
import argparse
from Bio import SeqIO
from pathlib import Path
from re import finditer, split
from effingpropa import parse_vcf, write_vcf
import pandas as pd
import glob
from collections import defaultdict

def get_indels(reference, sample):
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
        if indel["type"] == "insertion":
            variant = f'{reference.id}\t{indel["pos"] +1 - offset}\t.\t{indel["ref_base"]}\t{indel["alt_base"]}\t.\t.\tAC=1;AN=1\tGT\t1'
            vcf.append(variant)
            offset += indel["length"]
        else:
            variant = f'{reference.id}\t{indel["pos"] +1 -offset}\t.\t{indel["ref_base"]}\t{indel["alt_base"]}\t.\t.\tAC=1;AN=1\tGT\t1'
            vcf.append(variant)
    vcf = pd.DataFrame.from_records([sub.split("\t") for sub in vcf[1:]], columns = vcf[0].split(sep="\t"))
    return vcf

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_dir', metavar='muscle/', type=str,
    help='input muscle alignment')
    parser.add_argument('ref', metavar='MN908947.3', type=str,
    help='ref seq')
    parser.add_argument('outdir', metavar="", type = str,
    help="outputdirectory")
    parser.add_argument('--vcf_dir', metavar="", type = str,
    help="input vcf directory for merging snps and indels, only if vcf_dir specified")
    parser.add_argument('--write_indels', dest='write_indels', action='store_true')
    parser.set_defaults(write_indels=False)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    alignments=glob.glob(f'{args.input_dir}/*.muscle.aln')
    ref_aliases = ["NC_045512.2", "MN908947.3"]
    if args.ref in ref_aliases:
        ref = "NC_045512.2"
    else:
        ref = args.ref

    for alignment in alignments:
        for record in SeqIO.parse(alignment, "fasta"):
            if record.id == args.ref:
                reference = record
            else:
                sample = record
        indels = get_indels(reference, sample)
        if args.write_indels:
            indels.to_csv(Path.joinpath(outdir,f'{sample.id}.indels.tsv'), mode='w', index = False, sep = "\t")
        if args.vcf_dir:
            snps_header, snps = parse_vcf(Path.joinpath(Path(args.vcf_dir),f'{sample.id}.vcf'), split_info_cols=False)
            snps["POS"] = snps["POS"].astype(int)
            mnps_index = []
            mnps_pos = []
            for i in range(1, len(snps)):
                x = 0
                mnp_index = []
                mnp_pos = []
                if any((i + x) in sl for sl in mnps_index):
                    continue
                else:
                    while ((i + x) < len(snps)) and (snps.loc[i + x, "POS"] - snps.loc[i + x -1, "POS"] == 1):
                        if x == 0:
                            mnp_index.append(i + x -1)
                            mnp_pos.append(snps.loc[i + x -1, "POS"])
                            mnp_index.append(i + x)
                            mnp_pos.append(snps.loc[i + x, "POS"])
                        else:
                            mnp_index.append(i + x)
                            mnp_pos.append(snps.loc[i + x, "POS"])
                        x += 1
                if len(mnp_index) > 0:
                    mnps_index.append(mnp_index)
                    mnps_pos.append(mnp_pos)
            for index, pos in zip(mnps_index, mnps_pos):
                snps.loc[index[0],["ID"]] = "."
                snps.loc[index[0],["REF"]] = str(reference.seq[pos[0] - 1: pos[-1]])
                snps.loc[index[0],["ALT"]] = str(sample.seq[pos[0] - 1: pos[-1]])
                snps = snps.drop(index[1:])
                            
            fatovcf = pd.concat([indels,snps])
            fatovcf["POS"] = fatovcf["POS"].astype(int)
            fatovcf = fatovcf.sort_values(by = ["POS"], ascending = True)
            fatovcf["#CHROM"] = ref
            csv_vcf = fatovcf.copy()
            csv_vcf["#CHROM"] = sample.id
            csv_vcf.to_csv(Path.joinpath(outdir,f'{sample.id}.tsv'), mode='w', index = False, sep = "\t") #this file is only really necessary for comparison exercise
            write_vcf(snps_header, fatovcf,sample.id,outdir) #output the vcf header and body to file 
        else:
            indels["POS"] = indels["POS"].astype(int)
            indels = indels.sort_values(by = ["POS"], ascending = True)
            csv_vcf = indels
            csv_vcf["#CHROM"] = sample.id
            csv_vcf.to_csv(Path.joinpath(outdir,f'{sample.id}.tsv'), mode='w', index = False, sep = "\t") #this file is only really necessary for comparison exercise
            write_vcf([], indels,sample.id,outdir) #output the vcf header and body to file 

if __name__ == "__main__":
    main()