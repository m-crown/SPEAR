#!/usr/bin/env python

from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
from itertools import takewhile
import argparse
from pathlib import Path
import glob
import numpy as np
from summarise_snpeff import parse_vcf, write_vcf
import csv
import os
#take input directory as start

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('vcfs', metavar='path/to/vcfs', type = str,
        help='Input VCF file directory')
    parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
        help='Destination dir for summary csv files')
    args = parser.parse_args()
    
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    vcfs =glob.glob(f'{args.vcfs}/*.vcf')
    
    sample_summaries = []
    for vcf in vcfs:
        sample_name = Path(vcf).stem
        header, vcf , infocols = parse_vcf(vcf, split_info_cols = True)
        vcf["SPEAR"] = vcf["SPEAR"].str.split(",", expand = False)
        #REMOVE c. n. p.
        vcf = vcf.explode("SPEAR")
        vcf[["Gene", "HGVS.c" , "Annotation", "Variant", "Product" , "ID"]] = vcf["SPEAR"].str.split("|", expand = True)
        vcf.loc[vcf["Gene"].str.contains('-'), "Gene"] = "Intergenic_" + vcf.loc[vcf["Gene"].str.contains('-'), "Gene"]
        vcf.loc[vcf["Annotation"] == "synonymous_variant", "Variant"] = vcf.loc[vcf["Annotation"] == "synonymous_variant", "HGVS.c"]
        vcf["Variant"] = vcf["Annotation"] + "_" + vcf["Variant"]
        per_sample_output = vcf.copy()
        per_sample_output["Sample"] = sample_name
        cols = ["Sample", "Gene", "HGVS.c" , "Annotation", "Variant", "Product" , "ID"]
        per_sample_output = per_sample_output[cols]
        per_sample_output.to_csv(f'{args.output_dir}/{sample_name}.summary.csv', header = True, index = False)

        new_df = vcf[["Gene", "Variant"]].copy()
        new_df = new_df.groupby("Gene")["Variant"].apply(list).to_frame()
        new_df["Variant"] = new_df["Variant"].str.join("|")
        new_df = new_df.reset_index().rename(columns={new_df.index.name:'Gene'})
        new_df["cat"] =  new_df["Gene"] + ":" + new_df["Variant"]
        sample_variants = new_df["cat"].to_list()
        sample_variants = ";".join(sample_variants)
        sample_summary = [sample_name, sample_variants]
        sample_summaries.append(sample_summary)
    with open(f'{args.output_dir}/summary.csv', "a") as fp:
        wr = csv.writer(fp, delimiter=',')
        wr.writerows(sample_summaries)
if __name__ == "__main__":
    main()
    #GZIP OUTPUT FILE
    #sort alphanumerically before write
    #need to handle ambiguous 
    #need to handle multiple spear annotations