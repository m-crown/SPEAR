#!/usr/bin/env python

from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from summarise_snpeff import parse_vcf
import csv

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
        help='Destination dir for summary csv files')
    parser.add_argument('--vcfs', nargs='+', required=True,
        help='Input VCF files')
    args = parser.parse_args()
    
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    sample_summaries = []
    for vcf in args.vcfs:
        sample_name = Path(vcf).stem
        header, vcf , infocols = parse_vcf(vcf, split_info_cols = True)
        if len(vcf) != 0: #do not add summary if the vcf file is empty. 
            vcf["SPEAR"] = vcf["SPEAR"].str.split(",", expand = False)
            vcf = vcf.explode("SPEAR")
            vcf[["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]] = vcf["SPEAR"].str.split("|", expand = True)
            vcf.loc[vcf["Gene_Name"].str.contains('-'), "Gene_Name"] = "Intergenic_" + vcf.loc[vcf["Gene_Name"].str.contains('-'), "Gene_Name"]
            vcf.loc[vcf["Annotation"] == "synonymous_variant", "variant"] = vcf.loc[vcf["Annotation"] == "synonymous_variant", "HGVS.c"]
            per_sample_output = vcf.copy()
            per_sample_output["Sample"] = sample_name
            cols = ["Sample", "Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]
            per_sample_output = per_sample_output[cols]
            per_sample_output.to_csv(f'{args.output_dir}/{sample_name}.summary.csv', sep = '\t', header = True, index = False)
            vcf["variant"] = vcf["Annotation"] + "_" + vcf["variant"]
            new_df = vcf[["Gene_Name", "variant"]].copy()
            new_df = new_df.groupby("Gene_Name")["variant"].apply(list).to_frame()
            new_df["variant"] = new_df["variant"].str.join("|")
            new_df = new_df.reset_index().rename(columns={new_df.index.name:'Gene_Name'})
            new_df["cat"] =  new_df["Gene_Name"] + ":" + new_df["variant"]
            sample_variants = new_df["cat"].to_list()
            sample_variants = ";".join(sample_variants)
            sample_summary = [sample_name, sample_variants]
            sample_summaries.append(sample_summary)
        else:
            #touch file if empty, so snakemake doesnt fail. Implement a new version of this. Maybe check the vcf has one entry prior to entering pipeline.  
            Path(f'{args.output_dir}/{sample_name}.summary.csv').touch()
    sample_summaries.sort()
    with open(f'{args.output_dir}/summary.csv', "a") as fp:
        wr = csv.writer(fp, delimiter=',')
        wr.writerows(sample_summaries)

if __name__ == "__main__":
    main()
    #GZIP OUTPUT FILE
    #need to handle ambiguous 
    #need to handle multiple spear annotations