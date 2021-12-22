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
        help='Destination dir for summary tsv files')
    parser.add_argument('--vcfs', nargs='+', required=True,
        help='Input VCF files')
    args = parser.parse_args()
    
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    sample_summaries = []

    #["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "retained_escape", "ab_escape_fraction"] 
    #write these col headers to a summary file then also write the single sample csvs to this summary following creation

    for vcf in args.vcfs:
        sample_name = Path(vcf).stem
        header, vcf , infocols = parse_vcf(vcf, split_info_cols = True)
        if len(vcf) != 0: #do not add summary if the vcf file is empty.
            vcf["SPEAR"] = vcf["SPEAR"].str.split(",", expand = False)
            vcf = vcf.explode("SPEAR")
            vcf[["gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "BEC_EF", "BEC_RES"]] = vcf["SPEAR"].str.split("|", expand = True)
            vcf.loc[vcf["gene_name"].str.contains('-'), "gene_name"] = "Intergenic_" + vcf.loc[vcf["gene_name"].str.contains('-'), "gene_name"]
            vcf.loc[vcf["Annotation"] == "synonymous_variant", "variant"] = vcf.loc[vcf["Annotation"] == "synonymous_variant", "HGVS.c"]
            per_sample_output = vcf.copy()
            per_sample_output["sample_id"] = sample_name
            cols = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "BEC_EF", "BEC_RES"]
            per_sample_output[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "retained_escape", "ab_escape_fraction"]] = per_sample_output[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "BEC_EF", "BEC_RES"]].replace("[^-0-9a-zA-Z]+[~]+", "", regex = True)
            per_sample_output[["residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "retained_escape", "ab_escape_fraction"]] = per_sample_output[["residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "BEC_EF", "BEC_RES"]].replace("~",",", regex = True)
            per_sample_output = per_sample_output[cols]
            per_sample_output.columns = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape", "BEC_EF", "BEC_RES"] 
            per_sample_output.to_csv(f'{args.output_dir}/{sample_name}.summary.tsv', sep = '\t', header = True, index = False)
            vcf["variant"] = vcf["Annotation"] + "_" + vcf["variant"]
            new_df = vcf[["gene_name", "variant"]].copy()
            new_df = new_df.groupby("gene_name")["variant"].apply(list).to_frame()
            new_df["variant"] = new_df["variant"].str.join("|")
            new_df = new_df.reset_index().rename(columns={new_df.index.name:'gene_name'})
            new_df["cat"] =  new_df["gene_name"] + ":" + new_df["variant"]
            sample_variants = new_df["cat"].to_list()
            sample_variants = ";".join(sample_variants)
            sample_summary = [sample_name, sample_variants]
            sample_summaries.append(sample_summary)
        else:
            #touch file if empty  
            Path(f'{args.output_dir}/summary/{sample_name}.summary.tsv').touch()
    sample_summaries.sort()
    with open(f'{args.output_dir}/summary.csv', "a") as fp:
        wr = csv.writer(fp, delimiter=',')
        wr.writerows(sample_summaries)

if __name__ == "__main__":
    main()
    #GZIP OUTPUT FILE
    #need to handle ambiguous 
    #need to handle multiple spear annotations