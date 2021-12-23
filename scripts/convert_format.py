#!/usr/bin/env python

from pandas.io.formats.format import SeriesFormatter
from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from summarise_snpeff import parse_vcf
import csv

def main():

    def df_counts_to_string(summary_df):
        names = summary_df.index.to_list()
        counts = summary_df.to_list()
        list_items = [a + ":" + str(b) for a,b in zip(names,counts)]
        string_items = ",".join(c for c in list_items)
        return(string_items)

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
        help='Destination dir for summary tsv files')
    parser.add_argument('--vcfs', nargs='+', required=True,
        help='Input VCF files')
    args = parser.parse_args()
    
    Path(f'{args.output_dir}/per_sample_annotation').mkdir(parents=True, exist_ok=True)
    
    sample_summaries = []
    sample_summary_cols = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"]
    with open(f'{args.output_dir}/spear_annotation_summary.tsv', 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(sample_summary_cols)
    
    scores_columns = ["sample_id","total_variants","contact_type_variants", "region_variants", "domain_variants", "ACE2_contact_counts","ACE2_contact_score","trimer_contact_counts", "trimer_contact_score", "barns_class_variants", "bloom_ACE2" , "VDS", "serum_escape", "mAb_escape_all_classes", "class_1_mAb_escape", "class_2_mAb_escape", "class_3_mAb_escape", "class_4_mAb_escape", "BEC_sample_EF", "BEC_REF"]
    with open(f'{args.output_dir}/spear_score_summary.tsv', 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(scores_columns)
    
    for vcf in args.vcfs:
        sample_name = Path(Path(vcf).stem)
        sample_name = sample_name.stem #take off spear from input vcf - not very adaptable for other inputs but works for now 
        header, vcf , infocols = parse_vcf(vcf, split_info_cols = True)
        if len(vcf) != 0: #do not add summary if the vcf file is empty (but the empty file has to be created so need to handle). 
            vcf["SPEAR"] = vcf["SPEAR"].str.split(",", expand = False)
            vcf = vcf.explode("SPEAR")
            vcf[["gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"]] = vcf["SPEAR"].str.split("|", expand = True)
            vcf.loc[vcf["gene_name"].str.contains('-'), "gene_name"] = "Intergenic_" + vcf.loc[vcf["gene_name"].str.contains('-'), "gene_name"]
            vcf.loc[vcf["Annotation"] == "synonymous_variant", "variant"] = vcf.loc[vcf["Annotation"] == "synonymous_variant", "HGVS.c"]
            per_sample_output = vcf.copy()
            per_sample_output["sample_id"] = sample_name
            cols = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"]
            per_sample_output[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"]] = per_sample_output[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES","BEC_sample_EF"]].replace("[^-0-9a-zA-Z]+[~]+", "", regex = True)
            per_sample_output[["residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"]] = per_sample_output[["residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES","BEC_sample_EF"]].replace("~",",", regex = True)
            per_sample_output = per_sample_output[cols]
            per_sample_output.columns = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_EF", "BEC_RES", "BEC_sample_EF"] 
            per_sample_output.to_csv(f'{args.output_dir}/per_sample_annotation/{sample_name}.spear.annotation.summary.tsv', sep = '\t', header = True, index = False)
            per_sample_output.to_csv(f'{args.output_dir}/spear_annotation_summary.tsv', sep = '\t', mode='a', header=False, index = False)

            #now getting summary scores
            sample_variant_number = len(per_sample_output)
            consequence_type_counts = per_sample_output.groupby(["consequence_type"]).size()
            type_string = df_counts_to_string(consequence_type_counts)
            region_counts = per_sample_output["region"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).value_counts()
            region_string = df_counts_to_string(region_counts)
            domain_counts = per_sample_output["domain"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).value_counts()
            domain_string = df_counts_to_string(domain_counts)
            if per_sample_output["contact_type"].isin([""]).all():
                ace2_contacts_score = ""
                trimer_contacts_score = ""
                ace2_contacts_sum = ""
                trimer_contacts_sum = ""
            else:
                contacts_df = per_sample_output["contact_type"].str.split(" ").explode().replace("", np.nan).str.split(":", expand = True).reset_index(drop = True)
                contacts_df.rename(columns={0: "contact", 1: "contact_type"}, inplace = True)
                contacts_df["contact_type"] = contacts_df["contact_type"].str.split("+")
                contacts_df = contacts_df.explode(["contact_type"])
                contacts_df["contact_type"] = contacts_df["contact_type"].str.extract(r'^([a-z-A-Z]+)_*', expand = False)
                contacts_df["score"] = contacts_df["contact_type"].replace({"h-bond": 2, "contact": 1, "salt-bridge": 3})
                contacts_scores = contacts_df.groupby(["contact"]).sum()
                ace2_contacts_score = contacts_scores.loc[contacts_scores.index=="ACE2", "score"].values[0] if len(contacts_scores.loc[contacts_scores.index=="ACE2", "score"].values) == 1 else ""
                trimer_contacts_score = contacts_scores.loc[contacts_scores.index=="trimer", "score"].values[0] if len(contacts_scores.loc[contacts_scores.index=="trimer", "score"].values) == 1 else ""
                contacts_counts = contacts_df.groupby(["contact"]).count()
                ace2_contacts_sum = contacts_counts.loc[contacts_scores.index=="ACE2", "contact_type"].values[0] if len(contacts_counts.loc[contacts_scores.index=="ACE2", "contact_type"].values) == 1 else ""
                trimer_contacts_sum = contacts_counts.loc[contacts_scores.index=="trimer", "contact_type"].values[0] if len(contacts_counts.loc[contacts_scores.index=="trimer", "contact_type"].values) == 1 else ""
            if per_sample_output["barns_class"].isin([""]).all():
                barns_string = ""
            else:
                barns_counts = per_sample_output["barns_class"].explode(",").str.split("+").explode().replace(r'^\s*$', np.nan, regex=True).value_counts()
                barns_counts.index = "class_" + barns_counts.index
                barns_string = df_counts_to_string(barns_counts)
            
            if per_sample_output["bloom_ace2"].isin([""]).all():
                bloom_ace2_sum = ""
            else:
                bloom_ace2_sum = per_sample_output["bloom_ace2"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["VDS"].isin([""]).all():
                vds_sum = ""
            else:
                vds_sum = per_sample_output["VDS"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()

            if per_sample_output["serum_escape"].isin([""]).all():
                serum_escape_sum = ""
            else:
                serum_escape_sum = per_sample_output["serum_escape"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["mAb_escape"].isin([""]).all():
                mab_escape_all_sum = ""
            else:
                mab_escape_all_sum = per_sample_output["mAb_escape"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["mAb_escape_class_1"].isin([""]).all():
                mab_escape_class_1_sum = ""
            else:
                mab_escape_class_1_sum = per_sample_output["mAb_escape_class_1"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["mAb_escape_class_2"].isin([""]).all():
                mab_escape_class_2_sum = ""
            else:
                mab_escape_class_2_sum = per_sample_output["mAb_escape_class_2"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["mAb_escape_class_3"].isin([""]).all():
                mab_escape_class_3_sum = ""
            else:    
                mab_escape_class_3_sum = per_sample_output["mAb_escape_class_3"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["mAb_escape_class_4"].isin([""]).all():
                mab_escape_class_4_sum = ""
            else:
                mab_escape_class_4_sum = per_sample_output["mAb_escape_class_4"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            
            if per_sample_output["BEC_sample_EF"].isin([""]).all():
                bec_sample_ef_score = ""
            else:
                bec_sample_ef_score = per_sample_output["BEC_sample_EF"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["BEC_RES"].isin([""]).all():
                bec_res_score = ""
            else:
                bec_res_score = per_sample_output["BEC_RES"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).astype("float").sum()
            scores_list = [sample_name,sample_variant_number, type_string, region_string, domain_string,ace2_contacts_sum, ace2_contacts_score, trimer_contacts_sum,trimer_contacts_score, barns_string,bloom_ace2_sum , vds_sum, serum_escape_sum, mab_escape_all_sum, mab_escape_class_1_sum, mab_escape_class_2_sum, mab_escape_class_3_sum, mab_escape_class_4_sum, bec_sample_ef_score, bec_res_score]
            scores_df = pd.DataFrame([scores_list], columns=scores_columns)
            scores_df.to_csv(f'{args.output_dir}/spear_score_summary.tsv', sep = '\t', mode='a', header=False, index = False)

            #now gofasta style output 
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
            Path(f'{args.output_dir}/per_sample_annotation/{sample_name}.spear.annotation.summary.tsv').touch() #touch file if empty  
    sample_summaries.sort() #alphabetically ordered samples 
    with open(f'{args.output_dir}/spear_variant_summary.csv', "a") as fp:
        wr = csv.writer(fp, delimiter=',')
        wr.writerows(sample_summaries)

if __name__ == "__main__":
    main()