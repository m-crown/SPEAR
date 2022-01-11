#!/usr/bin/env python

from pandas.io.formats.format import SeriesFormatter
from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from summarise_snpeff import parse_vcf, write_vcf
import csv
import re
from bindingcalculator import BindingCalculator

def main():

    def df_counts_to_string(summary_df):
        names = summary_df.index.to_list()
        counts = summary_df.to_list()
        list_items = [a + ":" + str(b) for a,b in zip(names,counts)]
        string_items = ",".join(c for c in list_items)
        return(string_items)

    def get_contextual_bindingcalc_values(residues_list, respos, binding_calculator, option):
        res_ret_esc_df = binding_calculator.escape_per_site(residues_list)
        if option == "res_ret_esc":
            res_ret_esc = res_ret_esc_df.loc[res_ret_esc_df["site"] == respos, "retained_escape"].values[0]
            return(res_ret_esc)
        else:
            ab_escape_fraction = 1 - bindingcalc.binding_retained(residues_list)
            return(ab_escape_fraction)

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
        help='Destination dir for summary tsv files')
    parser.add_argument('data_dir', metavar='data/', type=str,
        help='Data dir for binding calculator data files')
    parser.add_argument('--vcfs', nargs='+', required=True,
        help='Input VCF files')
    args = parser.parse_args()
    
    Path(f'{args.output_dir}/per_sample_annotation').mkdir(parents=True, exist_ok=True)
    
    sample_summaries = []
    sample_summary_cols = ["sample_id", "POS", "REF", "ALT", "gene_name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4","BEC_RES", "BEC_EF", "BEC_EF_sample"]
    with open(f'{args.output_dir}/spear_annotation_summary.tsv', 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(sample_summary_cols)
    
    scores_columns = ["sample_id","total_variants","total_residue_variants", "consequence_type_variants", "region_variants", "domain_variants", "ACE2_contact_counts","ACE2_contact_score","trimer_contact_counts", "trimer_contact_score", "barns_class_variants", "bloom_ACE2" , "VDS", "serum_escape", "mAb_escape_all_classes", "mAb_escape_class_1", "mAb_escape_class_2", "mAb_escape_class_3", "mAb_escape_class_4", "BEC_RES", "BEC_EF", "BEC_EF_sample"]
    with open(f'{args.output_dir}/spear_score_summary.tsv', 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(scores_columns)
    
    for vcf_file in args.vcfs:
        sample_name = Path(Path(vcf_file).stem)
        sample_name = sample_name.stem #take off spear from input vcf - not very adaptable for other inputs but works for now 
        header, vcf , infocols = parse_vcf(vcf_file, split_info_cols = True)
        original_cols = vcf.columns.tolist()
        if len(vcf) != 0: #do not add summary if the vcf file is empty (but the empty file has to be created so need to handle).
            #need to have two things happen here, for intergenics need to produce an output with the SUM field, for non intergenics summary output is going to come from SPEAR (or from both?) 
            vcf = vcf.loc[vcf["ANN"] != "no_annotation"]
            total_variants = len(vcf)
            vcf["SUM"] = vcf["SUM"].str.split(",", expand = False) #split SUM field where multiple annotations remain (NSP11 RDRP overlap)
            vcf = vcf.explode("SUM") #explode on this list of SUM values
            vcf[["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues"]] = vcf["SUM"].str.split("|", expand = True)
            consequence_type_counts = vcf.groupby(["Annotation"]).size() #get the nucleotide variant conseqeunce type counts before exploding the vcf into per residue format. This will add up to the total number of SUM values in INFO field rather than total number of variants.
            vcf["SPEAR"] = vcf["SPEAR"].str.split(",", expand = False)
            vcf = vcf.explode("SPEAR")
            vcf[["residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES", "BEC_EF"]] = vcf["SPEAR"].str.split("|", expand = True)
            pattern = re.compile(r"[a-zA-Z]+([0-9]+)") #matches any point mutations or deletions , not insertions. 
            vcf["respos"] = vcf["residues"].str.extract(pattern).fillna(-1).astype("int")
            
            #replace the per sample bloom calculator scores with contextual per sample scores. 
            bindingcalc = BindingCalculator(csv_or_url = f'{args.data_dir}/escape_calculator_data.csv')
            respos_list = vcf.loc[(vcf["Gene_Name"] == "S") & (vcf["respos"] >= 331) & (vcf["respos"] <= 531), "respos"].values.tolist()
            vcf.loc[(vcf["Gene_Name"] == "S") & (vcf["respos"] >= 331) & (vcf["respos"] <= 531), ["BEC_RES"]] = vcf.loc[(vcf["Gene_Name"] == "S") & (vcf["respos"] >= 331) & (vcf["respos"] <= 531)].apply(lambda x: get_contextual_bindingcalc_values(respos_list, x["respos"], bindingcalc, "res_ret_esc"), axis=1)
            vcf.loc[(vcf["Gene_Name"] == "S") & (vcf["respos"] >= 331) & (vcf["respos"] <= 531), ["BEC_EF_sample"]] = vcf.loc[(vcf["Gene_Name"] == "S") & (vcf["respos"] >= 331) & (vcf["respos"] <= 531)].apply(lambda x: get_contextual_bindingcalc_values(respos_list, x["respos"], bindingcalc, "escape_fraction"), axis=1)
            vcf["BEC_EF_sample"] = vcf["BEC_EF_sample"].fillna("")

            vcf.loc[vcf["Gene_Name"].str.contains('-'), "Gene_Name"] = "Intergenic_" + vcf.loc[vcf["Gene_Name"].str.contains('-'), "Gene_Name"]
            vcf.loc[vcf["Annotation"] == "synonymous_variant", "variant"] = vcf.loc[vcf["Annotation"] == "synonymous_variant", "HGVS.c"]

            #writing the updated bloom score vcf
            final_vcf = vcf.copy()
            cols = ['residues', 'region', 'domain', 'contact_type', 'NAb', 'barns_class', 'bloom_ace2', 'VDS', 'serum_escape', 'mAb_escape', 'cm_mAb_escape', 'mAb_escape_class_1', 'mAb_escape_class_2', 'mAb_escape_class_3', 'mAb_escape_class_4', 'BEC_RES', 'BEC_EF', 'BEC_EF_sample']
            final_vcf["SPEAR"] = final_vcf[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
            all_cols = final_vcf.columns.tolist()
            final_vcf.drop([col for col in all_cols if col not in original_cols], axis = 1, inplace = True)
            cols = [e for e in final_vcf.columns.to_list() if e not in ("SUM", "SPEAR")]
            final_vcf = final_vcf.groupby(cols, as_index = False).agg({"SUM": set , "SPEAR": list})
            final_vcf["SUM"] = [','.join(map(str, l)) for l in final_vcf['SUM']]
            final_vcf["SPEAR"] = [','.join(map(str, l)) for l in final_vcf["SPEAR"].apply(lambda x: set(sorted(x, key = lambda y: re.search(r'^[a-zA-Z]+([0-9]+)|',y)[1])))] #sorting like this because the groupby list doesnt always put residues in correct order. use set around list to remove duplicate annotations on NSP11 and RDRP overlap.
            for col in infocols:
                final_vcf[col] = col + "=" + final_vcf[col]
            final_vcf['INFO'] = final_vcf[infocols].agg(';'.join, axis=1)
            final_vcf.drop(infocols, axis = 1, inplace = True)
            original_cols = [col for col in original_cols if col not in infocols]
            original_cols.insert(original_cols.index("FILTER") + 1, "INFO")
            final_vcf = final_vcf[original_cols]
            write_vcf(header,final_vcf,vcf_file)

            per_sample_output = vcf.copy()
            per_sample_output["sample_id"] = sample_name
            cols = ["sample_id", "POS", "REF", "ALT", "Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES","BEC_EF", "BEC_EF_sample"]
            per_sample_output = per_sample_output[cols]
            per_sample_output.columns = ["sample_id", "POS", "REF", "ALT", "Gene_Name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mAb_escape", "cm_mAb_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES", "BEC_EF", "BEC_EF_sample"] 
            per_sample_output.to_csv(f'{args.output_dir}/per_sample_annotation/{sample_name}.spear.annotation.summary.tsv', sep = '\t', header = True, index = False)
            per_sample_output.to_csv(f'{args.output_dir}/spear_annotation_summary.tsv', sep = '\t', mode='a', header=False, index = False)

            #now getting summary scores
            sample_residue_variant_number = len(per_sample_output)
            type_string = df_counts_to_string(consequence_type_counts)
            region_counts = per_sample_output["region"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).value_counts()
            region_string = df_counts_to_string(region_counts)
            domain_counts = per_sample_output["domain"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True).value_counts()
            domain_string = df_counts_to_string(domain_counts)
            if per_sample_output["contact_type"].isin([""]).all(): #if there are no contact types in SPEAR annotation
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
            
            #rethink this part at some point

            if per_sample_output["bloom_ace2"].isin([""]).all():
                bloom_ace2_sum = ""
            else:
                bloom_ace2_sum = per_sample_output["bloom_ace2"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["VDS"].isin([""]).all():
                vds_sum = ""
            else:
                vds_sum = per_sample_output["VDS"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()

            if per_sample_output["serum_escape"].isin([""]).all():
                serum_escape_sum = ""
            else:
                serum_escape_sum = per_sample_output["serum_escape"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["mAb_escape"].isin([""]).all():
                mab_escape_all_sum = ""
            else:
                mab_escape_all_sum = per_sample_output["mAb_escape"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["mAb_escape_class_1"].isin([""]).all():
                mab_escape_class_1_sum = ""
            else:
                mab_escape_class_1_sum = per_sample_output["mAb_escape_class_1"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["mAb_escape_class_2"].isin([""]).all():
                mab_escape_class_2_sum = ""
            else:
                mab_escape_class_2_sum = per_sample_output["mAb_escape_class_2"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["mAb_escape_class_3"].isin([""]).all():
                mab_escape_class_3_sum = ""
            else:    
                mab_escape_class_3_sum = per_sample_output["mAb_escape_class_3"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["mAb_escape_class_4"].isin([""]).all():
                mab_escape_class_4_sum = ""
            else:
                mab_escape_class_4_sum = per_sample_output["mAb_escape_class_4"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            
            if per_sample_output["BEC_EF"].isin([""]).all():
                bec_ef_score = ""
            else:
                bec_ef_score = per_sample_output["BEC_EF"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            if per_sample_output["BEC_EF_sample"].isin([""]).all():
                bec_ef_sample_score = ""
            else:
                bec_ef_sample_score = per_sample_output["BEC_EF_sample"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()
            if per_sample_output["BEC_RES"].isin([""]).all():
                bec_res_score = ""
            else:
                bec_res_score = per_sample_output["BEC_RES"].replace(r'^\s*$', np.nan, regex=True).astype("float").mean()

            scores_list = [sample_name,total_variants,sample_residue_variant_number, type_string, region_string, domain_string,ace2_contacts_sum, ace2_contacts_score, trimer_contacts_sum,trimer_contacts_score, barns_string,bloom_ace2_sum , vds_sum, serum_escape_sum, mab_escape_all_sum, mab_escape_class_1_sum, mab_escape_class_2_sum, mab_escape_class_3_sum, mab_escape_class_4_sum, bec_res_score, bec_ef_score, bec_ef_sample_score]
            scores_df = pd.DataFrame([scores_list], columns=scores_columns)
            scores_df.to_csv(f'{args.output_dir}/spear_score_summary.tsv', sep = '\t', mode='a', header=False, index = False)

            #now gofasta style output 
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
            Path(f'{args.output_dir}/per_sample_annotation/{sample_name}.spear.annotation.summary.tsv').touch() #touch file if empty  
    sample_summaries.sort() #alphabetically ordered samples 
    with open(f'{args.output_dir}/spear_variant_summary.csv', "a") as fp:
        wr = csv.writer(fp, delimiter=',')
        wr.writerows(sample_summaries)

if __name__ == "__main__":
    main()