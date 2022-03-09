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
from functools import reduce
from bindingcalculator import BindingCalculator
from itertools import takewhile

def get_contextual_bindingcalc_values(residues_list,binding_calculator, option):
    res_ret_esc_df = binding_calculator.escape_per_site(residues_list)
    if option == "res_ret_esc":
        res_ret_esc = res_ret_esc_df.loc[res_ret_esc_df["site"].isin(residues_list), "retained_escape"]
        return(res_ret_esc)
    else:
        ab_escape_fraction = 1 - binding_calculator.binding_retained(residues_list)
        return(ab_escape_fraction)

def summarise_score(summary_df, metric):
    #assumes grouping by sample_id and summarising for each sample
    summary_df_info = summary_df.groupby("sample_id").agg({metric: ['sum', 'min', 'max']})
    summary_df_info.columns = summary_df_info.columns.droplevel(0)
    summary_df_info = summary_df_info.reset_index()
    summary_df_info = summary_df_info.rename_axis(None, axis=1)

    summary_df_mins = pd.merge(left = summary_df, right = summary_df_info[["sample_id", "min"]], left_on = ["sample_id", metric], right_on = ["sample_id", "min"])
    summary_df_mins[metric + "_min"] = summary_df_mins["residues"] + ":" + summary_df_mins[metric].fillna("").astype(str)
    summary_df_mins = summary_df_mins[["sample_id",metric + "_min"]].groupby("sample_id").agg({metric + "_min" : lambda x : list(x)})
    summary_df_mins[metric + "_min"] = summary_df_mins[metric + "_min"].str.join(",")

    summary_df_max = pd.merge(left = summary_df, right = summary_df_info[["sample_id", "max"]], left_on = ["sample_id", metric], right_on = ["sample_id", "max"])
    summary_df_max[metric + "_max"] = summary_df_max["residues"] + ":" + summary_df_max[metric].fillna("").astype(str)
    summary_df_max = summary_df_max[["sample_id",metric + "_max"]].groupby("sample_id").agg({metric + "_max" : lambda x : list(x)})
    summary_df_max[metric + "_max"] = summary_df_max[metric + "_max"].str.join(",")

    summary_df_sum = summary_df.groupby("sample_id").agg({metric: sum})
    summary_df_sum.columns = [metric + "_sum"]

    summary_df_final = summary_df_sum.merge(summary_df_max,on='sample_id').merge(summary_df_mins,on='sample_id')

    return(summary_df_final)

def sample_header_format(item,sample,vcf,filtered,vcf_loc):
    #if filter and vcf comes from masked
    #if alignment or consensus comes from indels
    #\n are interpreted literally here
    if vcf:
        if item.startswith("##bcftools_mergeCommand=merge"):
            if filtered:
                item = re.sub(r'(?<=merged\.vcf )[a-zA-Z0-9_\. \/]+(?=\n)', vcf_loc, item)
            else:
                item = re.sub(r'(?<=merged\.vcf )[a-zA-Z0-9_\. \/]+(?=\n)', vcf_loc, item)
    else:
        if item.startswith("##reference="):
            item = re.sub(r'(?<=alignment\/)[a-zA-Z0-9_\.]+(?=\.fasta)', f'{sample}', item)
        if item.startswith("##reference="): 
            item = re.sub(r'(?<=alignment\/)[a-zA-Z0-9_\.]+(?=\.fasta)', f'{sample}', item)
            item = re.sub(r'(?<=fatovcf\/)[a-zA-Z0-9_\.]+(?=\.vcf)', f'{sample}', item)

        if item.startswith("##bcftools_mergeCommand=merge"):
            if filtered:
                item = re.sub(r'(?<=merged\.vcf )[a-zA-Z0-9_\. \/]+(?=\n)', vcf_loc, item)
            else:
                item = re.sub(r'(?<=merged\.vcf )[a-zA-Z0-9_\. \/]+(?=;)', vcf_loc, item)
    return(item)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', metavar='anno_concat.tsv', type=str,
        help='Concatenated SPEAR anno file')
    parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
        help='Destination dir for summary tsv files')
    parser.add_argument('data_dir', metavar='data/', type=str,
        help='Data dir for binding calculator data files')
    parser.add_argument('input_header', metavar='merged.vcf', type=str,
        help='Merged VCF file for header retrieval')
    parser.add_argument('sample_list', metavar='', nargs="+",
        help='list of inputs to detect no variant samples ')    
    parser.add_argument('--is_vcf_input', default="False", type=str,
        help = "Set input file type to VCF")
    parser.add_argument('--is_filtered', default="False", type=str,
        help = "Specify files come from filtered directory")

    args = parser.parse_args()

    Path(f'{args.output_dir}/per_sample_annotation').mkdir(parents=True, exist_ok=True)
    if args.is_vcf_input:
        if args.is_filtered:
            infiles = f'{args.output_dir}/intermediate_output/masked/*.masked.vcf'
        else:
            infiles = f'{args.output_dir}/intermediate_output/indels/*.indels.vcf'
    else:
        infiles = f'{args.output_dir}/intermediate_output/indels/*.indels.vcf'

    with open(args.input_header, 'r') as fobj:
        headiter = takewhile(lambda s: s.startswith('#'), fobj)
        merged_header = pd.Series(headiter)
    merged_header = merged_header.str.replace("\n", "")
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"]

    merged_header = merged_header[~merged_header.str.startswith('#CHROM')]
    input_file = pd.read_csv(args.input_vcf, sep = "\t", names = ["sample_id", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "end"])
    if len(input_file) > 0:
        input_file["sample_id"] = input_file["sample_id"].str.extract(r'final_vcfs\/([a-zA-Z0-9\._]+)\.spear\.vcf', expand = False)
        input_file[["AN", "AC", "ANN", "SUM", "SPEAR"]] = input_file["INFO"].str.split(';',expand=True)
        original_cols = input_file.columns.tolist()
        input_file[["AN", "AC", "ANN", "SUM", "SPEAR"]] = input_file[["AN", "AC", "ANN", "SUM", "SPEAR"]].replace("^[A-Z]+=", "", regex = True)
        input_file = input_file.loc[input_file["ANN"] != "no_annotation"]

        input_file["SUM"] = input_file["SUM"].str.split(",", expand = False)
        input_file = input_file.explode("SUM") #explode on this list of SUM values
        input_file[["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues"]] = input_file["SUM"].str.split("|", expand = True)
        total_variants = input_file.loc[input_file["Annotation"].isin(["intergenic", "non-synonymous"]) == False, ["sample_id", "POS"]].groupby("sample_id").nunique().reset_index()
        total_variants.columns = ["sample_id", "total_variants"]
        total_variants["total_variants"] = total_variants["total_variants"].astype("int")

        consequence_type_counts = input_file[["sample_id", "Annotation"]].groupby(["sample_id", "Annotation"], as_index = False)["Annotation"].value_counts()
        consequence_type_counts["consequence_count"] = consequence_type_counts["Annotation"] + ":" + consequence_type_counts["count"].astype(str)

        type_string = consequence_type_counts[["sample_id" ,"consequence_count"]].groupby("sample_id").agg({"consequence_count" : lambda x : list(x)})
        type_string["consequence_count"] = type_string["consequence_count"].str.join(",")
        type_string.columns = ["consequence_type_variants"]
        type_string.reset_index(drop = False)

        input_file["SPEAR"] = input_file["SPEAR"].str.split(",", expand = False)
        input_file = input_file.explode("SPEAR")
        input_file[["spear-product", "residues","region", "domain", "feature", "contact_type", "NAb", "barns_class", "bloom_ACE2", "VDS", "serum_escape", "mAb_escape_all_classes", "cm_mAb_escape_all_classes","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES", "BEC_EF"]] = input_file["SPEAR"].str.split("|", expand = True)
        input_file = input_file.loc[input_file["product"] == input_file["spear-product"]]
        pattern = re.compile(r"[a-zA-Z\*]+([0-9]+)") #matches any point mutations or deletions , not insertions.
        input_file["respos"] = input_file["residues"].str.extract(pattern).fillna(-1).astype("int")
        input_file["refres"] = input_file["residues"].str.extract(r"([a-zA-Z\*]+)[0-9]+[a-zA-Z\?\*]+")
        input_file["altres"] = input_file["residues"].str.extract(r"[a-zA-Z\*]+[0-9]+([a-zA-Z\?\*]+)")

        input_file.loc[input_file["Gene_Name"].str.contains('-'), "Gene_Name"] = "Intergenic_" + input_file.loc[input_file["Gene_Name"].str.contains('-'), "Gene_Name"]
        input_file.loc[input_file["Annotation"] == "synonymous_variant", "variant"] = input_file.loc[input_file["Annotation"] == "synonymous_variant", "HGVS.c"]

        bindingcalc = BindingCalculator(csv_or_url = f'{args.data_dir}/escape_calculator_data.csv')
        rbd_residues = input_file.loc[(input_file["Gene_Name"] == "S") & (input_file["respos"] >= 331) & (input_file["respos"] <= 531)]
        if len(rbd_residues) > 0:
            sample_ef = rbd_residues.groupby("sample_id").agg({"respos" : lambda x : get_contextual_bindingcalc_values(x, bindingcalc, "escape_fraction")}).reset_index()
            sample_ef.columns = ["sample_id", "BEC_EF_sample"]
            input_file = input_file.merge(sample_ef, on = "sample_id")
            input_file.loc[(input_file["Gene_Name"] != "S") | ((input_file["Gene_Name"] == "S") & ((input_file["respos"] < 331) | (input_file["respos"] > 531))), "BEC_EF_sample"] = ""

            sample_res = get_contextual_bindingcalc_values(rbd_residues,rbd_residues, bindingcalc, "res_ret_esc")
            print(sample_res)
            sample_res.columns = ["sample_id", "BEC_RES"]
            input_file = input_file.merge(sample_res, on = "sample_id")
            input_file.loc[(input_file["Gene_Name"] != "S") | ((input_file["Gene_Name"] == "S") & ((input_file["respos"] < 331) | (input_file["respos"] > 531))), "BEC_RES"] = ""
        else: 
            input_file["BEC_EF_sample"] = ""

        final_samples = input_file.copy()
        cols = ['spear-product', 'residues', 'region', 'domain', 'feature', 'contact_type', 'NAb', 'barns_class', 'bloom_ACE2', 'VDS', 'serum_escape', 'mAb_escape_all_classes', 'cm_mAb_escape_all_classes', 'mAb_escape_class_1', 'mAb_escape_class_2', 'mAb_escape_class_3', 'mAb_escape_class_4', 'BEC_RES', 'BEC_EF', 'BEC_EF_sample']
        final_samples["SPEAR"] = final_samples[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
        all_cols = final_samples.columns.tolist()


        final_samples.drop([col for col in all_cols if col not in original_cols], axis = 1, inplace = True)
        cols = [e for e in final_samples.columns.to_list() if e not in ("SUM", "SPEAR")]

        final_vcf = final_samples.groupby(cols, as_index = False).agg({"SUM": set , "SPEAR": list})
        final_vcf["SUM"] = [','.join(map(str, l)) for l in final_vcf['SUM']]
        final_vcf["SPEAR"] = [','.join(map(str, l)) for l in final_vcf["SPEAR"].apply(lambda x: set(sorted(x, key = lambda y: re.search(r'^[a-zA-Z\*]+([0-9]+)|',y)[1] if re.search(r'^[a-zA-Z\*]+([0-9]+)|',y)[1] else "")))] #sorting like this because the groupby list doesnt always put residues in correct order. use set around list to remove duplicate annotations on NSP11 and RDRP overlap.
        infocols = ["AC", "AN", "ANN", "SUM", "SPEAR"]
        for col in infocols:
            final_vcf[col] = col + "=" + final_vcf[col].fillna("").astype("str")

        final_vcf['INFO'] = final_vcf[infocols].agg(';'.join, axis=1)
        final_vcf.drop(infocols, axis = 1, inplace = True)
        original_cols = [col for col in original_cols if col not in infocols]
        final_vcf = final_vcf[original_cols]

        for sample in args.sample_list:
            sample_header = merged_header.copy().apply(lambda x : sample_header_format(x, sample, args.is_vcf_input, args.is_filtered, infiles))
            #sample_header.to_csv(f'output/{sample}.vcf', index = False)
            np.savetxt(f'{args.output_dir}/final_vcfs/{sample}.spear.vcf', sample_header.values, fmt = "%s")
            if sample in final_vcf["sample_id"]:
                sample_vcf = final_vcf.loc[final_vcf["sample_id"] == sample].copy()
                sample_vcf["sample_id"] = "NC_045512.2"
                sample_vcf.columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
                #append to header created above
                sample_vcf.to_csv(f'{args.output_dir}/final_vcfs/{sample}.spear.vcf', sep = "\t" ,  mode = 'a', index = False)
            

        cols = ["sample_id", "POS", "REF", "ALT", "Gene_Name", "HGVS.c", "Annotation", "variant", "spear-product", "protein_id", "residues","region", "domain", "feature", "contact_type", "NAb", "barns_class", "bloom_ACE2", "VDS", "serum_escape", "mAb_escape_all_classes", "cm_mAb_escape_all_classes","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES","BEC_EF", "BEC_EF_sample", "refres", "altres", "respos"]
        input_file = input_file[cols]
        input_file.columns = ["sample_id", "POS", "REF", "ALT", "Gene_Name", "HGVS.nt", "consequence_type", "HGVS", "description", "RefSeq_acc", "residues","region", "domain", "feature", "contact_type", "NAb", "barns_class", "bloom_ACE2", "VDS", "serum_escape", "mAb_escape_all_classes", "cm_mAb_escape_all_classes","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES", "BEC_EF", "BEC_EF_sample", "refres", "altres", "respos"] 
        input_file[[col for col in input_file.columns if col not in ["refres", "altres", "respos"]]].to_csv(f'{args.output_dir}/spear_annotation_summary.tsv', sep = "\t", index = False)
        for sample in args.sample_list:
            if sample in input_file["sample_id"].values:
                sample_summary = input_file.loc[input_file["sample_id"] == sample].copy()
                sample_summary.to_csv(f'{args.output_dir}/per_sample_annotation/{sample}.spear.annotation.summary.tsv', sep = "\t", index = False)
            else:
                Path(f'{args.output_dir}/per_sample_annotation/{sample}.spear.annotation.summary.tsv').touch() #touch file if empty    

        #now getting summary scores 
        #subset the dataframe to remove synonymous residue variants (or rather, keep anything that isnt synonymous)
        summary = input_file.loc[((input_file["refres"] != input_file["altres"]) & (input_file["residues"].isin([""]) == False)) | ((input_file["residues"].str.contains("[A-Z\*][0-9]+[A-Z\*\?]", regex = True) == False) & (input_file["residues"].isin([""]) == False))]
        summary.reset_index(inplace = True, drop = True)
        sample_residue_variant_number = summary.groupby("sample_id").size().to_frame("total_residue_variants").reset_index()
        sample_residue_variant_number["total_residue_variants"] = sample_residue_variant_number["total_residue_variants"].astype("int")

        region_counts = summary.loc[summary["region"] != "" , ["sample_id", "description", "region"]]
        region_counts["region"] = region_counts["region"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True)
        if region_counts["region"].isin([""]).all():
            region_counts["region_residues"] = ""
        else:
            region_counts = region_counts.groupby("sample_id", as_index = False).value_counts(["Gene_Name", "region"])
            region_counts["region_residues"] = region_counts["description"] + ":" + region_counts["region"] + ":" + region_counts["count"].astype(str)
            region_counts = region_counts[["sample_id" ,"region_residues"]].groupby("sample_id").agg({"region_residues" : lambda x : list(x)})
            region_counts["region_residues"] = region_counts["region_residues"].str.join(",")

        domain_counts = summary.loc[summary["domain"] != "" , ["sample_id", "description", "domain"]]
        domain_counts["domain"] = domain_counts["domain"].str.split(",").explode().replace(r'^\s*$', np.nan, regex=True)
        if domain_counts["domain"].isin([""]).all():
            domain_counts["domain_residues"] = ""
        else:
            domain_counts = domain_counts.groupby("sample_id", as_index = False).value_counts(["Gene_Name", "domain"])
            domain_counts["domain_residues"] = domain_counts["description"] + ":" + domain_counts["domain"] + ":" + domain_counts["count"].astype(str)
            domain_counts = domain_counts[["sample_id" ,"domain_residues"]].groupby("sample_id").agg({"domain_residues" : lambda x : list(x)})
            domain_counts["domain_residues"] = domain_counts["domain_residues"].str.join(",")

        contacts = summary[["sample_id", "description","contact_type"]].copy()
        contacts["contact_type"] = contacts["contact_type"].str.split(" ")
        contacts = contacts.explode("contact_type")
        if contacts["contact_type"].isin([""]).all():
            contacts[["contact", "contact_type"]] = np.nan
            contacts["score"] = 0
        else:
            contacts[["contact" , "contact_type"]] = contacts["contact_type"].replace("", np.nan).str.split(":", expand = True)
            contacts = contacts.reset_index(drop = True)
            contacts["contact_type"] = contacts["contact_type"].str.split("+")
            contacts = contacts.explode(["contact_type"])
            contacts["contact_type"] = contacts["contact_type"].str.extract(r'^([a-z-A-Z]+)_*', expand = False)
            contacts["score"] = contacts["contact_type"].replace({"h-bond": 2, "contact": 1, "salt-bridge": 3})

        #these dataframes can be used for more than just ACE2 and trimer contact scores in future
        contacts_scores = contacts[["sample_id", "description", "contact_type", "contact", "score"]].groupby(["sample_id","description", "contact_type", "contact"]).sum().reset_index()
        contacts_counts = contacts[["sample_id", "description", "contact_type", "contact"]].groupby(["sample_id","description", "contact"]).count().reset_index()

        ace2_contacts_score = contacts_scores.loc[contacts_scores["contact"]=="ACE2",["sample_id", "contact", "score"]].groupby(["sample_id", "contact"]).sum().astype("int")
        ace2_contacts_score.columns = ["ACE2_contact_score"]
        trimer_contacts_score = contacts_scores.loc[contacts_scores["contact"]=="trimer",["sample_id", "contact", "score"]].groupby(["sample_id", "contact"]).sum().astype("int")
        trimer_contacts_score.columns = ["trimer_contact_score"]

        ace2_contacts_sum = contacts_counts.loc[contacts_counts["contact"]=="ACE2",["sample_id", "contact_type"]]
        ace2_contacts_sum.columns = ["sample_id", "ACE2_contact_counts"]
        trimer_contacts_sum = contacts_counts.loc[contacts_counts["contact"]=="trimer",["sample_id", "contact_type"]]
        trimer_contacts_sum.columns = ["sample_id", "trimer_contact_counts"]

        barns = summary[["sample_id", "description","barns_class"]].copy()

        if barns["barns_class"].isin([""]).all():
            barns_counts_grouped = barns["sample_id"]
            barns_counts_grouped["barns_class_variants"] = ""
        else: 
            barns["barns_class"] = barns["barns_class"].str.split(",")
            barns = barns.explode(["barns_class"])
            barns["barns_class"] = barns["barns_class"].str.split("+")
            barns_counts = barns.explode(["barns_class"]).replace(r'^\s*$', np.nan, regex=True).dropna()
            barns_counts["barns_class_name"] = "class_" + barns_counts['barns_class']
            barns_counts_grouped = barns_counts[["sample_id", "barns_class_name", "barns_class"]].groupby(["sample_id", "barns_class_name"]).count().reset_index().sort_values(["sample_id", "barns_class"], ascending = [True, False])
            barns_counts_grouped["barns_class_variants"] = barns_counts_grouped["barns_class_name"] + ":" + barns_counts_grouped["barns_class"].astype("str")
            barns_counts_grouped = barns_counts_grouped[["sample_id", "barns_class_variants"]].groupby("sample_id").agg({"barns_class_variants" : lambda x : list(x)})
            barns_counts_grouped["barns_class_variants"] = barns_counts_grouped["barns_class_variants"].str.join(",")

        scores = ["bloom_ACE2", "VDS", "serum_escape", "mAb_escape_all_classes", "cm_mAb_escape_all_classes","mAb_escape_class_1", "mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4","BEC_EF_sample"]
        score_df_list = [total_variants, sample_residue_variant_number, type_string, region_counts,domain_counts, ace2_contacts_sum,ace2_contacts_score,trimer_contacts_sum,trimer_contacts_score, barns_counts_grouped]

        for score in scores:
            if score == "BEC_EF_sample":
                score_subset = summary.loc[summary[score].isin([""]) == False, ["sample_id", "description", "residues", score]].copy()
                score_subset[score] = score_subset[score].astype(float)
                score_subset_df_sum = score_subset.groupby("sample_id").mean()
                score_subset_df_sum.columns = [score]
                score_df_list.append(score_subset_df_sum)
            else: 
                score_subset = summary.loc[summary[score].isin([""]) == False, ["sample_id", "description", "residues", score]].copy()
                score_subset[score] = score_subset[score].astype(float)
                score_summary = summarise_score(score_subset, score)
                score_df_list.append(score_summary)

        scores_df = reduce(lambda left,right: pd.merge(left,right,on="sample_id", how = "outer"), score_df_list)
        scores_df = scores_df.sort_values(by = ["sample_id"], ascending = True).fillna("")
        scores_df.rename(columns={'mAb_escape_sum':'mAb_escape_all_classes_sum', 'mAb_escape_min':'mAb_escape_all_classes_min', 'mAb_escape_max':'mAb_escape_all_classes_max', 'cm_mAb_escape_sum':'cm_mAb_escape_all_classes_sum', 'cm_mAb_escape_min':'cm_mAb_escape_all_classes_min', 'cm_mAb_escape_max':'cm_mAb_escape_all_classes_max'}, inplace=True)
        scores_df.to_csv(f'{args.output_dir}/spear_score_summary.tsv' , sep = "\t", index = False)
    else:
        for sample in args.sample_list:
            Path(f'{args.output_dir}/per_sample_annotation/{sample}.spear.annotation.summary.tsv').touch() #touch file if empty
        Path(f'{args.output_dir}/spear_score_summary.tsv').touch()
        Path(f'{args.output_dir}/spear_annotation_summary.tsv').touch()
if __name__ == "__main__":
    main()