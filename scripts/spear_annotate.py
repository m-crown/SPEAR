#!/usr/bin/env python3

from typing_extensions import final

from pandas.io.pytables import dropna_doc
from Bio import SeqIO
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning) #ignoring the performance warning  PerformanceWarning: your performance may suffer as PyTables will pickle object types that it cannot map directly to c-types for spear_anno_file
from itertools import takewhile
import argparse
import numpy as np
import re
from shutil import copy
from bindingcalculator import BindingCalculator
from summarise_snpeff import parse_vcf , write_vcf

def natural_key(string_):
    """See https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
def annotate_residues(vcf, data_dir):
    def extract_df_scores(scores_matrix, respos, altres):
        return(scores_matrix.loc[respos, altres])

    def get_bindingcalc_values(respos, binding_calculator, option):
        res_ret_esc_df = binding_calculator.escape_per_site([respos])
        if option == "res_ret_esc":
            res_ret_esc = res_ret_esc_df.loc[res_ret_esc_df["site"] == respos, "retained_escape"].values[0]
            return(res_ret_esc)
        else:
            ab_escape_fraction = 1 - bindingcalc.binding_retained([respos])
            return(ab_escape_fraction)

    spear_anno_file = pd.read_hdf(f'{data_dir}/spear_data.h5', "spear_anno_file")
    spear_anno_file["mod_barns_class_mask_sum_gt0.75"] = spear_anno_file["mod_barns_class_mask_sum_gt0.75"].replace(["", np.nan], -1).astype(int)
    
    bloom_escape_all = pd.read_hdf(f'{data_dir}/spear_data.h5', "bloom_escape_all")
    bloom_escape_class1 = pd.read_hdf(f'{data_dir}/spear_data.h5', "bloom_escape_class1")
    bloom_escape_class2 = pd.read_hdf(f'{data_dir}/spear_data.h5', "bloom_escape_class2")
    bloom_escape_class3 = pd.read_hdf(f'{data_dir}/spear_data.h5', "bloom_escape_class3")
    bloom_escape_class4 = pd.read_hdf(f'{data_dir}/spear_data.h5', "bloom_escape_class4")
    greaney_serum_escape = pd.read_hdf(f'{data_dir}/spear_data.h5', "greaney_serum_escape")    
    
    bloom_ace2_file = pd.read_csv(f'{data_dir}/single_mut_effects.csv')
    vds = pd.read_csv(f'{data_dir}/vibentropy_occupancy_dmsdata.csv')

    bindingcalc = BindingCalculator(csv_or_url = f'{data_dir}/escape_calculator_data.csv')    
    vcf["SUM"] = vcf["SUM"].str.split(",", expand = False) #split SUM field where multiple annotations remain (NSP11 RDRP overlap)
    vcf = vcf.explode("SUM") #explode on this list of SUM values
    vcf[["gene_name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues"]] = vcf["SUM"].str.split("|", expand = True)
    vcf["residues"] = vcf["residues"].str.split('~')
    vcf = vcf.explode("residues")
    vcf.reset_index(inplace = True)
    
    pattern = re.compile(r"[a-zA-Z\*]+([0-9]+)") #matches any point mutations or deletions , not insertions. 
    vcf["respos"] = vcf["residues"].str.extract(pattern).fillna(-1).astype("int")
    vcf["refres"] = vcf["residues"].str.extract(r"([a-zA-Z\*])+[0-9]+")
    vcf = pd.merge(vcf,spear_anno_file[["ORF", "AA_coordinate","region", "domain", "feature", "contact_type", "NAb", "barns_class", "mod_barns_class_mask_sum_gt0.75", "product"]],left_on = ["product", "respos"], right_on = ["product","AA_coordinate"],how="left")
    vcf["mod_barns_class_mask_sum_gt0.75"] = vcf["mod_barns_class_mask_sum_gt0.75"].replace("", -1).fillna(-1).astype(int)
    bloom_ace2_file["ORF"] = "S"
    vcf = pd.merge(vcf,bloom_ace2_file[["ORF", "mutation", "bind_avg"]], left_on = ["gene_name", "residues"], right_on = ["ORF", "mutation"], how = "left")
    vcf.drop(["ORF_x", "AA_coordinate"], axis = 1, inplace = True)
    vds["ORF"] = "S"
    vds["mutation_residues"] = vds["wildtype"] + vds["site"].astype("str") + vds["mutation"]
    vcf = pd.merge(vcf, vds[["ORF", "mutation_residues", "mut_VDS"]], left_on = ["gene_name","residues"], right_on = ["ORF", "mutation_residues"], how = "left")
    vcf["alt_res"] = vcf["residues"].str.extract("[A-Z\*][0-9]+([A-Z\?\*])")
    
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) , "bloom_escape_all"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531)].apply(lambda x: extract_df_scores(bloom_escape_all, x["respos"], x["alt_res"]), axis=1)
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) , "serum_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531)].apply(lambda x: extract_df_scores(greaney_serum_escape, x["respos"], x["alt_res"]), axis=1)
    
    #make a copy of the barns classes to explode on
    vcf["mod_barns_class"] = vcf["mod_barns_class_mask_sum_gt0.75"].replace(-1, "").fillna("").astype("str")
    vcf["mod_barns_class"] = vcf["mod_barns_class"].fillna("")
    vcf["mod_barns_class"] = vcf["mod_barns_class"].astype("str")
    vcf.loc[vcf["mod_barns_class"] != "", "mod_barns_class"] = vcf.loc[vcf["mod_barns_class"] != "", "mod_barns_class"] + "*"
    vcf.loc[(vcf["barns_class"] != "") & (vcf["mod_barns_class"] != ""), "mod_barns_class"] = vcf.loc[(vcf["barns_class"] != "") & (vcf["mod_barns_class"] != ""),"barns_class"] + "+" + vcf.loc[(vcf["barns_class"] != "") & (vcf["mod_barns_class"] != ""),"mod_barns_class"]
    vcf.loc[vcf["mod_barns_class"] == "", "mod_barns_class"] = vcf.loc[vcf["mod_barns_class"] == "", "barns_class"]
    vcf.loc[(vcf["barns_class"] == "") & (vcf["mod_barns_class"] != ""), "mod_barns_class"] = vcf["mod_barns_class"]
    vcf["mod_barns_class"] = vcf["mod_barns_class"].str.split("+", expand = False)
    vcf.loc[vcf["mod_barns_class"].notna() , "mod_barns_class"] = vcf.loc[vcf["mod_barns_class"].notna() , "mod_barns_class"].apply(lambda x: sorted(x, key = natural_key))
    vcf.loc[vcf["mod_barns_class"].notna(), "mod_barns_class"] = vcf.loc[vcf["mod_barns_class"].notna(), "mod_barns_class"].str.join("+")
    
    vcf["barns_class"] = vcf["barns_class"].str.split("+")
    vcf = vcf.explode("barns_class")
    vcf["barns_class"] = vcf["barns_class"].replace("", -1).fillna(-1).astype("int")

    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 1), "mAb_class_1_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 1)].apply(lambda x: extract_df_scores(bloom_escape_class1, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 2), "mAb_class_2_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 2)].apply(lambda x: extract_df_scores(bloom_escape_class2, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 3), "mAb_class_3_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 3)].apply(lambda x: extract_df_scores(bloom_escape_class3, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 4), "mAb_class_4_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["barns_class"] == 4)].apply(lambda x: extract_df_scores(bloom_escape_class4, x["respos"], x["alt_res"]), axis=1).fillna("")

    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 1), "mAb_class_1_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 1)].apply(lambda x: extract_df_scores(bloom_escape_class1, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 2), "mAb_class_2_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 2)].apply(lambda x: extract_df_scores(bloom_escape_class2, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 3), "mAb_class_3_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 3)].apply(lambda x: extract_df_scores(bloom_escape_class3, x["respos"], x["alt_res"]), axis=1).fillna("")
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 4), "mAb_class_4_escape"] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["mod_barns_class_mask_sum_gt0.75"] == 4)].apply(lambda x: extract_df_scores(bloom_escape_class4, x["respos"], x["alt_res"]), axis=1)    
    
    vcf["mAb_class_1_escape"] = vcf["mAb_class_1_escape"].fillna("").astype("str")
    vcf["mAb_class_2_escape"] = vcf["mAb_class_2_escape"].fillna("").astype("str")
    vcf["mAb_class_3_escape"] = vcf["mAb_class_3_escape"].fillna("").astype("str")
    vcf["mAb_class_4_escape"] = vcf["mAb_class_4_escape"].fillna("").astype("str")

    vcf["mod_barns_class_mask_sum_gt0.75"] = vcf["mod_barns_class_mask_sum_gt0.75"].replace(-1, "").astype("str")

    cols = [e for e in vcf.columns.to_list() if e not in ("barns_class","mAb_class_1_escape", "mAb_class_2_escape", "mAb_class_3_escape", "mAb_class_4_escape")]
    vcf = vcf.groupby(cols, as_index = False, dropna = False).agg(set)

    vcf["mAb_class_1_escape"] = [''.join(map(str, l)) for l in vcf["mAb_class_1_escape"]]
    vcf["mAb_class_2_escape"] = [''.join(map(str, l)) for l in vcf["mAb_class_2_escape"]]
    vcf["mAb_class_3_escape"] = [''.join(map(str, l)) for l in vcf["mAb_class_3_escape"]]
    vcf["mAb_class_4_escape"] = [''.join(map(str, l)) for l in vcf["mAb_class_4_escape"]]

    vcf["mAb_class_1_escape"] = vcf["mAb_class_1_escape"].replace("", np.nan).astype("float")
    vcf["mAb_class_2_escape"] = vcf["mAb_class_2_escape"].replace("", np.nan).astype("float")
    vcf["mAb_class_3_escape"] = vcf["mAb_class_3_escape"].replace("", np.nan).astype("float")
    vcf["mAb_class_4_escape"] = vcf["mAb_class_4_escape"].replace("", np.nan).astype("float")
    
    
    vcf["barns_class"] = vcf["mod_barns_class"]
    vcf.drop(["mod_barns_class"], axis = 1 , inplace = True)

    vcf["cm_mab_escape"] = vcf[["mAb_class_1_escape","mAb_class_2_escape","mAb_class_3_escape","mAb_class_4_escape"]].mean(axis=1)
    bindingcalc_data = pd.read_csv(f'{data_dir}/escape_calculator_data.csv')
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) & (vcf["refres"] != vcf["alt_res"]), ["BEC_RES"]] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) & (vcf["refres"] != vcf["alt_res"]) & (vcf["respos"].isin(bindingcalc_data["site"].unique()))].apply(lambda x: get_bindingcalc_values(x["respos"], bindingcalc, "res_ret_esc"), axis=1)
    vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) & (vcf["refres"] != vcf["alt_res"]), ["BEC_EF"]] = vcf.loc[(vcf["gene_name"] == "S") & (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"]) == False) & (vcf["respos"] >= 331) & (vcf["respos"] <= 531) & (vcf["refres"] != vcf["alt_res"]) & (vcf["respos"].isin(bindingcalc_data["site"].unique()))].apply(lambda x: get_bindingcalc_values(x["respos"], bindingcalc, "escape_fraction"), axis=1)

    vcf = vcf.fillna("")
    vcf.loc[(vcf["refres"] == vcf["alt_res"]) | (vcf["alt_res"].isin([np.nan, "del", "fs", "*", "?"])), ['region', 'domain', 'feature', 'contact_type', 'NAb', 'barns_class', 'bind_avg', 'mut_VDS', 'serum_escape', 'bloom_escape_all', 'cm_mab_escape', 'mAb_class_1_escape', 'mAb_class_2_escape', 'mAb_class_3_escape', 'mAb_class_4_escape', 'BEC_RES', 'BEC_EF']] = "" #remove scores from residues with same and alt and ref residue.
    cols = ['product', 'residues', 'region', 'domain', 'feature', 'contact_type', 'NAb', 'barns_class', 'bind_avg', 'mut_VDS', 'serum_escape', 'bloom_escape_all', 'cm_mab_escape', 'mAb_class_1_escape', 'mAb_class_2_escape', 'mAb_class_3_escape', 'mAb_class_4_escape', 'BEC_RES', 'BEC_EF']
    vcf["SPEAR"] = vcf[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
    all_cols = vcf.columns.tolist()
    vcf.drop([col for col in all_cols if col not in ["original_index" , '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'AC','AN', 'ANN', 'SUM', "SPEAR", "problem_exc", "problem_filter" ]], axis = 1, inplace = True)
    return vcf

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('output_filename', metavar='spear_vcfs/merged.spear.vcf', type=str,
        help='Filename for SPEAR annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('vcf', metavar='path/to/vcf', type = str,
        help='Input VCF file')
    parser.add_argument('data_dir', metavar='path/to/data/', type = str,
        help ='Data files for peptide subpositions')
    args = parser.parse_args()

    header, vcf, infocols = parse_vcf(args.vcf)
    if "problem_exc" not in infocols:
        infocols.append("problem_exc")
        vcf["problem_exc"] = ""
    if "problem_filter" not in infocols: 
        infocols.append("problem_filter")
        vcf["problem_filter"] = ""
    infocols = ["AN", "AC", "problem_exc", "problem_filter", "ANN", "SUM"]
    
    if len(vcf) != 0: #do not add summary if the vcf file is empty.
        vcf = vcf.rename_axis('original_index').reset_index()
        df = vcf.iloc[:,:vcf.columns.get_loc("FORMAT")] 
        df = df.replace(np.nan, '', regex=True)
        samples = vcf.iloc[:,vcf.columns.get_loc("FORMAT"):] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
        samples["original_index"] = df["original_index"]
        header.append(f'##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'product | residue | region | domain | feature | contact_type | NAb | barns_class | bloom_ace2 | VDS | serum_escape | mAb_escape | cm_mAb_escape | mAb_escape_class_1 | mAb_escape_class_2 | mAb_escape_class_3 | mAb_escape_class_4 | BEC_RES | BEC_EF | BEC_EF_sample  \'">') #MAKE VARIANT HEADER HGVS
        df = annotate_residues(df.copy(), args.data_dir)

        cols = [e for e in df.columns.to_list() if e not in ("problem_exc", "problem_filter", "ANN", "SUM", "SPEAR")]
        df = df.groupby(cols, as_index = False).agg({
            "problem_exc" : set,
            "problem_filter" : set,
            "ANN": set,
            "SUM": set,
            "SPEAR": list})

        df["problem_exc"] = [','.join(map(str, l)) for l in df['problem_exc']]
        df["problem_filter"] = [','.join(map(str, l)) for l in df['problem_filter']]
        df["ANN"] = [','.join(map(str, l)) for l in df['ANN']]
        df["SUM"] = [','.join(map(str, l)) for l in df['SUM']]

        df["SPEAR"] = [','.join(map(str, l)) for l in df["SPEAR"].apply(lambda x: set(sorted(x, key = lambda y: re.search(r'^[a-zA-Z]+([0-9]+)|',y)[1] if re.search(r'^[a-zA-Z]+([0-9]+)|',y)[1] else "")))] #sorting like this because the groupby list doesnt always put residues in correct order. use set around list to remove duplicate annotations on NSP11 and RDRP overlap.
        infocols.append("SPEAR")
        for col in infocols:
            df[col] = col + "=" + df[col]
        df['INFO'] = df[infocols].agg(';'.join, axis=1)
        df.drop(infocols, axis = 1, inplace = True)

        cols = df.columns.to_list()
        cols.pop(cols.index("INFO"))
        cols.insert(cols.index("FILTER") + 1, "INFO")
        df = df[cols]
        #need to sanity check before concat shape
        vcf = pd.merge(left = df, right = samples, on = "original_index", how = "inner")
        vcf.drop(columns = ["original_index"], inplace = True, axis = 1)
        write_vcf(header,vcf,args.output_filename)
    else:
        copy(args.vcf, args.output_filename) #copy the empty vcf file so that snakemake sees command completion. 
  

if __name__ == "__main__":
    main()