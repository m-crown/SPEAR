#!/usr/bin/env python3

from typing_extensions import final
from Bio import SeqIO
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from itertools import takewhile
import argparse
import numpy as np
import re
from shutil import copy
from bindingcalculator import BindingCalculator

def convert_snpeff_annotation(vcf, gb_mapping, locus_tag_mapping, data_dir):

  def natural_key(string_):
    """See https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

  def filter_spear_summary(summaries):
    regexp = re.compile(r'ORF1(a|ab) polyprotein')
    if len(list(filter(None, list(regexp.search(summary) for summary in summaries)))) == len(summaries):
      summaries = [list(summaries)[0]]
    else:
      summaries = [summary for summary in summaries if not regexp.search(summary)]
    return summaries

  def annotate_s_residues(vcf, data_dir):
    s_residues = vcf.loc[vcf["product"] == "surface glycoprotein" , "residues"].str.split('~').to_list()
    
    spear_anno_file = pd.read_pickle(f'{data_dir}/spear.pkl')
    spear_anno_file["mod_barns_class_mask_sum_gt0.75"] = spear_anno_file["mod_barns_class_mask_sum_gt0.75"].astype("str")
    bloom_ace2_file = pd.read_csv(f'{data_dir}/single_mut_effects.csv').fillna("")
    bloom_escape_all = pd.read_csv(f'{data_dir}/Bloom_mAb_escape_all_class.csv', index_col = 0).fillna("")
    bloom_escape_class1 = pd.read_csv(f'{data_dir}/Bloom_mAb_escape_class_1.csv', index_col = 0)
    bloom_escape_class2 = pd.read_csv(f'{data_dir}/Bloom_mAb_escape_class_2.csv', index_col = 0)
    bloom_escape_class3 = pd.read_csv(f'{data_dir}/Bloom_mAb_escape_class_3.csv', index_col = 0)
    bloom_escape_class4 = pd.read_csv(f'{data_dir}/Bloom_mAb_escape_class_4.csv', index_col = 0)
    greaney_serum_escape = pd.read_csv(f'{data_dir}/Greaney_serum_escape.csv', index_col = 0).fillna("")
    vds = pd.read_csv(f'{data_dir}/vibentropy_occupancy_dmsdata.csv')
    bindingcalc = BindingCalculator(csv_or_url = f'{data_dir}/escape_calculator_data.csv')
    

    #these need to be sample specific and arent when calculated in here on merged sample
    #flat_s_respos = list(int(pattern.match(residue).group(1)) if pattern.match(residue) else -1 for residue in flat_s_residues)
    #flat_s_residues = [residue for variant in s_residues for residue in variant]
    #pattern = re.compile(r"[A-Z]+([0-9]+)([A-Z]+)")
    #s_binding_calc_pos = [pos for pos in flat_s_respos if (pos != -1) and (pos >= min(bindingcalc.sites)) and (pos <= max(bindingcalc.sites))] #s residues that are in the RBD as definded by bloom binding calc
    #sample_ef = 1 - bindingcalc.binding_retained(s_binding_calc_pos)
    #esc_per_site = bindingcalc.escape_per_site(s_binding_calc_pos)
    annotation = []
    for residues in s_residues:
      residues = [item for sublist in residues for item in sublist.split(",")]
      residue_anno = []
      for residue in residues:
        pattern = re.compile(r"[a-zA-Z]+([0-9]+)")
        respos = int(pattern.match(residue).group(1)) if pattern.match(residue) else -1 #it is an acknowledged limitation at this stage that this does not working for insertions
        if respos != -1:
          anno = spear_anno_file.loc[spear_anno_file["AA_coordinate"] == respos, ["region", "domain", "contact_type", "NAb", "barns_class", "mod_barns_class_mask_sum_gt0.75"]].values.tolist()
          anno = [item if len(item) == len(item) else [] for item in anno] #fill na's with empty lists
          anno = [item for sublist in anno for item in sublist]
          bloom_ace2 = bloom_ace2_file.loc[bloom_ace2_file["mutation"] == residue, "bind_avg"].values.tolist()
          bloom_ace2 = ["" if len(bloom_ace2) == 0 else str(bloom_ace2[0])]
          anno += bloom_ace2
        else:
          anno = ["","","","","","",""] #if residue doesnt match regex for annotation i.e. is an insertion       
        pattern = re.compile(r"[A-Z]+([0-9]+)([A-Z]+)")
        mutation_anno = []
        classes = list(int(x) if x != "" else -1 for x in anno[4].split("+"))
        mod_barns_class_mask_sum_gt = list(int(float(x)) if x != "" else -1 for x in anno[5].split("+"))
        search_classes = classes + mod_barns_class_mask_sum_gt
        str_classes = anno[4].split("+")
        str_mod_barns_class_mask_sum_gt = [str(int(float(s))) + "*" if s != "" else "" for s in anno[5].split("+")]
        str_mod_barns_class_mask_sum_gt = [value for value in str_mod_barns_class_mask_sum_gt if value != ""]
        final_class_combo = str_classes + str_mod_barns_class_mask_sum_gt
        final_class_list = sorted(final_class_combo, key = natural_key)
        final_class_list = '+'.join(final_class_list)
        final_class_list = "" if final_class_list == "+" else final_class_list
        anno[4] = final_class_list
        if pattern.match(residue): #cant do this annotation for del456 etc so only if it matches the pattern of a residue substitution
          respos = int(pattern.match(residue).group(1)) 
          altres = pattern.match(residue).group(2)
          if respos >= 14 and respos <= 913:
            vds_score = vds.loc[(vds["site"] == respos) & (vds["mutation"] == altres), "mut_VDS"].fillna("").to_list()
            mutation_anno += vds_score
            if respos >= 331 and respos <= 531: #data for these residues of spike only (RBD) - above this is calculated as max and min of calculator sites
              #ab_escape_fraction_sample = sample_ef needs to be set elsewhere to be sample specific 
              mab_escape = bloom_escape_all.loc[respos, altres] if bloom_escape_all.loc[respos, altres] != np.nan else ""
              serum_escape = greaney_serum_escape.loc[respos,altres] if greaney_serum_escape.loc[respos,altres] != np.nan else ""
              res_ret_esc_df = bindingcalc.escape_per_site([respos])
              res_ret_esc = res_ret_esc_df.loc[res_ret_esc_df["site"] == respos, "retained_escape"].values[0]
              ab_escape_fraction = 1 - bindingcalc.binding_retained([respos])
              cm_mab_escape = []
              mAb_class_1_escape = ""
              mAb_class_2_escape = ""
              mAb_class_3_escape = ""
              mAb_class_4_escape = ""
              for mab_class in search_classes:
                if mab_class == 1:
                  mAb_class_1_escape = bloom_escape_class1.loc[respos,altres]
                  cm_mab_escape.append(bloom_escape_class1.loc[respos,altres])  
                elif mab_class == 2:
                  mAb_class_2_escape = bloom_escape_class2.loc[respos,altres]
                  cm_mab_escape.append(bloom_escape_class2.loc[respos,altres])
                elif mab_class == 3:
                  mAb_class_3_escape = bloom_escape_class3.loc[respos,altres]
                  cm_mab_escape.append(bloom_escape_class3.loc[respos,altres])
                elif mab_class == 4:
                  mAb_class_4_escape = bloom_escape_class4.loc[respos,altres]
                  cm_mab_escape.append(bloom_escape_class4.loc[respos,altres])
              if len(cm_mab_escape) > 0 and sum(cm_mab_escape) > 0:
                cm_mab_escape = sum(cm_mab_escape)/len(cm_mab_escape)
              else:
                cm_mab_escape = ""
              mutation_anno += [str(serum_escape), str(mab_escape), str(cm_mab_escape),str(mAb_class_1_escape), str(mAb_class_2_escape), str(mAb_class_3_escape), str(mAb_class_4_escape), str(res_ret_esc), str(ab_escape_fraction)]
            else:
              mutation_anno += ["","","","","","","","","", ""] #if residue isnt in the RBD append 5 empty strings in replacement
          else:
            mutation_anno += ["","","","","","","","","","", ""] #if residue not in RBD and not in range for VDS score append 6 empty strings 
        else:
          mutation_anno += ["","","","","","","","","","", ""] #if residue not compatible with annotation i.e. not an AA substitution  
        anno += mutation_anno
        anno.pop(5) #remove internal antibody class from output
        #replace (4) with combined 4 and aterisk 5.  
        residue_anno.append(anno)  
      residue_anno = list(map(list, zip(*residue_anno)))
      residue_anno = ['~'.join([str(c) if c == c else "" for c in lst]) for lst in residue_anno]
      annotation.append(residue_anno)
    vcf.loc[vcf["product"] == "surface glycoprotein" , ["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS","serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4","BEC_RES","BEC_EF"]] = annotation
    vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS","serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4","BEC_RES", "BEC_EF"]] = vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS","serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES","BEC_EF"]].fillna("")
    vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS","serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4", "BEC_RES","BEC_EF"]] = vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS","serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4",  "BEC_RES","BEC_EF"]].replace("^[~]+[~]*|[~]+[~]*$", "", regex = True)
    return vcf
  #takes a input a dataframe row, splits the ann field into a new
  vcf["ANN2"] = vcf["ANN"].str.split(',') #put the split ANN field into a new column to preserve original snpeff value
  vcf = vcf.explode("ANN2")
  snpeff_anno_cols = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"]
  vcf[snpeff_anno_cols] = vcf["ANN2"].str.split('|', expand = True)
  vcf.loc[vcf["Feature_ID"] == "GU280_gp01","Feature_ID"]="YP_009724389.1"
  vcf.loc[vcf["Feature_ID"] == "GU280_gp01.2","Feature_ID"]="YP_009724389.1" #setting this id to be the same as orf1ab polyprotein for spear annotation purposes. 
  vcf.reset_index(inplace = True)
  vcf.loc[vcf["Annotation"] != "intergenic_region","variant"] = vcf["HGVS.p"]
  vcf.loc[vcf["Annotation"] == "intergenic_region", "variant"] = vcf["HGVS.c"]
  vcf["product"] = vcf["Feature_ID"].map(lambda x: gb_mapping.get(x))
  vcf.loc[vcf["product"].isna(), "product"] = vcf.loc[vcf["product"].isna(), "Gene_Name"]
  vcf["protein_id"] = vcf["product"].map(lambda x: locus_tag_mapping.get(x))
  vcf["protein_id"] = vcf["protein_id"].fillna("")

  #convert HGVS annotation to per residue for residue annotation. 
  #handle deletions and insertions residues first
  vcf[["start_res", "start_pos", "end_res", "end_pos", "change", "ins"]] = vcf["variant"].str.extract('([A-Z])([0-9]+)_*([A-Z])?([0-9]+)?([a-z]+)([A-Z]+)?')
  vcf["end_pos"] = vcf["end_pos"].fillna(vcf["start_pos"])
  vcf[["start_pos", "end_pos"]] = vcf[["start_pos", "end_pos"]].fillna(0).astype("int")
  vcf["ins_length"] = vcf['ins'].str.len().fillna(0).astype('int')
  range_df = pd.DataFrame(list(range(i, j+1)) for i, j in vcf[["start_pos", "end_pos"]].values).add_prefix("position_")
  vcf["length"] = range_df.count(axis = 1)
  range_df = "del" + range_df.fillna(0).astype('int').astype('str')
  range_df = range_df.replace("del0", "")
  cols = range_df.columns
  vcf = pd.concat([vcf,range_df], axis = 1)
  vcf["inserted_residues"] = vcf["start_res"] + vcf["start_pos"].astype('str') + vcf["ins"].fillna('').astype('str')
  vcf["range"] = vcf.apply(lambda x : list(range(x['ins_length'],x['length'])),1)
  vcf.loc[(vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]) == False) & (vcf["change"] == "delins"), "residues"] = vcf.loc[vcf["change"] == "delins","inserted_residues"] + "~" + vcf.apply(lambda x: '~'.join(x[b] for b in cols[x.range]),axis=1)
  vcf.loc[(vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]) == False) & (vcf["change"] == "del"), "residues"] = vcf.loc[vcf["change"] == "del"].apply(lambda x: '~'.join(x[b] for b in cols[x.range]),axis=1)
  vcf.loc[(vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]) == False) & (vcf["change"] == "ins"), "residues"] = vcf.loc[vcf["change"] == "ins", "start_res"] + vcf["start_pos"].astype("str") + "-" + vcf.loc[vcf["change"] == "ins","ins"] + "-" + vcf.loc[vcf["change"] == "ins","end_res"] + vcf.loc[vcf["change"] == "ins","end_pos"].astype('str')
  #frameshift and stop variants
  vcf.loc[vcf["Annotation"].astype("str") == "frameshift_variant", "residues"] = vcf.loc[vcf["Annotation"].astype("str") == "frameshift_variant", "variant"].str.replace("^p\.", '')
  vcf.loc[vcf["Annotation"].astype("str") == "stop_gained", "residues"] = vcf.loc[vcf["Annotation"].astype("str") == "stop_gained", "variant"].str.extract("^p\.([A-Z][0-9]+)\*") + "stop"
  #now handle missense variants and synonymous variants
  vcf[["start_res2", "start_pos2", "end_res2"]] = vcf["variant"].str.extract('([A-Z]+)([0-9]+)([A-Z]+)')
  vcf["end_pos2"] = vcf["end_pos"].fillna(0).astype("int")
  vcf["start_res2"] = vcf["start_res2"].fillna("").astype("str")
  vcf["end_res2"] = vcf["end_res2"].fillna("").astype("str")
  vcf["start_pos2"] = vcf["start_pos2"].fillna(0).astype("int")
  vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "end_pos2"] = vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "start_pos2"].fillna(0).astype("int") + vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "start_res2"].str.len() -1 
  vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "start_pos2"] = np.array([list(range(i, j+1)) for i, j in vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]),["start_pos2", "end_pos2"]].values],dtype = object)
  vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "start_res2"] = np.array([list(x) for x in vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "start_res2"]], dtype = object)
  vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "end_res2"] = np.array([list(x) for x in vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "end_res2"]], dtype = object)
  split_residues = vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]),["start_res2", "start_pos2", "end_res2"]].values.tolist()
  residues = []
  for item in split_residues:
    #allows iteration over lists to zip together, single sample 
    item[0] = item[0] if type(item[0]) is list else [item[0]]
    item[1] = item[1] if type(item[1]) is list else [item[1]]
    item[2] = item[2] if type(item[2]) is list else [item[2]]
    residue = '~'.join(list(map(''.join,list(x for x in zip(item[0], [str(i) for i in item[1]], item[2])))))
    residues.append(residue)
  vcf.loc[vcf["Annotation"].astype("str").isin(["missense_variant","synonymous_variant"]), "residues"] = residues
  vcf["residues"] = vcf["residues"].fillna("")
  cols = list(cols) + ["start_res", "start_pos", "end_res", "end_pos", "change", "ins", "ins_length", "length","inserted_residues","range", "start_res2", "start_pos2", "end_res2", "end_pos2"]
  vcf.drop(cols, axis = 1, inplace = True)
  #annotate residue specific information for each variant. 
  vcf = annotate_s_residues(vcf.copy(), data_dir)
  cols = ["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2", "VDS", "serum_escape", "mab_escape", "cm_mab_escape","mAb_escape_class_1","mAb_escape_class_2","mAb_escape_class_3","mAb_escape_class_4","BEC_RES", "BEC_EF"]
  vcf["SPEAR"] = vcf[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
  vcf.drop(list(set().union(snpeff_anno_cols, cols)), axis = 1, inplace = True)
  cols = [e for e in vcf.columns.to_list() if e not in ("SPEAR", "ANN2")]
  vcf = vcf.groupby(cols, as_index = False).agg(set)
  vcf.drop(["index", "ANN2"], axis = 1, inplace = True)
  #code block removes orf1ab from spear summaries where mat peptide products are described, otherwise retains the first annotation of position within polypeptide e.g. just orf1ab not orf1 
  vcf["SPEAR"] = vcf["SPEAR"].map(lambda a: filter_spear_summary(a))
  vcf["SPEAR"] = [','.join(map(str, l)) for l in vcf['SPEAR']]
  return vcf

def parse_vcf(vcf_file, split_info_cols = True):
  with open(vcf_file, 'r') as fobj:
    headiter = takewhile(lambda s: s.startswith('#'), fobj)
    header = list(headiter)

  header = [s.strip() for s in header]
  cols = header.pop(-1).split(sep="\t")
  df = pd.read_csv(vcf_file, comment="#", sep="\t", names=cols)
  if split_info_cols:
    info_list = df["INFO"].to_list()
    info_list = [dict(x.split("=") for x in item.split(";")) for item in info_list]
    info_df = pd.DataFrame(info_list)
    df = pd.concat([df, info_df], axis=1)
    cols = df.columns.to_list()
    infocols = info_df.columns.to_list()
    cols = [x for x in cols if x not in infocols] #would have used set subtraction here for a speedup but it loses the order of cols which needs to be retained
    cols[cols.index("INFO"):cols.index("INFO")] = infocols
    cols.pop(cols.index("INFO"))
    df = df[cols]
    return header,df, infocols
  else:
    return header, df

def write_vcf(header,body, output_filename):
  '''
  Function writes a vcf file. Two-step process writes header first, then appends a
  vcf file body parsed with parse_vcf using df.to_csv.
  '''
  with open(output_filename, 'w') as f:
      for item in header:
          f.write("%s\n" % item)
  body.to_csv(output_filename, mode='a', index = False, sep = "\t")

def main():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('output_filename', metavar='spear_vcfs/merged.spear.vcf', type=str,
      help='Filename for SPEAR annotated VCF') #ADD A DEFAULT FOR THIS 
  parser.add_argument('vcf', metavar='path/to/vcf', type = str,
      help='Input VCF file')
  parser.add_argument('--allow_seq_end', default=False, action='store_true',
      help = "Include variants in positions <= POS 55 and >= 29805") 
  parser.add_argument('data_dir', metavar='path/to/data/', type = str,
      help ='Data files for peptide subpositions')
  args = parser.parse_args()

  header, vcf, infocols = parse_vcf(args.vcf)
  if not args.allow_seq_end:
    vcf = vcf[vcf["POS"].between(56, 29804)].reset_index(drop = True)
  if len(vcf) != 0: #do not add summary if the vcf file is empty.
    df = vcf.iloc[:,:vcf.columns.get_loc("FORMAT")] # split vcf file columns up to ANN , could change this to LOC and up to format column to make more flexible ? 
    df = df.replace(np.nan, '', regex=True)
    samples = vcf.iloc[:,vcf.columns.get_loc("FORMAT"):] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
    header.append(f'##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'Gene | HGVS.c | Annotation | HGVS | Product | RefSeq_acc | residues | region | domain | contact_type | NAb | barns_class | bloom_ace2 | VDS | serum_escape |  mAb_escape | cm_mAb_escape | mAb_escape_class_1 | mAb_escape_class_2 | mAb_escape_class_3 | mAb_escape_class_4 | BEC_RES | BEC_EF  \'">') #MAKE VARIANT HEADER HGVS
    genbank = SeqIO.read(open(f'{args.data_dir}/NC_045512.2.gb',"r"), "genbank")
    genbank_mapping = {}
    locus_tag_mapping = {}
    for feature in genbank.features:
      if feature.type == "CDS" or feature.type == "mat_peptide":
        if feature.qualifiers["gene"][0] == "ORF1ab":
          genbank_mapping[feature.qualifiers["protein_id"][0]] = feature.qualifiers["product"][0]
          locus_tag_mapping[feature.qualifiers["product"][0]] = feature.qualifiers["protein_id"][0]
        else:
          genbank_mapping[feature.qualifiers["locus_tag"][0]] = feature.qualifiers["product"][0]
          locus_tag_mapping[feature.qualifiers["product"][0]] = feature.qualifiers["protein_id"][0]
    df = convert_snpeff_annotation(df.copy(), genbank_mapping, locus_tag_mapping, args.data_dir)
    infocols.append("SPEAR")
    df.loc[df["ANN"] == "", "ANN"] = "no_annotation"
    for col in infocols:
      df[col] = col + "=" + df[col]
    df['INFO'] = df[infocols].agg(';'.join, axis=1)
    df.drop(infocols, axis = 1, inplace = True)
    df['INFO'] = df['INFO'].str.replace('problem_exc=;','')
    df['INFO'] = df['INFO'].str.replace('problem_filter=;','')
    cols = df.columns.to_list()
    cols.pop(cols.index("INFO"))
    cols.insert(cols.index("FILTER") + 1, "INFO")
    df = df[cols]
    #sanity check before concat shape
    vcf = pd.concat([df, samples],axis=1)
    write_vcf(header,vcf,args.output_filename)
  else:
    copy(args.vcf, args.output_filename) #copy the empty vcf file so that snakemake sees command completion. 
  

if __name__ == "__main__":
    main()