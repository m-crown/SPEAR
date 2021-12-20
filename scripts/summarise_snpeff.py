#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
from itertools import takewhile
import argparse
import numpy as np
import re
from shutil import copy

def convert_snpeff_annotation(vcf, gb_mapping, locus_tag_mapping, data_dir):

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
    bloom_anno_file = pd.read_csv(f'{data_dir}/single_mut_effects.csv').fillna("")
    annotation = []
    for residues in s_residues:
      residues = [item for sublist in residues for item in sublist.split(",")]
      new_anno = []
      for residue in residues:
        pattern = re.compile(r"[a-zA-Z]+([0-9]+)")
        respos = int(pattern.match(residue).group(1)) if pattern.match(residue) else -1 #it is an acknowledged limitation at this stage that this does not working for insertions
        if respos != -1:
          anno = spear_anno_file.loc[spear_anno_file["AA_coordinate"] == respos, ["region", "domain", "contact_type", "NAb", "barns_class"]].values.tolist()
          anno = [item if len(item) == len(item) else [] for item in anno]
          anno = [item for sublist in anno for item in sublist]
          bloom_anno = bloom_anno_file.loc[bloom_anno_file["mutation"] == residue, "bind_avg"].values.tolist()
          bloom_anno = ["" if len(bloom_anno) == 0 else bloom_anno[0]]
          anno += bloom_anno
          new_anno.append(anno)
        else:
          new_anno.append(["","","","","",""])
            
      new_anno = list(map(list, zip(*new_anno)))
      new_anno = ['~'.join([str(c) if c == c else "" for c in lst]) for lst in new_anno]
      annotation.append(new_anno)
    vcf.loc[vcf["product"] == "surface glycoprotein" , ["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]] = annotation
    vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]] = vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]].fillna("")
    vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]] = vcf[["region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]].replace("[^0-9a-zA-Z]+[~]+", "", regex = True)
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
  vcf.loc[vcf["product"].isnull(), "product"] = vcf.loc[vcf["product"].isnull(), "Feature_ID"]
  vcf["protein_id"] = vcf["product"].map(lambda x: locus_tag_mapping.get(x))
  vcf.loc[vcf["protein_id"].isnull(), "protein_id"] = vcf.loc[vcf["protein_id"].isnull(), "Feature_ID"]

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
  #frameshift variants
  vcf.loc[vcf["Annotation"].astype("str") == "frameshift_variant", "residues"] = vcf.loc[vcf["Annotation"].astype("str") == "frameshift_variant", "variant"].str.replace("^p\.", '')
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
  cols = ["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues","region", "domain", "contact_type", "NAb", "barns_class", "bloom_ace2"]
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
  parser.add_argument('data_dir', metavar='path/to/data/', type = str, 
      help ='Data files for peptide subpositions')
  args = parser.parse_args()

  header, vcf, infocols = parse_vcf(args.vcf)
  if len(vcf) != 0: #do not add summary if the vcf file is empty.
    df = vcf.iloc[:,:vcf.columns.get_loc("FORMAT")] # split vcf file columns up to ANN , could change this to LOC and up to format column to make more flexible ? 
    df = df.replace(np.nan, '', regex=True)
    samples = vcf.iloc[:,vcf.columns.get_loc("FORMAT"):] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
    header.append(f'##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'Gene | HGVS.c | Annotation | HGVS | Product | RefSeq_acc | Residues | Bloom | BIS \'">') #MAKE VARIANT HEADER HGVS
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