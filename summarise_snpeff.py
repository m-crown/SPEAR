#!/usr/bin/env python

from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
from itertools import takewhile
import argparse
from pathlib import Path
import glob
import numpy as np
import re

def convert_snpeff_annotation(vcf, gb_mapping, locus_tag_mapping):

def convert_snpeff_annotation(vcf, gb_mapping):
  #takes a input a dataframe row, splits the ann field into a new
  vcf["ANN2"] = vcf["ANN"].str.split(',') #put the split ANN field into a new column to preserve original snpeff value
  vcf = vcf.explode("ANN2")
  snpeff_anno_cols = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"]
  vcf[snpeff_anno_cols] = vcf["ANN2"].str.split('|', expand = True)
  vcf.loc[vcf["Feature_ID"] == "GU280_gp01","Feature_ID"]="YP_009724389.1"
  vcf.loc[vcf["Feature_ID"] == "GU280_gp01.2","Feature_ID"]="YP_009725295.1"
  vcf.reset_index(inplace = True)
  vcf = vcf.drop(vcf[(vcf["Feature_ID"] == "YP_009724389.1") | (vcf["Feature_ID"] == "YP_009725295.1")].index)
  vcf.loc[vcf["Annotation"] != "intergenic_region","variant"] = vcf["HGVS.p"]
  vcf.loc[vcf["Annotation"] == "intergenic_region", "variant"] = vcf["HGVS.c"]
  vcf["product"] = vcf["Feature_ID"].apply(lambda x: gb_mapping.get(x))
  vcf.loc[vcf["product"].isnull(), "product"] = vcf.loc[vcf["product"].isnull(), "Feature_ID"]
  vcf["protein_id"] = vcf["product"].map(lambda x: locus_tag_mapping.get(x))
  vcf.loc[vcf["protein_id"].isnull(), "protein_id"] = vcf.loc[vcf["protein_id"].isnull(), "Feature_ID"]
  cols = ["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id"]
  vcf["SPEAR"] = vcf[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
  vcf.drop(list(set().union(snpeff_anno_cols, cols)), axis = 1, inplace = True)
  cols = [e for e in vcf.columns.to_list() if e not in ("SPEAR", "ANN2")]
  vcf = vcf.groupby(cols, as_index = False).agg(set)
  vcf.drop(["index", "ANN2"], axis = 1, inplace = True)
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
  parser.add_argument('vcf', metavar='path/to/vcfs', type = str,
      help='Input VCF file')
  parser.add_argument('data_dir', metavar='path/to/genbank+genpepts', type = str, 
      help ='Data files for peptide subpositions')
  args = parser.parse_args()

  header, vcf, infocols = parse_vcf(args.vcf)
  df = vcf.iloc[:,:vcf.columns.get_loc("FORMAT")] # split vcf file columns up to ANN , could change this to LOC and up to format column to make more flexible ? 
  df = df.replace(np.nan, '', regex=True)
  samples = vcf.iloc[:,vcf.columns.get_loc("FORMAT"):] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
  header.append(f'##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'Gene | HGVS.c | Annotation | Variant | Product | ID\'">') #MAKE VARIANT HEADER HGVS
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
  df = convert_snpeff_annotation(df.copy(), genbank_mapping, locus_tag_mapping)
  df["SPEAR"].to_csv("testing_spear_out.csv")
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
  vcf = pd.concat([df, samples],axis=1)
  write_vcf(header,vcf,args.output_filename)

if __name__ == "__main__":
    main()