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

def convert_snpeff_annotation(vcf, gb_mapping, locus_tag_mapping, respos_df):

  def filter_spear_summary(summaries):
    regexp = re.compile(r'ORF1(a|ab) polyprotein')
    if len(list(filter(None, list(regexp.search(summary) for summary in summaries)))) == len(summaries):
      summaries = [list(summaries)[0]]
    else:
      summaries = [summary for summary in summaries if not regexp.search(summary)]
    return summaries

  def get_ref_res(refres_info, respos, gene):
    if respos == respos:
        return(refres_info.loc[refres_info.index == int(respos), gene].values[0])
    else:
        return np.nan

  def populate_alt_res(variant, pos, respos):
    #extract the column residue position from HGVS.p
    #this is only required for variants matching delins, del, missense or synoynmous variants (e.g. not frameshift, stop or insertion - these are handled directly.)
    pos = int(pos)
    if re.match('p\.[A-Z][0-9]+_[A-Z][0-9]+del$', variant):
        max_res = int(re.match('p\.[A-Z][0-9]+_[A-Z]([0-9]+)+del$', variant)[1])
        #only give a deletion to positions which are within the deletion range
        if respos <= max_res:
            alt_res = "del" 
        else:
            alt_res = ""
    elif re.match('p\.[A-Z]+[0-9]+[A-Z\*]+|fs', variant):
        alt_residues = re.match('p\.[A-Z]+[0-9]+([A-Z\*]+|fs)', variant)[1]
        if alt_residues != "fs":
            alt_res_list = [res for res in alt_residues]
        else:
            alt_res = "fs"
        #try to get the alt_res from list and if position is out of range leave alt res empty (i.e. doesnt exist)
        try:
            alt_res = alt_res_list[pos]
        except:
            alt_res = ""           
    elif re.match('p\.[A-Z][0-9]+_[A-Z][0-9]+delins[A-Z\*]+?',variant):
        max_res = int(re.match('p\.[A-Z][0-9]+_[A-Z]([0-9]+)delins[A-Z\*]+', variant)[1])
        alt_residues = re.match('p\.[A-Z][0-9]+_[A-Z][0-9]+delins([A-Z\*]+)', variant)[1]
        alt_res_list = [res for res in alt_residues]
        if respos <= max_res:
            try:
                alt_res = alt_res_list[pos]
            except:
                alt_res = "del"
        else:
            alt_res = ""
    else:
        alt_res = ""      
    return alt_res


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
  vcf["protein_id"] = vcf["protein_id"]

  #convert HGVS annotation to per residue for residue annotation. 
  #handle deletions and insertions residues first
  vcf[["start_res", "start_pos", "end_res", "end_pos", "change", "ins"]] = vcf["variant"].str.extract('([A-Z])([0-9]+)_*([A-Z])?([0-9]+)?([a-z]+)([A-Z]+)?')
  vcf["end_pos"] = vcf["end_pos"].fillna(vcf["start_pos"])
  vcf[["start_pos", "end_pos"]] = vcf[["start_pos", "end_pos"]].fillna(0).astype("int")
  vcf["ins_length"] = vcf['ins'].str.len().fillna(0).astype('int')

  #Because multiple substitutions is incorrectly handled range_df wont necessarily expand to correct size, so need to set the minimum size of the dataframe to be the size of these multiple insertions, then expand to del/delins
  maxmultiaa = vcf.loc[vcf.variant.str.contains(r'p\.[A-Z]+[0-9]+[A-Z\*]{2,}') == True, "variant"].str.extract(r'p\.[A-Z]+[0-9]+([A-Z\*]{2,})', expand = False).str.len().max().astype("int")
  range_df= pd.DataFrame(list(range(i, j+1)) for i, j in vcf[["start_pos", "end_pos"]].values)
  additional_cols = [x for x in list(range(0, maxmultiaa)) if x not in range_df.columns.tolist()]
  range_df[additional_cols] = np.nan
  range_df = range_df.add_prefix("respos_")
  
  range_cols = range_df.columns.tolist()
  all_cols = []
  refres_cols = []
  altres_cols = []
  final_res_cols = []
  for respos in range_cols:
      refres_col = "ref_res" + re.match('respos(_[0-9]+)', respos)[1]
      altres_col = "alt_res" + re.match('respos(_[0-9]+)', respos)[1]
      final_res_col = "residue" + re.match('respos(_[0-9]+)', respos)[1]
      all_cols.append(refres_col)
      refres_cols.append(refres_col)
      altres_cols.append(altres_col)
      all_cols.append(respos)
      all_cols.append(altres_col)
      all_cols.append(final_res_col)
      final_res_cols.append(final_res_col)
  range_df[[refres_cols, altres_cols, final_res_cols]] = np.nan
  range_df = range_df[all_cols]
  range_df[["variant", "product"]] = vcf[["variant", "product"]] 
  residue_col_list = list(range(1, range_df.shape[1] - 2, 4))

  range_df.loc[range_df.variant.str.contains(r'p\.[A-Z]+[0-9]+[A-Z\*]+') == True, "respos_0"] = range_df.loc[range_df.variant.str.contains(r'p\.[A-Z]+[0-9]+[A-Z\*]+') == True, "variant"].str.extract(r'p\.[A-Z]+([0-9]+)[A-Z\*]+')[0]
  
  for x in range(1,len(residue_col_list)):
      range_df.iloc[:, residue_col_list[x]] = range_df.iloc[:, residue_col_list[x-1]].fillna(-2).astype(int) + 1

  for ref, pos, alt, final in zip(refres_cols, range_cols, altres_cols, final_res_cols):
      col_pos = re.search('([0-9]+)', pos)[1]
      range_df.loc[range_df["variant"].str.match(r'^p\.'), ref] = range_df.loc[range_df["variant"].str.match(r'^p\.')].apply(lambda x: get_ref_res(respos_df , x[pos], x["product"]), axis=1)
      range_df.loc[range_df["variant"].str.match(r'^p\.'), alt] = range_df.loc[range_df["variant"].str.match(r'^p\.')].apply(lambda x: populate_alt_res(x["variant"], col_pos, x[pos]), axis=1).fillna("")
      range_df.loc[range_df[alt].isin(["", np.nan]) == False , final] = range_df.loc[range_df[alt].isin(["", np.nan]) == False, ref].astype("str") + range_df.loc[range_df[alt].isin(["", np.nan]) == False , pos].astype("str") + range_df.loc[range_df[alt].isin(["", np.nan]) == False ,alt].astype("str")
  
  range_df["residues"] = np.nan
  range_df = range_df.fillna("")
  range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$', regex = True), "residues"] = range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$', regex = True), "variant"].str.extract(r'p\.([A-Z][0-9]+)_[A-Z][0-9]+ins[A-Z]+$', expand = False).astype("str") + "-" + range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$', regex = True), "variant"].str.extract(r'p\.[A-Z][0-9]+_[A-Z][0-9]+ins([A-Z]+)$', expand = False).astype("str") + "-" + range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$', regex = True), "variant"].str.extract(r'p\.[A-Z][0-9]+_([A-Z][0-9]+)ins[A-Z]+$', expand = False).astype("str")
  range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$', regex = True) == False, "residues"] = range_df[final_res_cols].agg('~'.join, axis=1)
  range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+fs$', regex = True), "residues"] = range_df.loc[range_df["variant"].str.contains('p\.[A-Z][0-9]+fs$', regex = True), "variant"].str.extract(r'p\.([A-Z][0-9]+fs)$', expand = False).astype("str")
  range_df["residues"] = range_df["residues"].str.replace('~+$', '', regex = True)
  vcf = pd.concat([vcf,range_df["residues"]], axis = 1)
  cols = ["start_res", "start_pos", "end_res", "end_pos", "change", "ins", "ins_length"]
  vcf.drop(cols, axis = 1, inplace = True)

  #annotate residue specific information for each variant.
  cols = ["Gene_Name", "HGVS.c", "Annotation", "variant", "product", "protein_id", "residues"]
  vcf["SUM"] = vcf[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
  vcf.drop(list(set().union(snpeff_anno_cols, cols)), axis = 1, inplace = True)
  cols = [e for e in vcf.columns.to_list() if e not in ("SUM", "ANN2")]
  vcf = vcf.groupby(cols, as_index = False).agg(set)
  vcf.drop(["index", "ANN2"], axis = 1, inplace = True)
  #code block removes orf1ab from spear summaries where mat peptide products are described, otherwise retains the first annotation of position within polypeptide e.g. just orf1ab not orf1 
  vcf["SUM"] = vcf["SUM"].map(lambda a: filter_spear_summary(a))
  vcf["SUM"] = [','.join(map(str, l)) for l in vcf['SUM']]
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
      help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
  parser.add_argument('vcf', metavar='path/to/vcf', type = str,
      help='Input VCF file')
  parser.add_argument('--allow_seq_end', default=False, action='store_true',
      help = "Include variants in positions <= POS 55 and >= 29805") 
  parser.add_argument('data_dir', metavar='path/to/data/', type = str,
      help ='Data files for peptide subpositions')
  args = parser.parse_args()

  respos_df = pd.read_pickle(f'{args.data_dir}/respos.pkl')
  header, vcf, infocols = parse_vcf(args.vcf)
  if not args.allow_seq_end:
    vcf = vcf[vcf["POS"].between(56, 29804)].reset_index(drop = True)
  if len(vcf) != 0: #do not add summary if the vcf file is empty.
    df = vcf.iloc[:,:vcf.columns.get_loc("FORMAT")] # split vcf file columns up to ANN , could change this to LOC and up to format column to make more flexible ? 
    df = df.replace(np.nan, '', regex=True)
    samples = vcf.iloc[:,vcf.columns.get_loc("FORMAT"):] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
    header.append(f'##INFO=<ID=SUM,Number=.,Type=String,Description="SnpEff Summary and per residue list: \'Gene | HGVS.c | Annotation | HGVS | Product | RefSeq_acc | residues \'">') #MAKE VARIANT HEADER HGVS
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
    df = convert_snpeff_annotation(df.copy(), genbank_mapping, locus_tag_mapping, respos_df)
    infocols.append("SUM")
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