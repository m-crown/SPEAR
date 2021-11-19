#!/usr/bin/env python

#running fatovcf
#for I in /home/covid19/Run_NXT59/results/ncovIllumina_Genotyping_alignSeqs/msa/*.muscle.aln; do J=$( basename $I .muscle.aln ); /home/covid19/faToVcf/faToVcf ${I} ${J}.vcf; done

#command to fix contig warning awk '/^#CHROM/ {printf("##contig=<ID=MN908947.3,length=29903,md5=d11d06b5d1eb1d85c69e341c3c026e08,URL=https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta>\n");} {print;}' test.vcf > test2.vcf

#snpeff command we went with
#for I in *.vcf; do J=$( basename $I .vcf ); snpEff -no-downstream -no-intron -no-upstream -no-utr  MN908947.3 ${I} > ${J}.ann.vcf; done


from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
from itertools import takewhile
import argparse
from pathlib import Path
import glob

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
  cols = ["Annotation", "variant", "product"]
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

def write_vcf(header,body, basename,desc,outdir):
  '''
  Function writes a vcf file. Two-step process writes header first, then appends a
  vcf file body parsed with parse_vcf using df.to_csv.
  '''
  with open(Path.joinpath(outdir,f'{basename}.{desc}.vcf'), 'w') as f:
      for item in header:
          f.write("%s\n" % item)
  body.to_csv(Path.joinpath(outdir,f'{basename}.{desc}.vcf'), mode='a', index = False, sep = "\t")

def main():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('output_dir', metavar='spear_vcfs/', type=str,
      help='Destination dir for SPEAR annotated VCFs')
  parser.add_argument('vcfs', metavar='path/to/vcfs', type = str,
      help='Input VCF file directory')
  parser.add_argument('data_dir', metavar='path/to/genbank+genpepts', type = str, 
      help ='Data files for peptide subpositions')
  args = parser.parse_args()

  vcfs =glob.glob(f'{args.vcfs}/*.vcf')
  outdir = Path(args.output_dir)
  #check if output directory exits, if not make it.
  outdir.mkdir(parents=True, exist_ok=True)
  print("Annotating VCFs")
  for vcf in vcfs:
    basename = Path(vcf).stem.split('.')[0]
    header, vcf, infocols = parse_vcf(vcf)
    df = vcf.iloc[:,:10] # split vcf file columns up to ANN 
    samples = vcf.iloc[:,10:] #split format and sample columns into separate dataframe to prevent fragmentation whilst annotating
    if df.empty:
      continue
    else:
      header.append(f'"##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'SnpEff Annotation | Position | Product\'">')
      genbank = SeqIO.read(open(f'{args.data_dir}/NC_045512.2.gb',"r"), "genbank")
      genbank_mapping = {}
      for feature in genbank.features:
        if feature.type == "CDS" or feature.type == "mat_peptide":
          if feature.qualifiers["gene"][0] == "ORF1ab":
            genbank_mapping[feature.qualifiers["protein_id"][0]] = feature.qualifiers["product"][0]
          else:
            genbank_mapping[feature.qualifiers["locus_tag"][0]] = feature.qualifiers["product"][0]
      df = convert_snpeff_annotation(df.copy(), genbank_mapping)
      infocols.append("SPEAR")
      for col in infocols:
        df[col] = col + "=" + df[col]
      df['INFO'] = df[infocols].agg(';'.join, axis=1)
      df.drop(infocols, axis = 1, inplace = True)
      cols = df.columns.to_list()
      cols.pop(cols.index("INFO"))
      cols.insert(cols.index("FILTER") + 1, "INFO")
      df = df[cols]
      vcf = pd.concat([df, samples],axis=1)
      csv_vcf = df.copy()
      csv_vcf["#CHROM"] = basename
      csv_vcf.to_csv(Path.joinpath(outdir,f'{basename}.spear.tsv'), mode='w', index = False, sep = "\t") #this file is only really necessary for comparison exercise
      write_vcf(header,vcf,basename,"spear",outdir)

if __name__ == "__main__":
    main()