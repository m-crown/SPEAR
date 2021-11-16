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
from go_propa import get_gb_coords
import glob

def convert_mp_position(row, coords):
  '''
  DEPRECATED
  Function converts a cds position to the position within a matured peptide (and returns the mat pepts name) for 
  variants within multi-product polypeptides.
  '''
  #converting the gene names to be lowercase when comparing to snpeff which uses lower - this means the gb_coords function remains same as go_propa.
  coords = {k.lower(): v for k, v in coords.items()}

  if row["Gene_Name"].lower() in coords:
    if "multi-product" in coords[row["Gene_Name"].lower()]:
      placeholder_names = []
      placeholder_pos = []
      for k,v in coords[row["Gene_Name"].lower()]["multi-product"].items():
        if v["sub-position"][0] <= int(row["pos"]) <= v["sub-position"][1]:
          placeholder_names.append(k)
          placeholder_pos.append(f'{row["ref_aa"]}{int(row["pos"]) - v["sub-position"][0]}{row["alt_aa"]}')
      if len(placeholder_names) > 1:
        name = "/".join(placeholder_names)
        if len(set(placeholder_pos)) == 1:
          pos = placeholder_pos[0]
      else:
        name = placeholder_names[0]
        pos = placeholder_pos[0]
      return name,pos
    else:
      return coords[row["Gene_Name"].lower()]["product"], row["oneletter"]
  else:
    return "intergenic", row["HGVS.c"]

def convert_snpeff_annotation(df_row, gb_mapping):
  snpeff_annotation = pd.DataFrame({"ANN" : [df_row["ANN"].split(',')]})
  snpeff_annotation = snpeff_annotation["ANN"].explode().copy().to_frame()
  snpeff_annotation.reset_index(inplace = True , drop=True)
  snpeff_annotation = snpeff_annotation["ANN"].str.split('|', expand = True)
  snpeff_annotation.columns = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"]
  snpeff_annotation.loc[snpeff_annotation["Feature_ID"] == "GU280_gp01","Feature_ID"]="YP_009724389.1"
  snpeff_annotation.loc[snpeff_annotation["Feature_ID"] == "GU280_gp01.2","Feature_ID"]="YP_009725295.1"
  snpeff_annotation = snpeff_annotation.drop(snpeff_annotation[(snpeff_annotation["Feature_ID"] == "YP_009724389.1") | (snpeff_annotation["Feature_ID"] == "YP_009725295.1")].index)
  snpeff_annotation.loc[snpeff_annotation["Annotation"].isin(["missense_variant", "synonymous_variant", "conservative_inframe_deletion", "disruptive_inframe_deletion"]),"variant"] = snpeff_annotation["HGVS.p"]
  snpeff_annotation.loc[snpeff_annotation["Annotation"].isin(["intergenic_region"]), "variant"] = snpeff_annotation["HGVS.c"]
  snpeff_annotation["product"] = snpeff_annotation["Feature_ID"].apply(lambda x: gb_mapping.get(x))
  snpeff_annotation.loc[snpeff_annotation["product"].isnull(), "product"] = snpeff_annotation.loc[snpeff_annotation["product"].isnull(), "Feature_ID"]
  cols = ["Annotation", "variant", "product"]
  snpeff_annotation["summary"] = snpeff_annotation[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
  summary = snpeff_annotation["summary"].to_list()
  summary = set(summary)
  summary = ','.join(summary)
  return (summary)

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
  print("Writing VCF file to: ", Path.joinpath(outdir,f'{basename}.{desc}.vcf'))
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
    header, df, infocols = parse_vcf(vcf)
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
      
      df["SPEAR"] = df.apply(lambda x : convert_snpeff_annotation(x, genbank_mapping), axis = 1)
      infocols.append("SPEAR")
      for col in infocols:
        df[col] = col + "=" + df[col]
      df['INFO'] = df[infocols].agg(';'.join, axis=1)
      df.drop(infocols, axis = 1, inplace = True)
      cols = df.columns.to_list()
      cols.pop(cols.index("INFO"))
      cols.insert(cols.index("FILTER") + 1, "INFO")
      df = df[cols]
      write_vcf(header,df,basename,"spear",outdir)

if __name__ == "__main__":
    main()