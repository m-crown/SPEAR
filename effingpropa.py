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

def write_vcf(header,body, basename,outdir):
  '''
  Function writes a vcf file. Two-step process writes header first, then appends a
  vcf file body parsed with parse_vcf using df.to_csv.
  '''
  print("Writing VCF file to: ", Path.joinpath(outdir,f'{basename}.spear.vcf'))
  with open(Path.joinpath(outdir,f'{basename}.spear.vcf'), 'w') as f:
      for item in header:
          f.write("%s\n" % item)
  body.to_csv(Path.joinpath(outdir,f'{basename}.spear.vcf'), mode='a', index = False, sep = "\t")

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

  coords = get_gb_coords(f'{args.data_dir}/NC_045512.2.gb', [f'{args.data_dir}/ORF1a.gp', f'{args.data_dir}/ORF1ab.gp'])

  for vcf in vcfs:
    basename = Path(vcf).stem.split('.')[0]
    header, df, infocols = parse_vcf(vcf)
    header.append(f'"##INFO=<ID=SPEAR,Number=.,Type=String,Description="SPEAR Tool Annotations: \'Position | Product | Peptide Position\'">')
    #find a way to get these header names automatically from the header list
    df[["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"]] = df['ANN'].str.split('|', expand=True)
    df["ref_aa"] = df["HGVS.p"].str[2:6].apply(seq1)
    df["pos"] = df["HGVS.p"].str.extract('(\d+)')
    df["alt_aa"] = df["HGVS.p"].str[-3:].apply(seq1)
    df["oneletter"] = df["ref_aa"] + df["pos"] + df["alt_aa"]

    df[["product", "peptide-position"]] = df.apply(lambda row: convert_mp_position(row,coords), axis=1, result_type='expand')
    df.loc[df["product"] != "intergenic", "SPEAR"] = df["oneletter"] + "|" + df["product"] + "|" + df["peptide-position"]
    df.loc[df["product"] == "intergenic", "SPEAR"] = df["HGVS.c"] + "|" + df["product"] + "|" + df["peptide-position"]
    infocols.append("SPEAR")
    for col in infocols:
      df[col] = col + "=" + df[col]
    
    df['INFO'] = df[infocols].agg(';'.join, axis=1)
    #find a way to get these header names automatically from the header list
    df = df.drop(["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO", "ref_aa", "pos", "alt_aa", "oneletter", "product", "peptide-position", "SPEAR", "AC", "AN", "ANN"], axis = 1)
    cols = df.columns.to_list()
    cols.pop(cols.index("INFO"))
    cols.insert(cols.index("FILTER") + 1, "INFO")

    df = df[cols]

    write_vcf(header,df,basename,outdir)
    
if __name__ == "__main__":
    main()