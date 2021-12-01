rule all:
   input: 
      expand(config["output_dir"] + '/final/{id}.spear.vcf', id=config["samples"]) if config["split_vcfs"] == True else config["output_dir"] + '/merged.spear.vcf'

rule split_vcfs:
   input:
      config["output_dir"] + "/merged.spear.vcf"
   output:
      config["output_dir"] + "/final/{id}.spear.vcf"
   shell:
      """for sample in `bcftools query -l {input}`; do bcftools view -Ov -s $sample -o {output} {input}; done"""

rule summarise_snpeff:
   input:
      config["output_dir"] + "/snpeff/merged.ann.vcf"
   output:
      config["output_dir"] + "/merged.spear.vcf"
   shell:
      "python3 summarise_snpeff.py {output} {input} {config[data]}"

rule snpeff:
   input:
      config["output_dir"] + "/indels/merged.vcf"
   output:
      config["output_dir"] + "/snpeff/merged.ann.vcf"
   shell:
      "java -Xmx8g -jar ~/snpEff/snpEff.jar -c ~/snpEff/snpEff.config -nodownload -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output}"

rule merge_vcfs:
   input:
      expand(config["output_dir"] + "/indels/{id}.indels.vcf.gz", id=config["samples"])
   output:
      config["output_dir"] + "/indels/merged.vcf"
   shell:
      "bcftools merge -m none -o {output} {input}"

rule index_vcfs:
   input:  
      config["output_dir"] + "/indels/{id}.indels.vcf"
   output:
      config["output_dir"] + "/indels/{id}.indels.vcf.gz"
   shell:
      "bgzip {input} && tabix -p vcf {output}"

rule get_indels:
   input:
      vcf_file = config["output_dir"] + "/masked/{id}.masked.vcf" if config["filter"] == True else config["output_dir"] + "/fatovcf/{id}.vcf",
      muscle_aln = config["output_dir"] + "/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}.muscle.aln"
   output:
      snps_indels = config["output_dir"] + "/indels/{id}.indels.vcf",
      snps_indels_tsv = config["output_dir"] + "/indels/{id}.indels.tsv"
   shell:
       "python3 {config[get_indels]} --vcf_file {input.vcf_file} --deletion_window 2 --write_indels {output.snps_indels_tsv} {input.muscle_aln} MN908947.3 {output.snps_indels}"

rule filter_problem_sites:
   input: 
      config["output_dir"] + "/masked/{id}.problem.vcf"
   output: 
      config["output_dir"] + "/masked/{id}.masked.vcf"
   shell:
      "java -Xmx8g -jar ~/snpEff/SnpSift.jar filter \"!( {config[filter_params]} )\" {input} > {output}"

rule mark_problem_sites:
   input:
      config["output_dir"] + "/fatovcf/{id}.vcf" 
   output:
      config["output_dir"] + "/masked/{id}.problem.vcf"
   shell:
      "~/GitHub/SPEAR/vcfanno_osx ~/GitHub/SPEAR/conf.toml {input} > {output}"

rule get_snps:
   input:
      config["output_dir"] + "/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}.muscle.aln"
   output:
      snps = config["output_dir"] + "/fatovcf/{id}.vcf"
   shell:
      "{config[fatovcf]} {config[exclude_ambiguous]} {input} {output.snps} ; ./update_vcf_header.sh {output.snps}; sed -i '' 's/MN908947.3/NC_045512.2/g' {output.snps}"


if config["align"]:
   rule muscle_alignment:
      input: config["input_dir"] + "/{id}.consensus.fa"
      output: 
         plus_ref =  config["output_dir"] + "/consensus/{id}.consensus_ref.fa",
         alignment = config["output_dir"] + "/muscle/{id}.muscle.aln"
      shell:
         "cat {config[reference_sequence]} {input} > {output.plus_ref} ; muscle -quiet -in {output.plus_ref} -out {output.alignment}"


# if config["input_dir"] == "/Users/matthewcrown/Desktop/muscle/":
#    print("OI")
# else:
#    print("OIOI")
#muscle -align sequences.fasta -output alignment.fasta -quiet  muscle.aln
