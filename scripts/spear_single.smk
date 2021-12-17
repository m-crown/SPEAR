rule all:
   input: 
      expand(config["output_dir"] + "/summary/{id}.summary.csv", id = config["samples"])

rule summarise_vcfs:
   input:
      expand(config["output_dir"] + "/final_vcfs/{id}.vcf" , id=config["samples"])
   output:
      expand(config["output_dir"] + "/summary/{id}.summary.csv", id = config["samples"])
   shell:
      """convert_format.py {config[output_dir]}/summary/ --vcf {input}"""

rule spear:
   input:
      config["output_dir"] + "/snpeff/{id}.ann.vcf"
   output:
      config["output_dir"] + "/final_vcfs/{id}.vcf" 
   shell:
      "summarise_snpeff.py {output} {input} {config[data_dir]}"

if config["vcf"] == True:
   rule snpeff:
      input:
         config["output_dir"] + "/masked/{id}.masked.vcf" if config["filter"] else config["input_dir"] + "/{id}" + config["extension"]
      output:
         config["output_dir"] + "/snpeff/{id}.ann.vcf" if config["filter"] else config["output_dir"] + "/snpeff/{id}.ann.vcf"
      shell:
         """
         java -Xmx8g -jar $CONDA_PREFIX/bin/snpEff.jar -hgvs1LetterAa -download -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output}
         """

if config["vcf"] != True:
   rule snpeff:
      input:
         config["output_dir"] + "/indels/{id}.indels.vcf"
      output:
         config["output_dir"] + "/snpeff/{id}.ann.vcf"
      shell:
         """
         java -Xmx8g -jar $CONDA_PREFIX/bin/snpEff.jar -hgvs1LetterAa -download -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output}
         """

rule get_indels:
   input:
      vcf_file = config["output_dir"] + "/masked/{id}.masked.vcf" if config["filter"] == True else config["output_dir"] + "/fatovcf/{id}.vcf",
      muscle_aln = config["output_dir"] + "/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
   output:
      snps_indels = config["output_dir"] + "/indels/{id}.indels.vcf",
      snps_indels_tsv = config["output_dir"] + "/indels/{id}.indels.tsv"
   shell:
      "get_indels.py --vcf {input.vcf_file} --window {config[del_window]} {config[allow_ambiguous]} --tsv {output.snps_indels_tsv} {input.muscle_aln} MN908947.3 {output.snps_indels}"

rule filter_problem_sites:
   input: 
      config["output_dir"] + "/masked/{id}.problem.vcf"
   output: 
      config["output_dir"] + "/masked/{id}.masked.vcf"
   shell:
      "java -jar $CONDA_PREFIX/SnpSift.jar filter \"!( {config[filter_params]} )\" {input} > {output}"

rule mark_problem_sites:
   input:
      config["input_dir"] + "/{id}" + config["extension"] if config["vcf"] == True else config["output_dir"] + "/fatovcf/{id}.vcf" 
   output:
      config["output_dir"] + "/masked/{id}.problem.vcf"
   shell:
      "vcfanno $CONDA_PREFIX/scripts/conf.toml {input} > {output}"

rule get_snps:
   input:
      config["output_dir"] + "/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
   output:
      snps = config["output_dir"] + "/fatovcf/{id}.vcf"
   shell:
      "faToVcf {config[exclude_ambiguous]} {input} {output.snps} ; update_vcf_header.sh {output.snps}; sed -i '' 's/MN908947.3/NC_045512.2/g' {output.snps}"

if config["align"]:
   rule muscle_alignment:
      input: config["input_dir"] + "/{id}" + config["extension"]
      output: 
         plus_ref =  config["output_dir"] + "/consensus/{id}.consensus_ref.fa",
         alignment = config["output_dir"] + "/muscle/{id}.muscle.aln"
      shell:
         "cat {config[reference_sequence]} {input} > {output.plus_ref} ; muscle -quiet -in {output.plus_ref} -out {output.alignment}"