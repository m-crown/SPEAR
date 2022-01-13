rule all:
   input: 
      expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"])

rule summarise_vcfs:
   input:
      expand(config["output_dir"] + "/final_vcfs/{id}.spear.vcf" , id=config["samples"])
   output:
      expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"])
   log: config["output_dir"] + "/intermediate_output/logs/summarise/summary.log"
   shell:
      """convert_format.py {config[output_dir]} {config[data_dir]} --vcf {input} 2> {log}"""

rule split_vcfs:
   input:
      config["output_dir"] + "/all_samples.spear.vcf"
   output:
      config["output_dir"] + "/final_vcfs/{id}.spear.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/split/{id}.split.log"
   shell:
      """
      bcftools view -Ov -c 1 -s {wildcards.id} -o {config[output_dir]}/final_vcfs/{wildcards.id}.spear.vcf {input} 2> {log}
      """

rule spear:
   input:
      config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf"
   output:
      config["output_dir"] + "/all_samples.spear.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/spear/spear.log"
   shell:
      "summarise_snpeff.py {output} {input} {config[data_dir]} ; spear_annotate.py {output} {output} {config[data_dir]} 2> {log}"

rule snpeff:
   input:
      config["output_dir"] + "/intermediate_output/merged.vcf" if config["filter"] else config["output_dir"] + "/intermediate_output/merged.vcf"
   output:
      config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf" if config["filter"] else config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/snpeff/snpeff.log"
   shell:
      """
      java -Xmx2g -jar $CONDA_PREFIX/snpEff/snpEff.jar -noShiftHgvs -hgvs1LetterAa -download -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output} 2> {log}
      """

if config["filter"] == True:
   rule merge_vcfs:
      input:
         expand(config["output_dir"] + "/intermediate_output/masked/{id}.masked.vcf", id=config["samples"]) if config["vcf"] else expand(config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf", id=config["samples"])
      output:
         config["output_dir"] + "/intermediate_output/merged.vcf"
      log: config["output_dir"] + "/intermediate_output/logs/merge/merge.log"
      shell:
         "bcftools merge --no-index -m none -o {output} {input} 2> {log}"

if config["filter"] != True:
   rule merge_vcfs:
      input:
         expand(config["input_dir"] + "/{id}" + config["extension"], id=config["samples"]) if config["vcf"] else expand(config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf", id=config["samples"])
      output:
         config["output_dir"] + "/intermediate_output/merged.vcf"
      log: config["output_dir"] + "/intermediate_output/logs/merge/merge.log"
      shell:
         "bcftools merge --no-index -m none -o {output} {input} 2> {log}"

rule get_indels:
   input:
      vcf_file = config["output_dir"] + "/intermediate_output/masked/{id}.masked.vcf" if config["filter"] == True else config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf",
      muscle_aln = config["output_dir"] + "/intermediate_output/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
   output:
      snps_indels = config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf",
      snps_indels_tsv = config["output_dir"] + "/intermediate_output/indels/{id}.indels.tsv"
   log: config["output_dir"] + "/intermediate_output/logs/indels/{id}.indels.log"
   shell:
      "get_indels.py --vcf {input.vcf_file} --window {config[del_window]} {config[allow_ambiguous]} --tsv {output.snps_indels_tsv} {input.muscle_aln} NC_045512.2 {output.snps_indels} 2> {log}"

rule filter_problem_sites:
   input: 
      config["output_dir"] + "/intermediate_output/masked/{id}.problem.vcf"
   output: 
      config["output_dir"] + "/intermediate_output/masked/{id}.masked.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/filter_problem_sites/{id}.filter_problem_sites.log"
   shell:
      "java -Xmx2g -jar $CONDA_PREFIX/snpEff/SnpSift.jar filter \"!( {config[filter_params]} )\" {input} > {output} 2> {log}"

rule mark_problem_sites:
   input:
      config["input_dir"] + "/{id}" + config["extension"] if config["vcf"] == True else config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf" 
   output:
      config["output_dir"] + "/intermediate_output/masked/{id}.problem.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/mark_problem_sites/{id}.mark_problem_sites.log"
   shell:
      "vcfanno --base-path $CONDA_PREFIX $CONDA_PREFIX/bin/conf.toml {input} > {output} 2> {log}"

rule get_snps:
   input:
      config["output_dir"] + "/intermediate_output/muscle/{id}.muscle.aln" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
   output:
      snps = config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/get_snps/{id}.get_snps.log"
   shell:
      "faToVcf {config[exclude_ambiguous]} {input} {output.snps} ; update_vcf_header.sh {output.snps} 2> {log}"

if config["align"]:
   rule muscle_alignment:
      input:
         config["input_dir"] + "/{id}" + config["extension"]
      output: 
         plus_ref =  config["output_dir"] + "/intermediate_output/consensus/{id}.consensus_ref.fa",
         alignment = config["output_dir"] + "/intermediate_output/muscle/{id}.muscle.aln" 
      log: config["output_dir"] + "/intermediate_output/logs/align/{id}.align.log"
      shell:
         "cat {config[reference_sequence]} {input} > {output.plus_ref} ; muscle -quiet -in {output.plus_ref} -out {output.alignment} 2> {log}"