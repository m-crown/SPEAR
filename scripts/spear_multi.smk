if config["report"] == False:
   rule all:
      input: 
         expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"])
else:
   rule all:
      input:
         samples = expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"]),
         report = config["output_dir"] + "/report/report.html"

rule produce_report:
   input:
      summary = config["output_dir"] + "/spear_score_summary.tsv",
      all_samples = config["output_dir"] + "/spear_annotation_summary.tsv",
      n_perc = config["output_dir"] + "/qc.csv" if config["vcf"] == False else config["output_dir"] + "/spear_score_summary.tsv"
   output:
      config["output_dir"] + "/report/report.html"
   log: config["output_dir"] + "/intermediate_output/logs/report/report.log"
   shell:
      """summary_report.py --n_perc {input.n_perc} {config[product_plots]} {input.summary} {input.all_samples} {config[baseline_scores]} {config[input_sample_num]} {config[qc_sample_num]} {config[images_dir]} {config[scripts_dir]} {config[data_dir]} {config[output_dir]}/report/ {config[baseline]} {config[global_n]} {config[s_n]} {config[s_contig]} {config[rbd_n]} 2> {log}"""

rule summarise_vcfs:
   input:
      expand(config["output_dir"] + "/final_vcfs/{id}.spear.vcf" , id=config["samples"])
   output:
      expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"]),
      config["output_dir"] + "/spear_annotation_summary.tsv",
      config["output_dir"] + "/spear_score_summary.tsv"
   log: config["output_dir"] + "/intermediate_output/logs/summarise/summary.log"
   shell:
      """for i in $(grep -lP '^NC_045512\.2' {config[output_dir]}/final_vcfs/*.spear.vcf); do BNAME=${{i##*/}}; BNAME2=${{BNAME%.spear.vcf}}; grep -v '^#' "$i" | sed "s|^NC_045512\.2|$BNAME2|g"; done > {config[output_dir]}/intermediate_output/anno_concat.tsv ; convert_format.py --is_vcf_input {config[vcf]} --is_filtered {config[filter]} {config[output_dir]}/intermediate_output/anno_concat.tsv {config[output_dir]} {config[data_dir]} {config[output_dir]}/all_samples.spear.vcf {config[samples]} 2> {log}"""
  

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

rule spear_annotate:
   input:
      config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf"
   output:
      config["output_dir"] + "/all_samples.spear.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/spear/spear.log"
   shell:
      "summarise_snpeff.py {output} {input} {config[data_dir]} ; spear_annotate.py {output} {output} {config[data_dir]} 2> {log}"

rule annotate_variants:
   input:
      config["output_dir"] + "/intermediate_output/masked/merged.masked.vcf" if config["filter"] else config["output_dir"] + "/intermediate_output/masked/merged.problem.vcf"
   output:
      config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/snpeff/snpeff.log"
   shell:
      """
      java -Xmx2g -jar $CONDA_PREFIX/snpEff/snpEff.jar -noShiftHgvs -hgvs1LetterAa -download -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output} 2> {log}
      """

if config["filter"]:
   rule filter_problem_sites:
      input: 
         config["output_dir"] + "/intermediate_output/masked/merged.problem.vcf"
      output: 
         config["output_dir"] + "/intermediate_output/masked/merged.masked.vcf"
      log: config["output_dir"] + "/intermediate_output/logs/filter_problem_sites/filter_problem_sites.log"
      shell:
         "java -Xmx2g -jar $CONDA_PREFIX/snpEff/SnpSift.jar filter \"!( {config[filter_params]} )\" {input} > {output} 2> {log}"

rule mark_problem_sites:
   input:
      config["output_dir"] + "/intermediate_output/merged.vcf" 
   output:
      config["output_dir"] + "/intermediate_output/masked/merged.problem.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/mark_problem_sites/mark_problem_sites.log"
   shell:
      "vcfanno --base-path $CONDA_PREFIX $CONDA_PREFIX/bin/conf.toml {input} > {output} 2> {log}"

rule merge_vcfs:
   input:
      expand(config["input_dir"] + "/{id}" + config["extension"], id=config["samples"]) if config["vcf"] else expand(config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf", id=config["samples"])
   output:
      config["output_dir"] + "/intermediate_output/merged.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/merge/merge.log"
   shell:
      "bcftools merge --no-index -m none -o {output} {input} 2> {log}"

rule merge_qc:
   input:
      expand(config["output_dir"] + "/intermediate_output/indels/{id}.nperc.csv", id=config["samples"])
   output:
      qc_file = config["output_dir"] + "/qc.csv"
   log: config["output_dir"] + "/intermediate_output/logs/qc/qc.log"
   shell:
      '''echo "sample_id,global_n,s_n,s_n_contig,rbd_n" > {config[output_dir]}/qc.csv ;  cat {input} >> {config[output_dir]}/qc.csv'''

if (config["align"]) & (config["aligner"] == "minimap2"):
   rule get_indels:
      input:
         vcf_file = config["output_dir"] + "/intermediate_output/masked/{id}.masked.vcf" if config["filter"] == True else config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf",
         aln = config["output_dir"] + "/intermediate_output/align/{id}.fasta" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
      output:
         snps_indels = config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf",
         sample_n_perc = config["output_dir"] + "/intermediate_output/indels/{id}.nperc.csv"
      log: config["output_dir"] + "/intermediate_output/logs/indels/{id}.indels.log"
      shell:
         "get_indels.py --vcf {input.vcf_file} --window {config[del_window]} {config[allow_ambiguous]} --nperc {output.sample_n_perc} {input.aln} NC_045512.2 {output.snps_indels} 2> {log}"
else:
   rule get_indels:
      input:
         vcf_file = config["output_dir"] + "/intermediate_output/masked/{id}.masked.vcf" if config["filter"] == True else config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf",
         aln = config["output_dir"] + "/intermediate_output/align/{id}.fasta" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
      output:
         snps_indels = config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf",
         sample_n_perc = config["output_dir"] + "/intermediate_output/indels/{id}.nperc.csv"
      log: config["output_dir"] + "/intermediate_output/logs/indels/{id}.indels.log"
      shell:
         "get_indels.py --vcf {input.vcf_file} --window {config[del_window]} {config[allow_ambiguous]} --nperc {output.sample_n_perc} {input.aln} NC_045512.2 {output.snps_indels} 2> {log}"

rule get_snps:
   input:
      config["output_dir"] + "/intermediate_output/align/{id}.fasta" if config["align"] == True else config["input_dir"] + "/{id}" + config["extension"]
   output:
      snps = config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/get_snps/{id}.get_snps.log"
   shell:
      "faToVcf {config[exclude_ambiguous]} {input} {output.snps} ; update_vcf_header.sh {output.snps} 2> {log}"

if config["align"]:
   if config["aligner"] == "minimap2":
      rule align:
         input:
            config["input_dir"] + "/{id}" + config["extension"]
         output: 
            alignment = config["output_dir"] + "/intermediate_output/align/{id}.fasta" 
         log: config["output_dir"] + "/intermediate_output/logs/vcf/{id}.vcf.log"
         shell:
            "minimap2 -ax asm20 --cs {config[reference_sequence]} {input} > {config[output_dir]}/intermediate_output/align/{wildcards.id}.sam 2> {log} ;  gofasta sam toPairAlign -r {config[reference_sequence]} -s {config[output_dir]}/intermediate_output/align/{wildcards.id}.sam --outpath {config[output_dir]}/intermediate_output/align 2> {log}"
   elif config["aligner"] == "muscle":
      rule align:
         input:
            config["input_dir"] + "/{id}" + config["extension"]
         output: 
            plus_ref = config["output_dir"] + "/intermediate_output/consensus/{id}.consensus_ref.fa",
            alignment = config["output_dir"] + "/intermediate_output/align/{id}.fasta" 
         log: config["output_dir"] + "/intermediate_output/logs/align/{id}.align.log"
         shell:
            "cat {config[reference_sequence]} {input} > {output.plus_ref} ; muscle -quiet -in {output.plus_ref} -out {output.alignment} 2> {log}"
