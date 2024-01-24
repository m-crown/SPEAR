problem_sites_in = config["output_dir"] + "/intermediate_output/merged.vcf" if not config["vcf"] else config["output_dir"] + "/input_files/input.vcf"
problem_sites_out = config["output_dir"] + "/intermediate_output/masked/merged.problem.vcf"
problem_sites_log = config["output_dir"] + "/intermediate_output/logs/mark_problem_sites/mark_problem_sites.log"
filter_sites_in = config["output_dir"] + "/intermediate_output/masked/merged.problem.vcf"
filter_sites_out = config["output_dir"] + "/intermediate_output/masked/merged.masked.vcf"
filter_sites_log = config["output_dir"] + "/intermediate_output/logs/filter_problem_sites/filter_problem_sites.log"
pangolin_input = config["output_dir"] + "/input_files/input.fasta.gz" if config["align"] else config["output_dir"] + "/intermediate_output/all_consensus.fa"
annotate_out = config["output_dir"] + "/intermediate_output/snpeff/merged.ann.vcf"
annotate_log = config["output_dir"] + "/intermediate_output/logs/snpeff/snpeff.log"
summarise_out = config["output_dir"] + "/intermediate_output/merged.summary.tsv"
spear_out = config["output_dir"] + "/all_samples.spear.vcf"
spear_log = config["output_dir"] + "/intermediate_output/logs/spear/spear.log"
summary_in = spear_out

if config["vcf"] == False:
   qc_file = config["output_dir"] + "/qc.csv"
   vcf_loc = config["output_dir"] + "/intermediate_output/indels/"
   vcf_ext = ".indels.vcf"
else:
   qc_file = config["output_dir"] + "/spear_score_summary.tsv"
   vcf_loc = config["output_dir"] + "/input_files/"
   vcf_ext = config["extension"]

if config["per_sample_outputs"] == "True":
   output = expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"])
else:
   output = config["output_dir"] + "/spear_annotation_summary.tsv"

if config["report"] == False:
   rule all:
      input: 
         outputs = output,
         lineages = config["output_dir"] + "/lineage_report.csv",
         qc = qc_file
else:
   rule all:
      input:
         outputs = output,
         report = config["output_dir"] + "/report/report.html",
         lineages = config["output_dir"] + "/lineage_report.csv",
         qc = qc_file

rule produce_report:
   input:
      summary = config["output_dir"] + "/spear_score_summary.tsv",
      all_samples = config["output_dir"] + "/spear_annotation_summary.tsv",
      n_perc = config["output_dir"] + "/qc.csv" if config["vcf"] == False else config["output_dir"] + "/spear_score_summary.tsv",
      lineage_report = config["output_dir"] + "/lineage_report.csv"
   output:
      config["output_dir"] + "/report/report.html"
   log: config["output_dir"] + "/intermediate_output/logs/report/report.log"
   shell:
      """summary_report.py --n_perc {input.n_perc} {config[product_plots]} {input.summary} {input.all_samples} {config[baseline_scores]} {config[input_sample_num]} {config[qc_sample_num]} {config[images_dir]} {config[scripts_dir]} {config[data_dir]} {config[output_dir]}/report/ {config[baseline]} {config[global_n]} {config[s_n]} {config[s_contig]} {config[rbd_n]} {input.lineage_report} {config[output_dir]}/pangolin_command.txt {config[spear_params]} 2> {log}"""

if config["per_sample_outputs"] == "True":
   rule summarise_vcfs:
      input:
         all_samples = expand(config["output_dir"] + "/final_vcfs/{id}.spear.vcf" , id=config["samples"]), 
         merged_samples = spear_out,
         lineages = config["output_dir"] + "/lineage_report.csv"
      output:
         config["output_dir"] + "/spear_annotation_summary.tsv",
         config["output_dir"] + "/spear_score_summary.tsv",
         expand(config["output_dir"] + "/per_sample_annotation/{id}.spear.annotation.summary.tsv", id = config["samples"])
      log: config["output_dir"] + "/intermediate_output/logs/summarise/summary.log"
      shell:
         """echo "{config[samples]}" > {config[output_dir]}/intermediate_output/sample_list.txt ; for i in $(find {config[output_dir]}/final_vcfs/ -type f -exec grep -lP '^NC_045512\.2' {{}} \;) ; do BNAME=${{i##*/}}; BNAME2=${{BNAME%.spear.vcf}}; grep -v '^#' "$i" | sed "s|^NC_045512\.2|$BNAME2|g"; done > {config[output_dir]}/intermediate_output/anno_concat.tsv 2> {log} ; convert_format.py --is_vcf_input {config[vcf]} --is_filtered {config[filter]} --per_sample_outputs {config[per_sample_outputs]} {config[output_dir]}/intermediate_output/anno_concat.tsv {config[output_dir]} {config[data_dir]} {input.merged_samples} {config[output_dir]}/intermediate_output/sample_list.txt 2> {log} """
else:
   rule summarise_vcfs:
      input:
         all_samples = expand(config["output_dir"] + "/final_vcfs/{id}.spear.vcf" , id=config["samples"]),
         merged_samples = spear_out,
         lineages = config["output_dir"] + "/lineage_report.csv"
      output:
         config["output_dir"] + "/spear_annotation_summary.tsv",
         config["output_dir"] + "/spear_score_summary.tsv"
      log: config["output_dir"] + "/intermediate_output/logs/summarise/summary.log"
      shell:
         """echo "{config[samples]}" > {config[output_dir]}/intermediate_output/sample_list.txt ; for i in $(find {config[output_dir]}/final_vcfs/ -type f -exec grep -lP '^NC_045512\.2' {{}} \;) ; do BNAME=${{i##*/}}; BNAME2=${{BNAME%.spear.vcf}}; grep -v '^#' "$i" | sed "s|^NC_045512\.2|$BNAME2|g"; done > {config[output_dir]}/intermediate_output/anno_concat.tsv 2> {log} ; convert_format.py --is_vcf_input {config[vcf]} --is_filtered {config[filter]} --per_sample_outputs {config[per_sample_outputs]} {config[output_dir]}/intermediate_output/anno_concat.tsv {config[output_dir]} {config[data_dir]} {input.merged_samples} {config[output_dir]}/intermediate_output/sample_list.txt 2> {log} """

if config["pangolin"] != "none":
   rule pangolin:
      input:
         pangolin_input
      output:
         config["output_dir"] + "/lineage_report.csv"
      log: config["output_dir"] + "/intermediate_output/logs/pangolin/pangolin.log"
      threads: workflow.cores
      shell:
         """run_pango.sh {input} {config[pangolin]} {config[output_dir]} {config[max_n]} {config[threads]} {log} 2> {log}"""

   if config["vcf"]:
      rule vcf_consensus:
         input:
            config["output_dir"] + "/final_vcfs/{id}.spear.vcf"
         output:
            config["output_dir"] + "/intermediate_output/vcf_consensus/{id}.fa"
         log: config["output_dir"] + "/intermediate_output/logs/vcf_consensus/{id}.vcf_cons.log"
         shell:
            """BNAME={wildcards.id} ; bgzip -k {config[output_dir]}/final_vcfs/{wildcards.id}.spear.vcf ; tabix {config[output_dir]}/final_vcfs/{wildcards.id}.spear.vcf.gz ; cat {config[reference_sequence]} | bcftools consensus {config[output_dir]}/final_vcfs/{wildcards.id}.spear.vcf.gz > {output} ; sed -i "s|NC_045512\.2|${{BNAME}}|g" {output}"""

      rule combine_vcf_consensus:
         input:
            expand(config["output_dir"] + "/intermediate_output/vcf_consensus/{id}.fa", id = config["samples"])
         output:
            config["output_dir"] + "/intermediate_output/all_consensus.fa"
         log: config["output_dir"] + "/intermediate_output/logs/vcf_consensus/combine.log"
         shell:
            """cat {input} > {output}"""

   elif config["align"] == False:
      rule combine_samples_for_pango:
         input:
            expand(config["output_dir"] + "/input_files/{id}" + config["extension"], id = config["samples"])
         output:
            config["output_dir"] + "/intermediate_output/consensus.fa"
         log: config["output_dir"] + "/intermediate_output/logs/pre_align_consensus/combine.log"
         shell:
            """for i in {input}; do BNAME=${{i##*/}}; BNAME2=${{BNAME%.*}}; sed -n '/>'"${{BNAME2}}"'/, />/{{ />NC_045512\.2/!p ;}}' ${{i}} ; done > {config[output_dir]}/intermediate_output/consensus.fa"""
   
else:
   rule pangolin_touch_file:
      input:
         pangolin_input
      output:
         report = config["output_dir"] + "/lineage_report.csv",
         command = config["output_dir"] + "/pangolin_command.txt"
      log: config["output_dir"] + "/intermediate_output/logs/pangolin/pangolin_touch.log"
      shell:
         """touch {output.report} ; echo "Pangolin not run\nNA" > {output.command}"""

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
      annotate_out
   output:
      summarise_out = summarise_out,
      spear_out = spear_out
   log: spear_log
   shell:
      "summarise_snpeff.py {output.summarise_out} {input} {config[data_dir]} ; spear_annotate.py {output.spear_out} {output.summarise_out} {config[data_dir]} 2> {log}"

rule annotate_variants:
   input:
      filter_sites_out if config["filter"] else problem_sites_out
   output:
      annotate_out
   log: annotate_log
   shell:
      """
      java -Xmx2g -jar $CONDA_PREFIX/snpEff/snpEff.jar -noShiftHgvs -hgvs1LetterAa -download -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 {input} > {output} 2> {log}
      """

if config["filter"]:
   rule filter_problem_sites:
      input: 
         filter_sites_in
      output: 
         filter_sites_out  
      log: filter_sites_log
      shell:
         "java -Xmx2g -jar $CONDA_PREFIX/snpEff/SnpSift.jar filter \"!( {config[filter_params]} )\" {input} > {output} 2> {log}"

rule mark_problem_sites:
   input:
      problem_sites_in 
   output:
      problem_sites_out
   log: problem_sites_log
   shell:
      "vcfanno --base-path $CONDA_PREFIX $CONDA_PREFIX/bin/conf.toml {input} > {output} 2> {log}"

rule merge_vcfs:
   input:
      expand(vcf_loc + "{id}" + vcf_ext, id=config["samples"])
   output:
      config["output_dir"] + "/intermediate_output/merged.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/merge/merge.log"
   shell:
      '''find {vcf_loc} -type f -name "*{vcf_ext}" > {config[output_dir]}/intermediate_output/merge_list.txt ; bcftools merge --no-index -m none -o {output} -l {config[output_dir]}/intermediate_output/merge_list.txt 2> {log}'''

rule get_indels:
   input:
      vcf_file = expand(config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf", id=config["samples"]),
      aln = expand(config["output_dir"] + "/intermediate_output/align/{id}.fasta" if config["align"] == True else config["output_dir"] + "/input_files/{id}" + config["extension"], id=config["samples"])
   output:
      snps_indels = expand(config["output_dir"] + "/intermediate_output/indels/{id}.indels.vcf", id=config["samples"]),
      qc_file = config["output_dir"] + "/qc.csv"
   params:
      vcf_dir = config["output_dir"] + "/intermediate_output/fatovcf/",
      aln_dir = config["output_dir"] + "/intermediate_output/align/" if config["align"] == True else config["output_dir"] + "/input_files/",
      out_dir = config["output_dir"] + "/intermediate_output/indels/",
      out_suffix = ".indels.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/indels/indels.log"
   threads: workflow.cores
   shell:
      "get_indels.py --threads {threads} --vcf-dir {params.vcf_dir} --alignments-dir {params.aln_dir} --window {config[del_window]} {config[allow_ambiguous]} --nperc {output.qc_file} --ref NC_045512.2 --out-dir {params.out_dir} --out-suffix {params.out_suffix} 2> {log}"

rule get_snps:
   input:
      config["output_dir"] + "/intermediate_output/align/{id}.fasta" if config["align"] == True else config["output_dir"] + "/input_files/{id}" + config["extension"]
   output:
      snps = config["output_dir"] + "/intermediate_output/fatovcf/{id}.vcf"
   log: config["output_dir"] + "/intermediate_output/logs/get_snps/{id}.get_snps.log"
   shell:
      "faToVcf {config[exclude_ambiguous]} {input} {output.snps} ; update_vcf_header.sh {output.snps} 2> {log}"

if config["align"]:
   if config["aligner"] == "minimap2":
      rule align:
         input:
            config["output_dir"] + "/input_files/input.fasta.gz"
         output:
            alignment = config["output_dir"] + "/intermediate_output/align/{id}.fasta"
         log: config["output_dir"] + "/intermediate_output/logs/align/{id}.align.log"
         shell:
            "seqkit grep {input} -p {wildcards.id} | minimap2 -ax asm20 --cs {config[reference_sequence]} - > {config[output_dir]}/intermediate_output/align/{wildcards.id}.sam 2> {log} ;  gofasta sam toPairAlign -r {config[reference_sequence]} -s {config[output_dir]}/intermediate_output/align/{wildcards.id}.sam --outpath {config[output_dir]}/intermediate_output/align 2> {log}"
   elif config["aligner"] == "muscle":
      rule align:
         input:
            config["output_dir"] + "/input_files/input.fasta.gz"
         output:
            plus_ref = config["output_dir"] + "/intermediate_output/consensus/{id}.consensus_ref.fa",
            alignment = config["output_dir"] + "/intermediate_output/align/{id}.fasta" 
         log: config["output_dir"] + "/intermediate_output/logs/align/{id}.align.log"
         shell:
            "seqkit grep -p {wildcards.id} {input} | seqkit seq {config[reference_sequence]} - > {output.plus_ref} ; muscle -quiet -in {output.plus_ref} -out {output.alignment} 2> {log}"
