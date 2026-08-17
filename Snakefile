# Snakefile
configfile: "src/config.yaml"

import os
import pandas as pd

#os.chdir(config["PTH"])
os.makedirs("log", exist_ok=True)

if not os.path.exists("analysis"):
  os.mkdir("analysis")

if not os.path.exists("alignment"):
  os.mkdir("alignment")

sample_names = pd.read_table(config["PTH"]+config["SAMPLE_LIST"],sep="\t",header=0)
pname= (list(set(sample_names.pname)))

rule all:
    input:
         expand("analysis/{sample}.del", sample=pname),
         "analysis/deletion_distribution.pdf"

rule align_samples:
    input:
        "raw_data/read1_{sample}.fastq.gz",
        "raw_data/read2_{sample}.fastq.gz"
    output:
        "alignment/{sample}.bam",
        "alignment/{sample}_sorted.bam.bai",
        "alignment/{sample}_sorted.bam",
        "alignment/{sample}_clean_sorted.bam",
        "alignment/{sample}_fix_clean_sorted.bam",
        "alignment/{sample}_rehead_fix_clean_sorted.bam",
        "alignment/{sample}_rehead_fix_clean_sorted.bam.bai",
        "alignment/{sample}.markdup_metrics",
    params:
        "resources/"+config["GENOME"]+".fa"
    resources:
        lsf_mem=4000,
        time="1:00",
        cpu=1
    conda:
        "lcprov2"
    shell:
        """
        bwa mem -t {resources.cpu} {params} {input[0]} {input[1]} | samtools view -b - > {output[0]}
        samtools sort -@ {resources.cpu} -o {output[2]} {output[0]} ; samtools index {output[2]}
        picard CleanSam -I {output[2]} -O {output[3]}              
        picard FixMateInformation -I {output[3]} -O {output[4]} --ASSUME_SORTED true --TMP_DIR tmp
        samtools addreplacerg -r "@RG\tID:ReadGroup1\tSM:{wildcards.sample}\tPL:Illumina\tLB:reporter.fa" -o {output[5]} {output[4]}
        #picard MarkDuplicates -I {output[5]} -O {output[6]} -M {output[7]} --REMOVE_DUPLICATES false --TMP_DIR tmp
        samtools index {output[6]}
        """

rule call_del:
     input:
         "alignment/{sample}_rehead_fix_clean_sorted.bam.bai"
     output:
         "analysis/{sample}.del"
     resources:
         lsf_mem=12000,
         time="5:00",
         cpu=8
     conda:
         "stingele",
     shell:
         """
         python bin/del_caller.py {input} {output[1]}
         """

rule plot:
      input:
          expand("analysis/{sample}.del", sample=pname)
      output:
          "analysis/deletion_distribution.pdf"
      resources:
          lsf_mem=4000,
          time="0:30",
          cpu=1
      conda:
          "lcprov2",
      shell:
          """
          Rscript --vanilla ../bin/plot.R 
          """

