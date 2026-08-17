This repo is for the characterisation of mutations falling over a short DNA construct. It was made by Craig Anderson at DKFZ in collaboration with Julian Stingele and Maximilian Donsbach at LMU.

The input is paired end Illumina sequencing data, where both the forward and the reverse was able to output the entirety of a DNA construct. We use this setup to establish if a cells with, or without HMCES, generate the same distribution of deletions as a result of abasic site processing near a region of microhomology, when compared to WT. Scripts collected in a Snakemake pipeline this in the following way:

1) Align fastq reads and QC reads to produce bam files, using 
2) Process bam files to characterise the position and size of deletions using pysam
3) Statistically analyse the difference between different cell lines and abasic site configurations in R.

I provide the scripts necessary to recapitulate these results. The raw data and processed alignments are on data dryad. Before running the Snakefile, please ensure you establish the two conda environments.

stingele_del
├── alignment #completely processed alignments available on dryad
│   ├── AAVS1_24_U_rehead_fix_clean_sorted.bam
│   ├── AAVS1_24_U_rehead_fix_clean_sorted.bam.bai
│   ├── AAVS1_31_U_rehead_fix_clean_sorted.bam
│   ├── AAVS1_31_U_rehead_fix_clean_sorted.bam.bai
│   ├── AAVS1_no_U_rehead_fix_clean_sorted.bam
│   ├── AAVS1_no_U_rehead_fix_clean_sorted.bam.bai
│   ├── HMCES_KO_24_U_rehead_fix_clean_sorted.bam
│   ├── HMCES_KO_24_U_rehead_fix_clean_sorted.bam.bai
│   ├── HMCES_KO_31_U_rehead_fix_clean_sorted.bam
│   ├── HMCES_KO_31_U_rehead_fix_clean_sorted.bam.bai
│   ├── HMCES_KO_no_U_rehead_fix_clean_sorted.bam
│   └── HMCES_KO_no_U_rehead_fix_clean_sorted.bam.bai
├── analysis #processed deletions- only single deletions, >1 nt and occurring on both reads are reported
│   ├── AAVS1_24_U.del
│   ├── AAVS1_31_U.del
│   ├── AAVS1_no_U.del
│   ├── deletion_distribution.pdf #plot of results
│   ├── HMCES_KO_24_U.del
│   ├── HMCES_KO_31_U.del
│   └── HMCES_KO_no_U.del
├── bin
│   ├── analysis.R #glm and other analyses
│   ├── del_caller.py #deletion characterisation- see separate write up in bin.
│   └── plot.R
├── log
├── raw_data #raw illumina reads available on dryad
│   ├── barcodes.txt
│   ├── demultiplex_summary.csv
│   ├── read1_AAVS1_24_U.fastq.gz
│   ├── read1_AAVS1_31_U.fastq.gz
│   ├── read1_AAVS1_no_U.fastq.gz
│   ├── read1_AAVS1_no_U_rep1.fastq.gz
│   ├── read1_HMCES_KO_24_U.fastq.gz
│   ├── read1_HMCES_KO_31_U.fastq.gz
│   ├── read1_HMCES_KO_no_U.fastq.gz
│   ├── read2_AAVS1_24_U.fastq.gz
│   ├── read2_AAVS1_31_U.fastq.gz
│   ├── read2_AAVS1_no_U.fastq.gz
│   ├── read2_HMCES_KO_24_U.fastq.gz
│   ├── read2_HMCES_KO_31_U.fastq.gz
│   └── read2_HMCES_KO_no_U.fastq.gz
├── resources #construct reference
│   ├── reporter.fa
│   ├── reporter.fa.amb
│   ├── reporter.fa.ann
│   ├── reporter.fa.bwt
│   ├── reporter.fa.fai
│   ├── reporter.fa.pac
│   ├── reporter.fa.sa
│   └── sample_names.txt
├── Snakefile #scripts run from snakemake
├── src #cluster organisation
│   ├── cluster.json
│   ├── config.yaml
│   ├── lcprov2.yaml
│   └── stingele.yaml
└── tmp

Please note: our assumption is that every read represents an individual mutagenic event, though without barcoding individual reads, we're unable to assert this.
