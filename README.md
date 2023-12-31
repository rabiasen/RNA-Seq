# RNA-seq Analysis Workflow: From FASTQ to Quantification
Rabia Şen

# Introduction
This guide includes my methodology for analyzing RNA-seq data. The materials and sources consulted to create this tutorial are listed in the resources section at the end.
Here is the diagram I prepared;
>![image](https://github.com/rabiasen/RNA-Seq/assets/58332251/099058c2-3f54-4b5c-bd8f-ee5ab6946869)
# Step 1: Trim adapters with  [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) 
Define placeholders for the forward, reverse adapter sequences
```
export FADAPTER=<FORWARD_ADAPTER_SEQUENCE>
export RADAPTER=<REVERSE_ADAPTER_SEQUENCE>
```
Define placeholders for input and output file names for Cutadapt
```
export INPUT_R1="<INPUT_FASTQ_R1>"
export INPUT_R2="<INPUT_FASTQ_R2>"
export OUTPUT_R1="<OUTPUT_TRIMMED_R1>"
export OUTPUT_R2="<OUTPUT_TRIMMED_R2>"
```
Run this for trimming adapters from paired-end data 
```
cutadapt -a $FADAPTER -A $RADAPTER -o $OUTPUT_R1 -p $OUTPUT_R2 $INPUT_R1 $INPUT_R2 -j 8 -m 32
```
# Step 2: Quality control of trimmed reads with [MultiQC](https://github.com/MultiQC/MultiQC)
```
multiqc ./ -o ./multiqc_report_pre_alignment
```
# Step 3: Align reads using [STAR](https://github.com/alexdobin/STAR)

```
STAR --genomeDir ${genomeIndexDir} \
     --runThreadN 8 \
     --readFilesIn $OUTPUT_R1 $OUTPUT_R2 \
     --outFileNamePrefix ${starOutputPrefix} \
     --outSAMtype BAM Unsorted \
     --outSAMmapqUnique 60 \
     --outSAMunmapped Within KeepPairs \
     --readFilesCommand zcat 
```
# Step 4: Quality control after alignment with MultiQC
```
multiqc ./ -o ./multiqc_report_post_alignment
```
# Step 5: Sort and index BAM files with [Samtools](https://github.com/samtools/samtools)
```
samtools sort -@ 8 ${starOutputPrefix}Aligned.out.bam > ${starOutputPrefix}Sorted.out.bam
samtools index ${starOutputPrefix}Sorted.out.bam
```

# Step 6: Quantify reads with [featureCounts](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)

```
featureCounts -t exon \
      -g gene_id \
      --primary \
      -a ${genomegtfFile} \
      -o ${featureCountsOutput} \
      ${starOutputPrefix}Sorted.out.bam
```

