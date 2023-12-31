# RNA-seq Analysis Workflow: From FASTQ to Quantification
Rabia Şen

## Introduction
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
Define placeholders for STAR alignment
```
export genomeIndexDir="<PATH_TO_GENOME_INDEX>"
export fastq1File=$OUTPUT_R1
export fastq2File=$OUTPUT_R2
export starOutputPrefix="<STAR_OUTPUT_PREFIX>"
```
```
STAR --runMode alignReads \
     --genomeDir ${genomeIndexDir} \
     --runThreadN 8 \
     --readFilesIn ${fastq1File} ${fastq2File} \
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
# Step 5: Sort and index BAM files with Samtools
Define placeholders for Samtools sorting and indexing
```
export inBam="${starOutputPrefix}Aligned.out.bam"
export outBam="${starOutputPrefix}Sorted.out.bam"
samtools sort -@ 8 ${inBam} > ${outBam}
samtools index ${outBam}
```

# Step 6: Quantify reads with featureCounts
Define placeholder for featureCounts
export gtfFile="<PATH_TO_GTF_FILE>"
export featureCountsOutput="<FEATURECOUNTS_OUTPUT_FILE>"

```
featureCounts -t exon \
      -g gene_id \
      --primary \
      -a ${gtfFile} \
      -o ${featureCountsOutput} \
      ${outBam}
```

