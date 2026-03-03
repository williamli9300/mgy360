# MGY360 Dry Lab Notebook

**Will Li | MGY360 H1 S | Winter 2026**

<br>

----
# Contents
- [Laboratory 1: Quality Control & Alignment](#lab1)
- Appendix

<br>

----
# Laboratory 1: Quality Control & Alignment <a name="lab1"></a>

## 1.1 Galaxy
In this lab, the **Galaxy** ([https://usegalaxy.org/](https://usegalaxy.org/)) platform was used.

## 1.2 Uploading Sequences
Two sequences prepared from *S. cerevisiae* samples and sequenced via paired-end Illumina NextSeq 2000 sequencing:

- `Andrews_008_H01_Will_9_S8_R1_5000.fastq` corresponds with the **forward** read, and
- `Andrews_008_H01_Will_9_S8_R2_5000.fastq` corresponds with the **reverse** read.

Sequences were uploaded to **Galaxy**, using the `S. cerevisae str. S288C (Saccharomyces_cerevisiae_S288C_SGD2010)` reference genome.

## 1.3 Initial Quality Inspection

The first 10 quality scores in both the **forward** *and* **reverse** reads were:

```
IIIIIIIIII
```

indicating that the first 10 nucleotide reads for each were assigned a Phred score of $Q30$, suggesting good quality reads.

## 1.4 Quality Inspection & `FastQC`

FastQC Read Quality Reports generated using Galaxy Version 0.74+galaxy1.

The FastQC report for the **forward**</abbr> read (<a href="#appx2">Appendix 2</a>)

## 1.5 Sequence Alignment with Reference Genome

<br>

----
# Appendices

## Appendix 1.2 FastQC Report for <abbr title="Andrews_008_H01_Will_9_S8_R1_5000.fastq">Forward</abbr> Read <a name="appx2"></a>

<iframe src="./fastqc/Andrews_008_H01_Will_9_S8_R1_5000_FastQC.html" title="description" width="500" height="300"></iframe>

## Appendix 1.3 FastQC Report for <abbr title="Andrews_008_H01_Will_9_S8_R2_5000.fastq">Reverse</abbr> Read <a name="appx3"></a>

<iframe src="./fastqc/Andrews_008_H01_Will_9_S8_R1_5000_FastQC.html" title="description" width="500" height="300"></iframe>