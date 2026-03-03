# MGY360 Dry Lab Notebook

**MGY360H1S** Whole-Genome Sequencing and Analysis (Winter 2026) | **Part 2:** Dry Lab & Bioinformatics

**Will Li** ([/williamli9300](https://github.com/williamli9300))

----
# Contents
- [**Laboratory 1:** Quality Control & Alignment](#lab1)
- [**Laboratory 2:** Viewing the Alignments](#lab2)
- [**Appendices**](#appx)

<br>

----
# Laboratory 1: Quality Control & Alignment <a name="lab1"></a>

## 1.1 Galaxy
In this lab, the **Galaxy**<sup>[1](#r1.1)</sup> ([https://usegalaxy.org/](https://usegalaxy.org/)) platform was used.

## 1.2 Uploading Sequences
Two sequences prepared from *S. cerevisiae* samples and sequenced via paired-end Illumina NextSeq 2000 sequencing (<a href="#appx1.2">Appendix 1.1</a>), downsampled to 5000 reads:

- `Andrews_008_H01_Will_9_S8_R1_5000.fastq` corresponds with the **forward** read, and
- `Andrews_008_H01_Will_9_S8_R2_5000.fastq` corresponds with the **reverse** read.

Sequences were uploaded to **Galaxy**, using the `S. cerevisae str. S288C (Saccharomyces_cerevisiae_S288C_SGD2010)` reference genome.

## 1.3 Initial Quality Inspection

The first 10 quality scores in both the **forward** *and* **reverse** reads were:

```
IIIIIIIIII
```

indicating that the first 10 nucleotide reads for each were assigned a Phred score of $Q40$, suggesting good quality reads.

## 1.4 Quality Inspection & `FastQC`

*FastQC Read Quality Reports<sup>[2](#r1.2)</sup> (<a href="#appx2">Appendix 1.2</a>) generated using **FastQC version 0.74+galaxy1**.*

### 1.4.1 FastQC Report for Forward Read

The FastQC report for the **forward** read returned a total of **755 kbp**, over **5000** total sequences of length **151**. The GC content was **40%**.

One [**$!$ Warning**] was thrown for **per sequence GC count**, indicating a deviation of between 15-30% from the theoretical distribution<sup>[3](#r1.3)</sup>. The peak of the actual distribution appeared to remain ±1.5 percent mean GC content from the theoretical distribution; the shape of the actual distribution curve appears slightly left skewed compared to the theoretical distribution.

Two [**$\times$ Error**]s were thrown, for **per base sequence content** and **adapter content**. The first indicates a >20% difference in A/T or G/C base ratios at at least one position<sup>[3](#r1.3)</sup>. Visual inspection of the Sequence Content vs. Position graph reveals highly variable sequence content at the first 10-15 bases. This is likely the result of potential nucleotide bias arising from random priming<sup>[4](#r1.4)</sup>. For the second, a visual inspection of the % Adapter vs. Position graph reveals cumulatively, around 25% of Nextera Transposase Sequence was observed by the last position of each read (with adapter observed starting around bases 45-49), as well as up to around 5% poly-guanine tracts. beginning at positions 95-99<sup>[3](#r1.3)</sup>. Both suggest the presence of inserts shorter than the read length, resulting in read-through into adapter sequences, or potentially the calling of calling of dark signal (i.e., no signal) as `G` even after termination of synthesis<sup>[5](#r1.5)</sup>.

### 1.4.2 FastQC Report for Reverse Read

The FastQC report for the Reverse read had the same one [**$!$ Warning**] and two [**$\times$ Error**]s as the report for the Forward read. These were assessed to be caused by the same factors as in the Forward read.

There was additionally one [**$!$ Warning**] for **per tile sequence quality**, where one standout tile around read 2202-2203 and around position 30-34 was indicated orange for poor quality<sup>[3](#r1.3)</sup>, and one other [**$!$ Warning**] for **overrepresented sequences**, which counted eleven poly-G tracts fifty `G`s long, likely arising for the same reason as mentioned in 1.4.1.

### 1.4.2 Characterization of Reads from FastQC Reports

Across both reports, quality per base was generally very good, with average quality always at 36 or higher or 38 or higher for the forward and reverse reads, respectively. In the forward read, the lower whisker began to fall into the orange zone (passable quality) around read 115-119, while in the reverse read the lower whisker always remained in the green zone. This terminal quality decay is expected and attributable to accumulation of error and noise or the presence of shorter-than-expected insert lengths. 

Both reports suggest, by visual inspection, a GC content of just under 40%, consistent with a known approximately 38% GC content in the reference genome<sup>[6](#r1.6)</sup>. 

Both sequences had an average quality score of 39. Given $Q = -10 log_{10}(P)$, we know that $P = 10^{-\frac{Q}{10}}$, so for a Phred score of $Q39$, we get a $P = 10^(-3.9) = 0.0126%$ expected rate of incorrect base calls. Given a read length of 151, we expect around 0.0196 incorrect base calls per read, and around 19026 incorrect base calls per million reads (126 incorrect base calls per million base calls). As expected, since both forward and reverse libraries were prepared from the same library, the quality and adapter content across the two sequence files are similar.

## 1.5 Sequence Alignment with Reference Genome

*Sequence Alignment performed using **Bowtie2**<sup>[7](#r1.7)</sup> **version 2.5.4+galaxy0**.*

Bowtie2 paired library alignment run with default settings, except maximum fragment length set to **800**, using the `sacCer3` reference genome index.

<details>
<summary><b>Bowtie2 Command Line</b> (<i>click to expand</i>)</summary>

```
set -o | grep -q pipefail && set -o pipefail;   ln -f -s '/corral4/main/objects/0/c/0/dataset_0c0c670a-7c7f-4e6b-822c-70792f762c3f.dat' input_f.fastq.gz &&  ln -f -s '/corral4/main/objects/e/0/1/dataset_e010227a-2b32-4ae9-9df1-3fbd6c1a31e4.dat' input_r.fastq.gz &&   THREADS=${GALAXY_SLOTS:-4} && if [ "$THREADS" -gt 1 ]; then (( THREADS-- )); fi &&   bowtie2  -p "$THREADS"  -x '/cvmfs/data.galaxyproject.org/managed/bowtie2_index/sacCer3/sacCer3'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz' -I 0 -X 800 --fr                     2> >(tee '/corral4/main/jobs/075/048/75048950/outputs/dataset_e7fc3948-60d2-4b1e-8730-c22aefb3182c.dat' >&2)  | samtools sort -l 0 -T "${TMPDIR:-.}" -O bam | samtools view --no-PG -O bam -@ ${GALAXY_SLOTS:-1} -o '/corral4/main/jobs/075/048/75048950/outputs/dataset_527fa994-cf90-41c0-8dd6-6cff1b68df19.dat'
```

</details>

Bowtie2 defaults to the `--sensitive` flag, meaning sensitive alignment in end-to-end alignment mode (`--end-to-end`). Local alignment allows alignment of substrings within a sequence to the reference genome, whereas end-to-end alignment prioritizes alignment of the entire sequence; this results in slightly different alignment score calculation. Sensitivity options (e.g. `--fast`, `--sensitive`) modulate `-D`, `-R`, `-N`, `-L`, and  `-i` flag arguments, which affect the number of allowed consecutive failed seed extension attempts, maximum number of allowed reseeding attempts, mismatches allowed during seed alignment, length of seed substrings, and interval between seed substrings, respectively<sup>[8](#r1.8)</sup>. Local alignment helpful when searching for local homology (e.g. conserved domains over several different sequences). More sensitive values increase runtime but increase sensitivity for alignment, which allows less perfect matches to be accepted and revealed (e.g., when matching a sequence with low, e.g. <70%, homology to any known reference sequence). Since our sequences are expected to match the reference sequence with high homology (genetically very close to the reference genome) across the entire genome (whole genome sequencing and alignment), more sensitive or local alignment options would not be necessary (and may not be helpful), regardless of how much compute power we have access to.

<details open>
<summary><b>Bowtie2 Standard Error</b> (<i>click to collapse</i>)</summary>

```
5000 reads; of these:
  5000 (100.00%) were paired; of these:
    1652 (33.04%) aligned concordantly 0 times
    2947 (58.94%) aligned concordantly exactly 1 time
    401 (8.02%) aligned concordantly >1 times
    ----
    1652 pairs aligned concordantly 0 times; of these:
      548 (33.17%) aligned discordantly 1 time
    ----
    1104 pairs aligned 0 times concordantly or discordantly; of these:
      2208 mates make up the pairs; of these:
        2013 (91.17%) aligned 0 times
        83 (3.76%) aligned exactly 1 time
        112 (5.07%) aligned >1 times
79.87% overall alignment rate
```

**66.96%** of reads were aligned at least once concordantly, and **10.96%** were aligned discordantly; a further **1.95%** of reads were mates that aligned at least once. The final overall alignment rate was **79.87%**. Concordant alignments satisfy both the upstream/downstream mate orientation requirements (`-fr`, `-rf`, `-ff`) as well as minimum and maximum fragment sizes (`-I`, `-X`); discordant alignments are alignments that do not meet all these requirements. An overall alignment of roughly ~80% is seems fairly acceptable for a first-pass alignment using downsampled data. Increasing sensitivity options may increase the overall alignment, but given the present experiment it is not likely to increase by much without changing the dataset (using the `--very-sensitive` flag only increased the overall alignment to 79.90%.)

</details>

## References

1. <a name="r1.1"></a> The Galaxy Community. The Galaxy platform for accessible, reproducible, and collaborative data analyses: 2024 update, *Nucleic Acids Res.* **52**, W83–W94 (2024). doi: [10.1093/nar/gkae410](https://doi.org/10.1093/nar/gkae410).
   
2. <a name="r1.2"></a> Andrews, S. FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. (2010) [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
   
3. <a name="r1.3"></a> Quality control: How do you read your FASTQC results? *CD Genomics Bioinformatics Analysis*. [https://bioinfo.cd-genomics.com/](https://bioinfo.cd-genomics.com/).
   
4. <a name="r1.4"></a> Hansen, K.D., Brenner, S.E., & Dudoit, S. Biases in Illumina transcriptome sequencing caused by random hexamer priming. *Nucleic Acids Res.* **38**, e131 (2010). doi: [10.1093/nar/gkq224](https://doi.org/10.1093/nar/gkq224).
   
5. <a name="r1.5"></a> Read Trimming. *Illumina DRAGEN Bio-IT Platform v4.0*. [https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/ReadTrimming.htm](https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/ReadTrimming.htm).

6. <a name="r1.6"></a> Wang, D. & Gao, F. Comprehensive analysis of replication origins in *Saccharomyces cerevisiae* genomes. *Front. Microbiol.* **10**, 2122 (2019). doi: [10.3389/fmicb.2019.02122](https://doi.org/10.3389/fmicb.2019.02122).

7. <a name="r1.7"></a> Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Bio.*, **10**, R25 (2009). doi: [10.1186/gb-2009-10-3-r25](https://doi.org/10.1186/gb-2009-10-3-r25).

8. <a name="r1.8"></a> Bowtie 2: Manual. *Bowtie (Sourceforge)*. [https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

<br>

----
# Laboratory 2: Viewing the Alignments <a name="lab2"></a>

Lorem ipsum dolor sit amet consectetuer adipiscing elit...

<br>

----
# Appendices <a name="appx"></a>

## Appendix 1.1 Sequence Files <a name="appx1.1"></a>

- `Andrews_008_H01_Will_9_S8_R1_5000.fastq`: [gzip](./fastq/Andrews_008_H01_Will_9_S8_R2_5000.fastq.gz)
- `Andrews_008_H01_Will_9_S8_R2_5000.fastq`: [gzip](./fastq/Andrews_008_H01_Will_9_S8_R2_5000.fastq.gz)

## Appendix 1.2 FastQC Reports Forward and Reverse Reads <a name="appx1.2"></a>

- `Andrews_008_H01_Will_9_S8_R1_5000_FastQC.html`: [html](./fastqc/Andrews_008_H01_Will_9_S8_R2_5000_FastQC.html)
- `Andrews_008_H01_Will_9_S8_R2_5000_FastQC.html`: [html](./fastqc/Andrews_008_H01_Will_9_S8_R2_5000_FastQC.html)