# Read subset from Riaz cohort Pt23

## Background

* The reads from the following samples were downsampled from the alignment file.
* Reads mapping ~1,500,000 up- and downstream of the KRAS gene were selected
  * chr12:23377886-27078289

| sequencing type | condition | SRA identifier |
|-----------------|-----------|----------------|
|      WES        |   tumor   |   SRR5134770   |
|      WES        |   normal  |   SRR5134767   |
|      RNA        |   tumor   |   SRR5088818   |

## Read extraction

* Select reads from region from a mapped bam file and sort the bam file by name (this step has to be executed for each sample)

```
samtools view -b -h \
    CI_testdata/reads/SRR5134767_1.fastq.gz.bam.markDup.bam \
    "chr12:23377886-27078289" \
    | samtools sort -n \
    > minireads/SRR5134767.mini.namesorted.bam

samtools view -b -h \
    CI_testdata/reads/SRR5134770_1.fastq.gz.bam.markDup.bam \
    "chr12:23377886-27078289" \
    | samtools sort -n \
    > minireads/SRR5134770.mini.namesorted.bam

samtools view -b -h \
    CI_testdata/reads/RNA_1_Aligned.out.bam \
    "chr12:23377886-27078289" \
    | samtools sort -n \
    > minireads/SRR5088818.mini.namesorted.bam
```

* Convert bam to fastq (also has to be done for each sample)

```
bedtools bamtofastq \
    -i minireads/SRR5134767.mini.namesorted.bam \
    -fq minireads/SRR5134767_mini_1.fastq \
    -fq2 minireads/SRR5134767_mini_2.fastq

bedtools bamtofastq \
    -i minireads/SRR5134770.mini.namesorted.bam \
    -fq minireads/SRR5134770_mini_1.fastq \
    -fq2 minireads/SRR5134770_mini_2.fastq

bedtools bamtofastq \
    -i minireads/SRR5088818.mini.namesorted.bam \
    -fq minireads/SRR5088818_mini_1.fastq \
    -fq2 minireads/SRR5088818_mini_2.fastq
```
