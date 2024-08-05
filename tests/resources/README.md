# Minimal reference generation as testdataset for neoantigen prediction pipeline

First step of ShrinkGenome is not required as it combines different regions from the genome (this could be interesting for fusion reference).

* Select a genomic region for the testdataset
  * For the given example the KRAS gene and ~1,500,000 bp up- and downstream were selected and written into `merged_intervals.tsv`

* Build minigenome

```
python ShrinkGenome/build_minigenome_from_intervals.py \
    --intervals merged_intervals.tsv \
    --genome CI_testdata/genome/GRCh38.primary_assembly.genome.fa \
    --spacer_len 0
```

* Translate coordinates

```
python ShrinkGenome/translate_fullgenome_to_minigenome_annot.py \
    --fullgenome_annot CI_testdata/genome/gencode.v46.basic.annotation.gtf \
    --translation_intervals minigenome.coord_translation.tsv \
    --output_gtf gencode.v46.basic.annotation.minigenome.fa
```

## Get subsetted reads

* Reads from given example selected from Riaz cohort (https://doi.org/10.1016/j.cell.2017.09.028)
  * WES tumor: SRR5134770
  * WES normal: SRR5134767
  * RNA tumor: SRR5088818

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

## Select the respective transcripts

* get the gene ids that fall into the specified region

```
bedtools intersect \
    -a gencode.v45.basic.annotation.gff3 \
    -b CI_testdata/merged_intervals.tsv \
    | less -S \
    | awk '$3 == "transcript" {print $9}' \
    | awk 'BEGIN{FS=";"}{print $4}' \
    | sed 's/transcript_id=//' > transcriptids.txt
```

* extract the transcript sequences that match the above gathered gene ids

```
for transcriptid in $(cat transcriptids.txt); 
do
    sed -n "/^>${transcriptid}|/,/^>/p" gencode.v45.transcripts.fa | head -n -1 \
    >> minitranscripts.fa;
done
```
