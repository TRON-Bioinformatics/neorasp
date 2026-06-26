Description on how to generate a minimal reference and testdataset for neoantigen prediction pipelines based on [ShrinkGenome](https://github.com/brianjohnhaas/ShrinkGenome).

First step of ShrinkGenome is not required as it combines different regions from the genome (this could be interesting for fusion reference).

- Select a genomic region for the testdataset

  - For the given example the KRAS gene and ~1,500,000 bp up- and downstream were selected and written into `merged_intervals.tsv`

- Build minigenome

```
python ShrinkGenome/build_minigenome_from_intervals.py \
    --intervals merged_intervals.tsv \
    --genome CI_testdata/genome/GRCh38.primary_assembly.genome.fa \
    --spacer_len 0
```

- Translate coordinates

```
python ShrinkGenome/translate_fullgenome_to_minigenome_annot.py \
    --fullgenome_annot CI_testdata/genome/gencode.v46.basic.annotation.gtf \
    --translation_intervals minigenome.coord_translation.tsv \
    --output_gtf gencode.v46.basic.annotation.minigenome.fa
```
