# Output

The output of the pipeline consists of two main files.

* `results/{sample}/fetchdata/sj_final.tsv`
* `results/{sample}/fetchdata/sj_final_neofox_annotation.tsv`

Within the files, each line describes a candidate junction transcript. The file `sj_final.tsv` contains all candidate junctions with annotated features and targeted re-quantification results. The file `sj_final_neofox_annotation.tsv` contains only junction candidates with a mutated peptide sequence in format
suitable for analysis with NeoFox.

## Column-description sj_final.tsv

* `junc_id`: Predicted splice junction in standardised format <chr>:<start>-<end>:<strand>