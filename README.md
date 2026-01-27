# NJRS
NJRS (Neighbour-Joining Read-support) constructs phylogenetic trees for user-defined groups from aligned read sequences. Supports low-frequency and homoplasy-associated SNPs.

## Configuration
Config file: `NJRS_config.yaml`
* Input and output paths
* Parameter settings

## Snakemake command
Run command (with recommended parallelization cap to memory-intensive tasks)
```
snakemake --snakefile NJRS_workflow.snakemake --configfile NJRS_config.yaml --cores 50 --use-conda --resources max_heavy_memory_slots=4
```

## Input structure
Input file structure **example** (paths set in config)
```
.
├── bams
│   ├── bam1.bam
│   └── bam2.bam
├── fasta
│   ├── assembly1.fasta
│   └── assembly2.fasta
├── fasta_groups.tsv
├── groups
│   ├── group1.list
│   └── group2.list
└── ref
    ├── FSC200.fasta
    └── FSC200.gff
```

- Bam files should be sorted
- Group-lists are single-column lists with read dataset name without extension on each line
- Fasta-sequences must have .fasta extension
- fasta_groups.tsv is a two-column file with (1) dataset file name without extension, (2) group identifier

## Parameters
### Candidate SNP selection for MSA (performed group-by-group)
| Parameter | Description |
| ------ | ------ |
| snp_min_count | Minimum count to consider, SNPs below this number are not considered |
| pos_min_cov | Minimum coverage to consider a position, positions below this number are considered unmapped |
| lowcount_thresh | SNP read count to classify as lowcount |
| lowfreq_frac_thresh| See above, determined by the read count fraction |
| fixed_min_count | Minimum SNP read count to classify as "fixed" (qualify for MSA) |
| fixed_min_frac | See above, determined by the read count fraction |
| ignore_same | When specified, will not output position to the MSA where all datasets within the group have the same nucleotide |
| allow_lowcov_fraction | Allow positions with up to this fraction of datasets in the group to have a low-coverage/missing-coverage classification |

### MSA position classification
At every position, each dataset nucleotide is classified as "certain" or "ambiguous" dependent on the read support:
- When coverage is greater than 20:
   - "Certain" if read count support fraction is greater or equal to 90%
- When coverage is lower than or equal to 20:
  - "Certain" if read count support fraction is greater or equal to 97%

| Parameter | Description |
| ------ | ------ |
| ambig_output_cutoff | Fraction of "certain"-classified dataset nucleotides at the position required to classify the position as "stable" |

### Misc Parameters
| Parameter | Description |
| ------ | ------ |
| include_leafs | When constructing the tree, include single-dataset SNPs |
| fit_fasta | Use fasta sequences and fit them to identified SNPs |
| clust_dist | In the summary-file, cluster SNPs within this distance to form "hotspots" |
| make_noncongruent_snp_plots | If specified, will make a tree figure and indicate non-congruent dataset SNPs |


## Short note on memory usage
Some steps may consume a lot of memory, dependent on the number of SNPs identified within a group and on the total number of datasets in the group: Memory increase with a greater number of SNPs identified in the group and a greater number of datasets. Apply flag `--resources max_heavy_memory_slots=X` to limit parallelization to X for these steps.

## TODO's
- Input validation checks, e.g. duplicate dataset names, sorted bams
- Performance and memory optimization. Probably use a database to store basecounts data.
