# Workflow


----------


This section describes the methodology used by SigProfilerSimulator to generate realistic simulations of somatic mutations.

----------

## Overview ##

SigProfilerSimulator takes real somatic mutations as input and produces simulated samples by randomly redistributing those mutations across the genome. The redistribution is performed in an unbiased manner, preserving three biological properties of the original data:

- **Sequence context** — each mutation is placed in a position with the same local nucleotide context as the original
- **Transcriptional strand bias** — the strand orientation of each mutation relative to the direction of transcription is maintained
- **Chromosomal mutation burden** — the number of mutations assigned to each chromosome reflects the same proportional distribution as the input

This approach ensures that simulated samples are realistic null hypothesis models of the original mutation landscape, suitable as background distributions for downstream statistical analyses.

## Simulation Procedure ##

For each simulation, SigProfilerSimulator performs the following steps:

1. **Input parsing** — the input file (VCF, MAF, simple text, or ICGC format) is parsed and mutations are catalogued by sample, chromosome, and mutational context.

2. **Context distribution** — the genomic distribution of available positions for each mutation context is computed from the reference genome. If a BED file or exome restriction is provided, only the targeted regions are considered.

3. **Random placement** — each mutation is randomly assigned to a new position selected from the pool of positions sharing its original context. The number of mutations per chromosome is preserved by sampling within each chromosome independently.

4. **Output generation** — simulated mutations are written to MAF or VCF files (one per simulation). Parallel execution across chromosomes accelerates this step for large genomes.

## Use With SigProfilerClusters ##

Simulated datasets produced by SigProfilerSimulator are directly used as the background model in [SigProfilerClusters][1]. For this use case, the `chrom_based=True` parameter must be set to ensure per-chromosome normalisation of mutation burden.

  [1]: https://sigprofilersuite.github.io/SigProfilerClusters/
