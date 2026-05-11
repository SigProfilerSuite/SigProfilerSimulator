# Using SigProfilerSimulator - Output


----------


This section describes the output files and directories produced by SigProfilerSimulator. All results are written under the project directory specified in `project_path`.

----------

## Output Overview ##

After a successful run, the project directory contains the following structure:

```
[project_path]/
├── input/
│   └── [input mutation files]
├── output/
│   └── simulations/
│       └── [project]_simulations_[genome]_[contexts]/
│           ├── 1.maf
│           ├── 2.maf
│           └── ...
└── logs/
    ├── SigProfilerSimulator_[project]_[genome]_[date].err
    └── SigProfilerSimulator_[project]_[genome]_[date].out
```

The output subdirectory name includes a suffix depending on the simulation scope:

| Condition | Output subdirectory suffix |
|-----------|---------------------------|
| Whole genome (default) | `[project]_simulations_[genome]_[contexts]/` |
| BED file restricted (`bed_file`) | `[project]_simulations_[genome]_[contexts]_BED/` |
| Exome restricted (`exome=True`) | `[project]_simulations_[genome]_[contexts]_exome/` |

## Simulation Files ##

### MAF format (default) ###

When `vcf=False` (default), one MAF file is produced per simulation, named `1.maf`, `2.maf`, etc. Each file contains all samples for that simulation iteration.

MAF files contain the following 17 columns:

| Column | Description |
|--------|-------------|
| `Hugo_symbol` | Gene symbol |
| `Entrez_gene_ID` | Entrez gene identifier |
| `Center` | Sequencing center |
| `Genome` | Reference genome used |
| `Chrom` | Chromosome |
| `Start_position` | Mutation start position (1-based) |
| `End_position` | Mutation end position (1-based) |
| `Strand` | Genomic strand |
| `Variant_Classification` | Mutation functional class |
| `Variant_Type` | Mutation type (SNP, INS, DEL) |
| `Reference_Allele` | Reference base(s) |
| `Tumor_Seq_Allele1` | First tumor allele |
| `Tumor_Seq_Allele2` | Second tumor allele |
| `dbSNP_RS` | dbSNP RS identifier |
| `dbSNP_Val_Status` | dbSNP validation status |
| `Tumor_Sample_Barcode` | Sample identifier |
| `matGenClass` | Mutational context classification |

### VCF format ###

When `vcf=True`, a separate subdirectory is created per sample within the output folder. Each subdirectory contains one VCF file per simulation:

```
[project]_simulations_[genome]_[contexts]/
├── [sample_1]/
│   ├── [sample_1]_1.vcf
│   ├── [sample_1]_2.vcf
│   └── ...
└── [sample_2]/
    ├── [sample_2]_1.vcf
    └── ...
```

## Log Files ##

Two log files are saved per run in the `logs/` subdirectory:

| File | Description |
|------|-------------|
| `SigProfilerSimulator_[project]_[genome]_[date].out` | Progress checkpoints and run parameters |
| `SigProfilerSimulator_[project]_[genome]_[date].err` | Error messages and warnings |

## Sequence Context Files ##

When `seqInfo=True`, additional files containing the sequence context of each simulated mutation are saved under:

```
[project_path]/output/vcf_files/simulations/[context]/
```

One file is saved per sample, chromosome, and simulation: `[sample]_[chrom]_seqinfo_[n].txt`.
