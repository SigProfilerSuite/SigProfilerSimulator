# Using SigProfilerSimulator


----------


This section describes SigProfilerSimulator's main function and all available parameters.

----------

## Function ##

The main function in SigProfilerSimulator is `SigProfilerSimulator`. It randomizes the position of each somatic mutation across the genome while preserving the sequence context, transcriptional strand bias, and chromosomal mutational burden of the input data.

### Input files ###

SigProfilerSimulator accepts four input file formats:

- **VCF** — one file per sample.
- **MAF** — standard Mutation Annotation Format file.
- **Simple text file** — tab-delimited plain text format as described in [SigProfilerMatrixGenerator][1].
- **ICGC Format** — ICGC simple somatic mutation format.

Input files must be placed in an `input/` subdirectory within the project folder.

### Running the function ###

First, start a Python interactive shell and import SigProfilerSimulator:

``` python
$ python
>>> from SigProfilerSimulator import SigProfilerSimulator as sigSim
```

Then call the function with the required parameters:

``` python
>>> sigSim.SigProfilerSimulator(project, project_path, genome, contexts)
```

### Required parameters ###

| Parameter | Variable Type | Parameter Description |
|-----------|---------------|-----------------------|
| `project` | String | Unique name for the given project |
| `project_path` | String | Path to the project directory. The `input/` subfolder containing the mutation files must exist within this path |
| `genome` | String | Reference genome to use. Must be installed using [SigProfilerMatrixGenerator][1]. Supported genomes: GRCh37, GRCh38, mm9, mm10, rn6, yeast |
| `contexts` | List of Strings | Mutational contexts to simulate. Must be provided as a list (e.g., `["96"]`, `["96", "ID"]`). See the full list of supported contexts below |

### Optional parameters ###

| Parameter | Variable Type | Parameter Description |
|-----------|---------------|-----------------------|
| `simulations` | Integer | Number of simulations to generate. Default: `1` |
| `exome` | Boolean | Restrict simulations to exome regions. Default: `None` (whole genome) |
| `chrom_based` | Boolean | Normalize mutation burden on a per-chromosome basis. Recommended when using the output as background model for [SigProfilerClusters][2]. Default: `False` |
| `gender` | String | Determines whether the Y chromosome is included. Accepted values: `"female"` (default, Y excluded), `"male"` (Y included) |
| `bed_file` | String | Path to a BED file to restrict simulations to user-defined genomic regions. Default: `None` |
| `vcf` | Boolean | Output simulated mutations as VCF files. When `False`, output is in MAF format. Default: `False` |
| `seqInfo` | Boolean | Save the sequence context information for each simulated mutation. Default: `False` |
| `seed_file` | String | Path to a file containing seeds for reproducible simulations. Default: `None` |
| `noisePoisson` | Boolean | Add Poisson-distributed noise to the simulated mutations. Default: `False` |
| `noiseUniform` | Float | Add uniform noise to the simulated mutations. Default: `0` |
| `spacing` | Integer | Minimum spacing (in bp) enforced between simulated mutations. Default: `1` |
| `cushion` | Integer | Cushion (in bp) around the edges of BED file regions within which mutations will not be placed. Default: `100` |
| `overlap` | Boolean | Allow simulated mutations to overlap. Default: `False` |
| `updating` | Boolean | Update mutation types during simulation. Default: `False` |
| `region` | String | Restrict simulations to a single chromosome (e.g., `"1"`). Default: `None` |
| `mask` | String | Path to a mask file to exclude specific genomic regions from simulations. Default: `None` |

### Supported contexts ###

| Mutation type | Accepted context values |
|---------------|------------------------|
| Single Base Substitutions (SBS) | `"6"`, `"24"`, `"96"`, `"288"`, `"384"`, `"1536"`, `"6144"` |
| Insertions and Deletions (ID) | `"ID"`, `"ID415"` |
| Double Base Substitutions (DBS) | `"DBS"`, `"DBS186"` |

  [1]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [2]: https://sigprofilersuite.github.io/SigProfilerClusters/
