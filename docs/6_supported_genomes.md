# Supported Genomes


----------


This section lists all reference genomes currently supported by SigProfilerSimulator. Reference genomes must be installed using [SigProfilerMatrixGenerator][1] before running simulations.

----------

## Available Reference Genomes ##

| Genome ID | Assembly | Species | Source | Last updated |
|-----------|----------|---------|--------|--------------|
| `GRCh38` | GRCh38.p12 (GCA_000001405.27) | *Homo sapiens* | ENSEMBL v93.38 | January 2018 |
| `GRCh37` | GRCh37.p13 (GCA_000001405.14) | *Homo sapiens* | ENSEMBL v93.37 | September 2013 |
| `mm10` | GRCm38.p6 (GCA_000001635.8) | *Mus musculus* | ENSEMBL v93.38 | March 2018 |
| `mm9` | GRCm37 (GCA_000001635.18) | *Mus musculus* | ENSEMBL release 67 | March 2012 |
| `rn6` | Rnor_6.0 (GCA_000001895.4) | *Rattus norvegicus* | ENSEMBL v96.6 | January 2017 |
| `yeast` | R64-2-1 | *Saccharomyces cerevisiae* S288C | NCBI | November 2014 |

## Installation ##

Install a reference genome from the command line:

```
$ SigProfilerMatrixGenerator install GRCh37
```

Or from a Python terminal:

``` python
$ python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37', rsync=False, bash=True)
```

Multiple genomes can be installed independently. For full installation details, refer to the [Installation][2] section.

  [1]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [2]: https://sigprofilersuite.github.io/SigProfilerSimulator/1_installation.html
