# Quick Start Example


----------


This section provides a minimal example to get started with SigProfilerSimulator. The following example generates 100 simulations from a VCF input using the GRCh37 reference genome and the SBS-96 context.

----------

## Prerequisites ##

This tutorial requires that you have completed all steps in the [installation guide][1], specifically:

- Installed SigProfilerSimulator
- Downloaded the **GRCh37** reference genome using [SigProfilerMatrixGenerator][2]

## Input data ##

SigProfilerSimulator accepts four input file formats: VCF, MAF, simple text file, and ICGC format. Input files must be placed in an `input/` subdirectory within the project folder:

```
path/to/project/
└── input/
    ├── sample1.vcf
    └── sample2.vcf
```

## Running SigProfilerSimulator ##

Start a Python interactive shell and import SigProfilerSimulator:

``` python
$ python
>>> from SigProfilerSimulator import SigProfilerSimulator as sigSim
```

Run the simulator on your data. **Note**: Update `"path/to/project/"` with the actual path to your project directory.

``` python
>>> sigSim.SigProfilerSimulator("my_project", "path/to/project/", "GRCh37",
                                 contexts=["96"], simulations=100, chrom_based=True)
```

After SigProfilerSimulator has finished, the simulated mutation files will be placed in the `output/` subdirectory of your project folder, organized by context and simulation number.

## Additional Information ##

In the above example, unspecified parameters use their default values. All function arguments are described in detail in the [Using the Tool][3] section. For the full list of supported reference genomes, refer to the [Supported Genomes][4] section.

  [1]: https://sigprofilersuite.github.io/SigProfilerSimulator/1_installation.html
  [2]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [3]: https://sigprofilersuite.github.io/SigProfilerSimulator/3_using_the_tool.html
  [4]: https://sigprofilersuite.github.io/SigProfilerSimulator/4_supported_genomes.html
