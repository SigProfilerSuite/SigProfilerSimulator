# SigProfilerSimulator

![Logo](assets/images/SigProfilerSimulator.png)

----------

**SigProfilerSimulator** is a Python framework for generating realistic simulations of mutational signatures in cancer genomes. The tool can simulate single base substitutions (SBS), double base substitutions (DBS), and insertions/deletions (ID) across complete genomes or user-defined genomic regions, using an unbiased random distribution methodology that preserves the sequence context, transcriptional strand bias, and chromosomal mutational burden of the input data.

Simulated datasets are widely used as background models for downstream statistical analyses and hypothesis testing, including clustered mutation detection with [SigProfilerClusters][1].

**SigProfilerSimulator** makes use of [SigProfilerMatrixGenerator][2] and [SigProfilerPlotting][3], enabling seamless integration with other tools in the SigProfiler suite.

The SigProfilerSimulator library is available on [GitHub](https://github.com/SigProfilerSuite/SigProfilerSimulator) and [PyPI](https://pypi.org/project/SigProfilerSimulator).

----------

### Citation

Bergstrom EN, Barnes M, Martincorena I, Alexandrov LB. Generating realistic null hypothesis of cancer mutational landscapes using SigProfilerSimulator. *BMC Bioinformatics*. 2020;21(1):438. [https://doi.org/10.1186/s12859-020-03772-3](https://doi.org/10.1186/s12859-020-03772-3)

### License

This software is copyrighted © 2020 by Erik Bergstrom, Alexandrov Lab. SigProfilerSimulator is distributed under the terms of the GNU General Public License. The software is provided in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

### Contact

For questions, support requests, or bug reports, please contact the SigProfilerSuite team via GitHub [issues](https://github.com/SigProfilerSuite/SigProfilerSimulator/issues) or by email at [contact@sigprofilersuite.org](mailto:contact@sigprofilersuite.org).

  [1]: https://sigprofilersuite.github.io/SigProfilerClusters/
  [2]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [3]: https://github.com/SigProfilerSuite/SigProfilerPlotting
