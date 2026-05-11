# Installation


----------


This section will help you set up the necessary software and packages required to run SigProfilerSimulator.

----------


## Prerequisites ##

- [Python][1] version >= 3.9
- [SigProfilerMatrixGenerator][2] with a downloaded reference genome
- Other dependencies are installed automatically during package installation

## Installation ##

SigProfilerSimulator can be executed on any Windows/macOS/Unix system. First follow the [SigProfilerMatrixGenerator][2] guide for installing `Python` and `pip`. Next, follow the instructions below for the latest stable release or the current GitHub version.

### Installation with `pip` ###

Install the latest `SigProfilerSimulator` PyPI version using `pip`:
```
$ pip install SigProfilerSimulator
```

To upgrade an existing installation to the most recent version:
```
$ pip install SigProfilerSimulator --upgrade
```

### Install specific GitHub Release ###

First, download the [zip file][3] or clone the GitHub repository:
```
$ git clone https://github.com/SigProfilerSuite/SigProfilerSimulator.git
```

Next, enter the downloaded directory and install the package:
```
$ cd SigProfilerSimulator
$ pip install .
```

## Download Reference Genome ##

SigProfilerSimulator requires a reference genome to perform simulations. To install the reference genome/s, use [SigProfilerMatrixGenerator][2].

The last PyPI [SigProfilerMatrixGenerator][2] version is installed with SigProfilerSimulator by default. Install your desired reference genome from the command line/terminal as follows.

### Installation from command line ###

```
$ SigProfilerMatrixGenerator install GRCh37
```

### Installation from Python terminal ###

``` python
$ python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37', rsync=False, bash=True)
```

If you have a firewall on your server, you may need to install `rsync` and use the `rsync=True` parameter. If bash is not available, use `bash=False`.

For a full list of supported reference genomes, refer to the [Supported Genomes][4] section.

  [1]: https://www.python.org/downloads
  [2]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/
  [3]: https://github.com/SigProfilerSuite/SigProfilerSimulator/releases
  [4]: https://sigprofilersuite.github.io/SigProfilerSimulator/6_supported_genomes.html
