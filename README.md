[![Build Status](https://travis-ci.org/RobersonLab/ngg2.svg?branch=master)](https://travis-ci.org/RobersonLab/ngg2)

# ngg2
Python package to identify NGGNGG Cas9 gRNA sites in any indexed FASTA file.

It is important to note that sites are identified using regular expressions (re). In standard mode, the sites are searched exhaustively.

In block scan mode, once characters are consumed in a match, they're gone. That means block scan only reports the first encountered gRNA site, *but* not second sites on the same strand that overlap it.

## Installation
is available via pip or GitHub download.
We **HIGHLY** recommend installing in a Python virtual environment.

```bash
pip install ngg2
```

Or user install

```bash
pip install --user ngg2
```

Or install from GitHub clone.

```bash
git clone https://github.com/RobersonLab/ngg2.git
git checkout vN.N.N # Choose highest version tag instead of vN.N.N

pip install -e .
```

## Depends on
* Python (tested on Python 2.7, 3.2-3.5)
* FASTA file of genome / contig of interest
* [pyfaidx](https://github.com/mdshw5/pyfaidx)
* [regex](https://pypi.python.org/pypi/regex)

## Usage

Find gRNA sites for in the first 10M bp of human chromosome 1.

```bash
ngg2 --outputFile myOutput.csv --region 1:1-10000000 human_genome.fa
```

Find all gRNA sites in a FASTA file.

```bash
ngg2 --outputFile myOutput.csv human_genome.fa
```

Find all gRNA sites in a FASTA file, but allow non-canonical (A, T, C) starting bases.

```bash
ngg2 --outputFile myOutput.csv --allowNoncanonical human_genome.fa
```

Find all gRNA sites in a FASTA file, using 10 processors

```bash
ngg2 --outputFile myOutput.csv --cores 10 human_genome.fa
```

Process in serial and write results to file, leaving out uniqueness of gRNA site

```bash
ngg2 --outputFile myOutput.csv --unbuffered human_genome.fa
```

Process in parallel, skipping uniqueness tests

```bash
ngg2 --outputFile myOutput.csv --cores 10 --skipUniqueScan human_genome.fa
```

Don't print log info

```bash
ngg2 --outputFile myOutput.csv --loglevel CRITICAL
```
