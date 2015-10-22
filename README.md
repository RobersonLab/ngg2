# ngg2
Python script to identify NGGNGG Cas9 gRNA sites in any indexed FASTA file.

It is important to note that sites are identified using regular expressions (re). In standard mode, the sites are searched exhaustively.

In block scan mode, once characters are consumed in a match. That means block scan only reports the first encountered gRNA site, *but* not second sites on the same strand that overlap it.

## Installation
The script is not prepackaged, so just make sure it's in your PATH and executable. Tested on Windows 7 and Ubuntu 64-bit. Should be cross-platform compatible

## Requirements
* Python (tested on Python 2.7)
* FASTA file of genome / contig of interest
* [pyfaidx](https://github.com/mdshw5/pyfaidx)
* [regex](https://pypi.python.org/pypi/regex)

## Usage

Find gRNA sites for in the first 10M bp of human chromosome 1.

```bash
ngg2.py --outputFile myOutput.csv --region 1:1-10000000 human_genome.fa
```

Find all gRNA sites in a FASTA file.

```bash
ngg2.py --outputFile myOutput.csv human_genome.fa
```

Find all gRNA sites in a FASTA file, but allow non-canonical (A, T, C) starting bases.

```bash
ngg2.py --outputFile myOutput.csv --allowNoncanonical human_genome.fa
```

Find all gRNA sites in a FASTA file, using 10 processors

```bash
ngg2.py --outputFile myOutput.csv --cores 10 human_genome.fa
```

Process in serial and write results to file, leaving out uniqueness of gRNA site

```bash
ngg2.py --outputFile myOutput.csv --unbuffered human_genome.fa
```

Process in parallel, skipping uniqueness tests

```bash
ngg2.py --outputFile myOutput.csv --cores 10 --skipUniqueScan human_genome.fa
```

Don't print log info

```bash
ngg2.py --outputFile myOutput.csv --quiet human_genome.fa
```
