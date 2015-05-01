# ngg2
Python script to identify NGGNGG Cas9 gRNA sites in any indexed FASTA file.

It is important to note that sites are identified using regular expressions (re). Once characters are consumed in a match, they aren't reused. As such, the tool reports the first encountered gRNA site, *but* not second sites on the same strand that overlap the first.

## Installation
The script is not prepackaged, so just make sure it's in your PATH and executable. Tested on Windows 7 and Ubuntu 64-bit. Should be cross-platform compatible

## Requirements
* Python (tested on Python 2.7)
* FASTA file of genome / contig of interest
* [pyfaidx](https://github.com/mdshw5/pyfaidx)

## Usage

Find gRNA sites for in the first 10M bp of human chromosome 1.

```bash
ngg2.py --outputFile myOutput.csv --region chr1:1-10000000 human_genome.fa
```

Find all gRNA sites in a FASTA file.

```bash
ngg2.py --outputFile myOutput.csv human_genome.fa
```

Find all gRNA sites in a FASTA file, but only keep canonical G starting sites.

```bash
ngg2.py --outputFile myOutput.csv --onlyGstart human_genome.fa
```

Find all gRNA sites in a FASTA file, only keep canonical sites, but allow N reference bases.

```bash
ngg2.py --outputFile myOutput.csv --onlyGstart --allowN human_genome.fa
```

Don't print log info

```bash
ngg2.py --outputFile myOutput.csv --quiet human_genome.fa
```
