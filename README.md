# ngg2
Python script to identify NGGNGG Cas9 gRNA sites in any indexed FASTA file.

## Installation
The script is not prepackaged, so just make sure it's in your PATH and executable. Tested on Windows 7 and Ubuntu 64-bit. Should be cross-platform compatible

## Requirements
* Python (tested on Python 2.7)
* FASTA file of genome / contig of interest
* samtools or other tool to index FASTA file, if it is not already indexed

## Usage

Find gRNA sites for in the first 10M bp of human chromosome 1.

```bash
ngg2.py --outputFile myOutput.csv --region chr1:1-10000000 human_genome.fa
```

Find all gRNA sites in a FASTA file.

```bash
ngg2.py --outputFile myOutput.csv --allContigs human_genome.fa
```

Find all gRNA sites in a FASTA file, but only keep canonical G starting sites.

```bash
ngg2.py --outputFile myOutput.csv --allContigs --onlyGstart human_genome.fa
```

Find all gRNA sites in a FASTA file, only keep canonical sites, but allow N reference bases.

```bash
ngg2.py --outputFile myOutput.csv --allContigs --onlyGstart --allowN human_genome.fa
```

Don't print log info

```bash
ngg2.py --outputFile myOutput.csv --allContigs --quiet human_genome.fa
```
