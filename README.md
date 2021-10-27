# PIPi
(P)estivirus (I)dentification (P)ipeline

## Dependencies
The following are required for PIPi to function

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bbduk.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
- [bwa](https://github.com/lh3/bwa)
- [samtools=1.13](http://www.htslib.org/)
- [pandas](https://pandas.pydata.org/)
- umpy](https://numpy.org/)
- atplotlib](https://matplotlib.org/)

## Installing dependencies
This can be easily achieved in a conda virtual environment
```
conda create -n PIPi fastqc bbmap bwa samtools=1.13 pandas numpy matplotlib
```

## Quick Usage
PIPi can be called as follows.
```
./PIPi.py -p <path/to/reads> -o <path/to/outdir> -t <num_threads> -r <path/to/reference.fasta>
```

## Options and Usage
```
Pestivirus Identification Pipeline

Required Arguments:
  -p PATH, --path PATH  Folder of fastq files from Ion torrent
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -r REFERENCE, --reference REFERENCE
                        Location of reference Pestivirus sequences

Optional Arguments:
  -t THREADS, --threads THREADS
                        Number of threads to use. Default 8
  -h, --help            show this help message and exit
```
