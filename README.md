# ctDNA_pypeliner

A pipeline to analyze low frequency mutations present in Circulating tumor DNA (ctDNA). Takes matched tumour normal fastq files and runs through analysis tools then creates a union from the results of each tool.

## Tools used

* Bioconductor deepSNV: https://bioconductor.org/packages/release/bioc/html/deepSNV.html
* LoLoPicker: https://github.com/jcarrotzhang/LoLoPicker
* VarScan: http://varscan.sourceforge.net/
* Mutationseq: https://github.com/shahcompbio/mutationseq
* Strelka2: https://github.com/Illumina/strelka

## Setup and Installation

Set up conda with the required packages.

First ensure you have the correct channels:

```
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/aroth85
conda config --add channels https://conda.anaconda.org/shahcompbio
conda config --add channels 'bioconda'
conda config --add channels 'r'
conda config --add channels 'conda-forge'
```

### From Source

Clone ctDNA_pypeliner:

```
git clone https://github.com/shahcompbio/ctDNA_pypeliner.git
cd ctDNA_pypeliner
```

Then create an environment with the required packages:

```
conda create --name ctDNApypeliner --file conda_packages.txt
```

Activate the environment:

```
source activate ctDNApypeliner
```

ctDNA_pypeliner uses LoLoPicker which is not available through conda.
To setup LoLoPicker:
```
git clone https://github.com/jcarrotzhang/LoLoPicker.git
cd LoLoPicker
python setup.py install
```

Add the ctDNA pipeline into the current site packages:
```
python setup.py install
```
## Usage

To run pipeline

```
ctdna_pypeliner /path/to/input.yaml /path/to/config.yaml --submit local
```

### Input
Pipeline accepts yaml as input. The yaml file contains input paths and metadata for each cell/patient. The format for each patient and their samples is as follows:
```
PATIENT_ID:
  normal:
    NORMAL_SAMPLE_ID:
      fastq1: /path/to/fastqfile/normal_sample_L001_R1_001.fastq.gz
      fastq2: /path/to/fastqfile/normal_sample_L001_R2_001.fastq.gz
      run: Run002
      sample_status: Normal
      type: saliva
  tumour:
    TUMOUR_SAMPLE1_ID:
      fastq1: /path/to/fastqfile/tumour_sample1_L001_R1_001.fastq.gz
      fastq2: /path/to/fastqfile/tumour_sample1_L001_R2_001.fastq.gz
      run: Run007
      sample_status: Treatment naive Primary Sx
      type: TNBC FFPE
    TUMOUR_SAMPLE2_ID:
      fastq1: /path/to/fastqfile/tumour_sample2_L001_R1_001.fastq.gz
      fastq2: /path/to/fastqfile/tumour_sample2_L001_R2_001.fastq.gz
      run: Run017
      sample_status: Post-Surgery/Post-Chemo
      type: plasma 18 month
    TUMOUR_SAMPLE3_ID:
      fastq1: /path/to/fastqfile/tumour_sample3_L001_R1_001.fastq.gz
      fastq2: /path/to/fastqfile/tumour_sample3_L001_R2_001.fastq.gz
      run: Run002
      sample_status: Pre-Surgery/Post-Chemo
      type: plasma baseline
```

## Config
Pipeline requires a config yaml with the following fields:
```
reference_genome: '/path/to/reference_genome/hg19.fa'
results_dir: '/path/to/results/'
bam_directory: '/path/to/bams/'
bed_file: '/path/to/beds/CG001v4.0.bed'
r_script_dir: '/path/to/ctDNA_pypeliner/ctDNA/r_scripts/'
```