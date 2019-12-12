# ctDNA_pypeliner

A pipeline to find high depth, low frequency mutations present in Circulating tumor DNA (ctDNA). Takes matched tumour normal fastq files and runs through somatic variant calling tools then creates an annotated union from the results of each tool. To reduce noise and false positive calls, only sites called by two or more tools are recorded in the output.

## Tools used

* Pypeliner (workflow manager): https://pypeliner.readthedocs.io
* Trim Galore (for UMI snipping): https://github.com/FelixKrueger/TrimGalore
* Burrows-Wheeler Aligner (for alignment): http://bio-bwa.sourceforge.net/
* samtools (for alignment): https://www.htslib.org/
* vcftools (for vcf merging and sorting): https://vcftools.github.io
* Annovar (for annotation): http://annovar.openbioinformatics.org
* bam-readcount (to fetch raw readcount): https://github.com/genome/bam-readcount

### For somatic variant calling

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

ctDNA_pypeliner uses LoLoPicker and Annovar which are not available through conda.
To setup LoLoPicker:
```
git clone https://github.com/jcarrotzhang/LoLoPicker.git
cd LoLoPicker
python setup.py install
```

To setup Annovar, follow the download procedures here: http://annovar.openbioinformatics.org/en/latest/user-guide/download/

Add the ctDNA pipeline into the current site packages:
```
python setup.py install
```
## Usage

To run pipeline

Locally
```
ctdna_pypeliner --input_yaml /path/to/input.yaml --config /path/to/config.yaml --map_q 25 --base_q 15 --p_cutoff 0.0005 --T_vaf_cutoff 0.004 --N_vaf_cutoff 0 --museq_cutoff 0.65 --tmpdir /path/to/tmp/ --pipelinedir /path/to/pipeline/ --submit local --maxjobs 4
```

On shahlab cluster
```
ctdna_pypeliner --input_yaml /path/to/input.yaml --config /path/to/config.yaml --map_q 25 --base_q 15 --p_cutoff 0.0005 --T_vaf_cutoff 0.004 --N_vaf_cutoff 0 --museq_cutoff 0.65 --tmpdir /path/to/tmp/ --pipelinedir /path/to/pipeline/ --submit asyncqsub --nativespec ' -V -hard -q shahlab.q -l h_vmem={mem}G -P shahlab_high -S /bin/bash' --maxjobs 128 --context_config /path/to/context.yaml
```

### Options
```
  --input_yaml INPUT_YAML     Input filename
  --config CONFIG             Configuration filename
  --map_q MAP_Q               Minimum mapping quality
  --base_q BASE_Q             Minimum base quality
  --p_threshold P_VALUE       Maximum p_value
  --museq_threshold MUSEQ     Minumum MutationSeq score
  --N_vaf_threshold N_VAF     Maximum normal variant allele frequency
  --T_vaf_threshold T_VAF     Minimum tumour variant allele frequency
  --umi_trim                  Set flag to true to trim trailing and leading UMI from reads
```
Pypeliner workflow arguments
```
  --tmpdir TMPDIR       location of temporary files
  --pipelinedir PIPELINEDIR
                        location of pipeline files
  --loglevel {CRITICAL,ERROR,WARNING,INFO,DEBUG}
                        logging level for console messages
  --submit SUBMIT       job submission system
  --submit_config SUBMIT_CONFIG
                        job submission system config file
  --nativespec NATIVESPEC
                        qsub native specification
  --storage STORAGE     file storage system
  --storage_config STORAGE_CONFIG
                        file storage system config file
  --maxjobs MAXJOBS     maximum number of parallel jobs
  --repopulate          recreate all temporaries
  --rerun               rerun the pipeline
  --nocleanup           do not automatically clean up temporaries
  --interactive         run in interactive mode
  --sentinel_only       no timestamp checks, sentinal only
  --context_config CONTEXT_CONFIG
                        container registry credentials and job context
                        overrides
  ```

### Input
Pipeline accepts yaml as input. The yaml file contains input paths and metadata for each cell/patient. Example input files can be found in the directory `test_files`. The format for each patient and their samples is as follows:
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

### Config
The config file provides supplemental information required by the pipeline. Example config files can be found in the directory `config`. The config file is a yaml with the following fields:
```
reference_genome: '/path/to/reference/GRCh37-lite.fa'
results_dir: '/path/to/results/'
bam_directory: '/path/to/bams/'
bed_file: '/path/to/beds/GRCh37.bed'
r_script_dir: '/path/to/ctDNA_pypeliner/ctDNA/r_scripts/'
museq_python: '/shahlab/pipelines/virtual_environments/museq_pipeline/bin/python'
museq_classify: '/path/to/mutationSeq_4.3.7_python2.7/classify.py'
museq_deep_model: '/path/to/mutationSeq_4.3.7_python2.7/model_deep_v0.2.npz'
museq_config: '/path/to/mutationSeq_4.3.7_python2.7/metadata.config'
snv_vcf_template: '/path/to/ctDNA_pypeliner/template_snv.vcf'
indel_vcf_template: '/path/to/ctDNA_pypeliner/template_indel.vcf'
annovar: '/path/to/annovar/'
annovar_humandb: '/path/to/annovar/humandb/'
```

## Output
Output files are created in the `results_dir` specified in the configuration file and are organized into subdirectories according to patient_id.
For each patient tumour sample, ctDNA_pypeliner will output 3 files for SNVs and 3 files for INDELS. One TSV, one TXT, and one VCF for each type of mutation. The TSV and VCF output files contain scoring information from each respective tool. Those scores are specified as follows:

* DeepSNV: The corrected p-value
* LoLoPicker: The corrected p-value
* MutationSeq: The probability score of mutation
* Varscan: The somatic p-value for Somatic/LOH events
* Strelka: Quality score for any somatic snv

Only VarScan and Strelka call INDEL mutations, so indel output files will only contain VarScan and Strelka scoring

### Tab-Separated Values (TSV) file
The output TSV file clearly states the score factors of each SNV calling tool, read depth, and allele counts/frequencies at each mutation site. This output is unannotated.

Below is a sample header from an SNV TSV file:
```
chr  pos ref alt count deepSNV LoLoPicker  VarScan MutationSeq Strelka N_coverage  N_A N_C N_G N_T N_N N_vaf T_coverage  T_A T_C T_G T_T T_N T_vaf
```

INDEL TSV file:
```
chr pos ref alt VarScan Strelka N_coverage  N_ref N_alt N_vaf T_coverage  T_ref T_alt T_vaf
```

### TXT file
The output TXT file contains annotation information of each mutation site. It does not contain scoring information from the tools, nor does it contain information about read depth and allele count/frequency.

Below is a sample header from an output TXT file:
```
Chr Start End Ref Alt Func.refGene  Gene.refGene  GeneDetail.refGene  ExonicFunc.refGene  AAChange.refGene  cytoBand
```

### VCF file
The output VCF file contains both annotation information as well as the information found in the TSV output file (scoring information, read depth, allele counts and frequencies)

Below is a sample meta information (header) from an output VCF file:
```
##fileformat=VCFv4.1
##source=ctDNA_pypeliner
##INFO=<ID=COUNT,Number=1,Type=Integer,Description="Number of consensus calls">
##INFO=<ID=DSNV,Number=1,Type=Float,Description="p-value from DeepSNV">
##INFO=<ID=LLP,Number=1,Type=Float,Description="p-value from LoLoPicker">
##INFO=<ID=VS,Number=1,Type=Float,Description="p-value from VarScan">
##INFO=<ID=MS,Number=1,Type=Float,Description="Probability score from MutationSeq">
##INFO=<ID=STR,Number=1,Type=Float,Description="Quality Score from Strelka">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AU,Number=1,Type=Integer,Description="Read count of nucleotide A">
##FORMAT=<ID=CU,Number=1,Type=Integer,Description="Read count of nucleotide C">
##FORMAT=<ID=GU,Number=1,Type=Integer,Description="Read count of nucleotide G">
##FORMAT=<ID=TU,Number=1,Type=Integer,Description="Read count of nucleotide T">
##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Read count of null">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,Description="Flag the start of ANNOVAR annotation for one alternative allele">
##INFO=<ID=Func.refGene,Number=.,Type=String,Description="Func.refGene annotation provided by ANNOVAR">
##INFO=<ID=Gene.refGene,Number=.,Type=String,Description="Gene.refGene annotation provided by ANNOVAR">
##INFO=<ID=GeneDetail.refGene,Number=.,Type=String,Description="GeneDetail.refGene annotation provided by ANNOVAR">
##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description="ExonicFunc.refGene annotation provided by ANNOVAR">
##INFO=<ID=AAChange.refGene,Number=.,Type=String,Description="AAChange.refGene annotation provided by ANNOVAR">
##INFO=<ID=cytoBand,Number=.,Type=String,Description="cytoBand annotation provided by ANNOVAR">
##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description="Flag the end of ANNOVAR annotation for one alternative allele">
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  NORMAL  TUMOR
```


