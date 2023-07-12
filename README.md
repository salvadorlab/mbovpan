[![image search api](https://user-images.githubusercontent.com/40697188/193487621-a4b91a1c-19b6-42df-9e63-7fcff0658be0.png)](https://github.com/noahaus)

<iframe src="https://giphy.com/embed/7SF5scGB2AFrgsXP63" width="480" height="272" frameBorder="0" class="giphy-embed" allowFullScreen></iframe><p><a href="https://giphy.com/gifs/pokemon-anime-7SF5scGB2AFrgsXP63">via GIPHY</a></p>

# mbovpan
Mbovpan is a nextflow bioinformatic pipeline for Mycobacterium bovis pangenome analysis. The goal of Mbovpan is to make the insights from M. bovis genomics easily accesible to reseachers.  

Mbovpan can be ran in three separate modes: SNP mode for only inferring Single Nucleotide Polymorphisms, PAN mode for assessing gene presence absences, or BOTH to do a full analysis. Below you can witness the main workflow encapsulated by mbovpan. 
![image](https://github.com/salvadorlab/mbovpan/assets/40697188/315e9533-1567-48c1-aa0c-f1b5c12e2589)


### Installation  

The first step would be to download the pipeline using the following git command.
```
$git clone https://github.com/salvadorlab/mbovpan.git
```
Successfully downloading will lead to a directory with the following layout 

```
.
├── auxilary
├── chrom_ranges.py
├── dev_nf_scripts
│   ├── abyss_assem.nf
│   ├── accessory_pca.nf
│   ├── bowtie_call.nf
│   ├── calling_snps.nf
│   ├── gene_prab.nf
│   ├── megahit_assem.nf
│   ├── minimap_call.nf
│   ├── smalt_call.nf
│   ├── spades_assem.nf
│   ├── spotyping.nf
│   └── tbprofiler.nf
├── envs
│   ├── consensus.yaml
│   ├── freebayes.yaml
│   ├── gene_prab.yaml
│   ├── iqtree.yaml
│   ├── longread.yaml
│   ├── mbovpan.yaml
│   ├── megahit.yaml
│   ├── multiqc.yaml
│   ├── picard.yaml
│   ├── prokka.yaml
│   ├── quast.yaml
│   ├── raven.yaml
│   ├── roary.yaml
│   ├── samtools.yaml
│   ├── scoary.yaml
│   ├── statistics.yaml
│   ├── tbprofile.yaml
│   ├── vcflib.yaml
│   └── vcflib.yaml~
├── mbovpan.nf
├── README.md
├── ref
│   ├── mbovAF212297_annotation.gb
│   ├── mbovAF212297_reference.fasta
│   ├── mbovAF212297_reference.fasta.amb
│   ├── mbovAF212297_reference.fasta.ann
│   ├── mbovAF212297_reference.fasta.bwt
│   ├── mbovAF212297_reference.fasta.dict
│   ├── mbovAF212297_reference.fasta.fai
│   ├── mbovAF212297_reference.fasta.pac
│   ├── mbovAF212297_reference.fasta.sa
│   ├── mbov_bowtie_index.1.bt2
│   ├── mbov_bowtie_index.2.bt2
│   ├── mbov_bowtie_index.3.bt2
│   ├── mbov_bowtie_index.4.bt2
│   ├── mbov_bowtie_index.rev.1.bt2
│   ├── mbov_bowtie_index.rev.2.bt2
│   ├── mbov_reference.fasta
│   ├── pe_ppe_regions.gff3
│   └── spacer.fasta
├── scripts
│   ├── accessory_pca.R
│   ├── gene_prab.R
│   ├── install_mbovpan.sh
│   ├── mbov_virulence.py
│   ├── pangenome_curve.R
│   ├── SpoTyping
│   │   ├── README.md
│   │   ├── ref
│   │   │   └── spacer.fasta
│   │   ├── SpoTyping_plot.r
│   │   ├── SpoTyping.py
│   │   ├── SpoTyping-README-md.pdf
│   │   └── Update-history.md
│   ├── SpoTyping_copy.py
│   └── statistics.py
└── seqs
    ├── SRR998656_1.fastq.gz
    ├── SRR998656_2.fastq.gz
    ├── SRR998657_1.fastq.gz
    ├── SRR998657_2.fastq.gz
    ├── SRR998658_1.fastq.gz
    ├── SRR998658_2.fastq.gz
    ├── SRR998659_1.fastq.gz
    └── SRR998659_2.fastq.gz

```  

We encourage the user to add the 'mbovpan' directory to the PATH variable. This makes running the pipeline much easier.

```
#if NOT in the directory, give the absolute path
export PATH=$PATH:/path/to/mbovpan

#if already in the directory
export PATH=$PATH:$(pwd)
```  

After downloading, the user will need to create the mbovpan environment that will make it possible to run the pipeline. You can install by using a provided script in 'scripts/install_mbovpan.sh'. This should take a few minutes. 

```
$conda create -n mbovpan
$bash path/to/install_mbovpan.sh
$conda activate mbovpan #replace 'conda' with 'source' based on conda version
(mbovpan)$ #ready for input 
```

### Quickstart

Once the conda environment is created, the user can execute a simple run of the pipeline by providing an input path, an output path, and specific parameters. 

```
(mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ 
```
In this command, **'nextflow run'** is the command used to look at and execute the pipeline instructions in **'mbovpan.nf'**. This file contains the general flow of the mbovpan pipeline

Mbovpan is downloaded with test data already included, and this is captured with the **'--input ./mbovpan/seqs/'** parameter. If paired end sequences are present in the input directory, they will be matched and ran through the pipeline. 

**'--run snp'** signifies what analysis mode mbovpan utilizes. **'snp'** mode maps the paired end files to the reference genome while, **'pan'** mode creates de novo genomes from scratch. if no option is supplied, the pipeline will run both. 

### Additional Usages

```
# boost the number of threads utilized
(mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --threads 16

# modulate the minimum quality and maximum SNP depth required
(mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --qual 20 --depth 20

# Using most of the parameters that mbovpan has to offer to decipher the pangenome
(mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --qual 100 --depth 5 --mapq 50 --threads 30 --run pan
```

### Help
```
    M B O V P A N (v0.1)    
=============================
A pangenomic pipeline for the analysis of
Mycobacterium bovis isolates 


usage: nextflow run mbovpan/mbovpan.nf [options] --input ./path/to/input --output ./path/to/output
  options:
    --run [all|snp|pan]: 
        Specifies in what mode to run mbovpan in [DEFAULT:all]
    --qual [INT]:
        The minimum QUAL score for a SNP to be considered [DEFAULT:150]
    --depth [INT]:
        The minimum DP score for a SNP to be considered [DEFAULT:10]
    --mapq [INT]:
        The minimum MQ score for a SNP to be considered [DEFAULT:55]
    --threads [INT]:
        How many threads to use for the programs [DEFAULT:(number of avail. threads)/2]
    --help
        Prints this help message
    --version
        Prints the current version 
=============================
```

### Authors
mbovpan is software initiated and maintained by Noah Legall, PhD under the supervision of Liliana C. M. Salvador, PhD.

### Maintainers 
Noah Legall 
    [![image search api](https://github.com/salvadorlab/mbovpan/assets/40697188/89705d13-4265-4cda-98d3-d21064344892)](https://github.com/noahaus)



### Citation 
The mbovpan manuscript is currently in preparation. 


