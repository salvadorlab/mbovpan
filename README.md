# mbovpan
Mbovpan is a nextflow bioinformatic pipeline for Mycobacterium bovis pangenome analysis. The goal of Mbovpan is to make the insights from M. bovis genomics easily accesible to reseachers.  

Mbovpan can be ran in three separate modes: SNP mode for only inferring Single Nucleotide Polymorphisms, PAN mode for assessing gene presence absences, or BOTH to do a full analysis. 

![image](https://user-images.githubusercontent.com/40697188/191386250-52b8a354-5611-44b5-8055-db29337cbe31.png)


### Quick Start 
```
$git clone https://github.com/salvadorlab/mbovpan.git
```

After downloading, the user will need to create the mbovpan environment that will make it possible to run the pipeline. The specifications in the environment yaml file will work across Linux and OS X operating systems. 

```
$conda env create -f mbovpan/envs/mbovpan.yaml 
$conda activate mbovpan 
(mbovpan)$ #ready for input
```

Once the conda environment is created, the user can execute a simple run of the pipeline by providing an input path, an output path, and specific parameters. 

```
(mbovpan)$nextflow run mbovpan/mbovpan.nf --input mbovpan/seqs/ --run snp --output path/to/output 
```

Mbovpan is downloaded with test data already included, therefore the above command should run with no problems and accept input from the user’s own curated dataset.


