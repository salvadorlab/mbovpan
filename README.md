# mbovpan
nextflow pipeline for Mycobacterium bovis pangenome analysis.

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


