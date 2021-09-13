//Author: Noah Austin Legall
// Note to self: do NOT run this program within the github folder
// test outside of the github folder. we do not want to make the github too full with run information 


/*
TODOS:
- finish full pipeline
- test on multiple 
- incorporate threading options
- call SNPs following vSNP pipeline


*/

log.info """ 
                   M B O V P A N (v0.1)    
             =============================
             A pangenomic pipeline for the analysis of
             Mycobacterium bovis isolates 

             Project : $workflow.projectDir
             Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
             Cmd line: $workflow.commandLine
             Manifest's pipeline version: $workflow.manifest.version
"""

// Define both the input and the output         
input = "$params.input"
output = "$params.output"



// Based on the input param, where to find the data.
// the 'seqs' folder contains a toy dataset to test/make sure the pipeline runs properly 
reads = Channel.fromFilePairs("$input/*_{1,2}.fastq.gz").ifEmpty { error "Cannot find any reads" }.view()


reads.into {
    reads_fastqc
    reads_trim
}

// Split fastqc operations into 'pre_fastqc' and 'post_fastqc'
// each process outputs information into the output directory
// also, each process has it's own environment to run (found in 'envs' directory). This will help users across different devices be able to have all dependencies
process pre_fastqc {

    publishDir = output

    conda "$workflow.projectDir/envs/fastqc.yaml"

    input:
    tuple sample_id, file(reads_file) from reads_fastqc

    output:
    file("pre_fastqc_${sample_id}_logs") into fastqc_ch1

    script:
    """
    mkdir  pre_fastqc_${sample_id}_logs
    fastqc -o  pre_fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

// This process runs the fastp command on our fastq files. 
// fastp is an 'all-in-one' trimming software that performs automatic adapter removal 
process fastp {

    publishDir = output

    conda "$workflow.projectDir/envs/fastp.yaml"

    input:
    tuple sample_id, file(reads_file) from reads_trim

    output:
    file("${sample_id}.trimmed_R*.fastq") into fastp_ch
 
    script:
    """
    fastp -q 30 --detect_adapter_for_pe -i  ${reads_file[0]} -I  ${reads_file[1]} -o  ${sample_id}.trimmed_R1.fastq -O  ${sample_id}.trimmed_R2.fastq
    """

}

// run this only after the read trimming process
// identical to the pre_fastqc process
process post_fastqc {

    publishDir = output

    conda "$workflow.projectDir/envs/fastqc.yaml"

    input:
    tuple file(trim1), file(trim2) from fastp_ch

    output:
    file("post_fastqc_${trim1.baseName}_logs") into fastqc_ch2

    script:
    """
    mkdir  post_fastqc_${trim1.baseName}_logs
    fastqc -o  post_fastqc_${trim1.baseName}_logs -f fastq -q ${trim1} ${trim2}
    """
}


 process multiqc {
    publishDir = output

   
    input:
    file(pre) from fastqc_ch1.collect()
    file(post) from fastqc_ch2.collect()

    output:
    file("mbovpan_seq_qual*")

    script:
    """
    multiqc -n mbovpan_seq_qual ${pre} ${post}
    """
} 



