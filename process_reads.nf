// Author: Noah Austin Legall
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

// Need to figure this out. 
// The path is not recognized if given straight to the program, yet works fine when its in relation to the $workflow.projectDir
// my first guess is that the files ALWAYS need a path
ref = "$workflow.projectDir/../mbovAF212297_reference.fasta"



reads = Channel.fromFilePairs("$input*_{1,2}.fastq.gz").ifEmpty { error "Cannot find any reads" }.view()


reads.into {
    reads_fastqc
    reads_trim
}


 // PART 1: Processing 

process pre_fastqc {

    publishDir = output

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

//change this around to get rid of the weird naming
process fastp {

    publishDir = output

    input:
    tuple sample_id, file(reads_file) from reads_trim

    output:
    file("${sample_id}.trimmed_R*.fastq") into fastp_ch
 
    script:
    """
    fastp -q 30 --detect_adapter_for_pe -i  ${reads_file[0]} -I  ${reads_file[1]} -o  ${sample_id}.trimmed_R1.fastq -O  ${sample_id}.trimmed_R2.fastq
    """

}

fastp_ch.into {
    fastp_reads1
    fastp_reads2
    fastp_reads3
}

process post_fastqc {

    publishDir = output

    input:
    tuple file(trim1), file(trim2) from fastp_reads1

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

// PART 2: Variant Calling 
/*
process setup {
    publishDir = output

    input:
    path(reference) from ref

    output:
    file("${reference}*")

    script:
    """
    echo 'performing setup for variant calling'
    samtools faidx ${reference}
    picard CreateSequenceDictionary REFERENCE=${reference} OUTPUT=${reference}.dict 
    bwa index ${reference} 
    echo 'setup complete'
    """
}

// now ready to parallelize
process freebayes_setup {
    publishDir = output 

    input:
    path(reference) from ref

    output:
    file("chrom_ranges.txt") into chrom_range

    script:
    """
    python $workflow.projectDir/chrom_ranges.py  ${reference}
    """

}

process read_map {
    publishDir = output 

    input:
    tuple file(trim1), file (trim2) from fastp_reads2 
    //path(reference) from ref

    output:
    file("${trim1.baseName}.bam")  into map_ch 

    script:
    """
    bwa mem -M -R "@RG\\tID:${trim1.baseName}\\tSM:${trim1.baseName}\\tPL:ILLUMINA\\tPI:250" ${ref} ${trim1} ${trim2} | samtools view -Sb | samtools sort -o ${trim1.baseName}.bam
    """
}

// picard MarkDuplicates INPUT={sorted_bamfile} OUTPUT={rmdup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv

// Important to have 'USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
process mark_dups {
    publishDir = output 

    input:
    file(bam) from map_ch

    output:
    file("${bam.baseName}.nodup.bam") into nodup_ch

    script:
    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam.baseName}.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    samtools index ${bam.baseName}.nodup.bam
    """
}

process freebayes {
    publishDir = output 

    cpus 2

    input:
    file(bam) from nodup_ch
    file(range) from chrom_range
    path(reference) from ref

    output:
    file("${bam.baseName}.vcf") into freebayes_ch 

    script:
    """
    freebayes-parallel ${range} ${task.cpus} -f ${reference} ${bam} > ${bam.baseName}.vcf
    """
}
*/

// vcf filtering + generate alignment? 
// pangenome steps (might need to separately create environment for pangenome software. inclusion in main environment might lead to conflicts)
// roary, quast OR busco

// PART 3: Pangenome 

process assembly {
    publishDir = output 

    cpus 3

    input:
    tuple file(trim1), file(trim2) from fastp_reads3

    output:
    file("${trim1.baseName}.scaffold.fasta") into assembly_ch

    script:
    """
    spades.py --threads ${task.cpus} --only-assembler --careful -1 ${trim1} -2 ${trim2} -o ${trim1.baseName}
    cd ${trim1.baseName}
    mv scaffolds.fasta  ../${trim1.baseName}.scaffold.fasta
    """
}

// slight problem with prokka downloaded from conda. 
// might require the use of pre-made environment. 
// tomorrow, take time to import this environment

process annotate {
    publishDir = output

    input:
    file(assembly) from assembly_ch

    output:
    file("${assembly.baseName}.annot.gff") into annotate_ch

    script:
    """
    prokka  --outdir ./${assembly.baseName} --prefix ${assembly.baseName}.annot ${assembly}
    """
}