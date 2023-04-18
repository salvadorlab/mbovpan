
// Where are we grabbing our sequence data from?     
input = ""

// Where should results be published to?
output = ""

// Are these paired-end short reads or longreads?
mode = "short"

// Are we computing SNPs or computing the pangenome? 
run_mode = "snp"

// How many threads will be available to run the pipeline. 
// Automatically uses all the cpus that are available 
// If not specified, use 50% of available resources 
threads = Math.floor(Runtime.getRuntime().availableProcessors()/2)

reads = ""

qual = 150

depth = 10 

mapq = 55

// Provide path to metadata for analysis
// isolate name as rows with metadata as columns
meta = ""

scoary_meta = ""

if(params.qual != null){
    qual = params.qual as Integer
    }
if(params.depth != null){
    depth = params.depth as Integer
    }
if(params.mapq != null){
    mapq = params.mapq as Integer
    }

if(params.meta != null){
    meta = params.meta 
    println "metadata loaded successfully"
}

if(params.scoary_meta != null){
    scoary_meta = params.scoary_meta 
    println "scoary metadata loaded successfully"
}

// record the path for the M. bovis reference genome
ref = "$workflow.projectDir/../ref/mbovAF212297_reference.fasta"
range = "$workflow.projectDir/../auxilary/chrom_ranges.txt" 

// are default parameters included?
if(params.input == null || params.output == null){
    println "necessary paramaters aren't supplied - supply values for input and/or output"
    exit 0
}

else {
    println "necessary paramaters supplied"
    input = params.input
    output = params.output
}

// record the mode the program should be ran in

    println "mbovpan will run in ${mode} read mode"


// what part of the pipeline should be ran?
if(params.run == "all" ){
    println "mbovpan will infer snps and pangenome"
    run_mode = "all"
}

else if(params.run == "snp"){
    println "mbovpan will only infer snps"
    run_mode = "snp"
}

else if(params.run == "pan"){
    println "mbovpan will only infer the pangenome"
    run_mode = "pan"
}

else {
    println "mbovpan will infer both snps and pangenome by default"
    run_mode = "all"
}

// how many threads will be utilized
if(params.threads != null){
    println "mbovpan will run using ${params.threads} threads"
    threads = params.threads
}
else{
    println "mbovpan will run using ${threads} threads by default"
}

reads = Channel.fromFilePairs("$input*{1,2}*.f*q*").ifEmpty { error "Cannot find the read files" }


reads.into {
    reads_process
    reads_trim
}

 // PART 1: Processing 


/* AUTOMATIC QC of read data */



//change this around to get rid of the weird naming
    process fastp {

    publishDir = output

    cpus threads

    input:
    tuple sample_id, file(reads_file) from reads_trim

    output:
    file("${sample_id}_trimmed_R*.fastq") into fastp_ch
 
    script:
    """
    fastp -w ${task.cpus} -q 30 --detect_adapter_for_pe -i  ${reads_file[0]} -I  ${reads_file[1]} -o  ${sample_id}_trimmed_R1.fastq -O  ${sample_id}_trimmed_R2.fastq
    """

    }

    fastp_ch.into {
        fastp_reads1
        fastp_reads2
        fastp_reads3
        fastp_reads4
    }



// MODE 1: Variant Calling 
bam = Channel.create()

if(run_mode == "snp" || run_mode == "all"){

    process read_map {
    publishDir = output 

    cpus threads
    
    conda "$workflow.projectDir/../envs/samtools.yaml"
   
    input:
    tuple file(trim1), file (trim2) from fastp_reads2 

    output:
    file("${trim1.baseName - ~/_trimmed_R*/}.bam")  into map_ch 

    script:
    """
    bwa mem -t ${task.cpus}-M -R "@RG\\tID:${trim1.baseName - ~/_trimmed_R*/}\\tSM:${trim1.baseName - ~/_trimmed_R*/}\\tPL:ILLUMINA\\tPI:250" ${ref} ${trim1} ${trim2} | samtools view -Sb | samtools sort -o ${trim1.baseName - ~/_trimmed_R*/}.bam
    """
    }

    bam = map_ch

  

// Important to have 'USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
    process mark_dups {
    publishDir = output 

    conda "$workflow.projectDir/../envs/picard.yaml"    

    input:
    file(bam) from bam

    output:
    file("${bam.baseName}.nodup.bam") into nodup_ch

    script:
    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam.baseName}.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    samtools index ${bam.baseName}.nodup.bam
    cp ${bam.baseName}.nodup.bam.bai $workflow.launchDir
    """
    }

    nodup_ch.into {
        nodup1_ch
        nodup2_ch
    }

    process freebayes {
    publishDir = output 

    cpus threads

    conda "$workflow.projectDir/../envs/freebayes.yaml"

    input:
    file(bam) from nodup1_ch
    //file(range) from chrom_range
    path(reference) from ref

    output:
    file("${bam.baseName - ~/.nodup/}.vcf") into freebayes_ch 

    script:
    """
    cp $workflow.launchDir/${bam.baseName - ~/.nodup/}.nodup.bam.bai ./
    freebayes-parallel ${range} ${task.cpus} -f ${ref} ${bam} > ${bam.baseName - ~/.nodup/}.vcf
    """
    }

//Start with a basic QUAL > 150, later add a parameter for changing this
    process vcf_filter {
    publishDir = output 

    conda "$workflow.projectDir/../envs/vcflib.yaml"

    input:
    file(vcf) from freebayes_ch

    output:
    file("${vcf.baseName}.filtered.vcf") into filter_ch

    script:
    """
    vcffilter -f "QUAL > ${qual}" ${vcf} | vcffilter -f "DP > ${depth}" | vcffilter -f "MQM > ${mapq}" |  vcffilter -f "TYPE = snp" | bedtools intersect -header -a - -b $workflow.projectDir/../ref/pe_ppe_regions.gff3 -v > ${vcf.baseName}.filtered.vcf
    """

    } 

}