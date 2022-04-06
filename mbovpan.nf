// Author: Noah Austin Legall
// Note to self: do NOT run this program within the github folder
// test outside of the github folder. we do not want to make the github too full with run information 

/*
Finish the pipeline
Extend Figure 2 with runtime reports
Testing in sapelo2!!!
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
=============================
"""

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

if(params.qual != null){
    qual = params.qual as Integer
    }
if(params.depth != null){
    depth = params.depth as Integer
    }
if(params.mapq != null){
    mapq = params.mapq as Integer
    }

// record the path for the M. bovis reference genome
ref = "$workflow.projectDir/ref/mbovAF212297_reference.fasta"
range = "$workflow.projectDir/auxilary/chrom_ranges.txt" 

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
    println "mbovpan will run using ${threads} threads"
    threads = params.threads
}
else{
    println "mbovpan will run using ${threads} threads by default"
}

reads = Channel.fromFilePairs("$input*{1,2}*.fastq*").ifEmpty { error "Cannot find the read files" }


reads.into {
    reads_process
    reads_trim
}

 // PART 1: Processing 

 log.info """ 
Summary of pipeline run

mode: $mode
run type: $run_mode
reference location: $ref
input: $input
output: $output
no. of threads: $threads
=====================================
"""

/* AUTOMATIC QC of read data */

    process pre_fastqc {

    publishDir = output

    input:
    tuple sample_id, file(reads_file) from reads_process

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


    process post_fastqc {

    publishDir = output

    input:
    tuple file(trim1), file(trim2) from fastp_reads1

    output:
    file("post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs") into fastqc_ch2

    script:
    """
    mkdir  post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs
    fastqc -o  post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs -f fastq -q ${trim1} ${trim2}
    """
    }

// MODE 1: Variant Calling 
bam = Channel.create()

if(run_mode == "snp" || run_mode == "all"){

    process read_map {
    publishDir = output 

    cpus threads
    
    conda "$workflow.projectDir/envs/samtools.yaml"
   
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

  

// picard MarkDuplicates INPUT={sorted_bamfile} OUTPUT={rmdup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv

// Important to have 'USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
    process mark_dups {
    publishDir = output 

    conda "$workflow.projectDir/envs/picard.yaml"    

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

    conda "$workflow.projectDir/envs/freebayes.yaml"

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

    conda "$workflow.projectDir/envs/vcflib.yaml"

    input:
    file(vcf) from freebayes_ch

    output:
    file("${vcf.baseName}.filtered.vcf") into filter_ch

    script:
    """
    vcffilter -f "QUAL > ${qual}" ${vcf} | vcffilter -f "DP > ${depth}" | vcffilter -f "MQM > ${mapq}" |  bedtools intersect -header -a - -b $workflow.projectDir/ref/pe_ppe_regions.gff3 -v > ${vcf.baseName}.filtered.vcf
    """

    } 

    filter_ch.into {
        filter1_ch
        filter2_ch
    }

    stats_ch = fastp_reads4.merge(nodup2_ch).merge(filter2_ch).view()

    process stats {
        publishDir = output

        conda "$workflow.projectDir/envs/statistics.yaml"

        input:
        file(nec_files) from stats_ch

        output:
        file("${nec_files[0].baseName}.stats") into output_stat_ch


        script:
        """
        python $workflow.projectDir/scripts/statistics.py ${nec_files[0]} ${nec_files[1]} ${nec_files[2]} ${nec_files[3]} > ${nec_files[0].baseName}.stats
        """

    }

    process psuedo_assembly {
        publishDir = output

        conda "$workflow.projectDir/envs/consensus.yaml"

        input:
        file(vcf) from filter1_ch 

        output:
        file("${vcf.baseName}.consensus.fasta") into fasta_ch

        script:
        """
        bgzip ${vcf} 
        bcftools index ${vcf}.gz
        cat ${ref} | vcf-consensus ${vcf.baseName}.vcf.gz > ${vcf.baseName}.dummy.fasta
        sed 's|LT708304.1 Mycobacterium bovis AF2122/97 genome assembly, chromosome: Mycobacterium_bovis_AF212297|${vcf.baseName}|g' ${vcf.baseName}.dummy.fasta > ${vcf.baseName}.consensus.fasta
        """
    }
    
    process iqtree {
        publishDir = output
        
        conda "$workflow.projectDir/envs/iqtree.yaml"
        
        cpus threads 

        memory '2 GB'
        
        //errorStrategy 'ignore'
        
        input:
        file(aln) from fasta_ch.collect()
        
        output:
        file("*") into phylo_ch 
        
        script:
        """
        cat *.fasta > mbovpan_align.fasta
        cat $workflow.projectDir/ref/mbov_reference.fasta >> mbovpan_align.fasta
        snp-sites -o mbovpan_align.snp_only.fasta mbovpan_align.fasta
        iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${task.cpus} -bb 1000 -pre mbovpan_align -o "LT708304.1"
        """
    }
}

// vcf filtering + generate alignment? 
// pangenome steps (might need to separately create environment for pangenome software. inclusion in main environment might lead to conflicts)
// roary, quast OR busco

// PART 3: Pangenome 
assembly = Channel.create()
if(run_mode == "pan" || run_mode == "all"){
    
    process assembly {
    publishDir = output 
    
    conda "$workflow.projectDir/envs/megahit.yaml"
    
    errorStrategy "ignore"

    cpus threads

    input:
    tuple file(trim1), file(trim2) from fastp_reads3

    output:
    file("${trim1.baseName}.scaffold.fasta") into shortassembly_ch

    script:
    """
    megahit -t ${task.cpus} -1 ${trim1} -2 ${trim2} -o ${trim1.baseName}
    cd ${trim1.baseName}
    mv final.contigs.fa  ../${trim1.baseName}.scaffold.fasta
    """
}
assembly_ch = shortassembly_ch

assembly_ch.into {
    assembly_ch1
    assembly_ch2
}

process quast {
    publishDir = output

    conda "$workflow.projectDir/envs/quast.yaml"
    
    cpus threads

    input:
    file(assemblies) from assembly_ch1.collect()
    
    output:
    file("*") into quast_ch
    
    script:
    """
    quast -o assembly_stats ${assemblies} -r ${ref} -t ${task.cpus}
    """
}


process annotate {
    publishDir = output
    
    cpus threads

    //conda "$workflow.projectDir/envs/prokka.yaml"
    //conda "/scratch/noahaus/aim_1/prokka"
    conda = 'bioconda::prokka'

    input:
    file(assembly) from assembly_ch2

    output:
    file("${assembly.baseName}.annot.gff") into annotate_ch

    script:
    """
    prokka  --outdir ./${assembly.baseName} --cpus ${task.cpus} --prefix ${assembly.baseName}.annot ${assembly} 
    cp ./${assembly.baseName}/${assembly.baseName}.annot.gff ./
    """
}

process roary {
    publishDir = output

    conda "$workflow.projectDir/envs/roary.yaml"

    cpus threads

    input:
    file(gff) from annotate_ch.collect()

    output:
    file("*") into roary_ch

    script:
    """
    roary -e --mafft -p ${task.cpus} $gff
    """
}

process pan_curve {
    publishDir = output
    
    conda 'r-ggplot2'

    errorStrategy 'ignore'
    
    input:
    file(input) from roary_ch.collect()
    
    output:
    file("pangenome_curve.png") into output_ch
    
    script:
    """
    Rscript $workflow.projectDir/scripts/pangenome_curve.R
    """
}   

}

process multiqc {
    publishDir = output
    
    conda "$workflow.projectDir/envs/multiqc.yaml"

    input:
    file(pre) from fastqc_ch1.collect().ifEmpty([])
    file(post) from fastqc_ch2.collect().ifEmpty([])
    //if( run_mode != "snp"){ file(quast) from quast_ch.collect().ifEmpty([]) }

    output:
    file("mbovpan_report*")

    script:
    """
    multiqc -n mbovpan_report .
    """
}



