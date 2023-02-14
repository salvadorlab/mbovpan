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
    file("${trim1.baseName - ~/_trimmed_R*/}.scaffold.fasta") into shortassembly_ch

    script:
    """
    megahit -t ${task.cpus} -1 ${trim1} -2 ${trim2} -o ${trim1.baseName}
    cd ${trim1.baseName}
    mv final.contigs.fa  ../${trim1.baseName - ~/_trimmed_R*/}.scaffold.fasta
    """
}