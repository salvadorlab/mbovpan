nextflow.enable.dsl=2

// Start off with a empty channel for prab input and sample metadata
input = ""
meta = ""

//pass in the paths for 'input' and 'metadata'
if(params.input != null){
    input = params.input
    }

if(params.meta != null){
    meta = params.meta
    }

params.help = false
if(params.help){
    println(
"""
    M B O V P A N (v1)    
=============================
hierarchical clustering of mbovpan output 

usage: nextflow run mbovpan/dev_nf_scripts [options] --input ./path/to/input --meta ./path/to/metadata
  options:
    --input 
        The gene_presence_absence.csv that comes as an output of mbovpan - should be located in the mbovpan_results/pangenome directory
    --meta
        A CSV file linking sequence names to external metadata. IDs should be labelled in the 'Name' column of the metadata. 
    --help
        Prints this help message
=============================
"""

    )
    exit(0)
}

//verify that it worked out 
println "${input}"
println "${meta}"


process gene_prab {
    publishDir = "./"
    
    conda "$workflow.projectDir/envs/gene_prab.yaml"

    //errorStrategy 'ignore'

    debug 'true'
    
    input:
    path x
    path y

    output:
    path 'gene_prab_figures.pdf'
    
    
    script:
    """
    sed 's/.annot//g' $x > prab.csv
    Rscript $workflow.projectDir/scripts/gene_prab.R prab.csv $y
    """
}

workflow {
  gene_prab(Channel.fromPath("${input}"), Channel.fromPath("${meta}"))
}