nextflow.enable.dsl=2

// Start off with a empty channel for prab input and sample metadata
pangenome_ref_input = ""


//pass in the paths for 'input' and 'metadata'
if(params.input != null){
    pangenome_ref_input = params.input
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
println "${prab_input}"


process filter_pan {
    publishDir = "./"
    
    conda "bioconda::blast=2.9.0"

    //errorStrategy 'ignore'

    debug 'true'
    
    input:
    path x

    output:
    path 'mb.out'
    
    
    script:
    """
    blastn -query $x -max_target_seqs 1 -db $workflow.projectDir/ref/mbovis_reference -out mb.out -outfmt "6 qseqid sseqid length qstart qend sstart send qlen slen"
    """
}

workflow {
  filter_pan(Channel.fromPath("${pangenome_ref_input}"))
}