nextflow.enable.dsl=2

// Start off with a empty channel
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

process accessory_pca {
  publishDir = "./"
    
  conda "r-ggplot2 r-dplyr"

  debug "true"
    
  input:
     path x
     path y
    
  output:
    path 'pca_figures.pdf'
    
    script:
      """
      sed 's/.annot//g' $x > prab.csv
      Rscript $workflow.projectDir/../scripts/accessory_pca.R prab.csv $y
      """
}

workflow {
  accessory_pca(Channel.fromPath( "${input}" ), Channel.fromPath( "${meta}" ))
}