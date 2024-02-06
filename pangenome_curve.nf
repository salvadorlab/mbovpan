nextflow.enable.dsl=2

// Start off with a empty channel
input = ""

if(params.input != null){
    input = params.input
    }

params.help = false
if(params.help){
    println(
"""
    M B O V P A N (v1)    
=============================
generate a pangenome curve to visualize gene conservation and expansion as a function of genomes

usage: nextflow run mbovpan/dev_nf_scripts [options] --input ./path/to/input --meta ./path/to/metadata
  options:
    --input 
        The gene_presence_absence.csv that comes as an output of mbovpan - should be located in the mbovpan_results/pangenome directory
    --help
        Prints this help message
=============================
"""

    )
    exit(0)
}


//verify that it worked out 
println "${input}"

process pangenome_curve {
  publishDir = "./"
    
  conda "r-ggplot2 r-dplyr"

  debug "true"
    
  input:
     path x

output:
    path 'pangenome_curve.pdf'
    
    script:
      """
      sed 's/.annot//g' $x > prab.csv
      Rscript $workflow.projectDir/../scripts/pangenome_curve.R prab.csv
      """
}

workflow {
  pangenome_curve(Channel.fromPath( "${input}" ))
}