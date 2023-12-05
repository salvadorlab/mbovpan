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

//verify that it worked out 
println "${input}"
println "${meta}"

process accessory_pca {
  publishDir = "./"

  debug true
    
  conda "r conda-forge::r-ggplot2 conda-forge::r-dplyr"

    //errorStrategy 'ignore'
    
  input:
     path x
     path y
    
  output:
    path 'pca_figures.pdf'
    
    script:
      """
      Rscript $workflow.projectDir/../scripts/accessory_pca.R $x $y
      """
}

workflow {
  accessory_pca(Channel.fromPath( "${input}" ), Channel.fromPath( "${meta}" ))
}