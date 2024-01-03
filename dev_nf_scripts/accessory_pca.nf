nextflow.enable.dsl=2
conda.enabled = true

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
    
  conda "r-ggplot2 r-dplyr"

    //errorStrategy 'ignore'
    
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