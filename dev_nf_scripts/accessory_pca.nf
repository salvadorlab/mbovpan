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
    
  conda "$workflow.projectDir/../envs/accessory_pca.yaml"

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