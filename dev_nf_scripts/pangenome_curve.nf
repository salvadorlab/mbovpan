nextflow.enable.dsl=2

// Start off with a empty channel
input = ""

if(params.input != null){
    input = params.input
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