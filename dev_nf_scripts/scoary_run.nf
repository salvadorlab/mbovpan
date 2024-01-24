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


process scoary {
    publishDir = "./"
    
    conda "$workflow.projectDir/../envs/scoary.yaml"


    debug 'true'

    input:
     path x
     path y

    output:
     path "*.csv"
    
    script:
    """
    sed 's/.annot//g' $x > prab.csv
    scoary -t $y -g prab.csv
    """
}

workflow {
  scoary(Channel.fromPath( "${input}" ), Channel.fromPath( "${meta}" ))
}
