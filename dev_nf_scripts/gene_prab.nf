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

//file with the virulence genes 
vir_genes = "$workflow.projectDir/../auxilary/M_bovis_virulence_genes.csv"

process gene_prab {
    publishDir = "./"

    debug true
    
    conda "$workflow.projectDir/../envs/gene_prab.yaml"

    //errorStrategy 'ignore'
    
    input:
    path x
    path y

    output:
    path 'gene_prab_figures.pdf'
    
    
    script:
    """
    python $workflow.projectDir/../scripts/mbov_virulence.py $x ${vir_genes}
    Rscript $workflow.projectDir/../scripts/gene_prab.R $y
    """
}

workflow {
  gene_prab(Channel.fromPath("${input}"), Channel.fromPath("${meta}"))
}