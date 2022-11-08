nextflow.enable.dsl=2

input = params.input 
meta = params.meta

vir_genes = "$workflow.projectDir/../auxilary/M_bovis_virulence_genes.csv"

process gene_prab {
     publishDir = "./"
    
    conda "$workflow.projectDir/../envs/gene_prab.yaml"

    //errorStrategy 'ignore'
    
    input:
    path 'gene_presence_absence.csv'
    
    output:
    file("mbov_virulent_prab.csv") 
    file("gene_prab_figures.pdf") 
    
    script:
    """
    python $workflow.projectDir/../scripts/mbov_virulence.py ${vir_genes}
    Rscript $workflow.projectDir/../scripts/gene_prab.R ${meta}
    """
}

workflow {
  def all_files = Channel.fromPath( "${input}/*" )
  gene_prab(all_files)
}