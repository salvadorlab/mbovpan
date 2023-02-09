nextflow.enable.dsl=2

input = "$workflow.projectDir/../auxilary/gene_presence_absence.csv" 
meta = "$workflow.projectDir/../auxilary/UK_meta.csv"


process accessory_pca {
    publishDir = "./"
    
    conda 'r conda-forge::r-ggplot2 conda-forge::r-dplyr'

    //errorStrategy 'ignore'
    
    input:
    path 'gene_presence_absence.csv'
    
    output:
    file("pca_figures.pdf")
    
    script:
    """
    Rscript $workflow.projectDir/../scripts/accessory_pca.R ${meta}
    """
}

workflow {
  def all_files = Channel.fromPath( "${input}" )
  accessory_pca(all_files)
}