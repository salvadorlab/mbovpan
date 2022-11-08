input = params.input 
meta = params.meta

process gene_prab {
     publishDir = output
    
    conda "$workflow.projectDir/envs/gene_prab.yaml"

    //errorStrategy 'ignore'
    
    input:
    path 'gene_presence_absence.csv'
    
    output:
    file("mbov_virulent_prab.csv") into geneprab_ch1
    file("gene_prab_figures.pdf") into geneprab_ch2
    
    script:
    """
    python $workflow.projectDir/scripts/mbov_virulence.py
    Rscript $workflow.projectDir/scripts/gene_prab.R ${meta}
    """
}

workflow {
  def all_files = Channel.fromPath( "${input}/*" )
  gene_prab(all_files)
}