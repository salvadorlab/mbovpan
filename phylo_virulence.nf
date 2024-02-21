nextflow.enable.dsl=2

// Start off with a empty channel for prab input and sample metadata
input = ""
meta = ""

blast_db = "$workflow.projectDir/ref/mb.out"

//pass in the paths for 'input' and 'metadata'
if(params.input != null){
    input = params.input
    }

if(params.meta != null){
    meta = params.meta
    }

params.help = false
if(params.help){
    println(
"""
    M B O V P A N (v1)    
=============================
plotting virulent m. bovis genes against the core genome phylogeny 

usage: nextflow run mbovpan/dev_nf_scripts [options] --input ./path/to/input --meta ./path/to/metadata
  options:
    --input 
        The gene_presence_absence.csv that comes as an output of mbovpan - should be located in the mbovpan_results/pangenome directory
    --meta
        A CSV file linking sequence names to external metadata. IDs should be labelled in the 'Name' column of the metadata. 
    --help
        Prints this help message
=============================
"""

    )
    exit(0)
}

//verify that it worked out 
println "${input}"
println "${meta}"

process phylo_vir {
    publishDir = "./"
    
    conda "$workflow.projectDir/envs/gene_prab.yaml"

    //errorStrategy 'ignore'

    debug 'true'
    
    input:
    path x //core phylo tree
    path y 

    output:
    path 'core_phylo.pdf'
    
    
    script:
    """
    #fix the names of the output
    sed 's/.annot//g' $x > mbovis_core.snp_only.nwk

    #use blast to find the best hit per query
    blastn -query output.fasta -max_target_seqs 1 -db mbovis_refernce -out mb.out -outfmt 6
    
    
    """
}

workflow {
  gene_prab(Channel.fromPath("${input}"), Channel.fromPath("${meta}"))
}