nextflow.enable.dsl=2

// Start off with a empty channel for prab input and sample metadata
pangenome_ref_input = ""
prab_input = ""

//pass in the paths for 'input' and 'metadata'
if(params.input != null){
    pangenome_ref_input = params.input
    }

if(params.prab != null){
    prab_input = params.prab
    }

params.help = false
if(params.help){
    println(
"""
    M B O V P A N (v1)    
=============================
hierarchical clustering of mbovpan output 

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
println "${pangenome_ref_input}"


process filter_pan {
    publishDir = "./"
    
    conda "bioconda::blast=2.9.0 pandas"

    //errorStrategy 'ignore'

    debug 'true'
    
    input:
    path x
    path y

    output:
    path 'mb.out.csv'
    path 'mbovis_filtered_cogs.csv'
    
    
    script:
    """
    # first determine the genes that are present and use blast to find their true annotation
    blastn -query $x -max_target_seqs 1 -db $workflow.projectDir/ref/mbovis_reference -out mb.out -outfmt "6 delim=, qseqid sseqid length qstart qend sstart send qlen slen"
    echo "qseqid,sseqid,length,qstart,qend,sstart,send,qlen,slen" > heading.txt 
    cat heading.txt mb.out >> mb.out.csv

    # once mb.out is created, we can use python to 1) create length %, 2) filter by 75% or higher, and 3) reduce repetitive blast results
    python $workflow.projectDir/scripts/blast_result_filter.py $y
    """
}

workflow {
  filter_pan(Channel.fromPath("${pangenome_ref_input}"),Channel.fromPath("${prab_input}"))
}