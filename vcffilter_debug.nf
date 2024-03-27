nextflow.enable.dsl=2

vcf = "ERR118935271.filtered.vcf"
ref = "mbov_reference.fasta"

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


    process psuedo_assembly {

        conda "./consensus.yaml"

        input:
        path x

        output:
        file("${vcf.baseName - ~/.filtered/}.consensus.fasta") into fasta_ch

        script:
        """
        bgzip $x
        bcftools index ${x}.gz
        cat ${ref} | bcftools norm --check-ref s --fasta-ref $ref -Ov ${x}.gz | vcf-consensus ${x.baseName}.vcf.gz > ${x.baseName - ~/.filtered/}.dummy.fasta
        sed 's|LT708304.1 Mycobacterium bovis AF2122/97 genome assembly, chromosome: Mycobacterium_bovis_AF212297|${x.baseName - ~/.filtered/}}|g' ${x.baseName - ~/.filtered/}.dummy.fasta > ${x.baseName - ~/.filtered/}.consensus.fasta
        """
    }

workflow {
  consensus(Channel.fromPath( "${vcf}" ))
}