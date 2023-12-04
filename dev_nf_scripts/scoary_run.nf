input_ch = Channel.fromPath("./*_roary.csv")
scoary_meta = "/work/n/noahlegall/mbovpan_testing/mbovpan/ref/scoary_examp.csv"

process scoary {
    publishDir = "/work/n/noahlegall/mbovpan_testing/mbovpan_results/pan_gwas"
    
    conda "$workflow.projectDir/../envs/scoary.yaml"

    errorStrategy 'ignore'
    
    input:
    file(input) from input_channel.collect()
    
    output:
    file("*.csv") into scoary_ch
    
    script:
    """
    sed 's/.annot//g' gene_presence_absence_roary.csv > prab.csv
    scoary -t ${scoary_meta} -g prab.csv
    """
}


process test {

    input:
    file(reads_file) from scoary_ch

    script:
    """
    cat ${reads_file}
    """
}