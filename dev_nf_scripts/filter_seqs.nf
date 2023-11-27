// can we create a tuple in real time of the 
spotyping = "$workflow.projectDir/scripts/SpoTyping/SpoTyping.py"
check = "$workflow.projectDir/scripts/lineage_check.py"

process spotyping {

    publishDir = "./mbovpan_results/spotyping"

    conda "$workflow.projectDir/../envs/spotyping.yaml"

    input:
    tuple sample_id, file(reads_file) from Channel.fromFilePairs("$input*{1,2}*.f*q*")

    output:
    tuple stdout result, file("${reads_file[0].baseName - ~/_1*/}.log") into spoligo_ch
 
    script:
    """
    python3 ${spotyping} ${reads_file[0]} ${reads_file[1]} -o ${reads_file[0].baseName - ~/_1*/}.log > stdout.out 
    python3 ${check} 
    """

    }

spoligo_ch.view()
