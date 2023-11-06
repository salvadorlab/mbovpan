threads = Math.floor(Runtime.getRuntime().availableProcessors())
gffs = Channel.fromPath("*.gff")

process panaroo {
    publishDir = "./"

    conda "$workflow.projectDir/envs/panaroo.yaml"

    cpus threads

    input:
    file(gff) from gffs.collect()

    output:
    file("*") into panaroo
    script:
    """
    panaroo -i *.gff -o ./ -t ${task.cpus} -a core --core_threshold 0.98 --clean-mode strict
    """
}