conda.enabled = true

conda {
    cacheDir = './mbovpan_conda_envs' // path where conda environments are stored
                    // probably smart to put it somewhere with ample storage
    useMamba = true
}

manifest {
    description = 'Pipeline does this and that'
    mainScript = 'mbovpan.nf'
    version = '1.0.0'
    name = 'mbovpan'
}