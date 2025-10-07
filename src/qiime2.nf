nextflow.enable.dsl=2

// Parâmetros
params.input    = "${params.input ?: 'data'}"      // Pasta com arquivos .fastq.gz
params.outdir   = "${params.outdir ?: 'results'}" // Pasta de saída

// Processo para listar arquivos R1 e R2
process LIST_FASTQ_FILES {
    cpus 1
    executor 'local'
    tag "Listando FASTQ"
    publishDir "${params.outdir}/qiime2", mode: 'copy'

    input:
    path input_dir

    output:
    path "R1.txt"
    path "R2.txt"

    script:
    """
    # Lista arquivos R1 e R2 separadamente
    find ${input_dir} -type f -name '*.fastq.gz' | grep '_R1' > R1.txt || true
    find ${input_dir} -type f -name '*.fastq.gz' | grep '_R2' > R2.txt || true
    """
}

// Workflow principal
workflow {
    // Passa o diretório inteiro como um único valor
    input_ch = Channel.value(file(params.input))

    LIST_FASTQ_FILES(input_ch)
}









