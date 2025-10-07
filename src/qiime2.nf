nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 1
    executor = 'local'
}

// Parâmetros
params.input    = "${params.input ?: 'data'}"         // Pasta com arquivos .fastq.gz
params.outdir   = "${params.outdir ?: 'results'}"      // Pasta de saída
params.marcador = "${params.marcador ?: '16S'}"        // Marcador alvo (ex: 16S, ITS, etc.)

// Processo para filtrar arquivos com base no marcador
process FILTRA_READS_QIIME2 {
    tag "$params.marcador"
    publishDir "${params.outdir}/qiime2", mode: 'copy'

    input:
    path input_dir

    output:
    path "arquivos_filtrados.txt"

    script:
    """
    # Busca arquivos .fastq.gz que contenham o marcador e padrão de leitura (R1 ou R2)
    find ${input_dir} -type f -name '*.fastq.gz' | \\
      grep '${params.marcador}-' | \\
      grep -E '_R[12]_.*\\.fastq\\.gz' > arquivos_filtrados.txt
    """
}

// Workflow principal
workflow {
    input_ch = Channel.fromPath("${params.input}")

    FILTRA_READS_QIIME2(input_ch)
}
