nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 4
    executor = 'local'
}

// Parâmetros
params.input     = "${params.input ?: 'data'}"         // Pasta com arquivos .fastq.gz
params.outdir    = "${params.outdir ?: 'results'}"      // Pasta de saída
params.threads   = "${params.threads ?: 4}"             // Threads para Trimmomatic
params.qual      = "${params.qual ?: 25}"               // Qualidade mínima
params.minlen    = "${params.minlen ?: 50}"             // Tamanho mínimo
params.window    = "${params.window ?: 4}"              // Tamanho da janela
params.leading   = "${params.leading ?: 3}"             // Corte inicial
params.trailing  = "${params.trailing ?: 3}"            // Corte final

// Processo para criar estrutura de saída
process CRIA_ESTRUTURA_TRIMMOMATIC {
    tag "estrutura"
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path "trimmomatic_reads"
    path "trimmomatic_temp"

    container:
    'ubuntu:latest'

    script:
    """
    mkdir -p trimmomatic_reads
    mkdir -p trimmomatic_temp
    """
}

// Processo para rodar Trimmomatic em pares de reads
process TRIMMOMATIC_PAIRED {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmomatic_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path estrutura_dir

    output:
    path "${sample_id}_R1_trimmed.fastq.gz"
    path "${sample_id}_R2_trimmed.fastq.gz"

    container:
    'gogeneticacr.azurecr.io/trimmomatic:0.39'

    script:
    """
    java -jar /data/trimmomatic.jar PE \\
      -threads ${params.threads} -phred33 \\
      ${read1} ${read2} \\
      ${sample_id}_R1_trimmed.fastq.gz ${estrutura_dir}/trimmomatic_temp/${sample_id}_R1_unpaired.fq.gz \\
      ${sample_id}_R2_trimmed.fastq.gz ${estrutura_dir}/trimmomatic_temp/${sample_id}_R2_unpaired.fq.gz \\
      ILLUMINACLIP:/data/adapters/NexteraPE-PE.fa:2:27:10:25:True \\
      LEADING:${params.leading} TRAILING:${params.trailing} \\
      SLIDINGWINDOW:${params.window}:${params.qual} MINLEN:${params.minlen}
    """
}

// Workflow principal
workflow {
    // Cria canal com pares de arquivos R1/R2
    reads_ch = Channel
      .fromFilePairs("${params.input}/*_R{1,2}_001.fastq.gz", flat: true)
      .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    // Cria estrutura de saída
    estrutura_ch = CRIA_ESTRUTURA_TRIMMOMATIC()

    // Executa Trimmomatic nos pares
    TRIMMOMATIC_PAIRED(reads_ch, estrutura_ch.out)
}
