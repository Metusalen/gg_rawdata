nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 4              // Número de threads por processo
    executor = 'local'    // Executor local (pode ser SLURM, AWS, etc.)
}

// Parâmetros com valores padrão
params.input     = "${params.input ?: 'data'}"         // Pasta com arquivos .fastq.gz
params.outdir    = "${params.outdir ?: 'results'}"      // Pasta de saída
params.database  = "${params.database ?: 'pluspfp-8'}"  // Nome do banco Kraken2
params.db_url    = "https://genome-idx.s3.amazonaws.com/kraken"  // Base de URLs
params.threads   = "${params.threads ?: 4}"             // Threads para Kraken2

// Processo para baixar e preparar o banco de dados do Kraken2
process DOWNLOAD_DB {
    tag "$params.database"
    publishDir "${params.outdir}/database_kraken2", mode: 'copy'

    output:
    path "${params.database}", emit: db_path

    container:
    'ubuntu:latest'

    script:
    """
    mkdir -p ${params.outdir}/database_kraken2
    wget ${params.db_url}/k2_${params.database}_20231009.tar.gz -O ${params.outdir}/database_kraken2/${params.database}.tar.gz
    tar -xvf ${params.outdir}/database_kraken2/${params.database}.tar.gz -C ${params.outdir}/database_kraken2/
    echo "${params.outdir}/database_kraken2/${params.database}" > ${params.database}
    """
}

// Processo para rodar Kraken2 em pares de reads
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path db_path

    output:
    path "out-kraken2-${sample_id}.txt"
    path "report-kraken2-${sample_id}.txt"

    container:
    'gogeneticacr.azurecr.io/kraken2:2.1.3'

    script:
    """
    kraken2 \\
      --db ${db_path} \\
      --paired ${read1} ${read2} \\
      --output out-kraken2-${sample_id}.txt \\
      --report report-kraken2-${sample_id}.txt \\
      --threads ${params.threads}
    """
}

// Workflow principal
workflow {
    // Cria canal com pares de arquivos R1/R2
    reads_ch = Channel
      .fromFilePairs("${params.input}/*_R{1,2}_001.fastq.gz", flat: true)
      .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    // Baixa e prepara o banco de dados
    db_ch = DOWNLOAD_DB()

    // Executa Kraken2 em cada par de reads
    KRAKEN2(reads_ch, db_ch.db_path)
}
