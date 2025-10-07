nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 1
    executor = 'local'
}

// Parâmetros
params.input       = "${params.input ?: 'data'}"         // Pasta com arquivos .fastq.gz
params.outdir      = "${params.outdir ?: 'results'}"      // Pasta de saída
params.num_reads   = "${params.num_reads ?: 1000}"        // Número de reads para subsample
params.padroes     = ["16S", "ITS", "18S", "12S", "GENOMA"] // Padrões esperados

// Processo para subsample com seqtk
process SEQTK_SUBSAMPLE {
    tag "$arquivo"

    publishDir "${params.outdir}/MICROBIOMAS/reads_sample/${padrao}", mode: 'copy'

    input:
    tuple val(padrao), val(arquivo), path arquivo_path

    output:
    path "${arquivo}_subsample.fastq.gz"

    container:
    'gogeneticacr.azurecr.io/seqtk:1.4'

    script:
    """
    seqtk sample -s100 ${arquivo_path} ${params.num_reads} > ${arquivo}_subsample.fastq.gz
    """
}

// Workflow principal
workflow {
    // Cria canal com todos os arquivos .fastq.gz
    all_reads_ch = Channel.fromPath("${params.input}/*.fastq.gz")
                          .map { file -> tuple(file.baseName, file) }

    // Filtra arquivos por padrão e cria canal com tupla (padrao, nome, caminho)
    filtered_ch = all_reads_ch
      .filter { name, file ->
          def match = params.padroes.find { padrao -> name.startsWith(padrao) }
          return match != null
      }
      .map { name, file ->
          def padrao = params.padroes.find { p -> name.startsWith(p) }
          tuple(padrao, name, file)
      }

    // Executa subsample com seqtk
    SEQTK_SUBSAMPLE(filtered_ch)
}
