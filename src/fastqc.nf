nextflow.enable.dsl=2

// Configurações básicas do processo
process {
    cpus = 2              // Número de threads por processo
    executor = 'local'    // Executor local (pode ser SLURM, AWS, etc.)
}

// Parâmetros com valores padrão
params.input   = "${params.input ?: 'data'}"       // Pasta com arquivos .fastq.gz
params.outdir  = "${params.outdir ?: 'results'}"   // Pasta de saída dos resultados
params.servico = "${params.servico ?: 'rawdata'}"  // Nome usado no relatório do MultiQC

// Processo para rodar FastQC em cada arquivo de leitura
process FASTQC {
    tag "$sample_id"  // Identificador da amostra para facilitar o rastreio
    publishDir "${params.outdir}/fastqc", mode: 'copy'  // Salva os resultados na pasta de saída

    input:
    tuple val(sample_id), path(reads)  // Recebe o nome da amostra e o arquivo .fastq.gz

    output:
    path "*.zip"   // Arquivo compactado com os resultados
    path "*.html"  // Relatório HTML do FastQC

    container:
    'gogeneticacr.azurecr.io/fastqc:latest'  // Imagem Docker personalizada da empresa

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}  // Executa FastQC com número de threads definido
    """
}

// Processo para rodar MultiQC nos resultados do FastQC
process MULTIQC {
    tag "multiqc"  // Tag fixa para identificar o processo
    publishDir "${params.outdir}/multiqc", mode: 'copy'  // Salva o relatório consolidado

    input:
    path fastqc_results  // Recebe os arquivos gerados pelo FastQC

    output:
    path "*_report.html"  // Relatório HTML interativo do MultiQC

    container:
    'gogeneticacr.azurecr.io/multiqc:1.14'  // Imagem Docker personalizada da empresa

    script:
    """
    multiqc -p -o . ${fastqc_results} --filename ${params.servico}_report.html --interactive
    """
}

// Workflow principal que conecta os processos
workflow {
    // Cria um canal com os arquivos .fastq.gz e seus nomes base
    reads_ch = Channel.fromPath("${params.input}/*.fastq.gz")
                      .map { file -> tuple(file.baseName, file) }

    // Executa FastQC em cada arquivo
    fastqc_out = FASTQC(reads_ch)

    // Executa MultiQC nos resultados do FastQC
    multiqc_out = MULTIQC(fastqc_out.out)
}
