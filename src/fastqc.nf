nextflow.enable.dsl=2

// -----------------------------
// Parâmetros com valores padrão
// -----------------------------
params.input   = "${params.input ?: 'data'}"
params.outdir  = "${params.outdir ?: 'results'}"
params.servico = "${params.servico ?: 'rawdata'}"
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")

    container 'gogeneticacr.azurecr.io/fastqc:latest'

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    """
}
process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus 2

    input:
    tuple val(sample_id), path(fastqc_zip), path(fastqc_html)

    output:
    path("${params.servico}_report.html")

    container 'gogeneticacr.azurecr.io/multiqc:1.14'

    script:
    """
    multiqc -o . . --filename ${params.servico}_report.html --interactive
    """
}
workflow {

    // Cria canal de arquivos FASTQ
    reads_ch = Channel
        .fromPath("${params.input}/*.fastq.gz")
        .map { file -> tuple(file.baseName, file) }

    // Roda FastQC
    fastqc_out = FASTQC(reads_ch)

    // Roda MultiQC
    MULTIQC(fastqc_out)
}

