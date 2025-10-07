nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 1
    executor = 'local'
}

// Parâmetros
params.input = "${params.input ?: 'data'}"         // Pasta raiz com subpastas de reads
params.outdir = "${params.outdir ?: 'results'}"     // Pasta de saída para estrutura QC

// Processo para identificar pastas com reads e criar estrutura de QC
process PREPARACAO_QC {
    tag "preparacao_qc"
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    path input_dir

    output:
    path "estrutura_qc.txt"

    container:
    'ubuntu:latest'

    script:
    """
    # Cria estrutura de pastas para QC com base nos padrões identificados
    mkdir -p ${params.outdir}/QC

    # Busca subpastas com padrão 'reads' ou 'NEXT_run'
    for sub in \$(find ${input_dir} -type d -regex '.*reads.*\\|.*NEXT_run.*'); do
        echo "Verificando pasta: \$sub"
        if ls \$sub/*.fastq.gz 1> /dev/null 2>&1; then
            echo "Pasta com reads encontrada: \$sub" >> estrutura_qc.txt

            # Identifica padrões nos nomes dos arquivos
            for file in \$(ls \$sub/*.fastq.gz); do
                base=\$(basename \$file)
                if [[ \$base == 16S* ]]; then
                    mkdir -p ${params.outdir}/QC/MICROBIOMAS/reads_sample/16S
                    mkdir -p ${params.outdir}/QC/MICROBIOMAS/resultados_Qiime2/16S
                elif [[ \$base == ITS* ]]; then
                    mkdir -p ${params.outdir}/QC/MICROBIOMAS/reads_sample/ITS
                    mkdir -p ${params.outdir}/QC/MICROBIOMAS/resultados_Qiime2/ITS
                elif [[ \$base == GENOMA* ]]; then
                    mkdir -p ${params.outdir}/QC/GENOMAS/reads
                    mkdir -p ${params.outdir}/QC/GENOMAS/resultados_Kraken2
                fi
            done
        fi
    done

    # Pastas comuns para todos os serviços
    mkdir -p ${params.outdir}/QC/seqkit
    mkdir -p ${params.outdir}/QC/fastqc
    mkdir -p ${params.outdir}/QC/multiqc

    echo "Estrutura de QC criada com sucesso." >> estrutura_qc.txt
    """
}

// Workflow principal
workflow {
    input_ch = Channel.fromPath("${params.input}")

    PREPARACAO_QC(input_ch)
}
