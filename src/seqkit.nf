nextflow.enable.dsl=2

// Configurações básicas
process {
    cpus = 4
    executor = 'local'
}

// Parâmetros
params.input     = "${params.input ?: 'data'}"         // Pasta com arquivos .fastq.gz
params.outdir    = "${params.outdir ?: 'results'}"      // Pasta de saída
params.threads   = "${params.threads ?: 4}"             // Threads para seqkit

// Processo para rodar seqkit stats
process SEQKIT_STATS {
    tag "seqkit"
    publishDir "${params.outdir}/seqkit", mode: 'copy'

    input:
    path reads_dir

    output:
    path "seqkit_stats.tsv"

    container:
    'gogeneticacr.azurecr.io/seqkit:v2.6.1'

    script:
    """
    seqkit stats -a -T -j ${params.threads} ${reads_dir}/*.fastq.gz > seqkit_stats.tsv
    """
}

// Processo para parsear o .tsv e gerar .csv com estatísticas combinadas
process PARSE_SEQKIT {
    tag "parse_seqkit"
    publishDir "${params.outdir}/seqkit", mode: 'copy'

    input:
    path stats_file

    output:
    path "seqkit_stats.csv"

    container:
    'ubuntu:latest'

    script:
    """
    # Renomeia colunas e arquivos, separa R1/R2, calcula estatísticas combinadas
    awk -F '\\t' '
    NR==1 {
        col_count = NF
        if (col_count == 16) {
            print "Arquivo,Numero_de_reads,Q30_R1,Q30_R2,Bases_Sequenciadas_(Mb)"
        }
    }
    NR>1 {
        file = \$1
        gsub(/.*\\//, "", file)
        gsub(/_S[0-9]+_L001_R[12]_001.fastq.gz/, "", file)
        gsub(/_R[12]_001.fastq.gz/, "", file)
        if (\$1 ~ /_R1/) {
            r1[file] = \$4
            q30_r1[file] = \$15
            bases_r1[file] = \$5
        } else if (\$1 ~ /_R2/) {
            r2[file] = \$4
            q30_r2[file] = \$15
            bases_r2[file] = \$5
        }
    }
    END {
        for (f in r1) {
            total_reads = r1[f] + r2[f]
            total_bases = (bases_r1[f] + bases_r2[f]) / 1000000
            printf "%s,%d,%s,%s,%.2f\\n", f, total_reads, q30_r1[f], q30_r2[f], total_bases
        }
    }
    ' ${stats_file} > seqkit_stats.csv
    """
}

// Workflow principal
workflow {
    reads_ch = Channel.fromPath("${params.input}")

    stats_out = SEQKIT_STATS(reads_ch)
    PARSE_SEQKIT(stats_out.out)
}
