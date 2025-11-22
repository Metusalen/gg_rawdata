nextflow.enable.dsl=2

// ------------------------------------------------------
// Parâmetros globais
// ------------------------------------------------------
params.input      = "${params.input ?: 'data'}"
params.outdir     = "${params.outdir ?: 'results'}"
params.servico    = "${params.servico ?: 'rawdata'}"
params.database   = "${params.database ?: 'pluspfp-8'}"
params.db_url     = "https://genome-idx.s3.amazonaws.com/kraken"
params.threads    = "${params.threads ?: 4}"
params.num_reads  = "${params.num_reads ?: 1000}"
params.padroes    = ["16S", "ITS", "18S", "12S", "GENOMA"]
params.qual       = "${params.qual ?: 25}"
params.minlen     = "${params.minlen ?: 50}"
params.window     = "${params.window ?: 4}"
params.leading    = "${params.leading ?: 3}"
params.trailing   = "${params.trailing ?: 3}"

// ------------------------------------------------------
// Process configurations
// ------------------------------------------------------
process {
    cpus = 4
    executor = 'local'
}

// ------------------------------------------------------
// Processos
// ------------------------------------------------------

// Estrutura QC
process PREPARACAO_QC {
    tag "preparacao_qc"
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    path input_dir

    output:
    path "estrutura_qc.txt"

    container 'ubuntu:latest'

    script:
    """
    mkdir -p ${params.outdir}/QC
    for sub in \$(find ${input_dir} -type d -regex '.*reads.*\\|.*NEXT_run.*'); do
        echo "Verificando pasta: \$sub"
        if ls \$sub/*.fastq.gz 1> /dev/null 2>&1; then
            echo "Pasta com reads encontrada: \$sub" >> estrutura_qc.txt
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
    mkdir -p ${params.outdir}/QC/seqkit
    mkdir -p ${params.outdir}/QC/fastqc
    mkdir -p ${params.outdir}/QC/multiqc
    echo "Estrutura de QC criada com sucesso." >> estrutura_qc.txt
    """
}

// Lista arquivos R1/R2
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
    find ${input_dir} -type f -name '*.fastq.gz' | grep '_R1' > R1.txt || true
    find ${input_dir} -type f -name '*.fastq.gz' | grep '_R2' > R2.txt || true
    """
}

// FastQC
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

// MultiQC
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

// Seqkit
process SEQKIT_STATS {
    tag "seqkit"
    publishDir "${params.outdir}/seqkit", mode: 'copy'

    input:
    path reads_dir

    output:
    path "seqkit_stats.tsv"

    container 'gogeneticacr.azurecr.io/seqkit:v2.6.1'

    script:
    """
    seqkit stats -a -T -j ${params.threads} ${reads_dir}/*.fastq.gz > seqkit_stats.tsv
    """
}

process PARSE_SEQKIT {
    tag "parse_seqkit"
    publishDir "${params.outdir}/seqkit", mode: 'copy'

    input:
    path stats_file

    output:
    path "seqkit_stats.csv"

    container 'ubuntu:latest'

    script:
    """
    awk -F '\\t' '
    NR==1 { col_count = NF; if (col_count == 16) { print "Arquivo,Numero_de_reads,Q30_R1,Q30_R2,Bases_Sequenciadas_(Mb)" } }
    NR>1 {
        file = \$1
        gsub(/.*\\//, "", file)
        gsub(/_S[0-9]+_L001_R[12]_001.fastq.gz/, "", file)
        gsub(/_R[12]_001.fastq.gz/, "", file)
        if (\$1 ~ /_R1/) { r1[file] = \$4; q30_r1[file] = \$15; bases_r1[file] = \$5 }
        else if (\$1 ~ /_R2/) { r2[file] = \$4; q30_r2[file] = \$15; bases_r2[file] = \$5 }
    }
    END { for (f in r1) { total_reads = r1[f] + r2[f]; total_bases = (bases_r1[f] + bases_r2[f]) / 1000000; printf "%s,%d,%s,%s,%.2f\\n", f, total_reads, q30_r1[f], q30_r2[f], total_bases } }
    ' ${stats_file} > seqkit_stats.csv
    """
}

// Subsample seqtk
process SEQTK_SUBSAMPLE {
    tag "$arquivo"
    publishDir "${params.outdir}/MICROBIOMAS/reads_sample/${padrao}", mode: 'copy'

    input:
    tuple val(padrao), val(arquivo), path arquivo_path

    output:
    path "${arquivo}_subsample.fastq.gz"

    container 'gogeneticacr.azurecr.io/seqtk:1.4'

    script:
    """
    seqtk sample -s100 ${arquivo_path} ${params.num_reads} > ${arquivo}_subsample.fastq.gz
    """
}

// Download Kraken2 DB
process DOWNLOAD_DB {
    tag "$params.database"
    publishDir "${params.outdir}/database_kraken2", mode: 'copy'

    output:
    path "${params.database}", emit: db_path

    container 'ubuntu:latest'

    script:
    """
    mkdir -p ${params.outdir}/database_kraken2
    wget ${params.db_url}/k2_${params.database}_20231009.tar.gz -O ${params.outdir}/database_kraken2/${params.database}.tar.gz
    tar -xvf ${params.outdir}/database_kraken2/${params.database}.tar.gz -C ${params.outdir}/database_kraken2/
    echo "${params.outdir}/database_kraken2/${params.database}" > ${params.database}
    """
}

// Kraken2
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path db_path

    output:
    path "out-kraken2-${sample_id}.txt"
    path "report-kraken2-${sample_id}.txt"

    container 'gogeneticacr.azurecr.io/kraken2:2.1.3'

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

// Estrutura Trimmomatic
process CRIA_ESTRUTURA_TRIMMOMATIC {
    tag "estrutura"
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path "trimmomatic_reads"
    path "trimmomatic_temp"

    container 'ubuntu:latest'

    script:
    """
    mkdir -p trimmomatic_reads
    mkdir -p trimmomatic_temp
    """
}

// Trimmomatic
process TRIMMOMATIC_PAIRED {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmomatic_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path estrutura_dir

    output:
    path "${sample_id}_R1_trimmed.fastq.gz"
    path "${sample_id}_R2_trimmed.fastq.gz"

    container 'gogeneticacr.azurecr.io/trimmomatic:0.39'

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

// ------------------------------------------------------
// Workflow principal
// ------------------------------------------------------
workflow {

    // 1. Estrutura de QC
    prep_qc_ch = Channel.fromPath("${params.input}")
    PREPARACAO_QC(prep_qc_ch)

    // 2. Lista arquivos R1/R2
    list_ch = Channel.value(file(params.input))
    LIST_FASTQ_FILES(list_ch)

    // 3. FastQC e MultiQC
    reads_ch_fastqc = Channel.fromPath("${params.input}/*.fastq.gz")
                         .map { file -> tuple(file.baseName, file) }
    fastqc_out = FASTQC(reads_ch_fastqc)
    MULTIQC(fastqc_out)

    // 4. Seqkit e Parse Seqkit
    reads_ch_seqkit = Channel.fromPath("${params.input}")
    stats_out = SEQKIT_STATS(reads_ch_seqkit)
    PARSE_SEQKIT(stats_out.out)

    // 5. Subsample seqtk
    all_reads_ch = Channel.fromPath("${params.input}/*.fastq.gz")
                         .map { file -> tuple(file.baseName, file) }
    filtered_ch = all_reads_ch
      .filter { name, file -> params.padroes.find { padrao -> name.startsWith(padrao) } != null }
      .map { name, file ->
          def padrao = params.padroes.find { p -> name.startsWith(p) }
          tuple(padrao, name, file)
      }
    SEQTK_SUBSAMPLE(filtered_ch)

    // 6. Download do banco de dados Kraken2
    db_ch = DOWNLOAD_DB()

    // 7. Kraken2
    reads_ch_kraken = Channel
      .fromFilePairs("${params.input}/*_R{1,2}_001.fastq.gz", flat: true)
      .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    KRAKEN2(reads_ch_kraken, db_ch.db_path)

    // 8. Trimmomatic: cria estrutura e executa
    estrutura_ch = CRIA_ESTRUTURA_TRIMMOMATIC()
    TRIMMOMATIC_PAIRED(reads_ch_kraken, estrutura_ch.out)
}
