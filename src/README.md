# Guia de Execução dos Pipelines Nextflow

Este repositório contém pipelines modulares para análise de dados de sequenciamento, todos escritos em Nextflow e integrados com ferramentas como FastQC, Kraken2, Trimmomatic, Qiime2, SeqKit e Seqtk. Abaixo estão os comandos para executar cada etapa, com explicações linha a linha.

## FastQC + MultiQC

```bash
nextflow run main.nf \                    # Executa o pipeline definido no arquivo main.nf
  --input data \                          # Pasta onde estão os arquivos .fastq.gz
  --outdir results \                      # Pasta onde os relatórios serão salvos
  --servico rawdata                       # Nome usado no relatório HTML do MultiQC (ex: rawdata_report.html)
  
## Kraken2

```bash
nextflow run main.nf \                    # Executa o pipeline Kraken2
  --input data \                          # Pasta com arquivos .fastq.gz em pares (R1/R2)
  --outdir results \                      # Pasta onde os relatórios serão salvos
  --database pluspfp-8 \                  # Nome do banco Kraken2 a ser baixado e usado
  --threads 4                             # Número de threads para acelerar a classificação taxonômica

  
