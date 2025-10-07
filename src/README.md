# Guia de Execução dos Pipelines Nextflow

Este repositório contém pipelines modulares para análise de dados de sequenciamento, todos escritos em Nextflow e integrados com ferramentas como FastQC, Kraken2, Trimmomatic, Qiime2, SeqKit e Seqtk. Abaixo estão os comandos para executar cada etapa, com explicações linha a linha.

## FastQC + MultiQC

```bash
nextflow run main.nf \                    # Executa o pipeline definido no arquivo main.nf
  --input data \                          # Pasta onde estão os arquivos .fastq.gz
  --outdir results \                      # Pasta onde os relatórios serão salvos
  --servico rawdata                       # Nome usado no relatório HTML do MultiQC (ex: rawdata_report.html)

Kraken2

nextflow run main.nf \                    # Executa o pipeline Kraken2
  --input data \                          # Pasta com arquivos .fastq.gz em pares (R1/R2)
  --outdir results \                      # Pasta onde os relatórios serão salvos
  --database pluspfp-8 \                  # Nome do banco Kraken2 a ser baixado e usado
  --threads 4                             # Número de threads para acelerar a classificação taxonômica

  
## Preparação de Pastas (estrutura de QC)

nextflow run main.nf \                    # Executa a etapa de preparação de estrutura de QC
  --input data \                          # Pasta raiz com subpastas contendo reads
  --outdir results                        # Pasta onde a estrutura QC será criada

## Qiime2 (filtragem por marcador)

nextflow run main.nf \                    # Executa a filtragem de arquivos por marcador
  --input data \                          # Pasta com arquivos .fastq.gz
  --outdir results \                      # Pasta onde os arquivos filtrados serão salvos
  --marcador 16S                          # Marcador alvo para filtragem (ex: 16S, ITS, GENOMA)

## SeqKit (estatísticas de qualidade)

nextflow run main.nf \                    # Executa o pipeline SeqKit
  --input data \                          # Pasta com arquivos .fastq.gz
  --outdir results \                      # Pasta onde os relatórios serão salvos
  --threads 4                             # Número de threads para acelerar o processamento

## Seqtk (subsample de reads)

nextflow run main.nf \                    # Executa o subsample com Seqtk
  --input data \                          # Pasta com arquivos .fastq.gz
  --outdir results \                      # Pasta onde os arquivos reduzidos serão salvos
  --num_reads 1000                        # Número de reads a serem extraídos de cada arquivo

## Trimmomatic (trimming de reads pareados)

nextflow run main.nf \                    # Executa o trimming com Trimmomatic
  --input data \                          # Pasta com arquivos .fastq.gz em pares (R1/R2)
  --outdir results \                      # Pasta onde os arquivos tratados serão salvos
  --threads 4 \                           # Número de threads para acelerar o trimming
  --qual 25 \                             # Qualidade mínima exigida na janela de corte
  --minlen 50 \                           # Tamanho mínimo de reads após o trimming
  --window 4 \                            # Tamanho da janela usada para avaliar qualidade
  --leading 3 \                           # Remove bases de baixa qualidade no início dos reads
  --trailing 3                            # Remove bases de baixa qualidade no final dos reads
