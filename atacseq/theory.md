---
title: Teoria
parent: ATAC-seq
nav_order: 1
---

# Teoria

## O que é ATAC-seq?

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) é uma técnica utilizada para identificar regiões acessíveis da cromatina em uma população celular.

A técnica utiliza uma transposase hiperativa (Tn5) que fragmenta preferencialmente regiões abertas da cromatina e simultaneamente adiciona adaptadores de sequenciamento.

As regiões acessíveis identificadas por ATAC-seq frequentemente correspondem a:

- enhancers
- promotores
- regiões regulatórias ativas
- sítios de ligação de fatores de transcrição

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_overview.png" width="700">
</div>

<p align="center">
<em>Figura 1. Visão geral da técnica de ATAC-seq. Template do Biorender </em>
</p>

---

## Perguntas biológicas comuns

ATAC-seq pode ser utilizado para investigar:

- regiões regulatórias ativas
- resposta epigenética a estímulos
- diferenças regulatórias entre tecidos
- fatores de transcrição potencialmente ativos

---

## Cromatina aberta e fechada

O DNA nuclear está organizado em estruturas chamadas nucleossomos.

Regiões altamente compactadas tendem a ser menos acessíveis à maquinaria transcricional.

Já regiões abertas da cromatina permitem:

- ligação de fatores de transcrição
- recrutamento de RNA polimerase
- ativação gênica

ATAC-seq explora justamente essa diferença de acessibilidade.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/chromatin_accessibility.png" width="700">
</div>

<p align="center">
<em>Figura 2. Regiões abertas e fechadas da cromatina. Fonte: https://www.nature.com/articles/s41576-018-0089-8 </em>
</p>

---

## Desenho experimental

Assim como em RNA-seq, o desenho experimental é fundamental para experimentos de ATAC-seq.

### Conceitos importantes

- réplicas biológicas
- qualidade nuclear
- número de células
- batch effects
- profundidade de sequenciamento

### Considerações importantes

ATAC-seq é particularmente sensível à:

- degradação celular
- lise excessiva
- contaminação mitocondrial
- baixa qualidade nuclear

---

## Visão geral do workflow de ATAC-seq

Um experimento típico de ATAC-seq envolve várias etapas computacionais.

### Fluxo geral

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figura 3. Visão geral de um workflow típico de ATAC-seq. </em>
</p>


---

## Reads e arquivos FASTQ

Assim como em RNA-seq, o sequenciamento gera arquivos FASTQ contendo reads de DNA.

Cada read representa um fragmento de DNA acessível identificado pela transposase Tn5.

### Particularidades de ATAC-seq

Em ATAC-seq:

- fragmentos curtos geralmente correspondem a regiões nucleosome-free
- fragmentos maiores podem refletir posicionamento nucleossomal
- reads mitocondriais frequentemente são abundantes

---

## Controle de qualidade

O controle de qualidade em ATAC-seq é essencial para avaliar:

- enriquecimento em regiões abertas
- qualidade da biblioteca
- sinal versus ruído

Ferramentas comuns incluem:

- FastQC
- MultiQC
- deepTools

---

## Métricas importantes em ATAC-seq

### Qualidade dos reads

Avalia:

- qualidade por base
- conteúdo GC
- adaptadores
- duplicação

---

### Taxa de alinhamento

Mede quantos reads alinham corretamente ao genoma de referência.

Baixas taxas podem indicar:

- contaminação
- baixa qualidade
- problemas na biblioteca

---

### Reads mitocondriais

ATAC-seq frequentemente gera muitos reads derivados de DNA mitocondrial.

Altas proporções de reads mitocondriais podem indicar:

- lise excessiva
- baixa qualidade nuclear

---

### Distribuição do tamanho dos fragmentos

Uma das métricas mais características de ATAC-seq.

Fragmentos pequenos geralmente correspondem a regiões abertas livres de nucleossomos.

Fragmentos maiores podem refletir:

- mono-nucleossomos
- di-nucleossomos
- organização da cromatina

<div align="center">
<img src="/genomics-course-usp2026/assets/images/fragment_distribution.png" width="700">
</div>

<p align="center">
<em>Figura 4. Distribuição de fragmentos em ATAC-seq. Fonte: https://www.activemotif.com/blog-library-qc </em>
</p>

---

## Alinhamento

Após a etapa de QC, os reads são alinhados ao genoma de referência.

Ferramentas comuns incluem:

| Ferramenta | Características |
|---|---|
| Bowtie2 | amplamente utilizado em ATAC-seq |
| BWA | eficiente para reads curtos |

O objetivo é determinar a posição genômica dos fragmentos acessíveis.

### Arquivos BAM

Após o alinhamento, os reads são armazenados em arquivos BAM contendo:

- posição genômica
- orientação
- qualidade do alinhamento
- informações de pareamento

---

## Remoção de duplicatas

Durante PCR, fragmentos podem ser amplificados excessivamente. Como estamos tratando de DNA e não mais de RNA, em que uma mesma molécula de RNA biologicamente tem muitas cópias, precisamos limpar duplicatas geradas pela técnica.

Reads duplicados podem inflar artificialmente o sinal.

Ferramentas comuns incluem:

- samtools
- Picard

---

## Peak calling

Uma das etapas centrais de ATAC-seq.

O objetivo do peak calling é identificar regiões do genoma com enriquecimento significativo de reads.

Essas regiões representam potenciais elementos regulatórios acessíveis.

### Ferramentas comuns

| Ferramenta | Características |
|---|---|
| MACS2/MACS3 | padrão mais utilizado |
| Genrich | otimizado para ATAC-seq |

### Conceito de peak

Um peak representa:

```text
alta densidade local de fragmentos
→ maior acessibilidade cromatínica
```

<div align="center">
<img src="/genomics-course-usp2026/assets/images/peak_calling.png" width="700">
</div>

<p align="center">
<em>Figura 4. Exemplo de identificação de peaks em ATAC-seq. Fonte: https://www.nature.com/articles/s41467-025-67491-0</em>
</p>

---

## Quantificação de acessibilidade

Após identificar os peaks, é possível quantificar acessibilidade da cromatina entre amostras.

O resultado geralmente é uma matriz:

| Peak | Sample_1 | Sample_2 |
|---|---|---|
| chr1:1-100 | 120 | 340 |
| chr2:200-300 | 540 | 210 |

Nessa matriz:

- linhas representam regiões acessíveis
- colunas representam amostras
- valores representam abundância de fragmentos

---

## Normalização

Assim como em RNA-seq, contagens brutas de ATAC-seq não são diretamente comparáveis entre amostras.

Diferenças técnicas podem surgir devido a:

- profundidade de sequenciamento
- eficiência da transposição
- composição da biblioteca
- proporção de reads mitocondriais
- número total de peaks detectados

A normalização busca corrigir essas diferenças para permitir comparações biológicas mais confiáveis.

### Abordagens comuns

| Método | Objetivo |
|---|---|
| CPM | corrigir profundidade de sequenciamento |
| TMM | normalização robusta entre amostras |
| DESeq2 size factors | modelagem estatística de contagens |

Após normalização, as matrizes podem ser utilizadas em análises downstream como:

- PCA
- clustering
- acessibilidade diferencial
- heatmaps

---

## Acessibilidade diferencial

ATAC-seq pode ser utilizado para identificar regiões diferencialmente acessíveis entre condições.

Essas análises permitem investigar:

- ativação regulatória
- remodelamento cromatínico
- mudanças epigenéticas

Ferramentas comuns incluem:

- DESeq2
- edgeR
- DiffBind

---

## Limitações de ATAC-seq

Embora extremamente poderosa, a técnica possui limitações importantes.

### Algumas limitações incluem:

- não mede ligação direta de proteínas
- acessibilidade não implica atividade funcional
- regiões repetitivas dificultam alinhamento
- sensível à qualidade celular
- resolução limitada em populações heterogêneas