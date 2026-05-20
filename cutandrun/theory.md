---
title: Teoria
parent: CUT&RUN
nav_order: 1
---
---
title: Teoria
parent: Interações proteína-DNA (CUT&RUN)
nav_order: 1
---

# Teoria

## O que é CUT&RUN?

CUT&RUN (Cleavage Under Targets and Release Using Nuclease) é uma técnica utilizada para mapear interações entre proteínas e DNA no genoma.

A abordagem utiliza anticorpos específicos para direcionar uma nuclease às regiões do DNA associadas a proteínas de interesse, como:

- fatores de transcrição
- histonas modificadas
- proteínas regulatórias

Após a clivagem do DNA nas regiões-alvo, os fragmentos liberados são sequenciados e analisados computacionalmente.

CUT&RUN permite identificar:

- sítios de ligação de fatores de transcrição
- marcas epigenéticas
- regiões regulatórias ativas ou reprimidas

<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutandrun_overview.png" width="500">
</div>

<p align="center">
<em>Figura 1. Visão geral da técnica de CUT&RUN. Fonte: https://elifesciences.org/articles/21856 </em>
</p>

---

## Perguntas biológicas comuns

CUT&RUN pode ser utilizado para investigar:

- sítios de ligação de fatores de transcrição
- distribuição de modificações de histonas
- elementos regulatórios ativos
- mecanismos de regulação gênica
- remodelamento epigenético

---

## Proteínas regulatórias e cromatina

A expressão gênica é controlada por proteínas capazes de interagir com regiões específicas do DNA.

Essas proteínas incluem:

- fatores de transcrição
- co-fatores
- proteínas remodeladoras
- histonas modificadas

CUT&RUN permite mapear essas interações diretamente no genoma.

---

## Modificações de histonas

As histonas podem sofrer diferentes modificações químicas associadas a estados regulatórios distintos.

### Exemplos comuns

| Marca | Associação biológica |
|---|---|
| H3K27ac | enhancers ativos |
| H3K4me3 | promotores ativos |
| H3K27me3 | repressão gênica |
| H3K4me1 | enhancers potenciais |

---

## Desenho experimental

Assim como em RNA-seq e ATAC-seq, o desenho experimental é fundamental.

### Conceitos importantes

- réplicas biológicas
- qualidade do anticorpo
- controles negativos
- profundidade de sequenciamento
- batch effects

### Controles importantes

Experimentos de CUT&RUN frequentemente incluem:

- IgG control
- no-antibody control
- input DNA (menos comum)

---

## Visão geral do workflow de CUT&RUN

Um experimento típico de CUT&RUN envolve várias etapas computacionais.

### Fluxo geral

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figura 2. Visão geral de um workflow típico de CUT&RUN. </em>
</p>

---

## Reads e arquivos FASTQ

O sequenciamento gera arquivos FASTQ contendo reads derivados das regiões ligadas à proteína de interesse.

Em CUT&RUN:

- os fragmentos costumam ser curtos
- o background geralmente é baixo
- o enriquecimento tende a ser altamente específico

---

## Controle de qualidade

O controle de qualidade em CUT&RUN avalia:

- qualidade dos reads
- eficiência do enriquecimento
- background experimental
- complexidade da biblioteca

Ferramentas comuns incluem:

- FastQC
- MultiQC
- deepTools

---

## Métricas importantes

### Qualidade dos reads

Avalia:

- qualidade por base
- conteúdo GC
- presença de adaptadores
- duplicação

---

### Taxa de alinhamento

Indica quantos reads alinham corretamente ao genoma de referência.

Baixas taxas podem indicar:

- contaminação
- baixa qualidade da biblioteca
- fragmentos muito curtos

---

### Distribuição dos fragmentos

CUT&RUN geralmente produz fragmentos menores e mais específicos do que ChIP-seq.

Distribuições anormais podem sugerir:

- digestão excessiva
- baixa eficiência do experimento
- degradação do DNA

---

## Alinhamento

Após QC, os reads são alinhados ao genoma de referência.

Ferramentas comuns incluem:

| Ferramenta | Características |
|---|---|
| Bowtie2 | amplamente utilizado |
| BWA | eficiente para reads curtos |

O objetivo é identificar onde os fragmentos se originaram no genoma.

---

## Arquivos BAM

Após alinhamento, os reads são armazenados em arquivos BAM contendo:

- posição genômica
- orientação
- qualidade do alinhamento
- informações de pareamento

Esses arquivos podem ser visualizados em genome browsers como:

- IGV
- UCSC Genome Browser

---

## Peak calling

Peak calling é a etapa utilizada para identificar regiões enriquecidas de reads.

Essas regiões representam potenciais sítios de ligação da proteína analisada.

### Ferramentas comuns

| Ferramenta | Características |
|---|---|
| MACS2/MACS3 | amplamente utilizado |
| SEACR | otimizado para CUT&RUN |

---

## Conceito de peak

Um peak representa:

```text
alta densidade local de fragmentos
→ possível interação proteína-DNA
```

<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutandrun_peaks.png" width="700">
</div>

<p align="center">
<em>Figura 3. Exemplo de peaks em CUT&RUN. fonte: https://genome.cshlp.org/content/30/1/35</em>
</p>

---

## Narrow peaks e broad peaks

Diferentes proteínas geram padrões distintos de enriquecimento.

### Narrow peaks

Associados a:

- fatores de transcrição
- ligação localizada

### Broad peaks

Associados a:

- modificações de histonas
- domínios cromatínicos amplos

---

## Quantificação de sinal

Após identificar os picos, é possível quantificar intensidade de sinal entre amostras.

O resultado geralmente é uma matriz:

| Peak | Sample_1 | Sample_2 |
|---|---|---|
| chr1:1-100 | 120 | 340 |
| chr2:200-300 | 540 | 210 |

Nessa matriz:

- linhas representam regiões enriquecidas
- colunas representam amostras
- valores representam abundância de fragmentos

---

## Normalização

Em CUT&RUN, diferenças técnicas entre bibliotecas podem influenciar a intensidade do sinal observado entre amostras.

Essas diferenças incluem:

- profundidade de sequenciamento
- eficiência de digestão
- enriquecimento do anticorpo
- complexidade da biblioteca
- background experimental

A normalização busca reduzir esses efeitos técnicos antes de comparações biológicas.

### Abordagens comuns

| Método | Objetivo |
|---|---|
| CPM | corrigir profundidade de sequenciamento |
| RPGC | normalização por cobertura genômica |
| DESeq2 size factors | comparação robusta entre amostras |

Após normalização, os dados podem ser utilizados para:

- comparação de enriquecimento
- geração de heatmaps
- metaplots
- análise diferencial de peaks

---

## Enriquecimento diferencial

CUT&RUN pode ser utilizado para comparar enriquecimento entre condições biológicas.

Essas análises permitem investigar:

- ganho ou perda de ligação
- remodelamento epigenético
- mudanças regulatórias

Ferramentas comuns incluem:

- DiffBind
- DESeq2
- edgeR

---

## Motif analysis

Após identificar peaks, é possível buscar sequências enriquecidas de DNA chamadas motifs.

Motifs podem indicar:

- fatores de transcrição ativos
- programas regulatórios
- redes gênicas potenciais

Ferramentas comuns incluem:

- HOMER
- MEME
- chromVAR

---

## Integração com RNA-seq e ATAC-seq

CUT&RUN frequentemente é integrado com outras abordagens genômicas.

### Exemplos

| Técnica | Informação |
|---|---|
| RNA-seq | expressão gênica |
| ATAC-seq | acessibilidade cromatínica |
| CUT&RUN | ligação proteína-DNA |

A integração dessas camadas permite investigar mecanismos regulatórios de forma mais completa.

---

## Limitações de CUT&RUN

Embora extremamente poderosa, a técnica possui limitações importantes.

### Algumas limitações incluem:

- dependência da qualidade do anticorpo
- resolução depende da proteína analisada
- background experimental ainda pode ocorrer
- proteínas pouco abundantes podem gerar sinal fraco
- interpretação funcional nem sempre é direta