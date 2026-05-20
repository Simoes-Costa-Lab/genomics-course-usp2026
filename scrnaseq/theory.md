---
title: Teoria
parent: Single‑cell RNA-seq
nav_order: 1
---


# Teoria

## O que é scRNA-seq?

Single-cell RNA-seq (scRNA-seq) é uma técnica utilizada para medir expressão gênica individualmente em milhares de células simultaneamente.

Diferentemente do RNA-seq convencional (bulk RNA-seq), em que o sinal observado representa a média de uma população celular, scRNA-seq permite investigar heterogeneidade celular.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrnaseq_overview.png" width="700">
</div>

<p align="center">
<em>Figura 1. Visão geral da técnica de single-cell RNA-seq. Fonte: https://doi.org/10.1371/journal.ppat.1011898 </em>
</p>

---

## Perguntas biológicas comuns

scRNA-seq pode ser utilizado para investigar:

- composição celular de tecidos
- desenvolvimento embrionário
- heterogeneidade tumoral
- resposta celular a estímulos
- trajetórias de diferenciação
- estados celulares transitórios

---

## Bulk RNA-seq vs scRNA-seq

No bulk RNA-seq, o sinal representa a média de milhares ou milhões de células.

Em contraste, scRNA-seq mede expressão gênica célula por célula.

### Comparação

| Bulk RNA-seq | scRNA-seq |
|---|---|
| média populacional | resolução celular |
| menor ruído técnico | maior variabilidade |
| menor custo | maior custo computacional |
| expressão média | heterogeneidade celular |

<div align="center">
<img src="/genomics-course-usp2026/assets/images/bulk_vs_singlecell.png" width="700">
</div>

<p align="center">
<em>Figura 2. Comparação entre bulk RNA-seq e scRNA-seq. Fonte: https://www.completegenomics.com/methods/single-cell-rna-sequencing/ </em>
</p>

---

## Isolamento celular e barcoding

Em plataformas como 10x Genomics, células individuais são encapsuladas em droplets contendo beads com barcodes moleculares.

Cada molécula de RNA recebe:

- barcode celular
- UMI (Unique Molecular Identifier)

### Barcodes celulares

Permitem identificar de qual célula cada read se originou.

### UMIs

Permitem identificar moléculas únicas e reduzir viés de PCR.

---

## Visão geral do workflow de scRNA-seq

### Fluxo geral

<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrnaseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figura 2. Visão geral de um workflow típico de single cell RNA-Seq. </em>
</p>

---

## Matriz célula-gene

Após alinhamento e quantificação, os dados são organizados em uma matriz:

| Gene | Cell_1 | Cell_2 | ... | Cell_n |
|---|---|---|---|---|
| SOX10 | 4 | 0 | ... | 2 |
| PAX7 | 12 | 3 | ... | 0 |
| TFAP2A | 0 | 8 | ... | 23 |

Nessa matriz:

- linhas representam genes
- colunas representam células
- valores representam abundância de RNA

---

## Sparsity e dropouts

Matrizes de scRNA-seq costumam ser extremamente esparsas.

Muitos valores aparecem como zero devido a:

- baixa captura de RNA
- baixa profundidade
- expressão muito baixa
- limitações técnicas

Esses zeros são frequentemente chamados de dropouts.

---

## Controle de qualidade

O controle de qualidade em scRNA-seq é realizado no nível de células individuais.

O objetivo é remover células de baixa qualidade antes da análise downstream.

### Métricas comuns

| Métrica | Interpretação |
|---|---|
| número de genes detectados | complexidade celular |
| número total de UMIs | profundidade |
| porcentagem mitocondrial | integridade celular |

---

## Células de baixa qualidade

Células problemáticas frequentemente apresentam:

- poucos genes detectados
- baixa contagem total
- alta porcentagem de genes mitocondriais

Essas células geralmente correspondem a:

- células mortas
- debris
- droplets vazios

---

## Doublets

Em alguns casos, duas células podem ser encapsuladas no mesmo droplet.

Esses eventos são chamados de doublets.

Doublets podem gerar perfis artificiais misturando dois tipos celulares.

Ferramentas comuns incluem:

- DoubletFinder
- Scrublet

---

## Normalização

Células individuais possuem diferentes profundidades de sequenciamento.

A normalização busca corrigir essas diferenças técnicas antes da comparação entre células.

### Objetivos da normalização

- corrigir diferenças de profundidade
- estabilizar variância
- permitir comparação entre células

### Métodos comuns

| Método | Ferramenta |
|---|---|
| LogNormalize | Seurat |
| SCTransform | Seurat |
| CPM normalization | abordagens clássicas |

---

## Seleção de genes variáveis

Nem todos os genes são igualmente informativos.

Genes altamente variáveis ajudam a identificar:

- diferenças celulares
- estados biológicos
- populações distintas

Esses genes são frequentemente utilizados em análises downstream.

---

## Redução de dimensionalidade

Experimentos de scRNA-seq possuem milhares de genes por célula.

Técnicas de redução de dimensionalidade ajudam a resumir os dados.

### Métodos comuns

| Método | Objetivo |
|---|---|
| PCA | redução linear |
| UMAP | visualização |
| t-SNE | visualização local |

<div align="center">
<img src="/genomics-course-usp2026/assets/images/umap_example.png" width="700">
</div>

<p align="center">
<em>Figura 3. Exemplo de projeção UMAP em scRNA-seq. Fonte: https://www.nature.com/articles/s41586-023-05869-0 </em>
</p>

---

## Clustering celular

Após a redução de dimensionalidade, células semelhantes podem ser agrupadas.

Esses agrupamentos frequentemente correspondem a:

- tipos celulares
- estados celulares
- populações biológicas distintas

Ferramentas comuns incluem:

- Seurat
- Scanpy

---

## Anotação celular

Após identificar os clusters, é necessário interpretar biologicamente cada população.

Isso geralmente é feito utilizando:

- genes marcadores conhecidos
- bancos de referência
- literatura

---

## Expressão diferencial

scRNA-seq também permite comparar expressão gênica entre:

- clusters
- condições
- estados celulares

---

## Limitações de scRNA-seq

### Algumas limitações incluem:

- alta sparsity
- dropouts técnicos
- maior ruído experimental
- custo computacional elevado
- perda de informação espacial
- dissociação celular pode alterar expressão gênica