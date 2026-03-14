---
title: Sessão prática
parent: CUT&RUN
nav_order: 2
---

# Sessão prática

Nesta sessão iremos explorar **dados reais** utilizando ferramentas computacionais para análise genômica.

## Objetivos da prática

- carregar e explorar um conjunto de dados
- realizar análises exploratórias
- visualizar padrões nos dados
- interpretar resultados biológicos

## Conjunto de dados

O dataset utilizado neste exercício inclui:

- matriz de dados
- metadados das amostras

## Passo 1 — Carregar os dados

Exemplo em R:

```r
data <- read.csv("dataset.csv", row.names = 1)
head(data)
```

## Passo 2 — Explorar os dados

```r
summary(data)
dim(data)
```

## Passo 3 — Visualização

Visualizações comuns incluem:

- PCA
- heatmaps
- gráficos de distribuição

## Interpretação

Perguntas para considerar:

- Existem padrões claros entre as amostras?
- Amostras semelhantes agrupam juntas?
- Existem possíveis efeitos técnicos?
