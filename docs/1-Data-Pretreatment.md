# 数据预处理

## 物种


``` r
species = 'human'
```

## Counts 数据预处理


``` r
# 加载数据
data<-data.table::fread("raw_counts.txt", data.table = F)
head(data, 3)
```

```
##                ID a-1 a-2 a-3 b-1 b-2 b-3 c-1 c-2 c-3 d-1 d-2 d-3
## 1 ENSG00000282222   0   0   0   0   0   0   0   0   0   0   0   0
## 2 ENSG00000282221   2   0   0   0   0   0   0   0   0   0   0   0
## 3 ENSG00000308368   4   3   0   2   0   1   6   3   9   5   4   4
```

``` r
dim(data)
```

```
## [1] 80170    13
```

``` r
data <- data[!duplicated(data$ID),]

# ENSEMBLE ID to SYMBOL
library(AnnoProbe)
```

```
## AnnoProbe v 0.1.7  welcome to use AnnoProbe!
## If you use AnnoProbe in published research, please acknowledgements:
## We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.
```

``` r
ids <- annoGene(data$ID, ID_type = 'ENSEMBL',species = species)
```

```
## Warning in annoGene(data$ID, ID_type = "ENSEMBL", species = species): 27.64% of
## input IDs are fail to annotate...
```

``` r
ids=ids[!duplicated(ids$ENSEMBL),] # 重复ID直接删除
ids=ids[!duplicated(ids$SYMBOL),] # 重复ID直接删除

data = data[data$ID %in% ids$ENSEMBL,]
pos = match(data$ID, ids$ENSEMBL)
data = cbind(data, ids[pos, ])
rownames(data) = NULL
```


## 构建样本分组信息


``` r
library(stringr)
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ purrr     1.0.2
## ✔ forcats   1.0.0     ✔ readr     2.1.5
## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
count <- data %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  select(-c("ID", "biotypes", "ENSEMBL", "chr", "start", "end"))

## 也可以手动创建一个分组信息
group_list <- str_split(colnames(count),'[-_]',simplify = T)[,1]
group <- data.frame(sample = colnames(count), 
                    group = group_list)

table(group)
```

```
##       group
## sample a b c d
##    a-1 1 0 0 0
##    a-2 1 0 0 0
##    a-3 1 0 0 0
##    b-1 0 1 0 0
##    b-2 0 1 0 0
##    b-3 0 1 0 0
##    c-1 0 0 1 0
##    c-2 0 0 1 0
##    c-3 0 0 1 0
##    d-1 0 0 0 1
##    d-2 0 0 0 1
##    d-3 0 0 0 1
```

## 构建组间比较信息


``` r
group_compare <- list(
  a_vs_b = c("a", "b") %>% factor(levels = c("a", "b")),
  c_vs_d = c("c", "d") %>% factor(levels = c("c", "d")),
  a_vs_c = c("a", "c") %>% factor(levels = c("a", "c"))
)
```

## 保存数据


``` r
dir.create("Expression_Annotation/1_Expression_result", recursive = T)
```

```
## Warning in dir.create("Expression_Annotation/1_Expression_result", recursive =
## T): 'Expression_Annotation/1_Expression_result' already exists
```

``` r
library(writexl)
write_xlsx(data, "Expression_Annotation/1_Expression_result/All_gene_counts.xlsx")
```


``` r
dir.create("Rdata", recursive = T)
```

```
## Warning in dir.create("Rdata", recursive = T): 'Rdata' already exists
```

``` r
save(species, count, data, group, group_compare, file = "Rdata/data.RData")
```


