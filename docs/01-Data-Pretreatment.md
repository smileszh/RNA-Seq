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
##                ID a-1 a-2 a-3 b-1 b-2 b-3 c-1 c-2 c-3 d-1 d-2 d-3
## 1 ENSG00000282222   0   0   0   0   0   0   0   0   0   0   0   0
## 2 ENSG00000282221   2   0   0   0   0   0   0   0   0   0   0   0
## 3 ENSG00000308368   4   3   0   2   0   1   6   3   9   5   4   4
dim(data)
## [1] 80170    13

data <- data[!duplicated(data$ID),]

# ENSEMBLE ID to SYMBOL
library(AnnoProbe)
ids <- annoGene(data$ID, ID_type = 'ENSEMBL',species = species)
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
count <- data %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  select(-c("ID", "biotypes", "ENSEMBL", "chr", "start", "end"))

## 也可以手动创建一个分组信息
group_list <- str_split(colnames(count),'[-_]',simplify = T)[,1]
group <- data.frame(sample = colnames(count), 
                    group = group_list)

table(group)
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
library(writexl)
write_xlsx(data, "Expression_Annotation/1_Expression_result/All_gene_counts.xlsx")
```


``` r
dir.create("Rdata", recursive = T)
save(species, count, data, group, group_compare, file = "Rdata/data.RData")
```


