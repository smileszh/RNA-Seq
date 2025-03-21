# (PART) counts 数据预处理 {-}



# 数据预处理

## 加载数据

这一步是加载我们的原始数据，分为两个部分，一个是实验设计文件，一个是原始的表达矩阵。
可以根据需求替换成自己的数据。但是需要注意的是，实验设计文件中的分组信息需要和表达矩阵中的样本名一一对应。
设计文件中，有sample 和 group 两列，分别代表样本名和分组信息。
表达矩阵中，第一列是基因名，后面的列是样本名，且样本名顺序和设计文件中的样本名顺序一致。

``` r
# 加载包
library(data.table)
library(tidyverse)

# 加载分组设计文件
design <- fread("design.csv") %>%
  mutate(group = factor(group, levels = c("s_MiaFRFP", "s_Part1Lag")))
group_list <- design$group
head(design)
##            run      group       sample
##         <char>     <fctr>       <char>
## 1: SRR19913146  s_MiaFRFP  s_MiaFRFP_1
## 2: SRR19913145  s_MiaFRFP  s_MiaFRFP_2
## 3: SRR19913144  s_MiaFRFP  s_MiaFRFP_3
## 4: SRR19913143 s_Part1Lag s_Part1Lag_1
## 5: SRR19913142 s_Part1Lag s_Part1Lag_2
## 6: SRR19913141 s_Part1Lag s_Part1Lag_3
head(group_list)
## [1] s_MiaFRFP  s_MiaFRFP  s_MiaFRFP  s_Part1Lag s_Part1Lag s_Part1Lag
## Levels: s_MiaFRFP s_Part1Lag

# 加载表达矩阵，去除重复基因，并将样本名顺序和设计文件中的样本名顺序一致
counts <- fread("counts.csv") %>%
  distinct(name, .keep_all = T) %>%
  select(name, design$sample, everything())
head(counts)
##               name s_MiaFRFP_1 s_MiaFRFP_2 s_MiaFRFP_3 s_Part1Lag_1
##             <char>       <int>       <int>       <int>        <int>
## 1: ENSG00000228037           2           0           0            2
## 2: ENSG00000142611        2088        1753        2396          358
## 3: ENSG00000284616           0           0           0            0
## 4: ENSG00000157911        5039        3794        4870         2294
## 5: ENSG00000260972           0           0           0            0
## 6: ENSG00000224340           0           0           1            1
##    s_Part1Lag_2 s_Part1Lag_3
##           <int>        <int>
## 1:            9            0
## 2:          272          264
## 3:            0            2
## 4:         1996         2352
## 5:            0            0
## 6:            0            2
```


## 过滤低表达数据
为了避免低表达基因对后续分析的影响，我们需要过滤掉低表达基因。


``` r
# 将基因 ID 作为行名
mat <- dplyr::select(counts, design$sample)
rownames(mat) <- counts$name
# 筛选出至少在两个样本中表达的基因
keep_feature <- rowSums (mat > 1) >= 2 ;table(keep_feature)
## keep_feature
## FALSE  TRUE 
## 32850 30290
ensembl_matrix <- mat[keep_feature, ]  
rownames(ensembl_matrix) <- rownames(mat)[keep_feature]
ensembl_matrix <- rownames_to_column(ensembl_matrix, var = "name")
```


## 转换基因 ID
一般情况下，我们会使用 SYMBOL ID 进行分析，因为SYMBOL ID 的命名与基因的功能更加接近，
所以需要将 ENSEMBL ID 转换为 SYMBOL ID，也方便后续我们的分析。


``` r
# 加载包
library(biomaRt)
# 连接ensembl数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# atr <- listAttributes(ensembl) # 查看数据库中有哪些信息
# 获得数据库中关注的 ID 信息
ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","entrezgene_id"),
             mart = ensembl)
# 修改列名
colnames(ids) <- c("ENSEMBL", "SYMBOL", "ENTREZ")
# 去除重复 ID
ids=ids[!duplicated(ids$ENSEMBL),]
# 连接基因 counts 矩阵和 ID 信息
All_gene_counts <- merge(ensembl_matrix, ids, by.x = "name", by.y = "ENSEMBL", all.x = TRUE) %>%
  dplyr::select(name, SYMBOL, ENTREZ, everything())
# 输出注释结果
dir.create("Expression_result", showWarnings = F)
write.csv(All_gene_counts, "Expression_result/All_gene_counts.csv", row.names = F)
```

## 保存变量

``` r
# 以 SYMBOL 为行名的矩阵
symbol_matrix <- All_gene_counts %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = T) %>%
  column_to_rownames("SYMBOL") %>%
  dplyr::select(design$sample)
save(symbol_matrix, design, group_list,file = 'Rdata/symbol_matrix.Rdata')
```
