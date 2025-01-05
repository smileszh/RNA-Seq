# 基因表达定量


## 基因表达定量



### FPKM / RPKM

转录本的片段数量受测序数据量（或比对到的reads数量）、转录本长度、转录本表达水平的影响。为了更准确地揭示每个转录本的表达水平，需要通过其转录本的长度对比对到的reads数量进行标准化。

FPKM（每千碱基转录本每百万比对片段数）被应用于测量基因或转录本的表达水平。FPKM的计算公式如下所示。

$$
\text{FPKM} = \frac{\text{number of fragments}}{\text{length of transcript (in kilobases)} \times \text{total number of mapped reads (in millions)}}
$$


``` r
# 下载并加载必要的包
library(tidyverse)
library(writexl)

counts <- count %>% as.matrix()
# 计算每个样本中的总的测序读数
total_count <- colSums(count)

# 计算基因长度（以kb为单位）
gene_positions <- data %>%
  select(SYMBOL, start, end) %>% 
  distinct()
gene_positions$Length_kb <- (gene_positions$end - gene_positions$start + 1) / 1000

# 创建 DGEList 对象
library(edgeR)
dge <- DGEList(counts = count)

# 计算FPKM
# edgeR 使用 RPKM 函数来计算 RPKM 和 FPKM
# 传入基因长度（以kb为单位）和总的测序读数
fpkm <- rpkm(dge, gene.length = gene_positions$Length_kb * 1000) %>%
  as.data.frame()


## 保存结果
dir.create("Expression_Annotation/1_Expression_result", recursive = T)
All_gene_fpkm <- fpkm %>%
  rename_with(~ paste0(.x, ".FPKM")) %>%
  rownames_to_column(var = "SYMBOL") %>%
  left_join(data %>% select(ID, SYMBOL, everything()), by = "SYMBOL")

writexl::write_xlsx(All_gene_fpkm, "Expression_Annotation/1_Expression_result/All_gene_fpkm.xlsx")
```



``` r
save(species, count, fpkm, data, group, group_compare, file = "Rdata/data.RData")
```

