## 差异表达基因筛选





``` r
load("Rdata/comparison_results_DESeq2.Rdata")
load("Rdata/comparison_results_edgeR.Rdata")
```


``` r
# 可以自定义设定筛选差异基因的参数
FDR_cutoff <- 0.05
log2FC_cutoff <- 1
```


``` r
# 加载必要的包
library(tidyverse)

# 定义筛选差异基因的函数
filter_DEGs <- function(df, FDR_cutoff = 0.05, log2FC_cutoff = 1) {
  df <- df %>%
    mutate(regulated = case_when(
      FDR < FDR_cutoff & log2FoldChange > log2FC_cutoff ~ "up",
      FDR < FDR_cutoff & log2FoldChange < -log2FC_cutoff ~ "down",
      TRUE ~ "normal"
    ))
  return(df)
}

# 对列表中的每个数据框进行处理
cr_DESeq2 <- setNames(lapply(comparison_results_DESeq2, filter_DEGs),
                      names(comparison_results_DESeq2))
cr_edgeR <- setNames(lapply(comparison_results_edgeR, filter_DEGs),
                     names(comparison_results_edgeR))

# 导出处理后的数据框
for (i in names(cr_DESeq2)) {
  write.csv(cr_DESeq2[[i]], 
            paste0("DEG_Analysis/1_DEG_stat/DESeq2/", i, "_regulated_DESeq2.csv"),
            row.names = FALSE)
}

for (i in names(cr_edgeR)) {
  write.csv(cr_edgeR[[i]],
            paste0("DEG_Analysis/1_DEG_stat/edgeR/", i, "_regulated_edgeR.csv"), 
            row.names = FALSE)
}
```


``` r
save(cr_DESeq2, cr_edgeR, file = "Rdata/Regulated_DEG.Rdata")
```
