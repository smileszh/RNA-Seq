# 差异表达分析


## 差异表达分析

### DESeq2 




``` r
# 安装加载必要的包
library(DESeq2)
library(tidyverse)

# 定义一个函数来进行 DESeq2 差异基因表达分析
run_deseq2 <- function(count, group_df, group1, group2) {
  samples_group1 <- group_df %>% filter(group == group1) %>% pull(sample)
  samples_group2 <- group_df %>% filter(group == group2) %>% pull(sample)
  
  selected_samples <- c(samples_group1, samples_group2)
  selected_counts <- count[, selected_samples, drop = FALSE]
  selected_group <- group_df %>% filter(sample %in% selected_samples)
  
  # 创建 DESeq2 数据集
  dds <- DESeqDataSetFromMatrix(countData = selected_counts, 
                                colData = selected_group, 
                                design = ~ group)
  
  # 去除低表达基因
  keep <- rowSums(counts(dds)) > 3
  # 或者更严格筛选至少3个样本中计数大于或等于10的基因
  # keep <- rowSums(counts(dds) >= 10) >= 3

  dds <- dds[keep, ]
  
  # 运行 DESeq2
  dds <- DESeq(dds)
  
  # 获取结果
  res <- results(dds) %>%
    as.data.frame() %>%
    rownames_to_column(var = "SYMBOL")
  
  # 得到标准化后的count
  normed = counts(dds, normalized=TRUE) %>%
    round(1)
  
  # 整理结果
  
  ## 创建 foldChange 列.
  res$foldChange = 2 ^ res$log2FoldChange
  
  ## 创建 FDR 列，原因[在此](https://www.biostars.org/p/462897/).
  res =  dplyr::rename(res, FDR=padj)
  res$FDR[is.na(res$FDR)] <- 1
  
  ## 创建 padj 列.
  res$padj = p.adjust(res$pvalue, method="hochberg")
  

  ## 合并结果和标准化后的 count 
  total <- bind_cols(res, normed) %>%
    arrange(FDR)
  
  # 创建其他想呈现的信息

  ## 计算排序表的错误发现计数。
  total$falsePos = 1:nrow(total) * total$FDR

  ## 创建各族均值
  total$baseMeanA = 1
  total$baseMeanB = 1
  
  total$baseMeanA = rowMeans(total[, samples_group1])
  total$baseMeanB = rowMeans(total[, samples_group2])

  
  # 整理结果，美化输出
  # total$foldChange = round(total$foldChange, 3)
  # total$log2FoldChange = round(total$log2FoldChange, 1)
  # total$baseMean  = round(total$baseMean, 1)
  # total$baseMeanA = round(total$baseMeanA, 1)
  # total$baseMeanB =  round(total$baseMeanB, 1)
  # total$lfcSE = round(total$lfcSE, 2)
  # total$stat = round(total$stat, 2)
  # total$FDR = round(total$FDR, 4)
  # total$falsePos = round(total$falsePos, 0)

  new_cols = c(
             "SYMBOL","baseMean","baseMeanA","baseMeanB",
             "foldChange", "log2FoldChange",
             "lfcSE","stat","pvalue","padj", 
             "FDR","falsePos", samples_group1, samples_group2)

  total = total[, new_cols]
  

  return(total)
}
```


``` r
# 进行两两比对并存储结果
dir.create("DEG_Analysis/1_DEG_stat/DESeq2", recursive = T)

comparison_results_DESeq2 <- lapply(names(group_compare), function(comp_name) {
  # 获取因子的水平
  # comp_name <- "a_vs_b"
  groups <- levels(group_compare[[comp_name]])
  group1 <- groups[1]
  group2 <- groups[2]
  
  # 打印调试信息
  print(paste("Comparing", group1, "vs", group2))
  print(head(group))
  
  # 运行 DESeq2
  result <- run_deseq2(count, group, group1, group2)

  
  # 保存结果
  ## 定义结果文件名
  result_file <- paste0(comp_name, "_DEG_DESeq2.csv")
  
  ## 将结果保存为 CSV 文件
  write.csv(result, paste0("DEG_Analysis/1_DEG_stat/DESeq2/", result_file) , row.names = FALSE)
  
  return(result)
})
## [1] "Comparing a vs b"
##   sample group
## 1    a-1     a
## 2    a-2     a
## 3    a-3     a
## 4    b-1     b
## 5    b-2     b
## 6    b-3     b
## [1] "Comparing c vs d"
##   sample group
## 1    a-1     a
## 2    a-2     a
## 3    a-3     a
## 4    b-1     b
## 5    b-2     b
## 6    b-3     b
## [1] "Comparing a vs c"
##   sample group
## 1    a-1     a
## 2    a-2     a
## 3    a-3     a
## 4    b-1     b
## 5    b-2     b
## 6    b-3     b

names(comparison_results_DESeq2) <- names(group_compare)
```



``` r
save(comparison_results_DESeq2, file = "Rdata/comparison_results_DESeq2.Rdata")
```


