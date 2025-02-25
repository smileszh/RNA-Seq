### edgeR




``` r
# 安装加载必要的包
library(edgeR)
```

```
## Loading required package: limma
```

``` r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
# 定义一个函数来进行 edgeR 差异基因表达分析
run_edger <- function(count, group_df, group1, group2) {
  samples_group1 <- group_df %>% filter(group == group1) %>% pull(sample)
  samples_group2 <- group_df %>% filter(group == group2) %>% pull(sample)
  
  selected_samples <- c(samples_group1, samples_group2)
  selected_counts <- count[, selected_samples, drop = FALSE]
  selected_group <- group_df %>% filter(sample %in% selected_samples)
  
  # 创建 DGEList 对象
  group <- factor(selected_group$group)
  deg <- DGEList(counts = selected_counts, group = group)
  
  # 过滤低表达基因
  keep <- filterByExpr(deg)
  deg <- deg[keep, , keep.lib.sizes = FALSE]
  
  # 规范化
  deg <- calcNormFactors(deg)
  
  # 设计矩阵
  design <- model.matrix(~ group)
  
  # 估计离散度
  dge <- estimateDisp(deg, design)
  
  # 进行广义线性模型拟合
  fit <- glmFit(dge, design)
  
  # 进行差异表达检验
  qlf <- glmLRT(fit, coef = 2)
  
  # 获取结果
  etp = topTags(qlf, n=Inf)
  
  ## 获取数据 scale
  scale = dge$samples$lib.size * dge$samples$norm.factors

  ## 获取标准化 counts
  normed = round(t(t(selected_counts)/scale) * mean(scale))
  
  ## 获取差异分析结果
  dat = etp$table %>%
    as.data.frame() %>%
    rownames_to_column(var = "SYMBOL")
  row.names(dat) = dat[,1]
  
  # 整理结果
  ## Create column placeholders.
  dat$baseMean = 1
  dat$baseMeanA = 1
  dat$baseMeanB = 1
  dat$foldChange = 2 ^ dat$logFC
  dat$log2FoldChange=dat$logFC
  dat$pvalue = dat$PValue
  dat$falsePos = 1
  
  # Compute the adjusted p-value
  dat$padj = p.adjust(dat$pvalue, method="hochberg")

  # Create a merged output that contains the normalized counts.
  total <- merge(dat, normed, by='row.names')

  # Get rid of extra column it gained
  total = total[, 2:ncol(total)]

  # Sort again by P-value.
  total = arrange(total, pvalue)
  
  # Compute the false discovery counts on the sorted table.
  total$falsePos = 1:nrow(total) * total$FDR
  
  # Create the individual baseMean columns.
  total$baseMeanA = rowMeans(total[, samples_group1])
  total$baseMeanB = rowMeans(total[, samples_group2])
  total$baseMean = total$baseMeanA + total$baseMeanB
  
  # Round the numbers to make them look better
  # total$foldChange = round(total$foldChange, 3)
  # total$FDR = round(total$FDR, 4)
  # total$padj = round(total$padj, 4)
  # total$logCPM = round(total$logCPM, 1)
  # total$log2FoldChange = round(total$log2FoldChange, 1)
  # total$baseMean = round(total$baseMean, 1)
  # total$baseMeanA = round(total$baseMeanA, 1)
  # total$baseMeanB = round(total$baseMeanB, 1)
  # total$falsePos = round(total$falsePos, 0)

  # Reorganize columns names to make more sense.
  new_cols = c( "SYMBOL","baseMean","baseMeanA","baseMeanB",
                "logCPM","foldChange", "log2FoldChange",
                "pvalue","padj", 
                "FDR","falsePos", samples_group1, samples_group2)

  total = total[, new_cols]


  # Reformat these columns as string.
  total$padj = formatC(total$padj, format = "e", digits = 1)
  total$pvalue = formatC(total$pvalue, format = "e", digits = 1)
  
  return(total)
}
```



``` r
# 进行两两比对并存储结果
dir.create("DEG_Analysis/1_DEG_stat/edgeR", recursive = TRUE)
```

```
## Warning in dir.create("DEG_Analysis/1_DEG_stat/edgeR", recursive = TRUE):
## 'DEG_Analysis/1_DEG_stat/edgeR' already exists
```

``` r
comparison_results_edgeR <- lapply(names(group_compare), function(comp_name) {
  # 获取因子的水平
  groups <- levels(group_compare[[comp_name]])
  group1 <- groups[1]
  group2 <- groups[2]
  
  # 打印调试信息
  print(paste("Comparing", group1, "vs", group2))
  print(head(group))
  
  # 运行 edgeR
  result <- run_edger(count, group, group1, group2)
  
  # 保存结果
  result_file <- paste0(comp_name, "_DEG_edgeR.csv")
  write.csv(result, paste0("DEG_Analysis/1_DEG_stat/edgeR/", result_file), row.names = FALSE)
  
  return(result)
})
```

```
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
```

``` r
names(comparison_results_edgeR) <- names(group_compare)
```



``` r
save(comparison_results_edgeR, file = "Rdata/comparison_results_edgeR.Rdata")
```
