


# DESeq2

## 加载数据

``` r
# 加载包
bio.p <- c("DESeq2", "tibble", "dplyr", "tools")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
  )

load(file = 'Rdata/symbol_matrix.Rdata')
group <- group_list # 分组标准名字
count_data <- symbol_matrix # 计数矩阵标准名字

# title
title <- paste0(group[1], "_vs_", group[2])

# 输出文件夹
OUT_FOLDER <- "DEG_Analysis"
dir.create(OUT_FOLDER, showWarnings = FALSE)
# 输出文件名
OUTPUT_FILE <- paste0(title, "_all_DESeq2.csv")

# 创建一个样本分组信息
df = data.frame(
  sample = colnames(symbol_matrix),
  group = group
)
# 样本名
sample_names <- df$sample
col_names_A <- filter(df, group == levels(df$group)[1])$sample
col_names_B <- filter(df, group == levels(df$group)[2])$sample
```


## DESeq2 分析

``` r
# 创建 DESeq2 数据集
dds = DESeqDataSetFromMatrix(
  countData = count_data,
  colData = df,
  design = ~ group
)

# 过滤低表达基因
keep <- rowSums(counts(dds)) > 3
# 严格过滤：至少 3 个样本有 10 个以上的计数
# keep <- rowSums(counts(dds) >= 10) >= 3
# 过滤
dds <- dds[keep, ]

# 运行 DESeq2
dds = DESeq(dds)

# 得到归一化的 counts 矩阵
normed = counts(dds, normalized=TRUE)
normed = round(normed, 1)

# 提取差异基因
res = results(dds)
```

## 整理分析结果，规范输出

``` r
# 差异分析结果
data = as.data.frame(res)

# 添加想要的信息
data$foldChange = 2 ^ data$log2FoldChange
data =  dplyr::rename(data, PValue=pvalue, FDR=padj)
data$FDR[is.na(data$FDR)] <- 1
data$PAdj = p.adjust(data$PValue, method="hochberg")
data$baseMeanA = 1
data$baseMeanB = 1

# 合并差异分析结果和标准化的 counts
total <- bind_cols(data, normed) %>%
  rownames_to_column("SYMBOL")

# 按照 P-value排序
total = arrange(total, FDR)

# 计算 false discovery counts
total$falsePos = 1:nrow(total) * total$FDR

# 创建各组平均值
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])
total$baseMean = total$baseMeanA + total$baseMeanB

# # 将数字四舍五入，使其更美观
# total$foldChange = round(total$foldChange, 3)
# total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
# total$FDR = round(total$FDR, 4)
total$falsePos = round(total$falsePos, 0)
total$PAdj = formatC(total$PAdj, format = "e", digits = 3)
total$PValue = formatC(total$PValue, format = "e", digits = 3)

# 重新组织列名
new_cols = c("SYMBOL", 
             "baseMean","baseMeanA","baseMeanB",
             "foldChange", "log2FoldChange",
             "lfcSE","stat","PValue","PAdj", 
             "FDR","falsePos", col_names_A, col_names_B)
total = total[, new_cols]
```

## 保存结果
保存结果

``` r
write.csv(total, file=paste0(OUT_FOLDER, "/", OUTPUT_FILE), row.names = F, quote = F)
```

保存变量

``` r
save(total, file = "Rdata/DEG_DESeq2.Rdata") 
```
