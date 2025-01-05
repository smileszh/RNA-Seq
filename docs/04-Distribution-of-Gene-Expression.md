## 基因表达分布





``` r
# 下载并加载必要的包

library(tidyverse)
library(reshape2)
fpkm_df <- fpkm[rowSums(fpkm >1) > 1,]
fpkm_long <- log10(fpkm_df+1) %>%
  melt(variable.name = "variable", value.name = "value")
```

### 箱线图

``` r
ALL_box <- ggplot(fpkm_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=13, angle = 45, vjust = 1, hjust=1, color = "black"),
    axis.text.y = element_text(size=13, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) +
  labs(x = "Samples", y = "log2(FPKM+1)")

dir.create("Expression_Annotation/2_Expression_distribution")
ggsave(plot = ALL_box, filename = "./Expression_Annotation/2_Expression_distribution/all.fpkm_box.pdf", width = 7, height = 4)
```

### 密度图

``` r
ALL_density <- ggplot(fpkm_long) +
  geom_density(aes(x = value, fill = variable), alpha = 0.6) +
  scale_fill_brewer(palette = "Set3") +  # 使用颜色调色板
  theme_bw() +  # 使用简洁主题
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中且加粗
    axis.text = element_text(size = 12),  # 轴文字字体大小
    legend.position = "right",  # 图例位置
    legend.title = element_text(size = 14),  # 图例标题字体大小
    legend.text = element_text(size = 12), # 图例文字字体大小
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()) +   # 去掉次要网格线
  scale_x_continuous(limits = c(-1, NA)) +  # 设置X轴最小值为-5，最大值自动
  labs(x = "log2(FPKM+1)", y = "density", fill = "Sample")  # 修改图例标题

dir.create("Expression_Annotation/2_Expression_distribution")
ggsave(plot = ALL_density, filename = "./Expression_Annotation/2_Expression_distribution/all.fpkm_density.pdf", width = 6, height = 4)
```


``` r
# 获取样本名列表
sample_names <- unique(fpkm_long$variable)

# 循环生成每个样本的密度图并保存
for (sample in sample_names) {
  # 过滤数据，仅保留当前样本的数据
  sample_data <- fpkm_long[fpkm_long$variable == sample, ]
  
  # 生成密度图
  sample_density <- ggplot(sample_data) +
    geom_density(aes(x = value, fill = variable), alpha = 0.6) +
    scale_fill_brewer(palette = "Set3") +  # 使用颜色调色板
    theme_bw() +  # 使用简洁主题
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中且加粗
      axis.text = element_text(size = 12),  # 轴文字字体大小
      legend.position = "right",  # 图例位置
      legend.title = element_text(size = 14),  # 图例标题字体大小
      legend.text = element_text(size = 12), # 图例文字字体大小
      panel.grid.major = element_blank(),  # 去掉主要网格线
      panel.grid.minor = element_blank()) +   # 去掉次要网格线
    scale_x_continuous(limits = c(-1, NA)) +  # 设置X轴最小值为-5，最大值自动
    labs(title = paste("Density Plot of", sample),
         x = "log2(FPKM+1)", 
         y = "density", 
         fill = "Sample")  # 修改图例标题
  
  # 保存密度图
  ggsave(plot = sample_density, 
         filename = paste0("./Expression_Annotation/2_Expression_distribution/", sample, "_fpkm_density.pdf"), 
         width = 6, 
         height = 4)
}
```




``` r
save(species, count, fpkm, fpkm_long, data, group, group_compare, file = "Rdata/data.Rdata")
```
