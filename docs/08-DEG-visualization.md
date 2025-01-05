## 差异基因可视化


### 火山图




``` r
# 加载必要的包
library(tidyverse)
library(ggrepel)



# 定义画图函数
 vol <- function(data){
   
   # 统计上下调基因的数量
   regulated_counts <- table(data$regulated)
   up_count <- regulated_counts["up"]
   down_count <- regulated_counts["down"]
   
   # 画图
   p <- ggplot(data, aes(x = log2FoldChange, y = -log10(FDR))) +
    geom_point(aes(color = regulated), size = 1, alpha = 1) +
    scale_color_manual(values = c('up' = '#F07C79', 'down' = '#4588C8', 'normal' = 'gray')) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size=1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      axis.text = element_text(size = 16, color = "black", face = "bold"),
      panel.border = element_rect(color = "black", size = 1.5)
    )
      
   # 计算图的范围
   x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range
   y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range

   # 定义注释位置
   x_left <- x_range[1] + (x_range[2] - x_range[1]) * 0.05
   x_right <- x_range[2] - (x_range[2] - x_range[1]) * 0.05
   y_top <- y_range[2] - (y_range[2] - y_range[1]) * 0.05
  
   # 添加注释
  vol_p <- p +
    annotate("text", x = x_right, y = y_top, label = paste0(up_count), size = 7, hjust = 1, 
           color = '#F07C79') +
    annotate("text", x = x_left, y = y_top, label = paste0(down_count), size = 7, hjust = 0,
           color = '#4588C8')
      
  return(vol_p)
   }

# 对列表中的每个数据框进行处理
vol_DESeq2 <- setNames(lapply(cr_DESeq2, vol),
                      names(cr_DESeq2))
vol_edgeR <- setNames(lapply(cr_edgeR, vol),
                     names(cr_edgeR))

# 创建图像保存目录
dir.create("DEG_Analysis/2_DEG_visualization/DESeq2", recursive = T)
dir.create("DEG_Analysis/2_DEG_visualization/edgeR", recursive = T)

for (i in names(vol_DESeq2)) {
  ggsave(paste0("DEG_Analysis/2_DEG_visualization/DESeq2/", i, "_volcano_DESeq2.pdf"),
            vol_DESeq2[[i]]
         )
}

for (i in names(vol_edgeR)) {
 ggsave(paste0("DEG_Analysis/2_DEG_visualization/edgeR/", i, "_volcano_edgeR.pdf"), 
        vol_edgeR[[i]]
        )
}

```

### 热图




``` r
# 加载必要的包
library(pheatmap)
library(tidyverse)

# 定义画图函数
 phe <- function(DEG_exp){
   exp <- DEG_exp %>%
     filter(regulated == "up" | regulated == "down") %>% 
     column_to_rownames(var = "SYMBOL") %>%
     select(any_of(group$sample))
   
   phe_p <- 
     pheatmap(exp,
         show_colnames =T,
         show_rownames = F,
         scale = "row",
         cluster_cols = T
         )
   return(phe_p)
 }
 
 # 对列表中的每个数据框进行处理
phe_DESeq2 <- setNames(lapply(cr_DESeq2, phe),
                      names(cr_DESeq2))
```

<img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-1.png" width="672" /><img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-2.png" width="672" /><img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-3.png" width="672" />

``` r
phe_edgeR <- setNames(lapply(cr_edgeR, phe),
                     names(cr_edgeR))
```

<img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-4.png" width="672" /><img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-5.png" width="672" /><img src="08-DEG-visualization_files/figure-html/unnamed-chunk-4-6.png" width="672" />

``` r
 
# 创建图像保存目录
dir.create("DEG_Analysis/2_DEG_visualization/DESeq2", recursive = T)
dir.create("DEG_Analysis/2_DEG_visualization/edgeR", recursive = T)

for (i in names(phe_DESeq2)) {
  ggsave(paste0("DEG_Analysis/2_DEG_visualization/DESeq2/", i, "_pheatmap_DESeq2.pdf"),
            phe_DESeq2[[i]]
         )
}

for (i in names(phe_edgeR)) {
 ggsave(paste0("DEG_Analysis/2_DEG_visualization/edgeR/", i, "_pheatmap_edgeR.pdf"), 
        phe_edgeR[[i]]
        )
}
```
