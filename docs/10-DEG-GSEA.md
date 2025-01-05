## GSEA 分析





``` r
# 加载必要的包
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# 确定物种
if (species == "human") {
    org_db <- "org.Hs.eg.db"
    kegg_prefix <- "hsa"
  } else if (species == "mouse") {
    org_db <- "org.Mm.eg.db"
    kegg_prefix <- "mmu"
  } else {
    stop("Unsupported species. Please use 'human' or 'mouse'.")
  }
```

## KEGG数据库

``` r

# 定义GSEA函数
GSE_kk <- function(DEG_exp) {
  gene_df <- bitr(DEG_exp$SYMBOL, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org_db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
    
  DEG_exp <- left_join(DEG_exp, gene_df, by = "SYMBOL") %>%
    na.omit()
  
  geneList <- DEG_exp$log2FoldChange
  names(geneList) <- DEG_exp$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  gse_kk <- gseKEGG(geneList = geneList,
                    organism = kegg_prefix,
                    minGSSize = 10,
                    pvalueCutoff = 0.99,
                    verbose = FALSE)
  
  gse_kk <- DOSE::setReadable(gse_kk, OrgDb = org_db, keyType = 'ENTREZID')
  
  # 可视化
  up_k <- gse_kk[head(order(gse_kk$enrichmentScore, decreasing = TRUE)), ]
  up_k$group <- 1
  down_k <- gse_kk[tail(order(gse_kk$enrichmentScore, decreasing = TRUE)), ]
  down_k$group <- -1

  dat <- rbind(up_k, down_k)
  dat$pvalue <- -log10(dat$pvalue)
  dat$pvalue <- dat$pvalue * dat$group 
  dat <- dat[order(dat$pvalue, decreasing = FALSE), ]

  gse_p <- ggplot(dat, aes(x = reorder(Description, order(pvalue, decreasing = FALSE)),
                           y = pvalue, fill = group)) + 
    geom_bar(stat = "identity") + 
    scale_fill_gradient(low = "#34bfb5", high = "#ff6633", guide = FALSE) + 
    scale_x_discrete(name = "Pathway names") +
    scale_y_continuous(name = "log10P-value") +
    coord_flip() + 
    ggstatsplot::theme_ggstatsplot() +
    theme(plot.title = element_text(size = 15, hjust = 0.5),  
          axis.text = element_text(size = 12, face = 'bold'),
          panel.grid = element_blank()) +
    ggtitle("Pathway Enrichment")
  
  return(list(gse_kk = gse_kk, gse_p = gse_p))
}

# 对列表中的每个数据框进行处理
GSE_kk_DESeq2 <- setNames(lapply(cr_DESeq2, GSE_kk), names(cr_DESeq2))
GSE_kk_edgeR <- setNames(lapply(cr_edgeR, GSE_kk), names(cr_edgeR))

# 创建文件保存目录
dir.create("DEG_Analysis/4_GSEA/DESeq2", recursive = TRUE)
dir.create("DEG_Analysis/4_GSEA/edgeR", recursive = TRUE)

# 保存结果
for (i in names(GSE_kk_DESeq2)) {
  ggsave(paste0("DEG_Analysis/4_GSEA/DESeq2/", i, "_GSE_kk_DESeq2.pdf"), 
         GSE_kk_DESeq2[[i]]$gse_p)
  write_csv(GSE_kk_DESeq2[[i]]$gse_kk@result, 
            paste0("DEG_Analysis/4_GSEA/DESeq2/", i, "_GSE_kk_DESeq2.csv"))
}

for (i in names(GSE_kk_edgeR)) {
  ggsave(paste0("DEG_Analysis/4_GSEA/edgeR/", i, "_GSE_kk_edgeR.pdf"), 
         GSE_kk_edgeR[[i]]$gse_p)
  write_csv(GSE_kk_edgeR[[i]]$gse_kk@result, 
            paste0("DEG_Analysis/4_GSEA/edgeR/", i, "_GSE_kk_edgeR.csv"))
}

```
