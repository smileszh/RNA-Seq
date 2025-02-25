# 差异基因富集分析





``` r
# 加载必要的包
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

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

## GO 富集分析

``` r
GO <- function(DEG_exp){
  geneList <- DEG_exp$SYMBOL %>% unique()
  
  exp <- DEG_exp %>%
    filter(regulated == "up" | regulated == "down")
  gene <- exp$SYMBOL %>% unique()
  
  ego <- enrichGO(gene          = gene,
                  universe      = geneList,
                  OrgDb         = org_db,
                  keyType       = 'SYMBOL',
                  ont           = "all",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99,
                  readable      = TRUE)
  
  go_bar <- barplot(ego, split="ONTOLOGY", font.size = 10) + 
      facet_grid(ONTOLOGY~., scales="free") +
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  
  return(list(ego = ego, go_bar = go_bar))
}


# 对列表中的每个数据框进行处理
GO_DESeq2 <- setNames(lapply(cr_DESeq2, GO),
                      names(cr_DESeq2))
GO_edgeR <- setNames(lapply(cr_edgeR, GO),
                     names(cr_edgeR))

# 创建文件保存目录
dir.create("DEG_Analysis/3_DEG_ORA/DESeq2", recursive = T)
dir.create("DEG_Analysis/3_DEG_ORA/edgeR", recursive = T)

for (i in names(GO_DESeq2)) {
  ggsave(paste0("DEG_Analysis/3_DEG_ORA/DESeq2/", i, "_GO_DESeq2.pdf"),
            GO_DESeq2[[i]]$go_bar
         )
  write_csv(GO_DESeq2[[i]]$ego@result, paste0("DEG_Analysis/3_DEG_ORA/DESeq2/", i, "_GO_DESeq2.csv"))
  
}

for (i in names(GO_edgeR)) {
 ggsave(paste0("DEG_Analysis/3_DEG_ORA/edgeR/", i, "_GO_edgeR.pdf"), 
        GO_edgeR[[i]]$go_bar
        )
  write_csv(GO_edgeR[[i]]$ego@result, paste0("DEG_Analysis/3_DEG_ORA/edgeR/", i, "_GO_edgeR.csv"))
}

```


## KEGG 富集分析


``` r
KEGG <- function(DEG_exp){
  
  exp <- DEG_exp %>%
    filter(regulated == "up" | regulated == "down")
  gene <- exp$SYMBOL %>% unique()
  
  gene_df <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org_db)
  
  kk <- enrichKEGG(gene        = gene_df$ENTREZID,
                  organism     = kegg_prefix,,
                  pvalueCutoff = 0.999,
                  qvalueCutoff =0.999)
  kk <- DOSE::setReadable(kk, 
                          OrgDb = org_db,
                          keyType='ENTREZID') #按需替换
  
  kegg_bar <- barplot(kk, showCategory=20)
  
  return(list(kk = kk, kegg_bar = kegg_bar))
}


# 对列表中的每个数据框进行处理
KEGG_DESeq2 <- setNames(lapply(cr_DESeq2, KEGG),
                      names(cr_DESeq2))
KEGG_edgeR <- setNames(lapply(cr_edgeR, KEGG),
                     names(cr_edgeR))

# 创建文件保存目录
dir.create("DEG_Analysis/3_DEG_ORA/DESeq2", recursive = T)
dir.create("DEG_Analysis/3_DEG_ORA/edgeR", recursive = T)

for (i in names(KEGG_DESeq2)) {
  ggsave(paste0("DEG_Analysis/3_DEG_ORA/DESeq2/", i, "_KEGG_DESeq2.pdf"),
            KEGG_DESeq2[[i]]$kegg_bar
         )
  write_csv(KEGG_DESeq2[[i]]$kk@result, paste0("DEG_Analysis/3_DEG_ORA/DESeq2/", i, "_KEGG_DESeq2.csv"))
  
}

for (i in names(KEGG_edgeR)) {
 ggsave(paste0("DEG_Analysis/3_DEG_ORA/edgeR/", i, "_KEGG_edgeR.pdf"), 
        KEGG_edgeR[[i]]$kegg_bar
        )
  write_csv(KEGG_edgeR[[i]]$kk@result, paste0("DEG_Analysis/3_DEG_ORA/edgeR/", i, "_KEGG_edgeR.csv"))
}

```
