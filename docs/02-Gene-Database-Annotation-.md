# 基因 ID 注释





``` r
# 安装并加载必要的包
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

get_annotations <- function(gene_list, species = "human") {
  # 设置物种相关信息
  if (species == "human") {
    ensembl_dataset <- "hsapiens_gene_ensembl"
    org_db <- "org.Hs.eg.db"
    kegg_prefix <- "hsa"
  } else if (species == "mouse") {
    ensembl_dataset <- "mmusculus_gene_ensembl"
    org_db <- "org.Mm.eg.db"
    kegg_prefix <- "mmu"
  } else {
    stop("Unsupported species. Please use 'human' or 'mouse'.")
  }

# 使用biomaRt获取基因基本注释信息
ensembl <- useMart("ensembl", dataset = ensembl_dataset)
attributes <- c(
  'hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id',
  'description'
                )
gene_annotations <- getBM(attributes = attributes,
                            filters = 'hgnc_symbol', 
                            values = gene_list, 
                            mart = ensembl)
}


# 定义基因列表
gene_list <- data$SYMBOL

# 获取注释信息
annotations <- get_annotations(gene_list, species)

# 保存注释结果
dir.create("Expression_Annotation/1_Expression_result/", recursive = T)
writexl::write_xlsx(annotations, 
                   "Expression_Annotation/1_Expression_result/All_gene_annotations.xlsx")
```

