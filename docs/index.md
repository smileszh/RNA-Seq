--- 
title: "RNA-Seq"
author:
  - 小苏
date: "2025-03-25"
documentclass: ctexbook
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
colorlinks: yes
lot: yes
lof: yes
geometry: [b5paper, tmargin=2.5cm, bmargin=2.5cm, lmargin=3.5cm, rmargin=2.5cm]
site: bookdown::bookdown_site

description: "This is an online note of the scRNA-Seq analysis."
github-repo: smileszh/RNA-Seq
always_allow_html: true
#cover-image: images/*.jpg

---





# 前言 {-}

*我不是代码的创作者，我只是代码的搬运工！*

RNA-Seq 数据分析的核心在于比较不同样本中各基因的表达量。经典的 RNA-Seq 数据分析流程可以概括为以下几个步骤：

**RNA 提取与片段化**：从样本中提取 RNA，并将其打碎成片段，这些片段被称为 reads。

**序列比对与计数**：检测每个 Read 的序列，将其比对到基因组上，从而获得每个基因的 Reads 数目。结果通常以 Counts 矩阵的形式呈现，其中每一列代表一个样本，每一行代表一个基因。

**归一化**：对每个基因的 Reads 数目进行归一化处理，以校正不同样本间总 reads 数目和基因长度的差异。

**差异表达分析**：比较同一基因在不同样本中的表达量，识别出差异表达基因。

**生物学意义分析**：进一步分析差异表达基因的生物学意义，揭示其在生物过程中的作用。

在实际操作中，通常前两个步骤由专业公司完成，研究者拿到的是一个 counts 矩阵。因此，研究者主要负责后续的数据分析和结果解读。
