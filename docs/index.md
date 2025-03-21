--- 
title: "RNA-Seq"
author:
  - 苏总华
date: "2025-03-22"
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

  RNA-Seq Pipeline
  
  我不是代码的创作者，我只是代码的搬运工。




RNA-Seq 数据分析从本质上来说就是对每一个基因在不同样本中的表达量进行比较。粗略来说经典的
RNA-Seq 数据分析流程包括以下几个步骤：

- 第一步：提取出样本中的 RNA，将其打碎，每一个碎片我们称为 Read。  
- 第二步：检测每一个 Read 的序列，将其比对到基因组上，得到每一个基因的 Reads 数目，通常呈现形式是一个行列表，每一列代表一个样本，每一行代表一个基因，我们称其为
Counts 矩阵。
- 第三步：对每一个基因的 Reads 数目进行归一化（因为不同样本间总 Reads 数目和每个基因长度不一致）。  
- 第四步：对同一个基因在不同样本中的表达量进行比较，找到差异表达基因。
- 第五步：分析差异表达基因的生物学意义。  


在 RNA-Seq 分析中，通常第一步和第二步是由公司完成的，我们拿到的数据是一个 Counts 矩阵。
所以我们只需要对后续进行分析和解读就可以了。
