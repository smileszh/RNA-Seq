# RNA-Seq

## 分析前准备

如不想进行代码改动，则有以下要求：

**counts 矩阵（counts.csv）**

> - 第一列为基因 ENSEMBL ID，命名为 "name"。
> - 其余列为样本 counts 值，一列一个样本。
> - 样本从左到右依次为对照组和处理组。
> - 不要有其他列。
> - CSV 格式

<img src="./readme.assets/Screenshot 2025-03-22 at 18.32.24.png" alt="Screenshot 2025-03-22 at 18.32.24" style="zoom:20%;" />

**分组设计（design.csv）**

> - 一列为 ”sample“，一列为 “group”，可有其他列。
> - sample 列样本名和 counts 矩阵样本名一致。
> - Group 列从上到下依次为对照组和处理组
> - CSV 格式

<img src="./readme.assets/Screenshot 2025-03-22 at 18.34.02.png" alt="Screenshot 2025-03-22 at 18.34.02" style="zoom:20%;" />



## 章节概述

### 数据预处理

1. 将 counts 矩阵转换为 SYMBOL ID 为行名的矩阵。

###  样本表达分布

1. 查看是否有离群样本。
2. 查看组内重复性和组间差异性。

### 差异基因分析

1. 归一化counts 值，计算差异表达基因。

### 差异基因可视化

1. 差异基因可视化

### 富集分析

1. 差异表达 ORA 富集分析 （GO/KEGG）。
2. GSEA 分析



