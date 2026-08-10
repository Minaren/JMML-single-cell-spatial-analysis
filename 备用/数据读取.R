 ####数据读取方式####
 #1.cell ranger生成的raw count(包含barcode,gene,matrix)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
#自动读取cellranger(LINUX)输出的feature barcode matric
rm(list = ls())
pbmc.data <- Read10X(data.dir = "JMML_CB35") 
#load the PBMC datasets
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,min.features=200
)
#initialize the seurat object with the raw(non-normalized data)


#2.仅有一个稀疏矩阵时的读取方法:Count matrix导入
matrix_data <- read.table("single_cell_datamatrix.txt", sep="\t", header=T, row.names=1)
dim(matrix_data)
seurat_obj <- CreateSeuratObject(counts = matrix_data)

#读取RDS文件
rm(list = ls())
pbmc <- readRDS("panc8.rds")
saveRDS(pbmc,"pbmc.rds")

#seurat功能思维导图
library(tidyverse)
str(pbmc)
library(mindr)
(out <- capture.output(str(pbmc)))
out2 <- paste(out, collapse="\n")
mm(gsub("\\.\\.@","# ",gsub("\\.\\. ","#",out2)),type ="text",root= "Seurat")