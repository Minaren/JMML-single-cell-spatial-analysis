HSC_deg<- read.csv("D:/aJMML/scRNA1/HSC.csv",header = TRUE,,sep = ",")
head(HSC_deg)
class(HSC_deg)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
row.names(HSC_deg) <- HSC_deg[, 1]
EnhancedVolcano(HSC_deg,
                lab = rownames(HSC_deg),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'HSC_deg')
ggsave("HSC_deg.pdf",width = 9, height = 8)



#MPP差异基因
MPP_deg<- read.csv("D:/aJMML/scRNA1/MPP.csv",header = TRUE,,sep = ",")
head(MPP_deg)
class(MPP_deg)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
row.names(MPP_deg) <- MPP_deg[, 1]
EnhancedVolcano(MPP_deg,
                lab = rownames(MPP_deg),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'MPP_deg')
ggsave("MPP_deg.pdf",width = 9, height = 8)
