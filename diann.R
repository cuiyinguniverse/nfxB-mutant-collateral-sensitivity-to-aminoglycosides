#安装diann
install.packages("devtools")
library(usethis)
library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)

getwd()
setwd("G:/D2_20231117134316")
getwd()
df <- diann_load("./result.tsv")#输入搜库后的结果文件.tsv
precursors <- diann_matrix(df, pg.q = 0.01)
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = 0.01)
peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Stripped.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
unique.genes <- diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01)
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

#给第一列加上列名description
protein.groups
locus <- rownames(protein.groups)
locus
protein.groups <- cbind(locus,protein.groups)
#将列名改成1-dim[],便于后续生成列表
protein.groups
dim(protein.groups)
dim(protein.groups)[1]
rownames(protein.groups) <- seq(1:dim(protein.groups)[1])
#生成结果
write.table(protein.groups,file = "./result.csv", sep =",", row.names =F, col.names =TRUE, quote =TRUE)
