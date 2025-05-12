getwd()
setwd("E:/")
getwd()

install.packages('BiocManager') 
BiocManager::install('clusterProfiler')
library(clusterProfiler)
library(eoffice)

go_bg <- read.csv("./gene_ontology_csv.csv")
gene_change <- read.csv("./nfxB_Difference_change.csv")

genes <- gene_change$locus

go_rich <- enricher(gene = genes,  #待富集的基因列表
                    TERM2GENE = go_bg[c('ID', 'locus')],  #背景基因集
                    TERM2NAME = go_bg[c('ID', 'GO_Term')], 
                    pAdjustMethod = 'BH',  #指定 p 值校正方法
                    pvalueCutoff = 0.05,  #指定 p 值阈值 
                    qvalueCutoff = 0.05)  #指定 q 值阈值

write.csv(go_rich,"./nfxB_go_rich.csv")

#再把 GO Ontology 信息添加在上述 GO 富集结果中
tmp <- read.csv('./nfxB_go_rich.csv')
gene_GO <- go_bg[!duplicated(go_bg$ID), ]
tmp2 <- merge(tmp, gene_GO[c('ID', 'ONTOLOGY')], by = 'ID')
tmp3 <- tmp2[c(10, 1:9)]
tmp4 <- tmp3[order(tmp$pvalue), ]
write.csv(tmp4,"./nfxB_GO_cluster_anno.csv")

#绘图
b1 <- barplot(go_rich)  #富集柱形图
topptx(b1,filename = "nfxB_Gorich1.pptx",width = 6,height = 4)
d1 <- dotplot(go_rich)  #富集气泡图
topptx(d1,filename = "nfxB_Gorich2.pptx",width = 6,height = 4)



###############################################################################
KEGG <- read.csv("./pathways.csv")
gene_change <- read.csv("./nfxB_Difference_change.csv")

genes <- gene_change$locus

kegg_rich <- enricher(gene = genes,  #待富集的基因列表
                      TERM2GENE = KEGG[c('Pathway', 'locus')],  #背景基因集
                      TERM2NAME = KEGG[c('Pathway', 'Description')], 
                      pAdjustMethod = 'BH',  #指定 p 值校正方法
                      pvalueCutoff = 0.5,  #指定 p 值阈值
                      qvalueCutoff = 0.5)  #指定 q 值阈值
#绘图
d2 <- dotplot(kegg_rich)

topptx(d2,filename = "nfxB_keggrich.pptx",width = 6,height = 4)
