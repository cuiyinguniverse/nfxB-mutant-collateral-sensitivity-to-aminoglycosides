if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

setwd("F:/work")

library(grid)
library(ComplexHeatmap)
rm(list=ls())


df <- read.csv("./other_efflux_pump_spectrum.csv")
head(df)
substrate <- as.matrix(df) 
#head(substrate)
rownames(substrate) <- c('A_baumannii_AdeABC','Enterobacteriaceae_AcrAB_TolC','Enterobacteriaceae_OqxAB_TolC','E_coli_AcrAB-TolC','K_pneumoniae_AcrAB-TolC','K_pneumoniae_AcrAB-KocC','K_pneumoniae_KexD','K_pneumoniae_OqxAB') 
#head(substrate)
#Heatmap(substrate)
substrate <- t(substrate)


tiff("other_RND_efflux_pumps_substrate_spectrum.tiff", width = 6, height = 2.2, units = "in", res = 300)  # 设置宽高为英寸，分辨率为300 DPI
pdf("other_RND_efflux_pumps_substrate_spectrum.pdf", width = 6, height = 2.2)
svg("other_RND_efflux_pumps_substrate_spectrum.svg", width = 6, height = 2.2)


Heatmap(substrate,
        name = "Efflux-pump action",
        #cell_fun = cell_fun,
        col = c('#9ecae1',"#d73027"),
        show_row_names = TRUE, 
        row_names_side = "left",
        show_column_names = TRUE, 
        column_names_side = "top",
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE, 
        width = ncol(substrate)*unit(10, "mm"), 
        height = nrow(substrate)*unit(10, "mm"),
        column_split = rep(c("A","B","C",'D','a','b','c','d'),each=1),
        row_split = rep(c("A","B","C","D",'E','F','G','a','b','c','d','e','f','g'),each=1), 
        row_title = NULL,
        column_title = NULL,
        column_names_rot = 60,
        column_names_gp = gpar(fontface = "italic"),
        #column_names_gp = grid::gpar(fontsize = 12),
        #row_names_gp = grid::gpar(fontsize = 10),
        heatmap_legend_param = list(
            at = c(1,0),
            labels = c("Excreted", "Non-excreted"),
            title = "Efflux-pump action", 
            legend_gp = gpar(fill = 1:2))
        )


#MIC of delta_nfxB_mexC_mexD_oprJ in Tob or Gen or Cip

df1 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(8,4,0.25),
  delta_n = c(2,1,4),
  delta_nC = c(8,4,0.25),
  delta_nD = c(8,4,0.25),
  delta_nJ = c(2,1,1),
  delta_nCJ = c(8,4,0.25),
  delta_nDJ = c(8,4,0.25)
)


df1 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(8,4,0.25),
  delta_n = c(2,1,4)
)


rownames(df1) <- df1$Antibiotics

df1 <- df1[ , -1]

df_log2 <- df1

df_log2[,] <- lapply(df_log2, log2)

df1 <- df_log2

df_adjusted <- df1
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df1 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  expression(Delta * "nfxB" * Delta * "mexC"),
  expression(Delta * "nfxB" * Delta * "mexD"),
  expression(Delta * "nfxB" * Delta * "oprJ"),
  expression(Delta * "nfxB" * Delta * "mexC" * Delta * "oprJ"),
  expression(Delta * "nfxB" * Delta * "mexD" * Delta * "oprJ")
)

col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df1[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df1, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(15, "mm"), 
        height = nrow(df1)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B"),each=1),
        heatmap_legend_param = list(
          at = c(-2,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))


#MIC of delta_nfxB_mexC_mexD_oprJ in Tob or Gen or Cip

################################################################################
#MG1655_over expression mexCD-oprJ

df2 <- data.frame(
  Antibiotics = c("Aztreonam","Colistin","Ampicillin","Trimethoprim","Carbenicillin","Chloramphenicol","Ciprofloxacin","Ceftazidime","Polymixin B","Nitrofurantoin","Piperacillin tazobactam","Meropenem","Gentamicin","Tobramycin"),
  MG1655 = c(16,0.25,8,0.5,16,16,0.0625,0.5,0.25,32,4,0.0625,8,8),
  empty_vector = c(16,0.125,8,0.5,16,8,0.0625,0.25,0.125,16,2,0.0625,4,4),
  mexCD_oprJ = c(1,0.0625,8,0.5,16,8,0.25,0.5,0.125,16,2,0.0625,4,4)
)

rownames(df2) <- df2$Antibiotics

df2 <- df2[ , -1]

df_log2 <- df2

df_log2[,] <- lapply(df_log2, log2)

df2 <- df_log2

df_adjusted <- df2
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df2 <- as.matrix(df_adjusted)
df2 <- t(df2)
#col_labels <- c(
#  "PAO1",
#  expression(Delta * "nfxB"),
#  expression(Delta * "nfxB" * Delta * "mexC"),
#  expression(Delta * "nfxB" * Delta * "mexD"),
#  expression(Delta * "nfxB" * Delta * "oprJ"),
#  expression(Delta * "nfxB" * Delta * "mexC" * Delta * "oprJ"),
#  expression(Delta * "nfxB" * Delta * "mexD" * Delta * "oprJ")
#)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df2[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df2, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,
        row_names_gp = gpar(fontface = "italic", rot = 45),
        row_title = "Escherichia coli strains",
        row_title_gp = gpar(fontface = "italic"),
        row_title_side = "left",
        column_title = "Antibiotics",
        #column_title_gp = gpar(fontface = "italic"),
        column_title_side = "bottom",
        column_names_rot = 60,
        width = ncol(df2)*unit(8, "mm"), 
        height = nrow(df2)*unit(10, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C",'D','E','F','G','a','b','c','d','e','f','g'),each=1),
        heatmap_legend_param = list(
          at = c(-4,0,2),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))
#MIC of MG1655_over expression mexCD-oprJ in many antibiotics



###########################################################################

#Pf-5 over expression mexCD-oprJ


df3 <- data.frame(
  Antibiotics = c("Aztreonam","Meropenem","Piperacillin tazobactam","Nitrofurantoin","Amikacin","Ceftazidime","Carbenicillin","Trimethoprim","Ampicillin","Polymixin B","Colistin","Ciprofloxacin","Gentamicin","Tobramycin"),
  Pf_5 = c(128,2,1,0.03125,0.5,16,0.03125,128,256,2,16,0.0009765625,0.5,0.5),
  empty_vector = c(128,2,0.5,0.03125,0.5,2,0.03125,128,128,2,2,0.0009765625,0.25,0.25),
  mexCD_oprJ = c(32,2,0.5,0.03125,0.5,2,0.03125,128,128,2,2,0.00390625,0.25,0.25)
)

rownames(df3) <- df3$Antibiotics

df3 <- df3[ , -1]

df_log2 <- df3

df_log2[,] <- lapply(df_log2, log2)

df3 <- df_log2

df_adjusted <- df3
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df3 <- as.matrix(df_adjusted)
df3 <- t(df3)
#col_labels <- c(
#  "PAO1",
#  expression(Delta * "nfxB"),
#  expression(Delta * "nfxB" * Delta * "mexC"),
#  expression(Delta * "nfxB" * Delta * "mexD"),
#  expression(Delta * "nfxB" * Delta * "oprJ"),
#  expression(Delta * "nfxB" * Delta * "mexC" * Delta * "oprJ"),
#  expression(Delta * "nfxB" * Delta * "mexD" * Delta * "oprJ")
#)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df3[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df3, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,
        row_names_gp = gpar(fontface = "italic", rot = 45),
        row_title = "Pseudomonas protegens\nstrains",
        row_title_gp = gpar(fontface = "italic"),
        row_title_side = "left",
        column_title = "Antibiotics",
        #column_title_gp = gpar(fontface = "italic"),
        column_title_side = "bottom",
        column_names_rot = 60,
        width = ncol(df3)*unit(8, "mm"), 
        height = nrow(df3)*unit(10, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C",'D','E','F','G','a','b','c','d','e','f','g'),each=1),
        heatmap_legend_param = list(
          at = c(-3,0,2),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of Pf-5 over expression mexCD-oprJ in many antibiotics

##############################################################################




df4 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(0,0,0),
  delta_n = c(-2,-2,3),
  empty_vector = c(-1,-1,0),
  mexC = c(-1,-1,0),
  mexD = c(-3,-2,0),
  oprJ = c(-1,-1,0),
  mexCD = c(-3,-2,2),
  mexC_oprJ = c(-2,-1,-1),
  mexD_oprJ = c(-2,-1,0),
  mexCD_oprJ = c(-2,-1,2)
)

rownames(df4) <- df4$Antibiotics

df4 <- df4[ , -1]

#

df_adjusted <- df4
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df4 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  "empty_vector",
  "mexC",
  "mexD",
  "oprJ",
  "mexCD",
  "mexC_oprJ", 
  "mexD_oprJ", 
  "mexCD_oprJ")



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df4[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df4, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(15, "mm"), 
        height = nrow(df1)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C",'D','E','F','G','H','I','J'),each=1),
        heatmap_legend_param = list(
          at = c(-3,0,3),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of PAO1 over expression mexCD_oprJ in Gen or Tob or Cip 
###################################################################



df5 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(0,0,0),
  empty_vector = c(-1,-1,0),
  mexD = c(-3,-2,0),
  mexD_N256A = c(-2,-2,0),
  mexD_D258A = c(-1,0,0),
  mexD_S260A = c(-2,-2,0),
  mexD_D656A = c(-1,-1,0),
  mexD_D657A = c(-1,-1,0)
)

rownames(df5) <- df5$Antibiotics

df5 <- df5[ , -1]

#

df_adjusted <- df5
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df5 <- as.matrix(df_adjusted)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df5[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df5, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(15, "mm"), 
        height = nrow(df1)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C",'D','E','F','G','H'),each=1),
        heatmap_legend_param = list(
          at = c(-3,0),
          labels = c("Susceptible",
                     "Intermediate"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of mexD point mutation in MexCD Interaction interface
##################################################



df6 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin"),
  PAO1 = c(16,16),
  delta_n = c(4,4),
  delta_x = c(8,8),
  delta_y = c(4,8),
  delta_xy = c(8,8),
  delta_nx = c(0.5,1),
  delta_ny = c(1,2),
  delta_nxy = c(2,2)
)

rownames(df6) <- df6$Antibiotics

df6 <- df6[ , -1]

df_log2 <- df6

df_log2[,] <- lapply(df_log2, log2)

df6 <- df_log2

df_adjusted <- df6
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df6 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  expression(Delta * "mexX"),
  expression(Delta * "mexY"),
  expression(Delta * "mexXY"),
  expression(Delta * "nfxB" * Delta * "mexX"),
  expression(Delta * "nfxB" * Delta * "mexY"),
  expression(Delta * "nfxB" * Delta * "mexXY")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df6[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df6, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(15, "mm"), 
        height = nrow(df1)*unit(15, "mm"),
        row_split = rep(c("A","B"),each=1), 
        column_split = rep(c("A","B","C",'D','a','b','c','d'),each=1),
        heatmap_legend_param = list(
          at = c(-5,0),
          labels = c("Susceptible",
                     "Intermediate"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

# MIC of delta_nfxB_mexX_mexY in PAO1
################################################################################



df7 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin"),
  PAO1 = c(8,8),
  delta_n = c(2,2),
  delta_na = c(2,2),
  delta_nb = c(2,2),
  delta_nab = c(2,2),
  delta_nxy = c(1,1),
  delta_nxyab = c(1,0.5)
)

rownames(df7) <- df7$Antibiotics

df7 <- df7[ , -1]

df_log2 <- df7

df_log2[,] <- lapply(df_log2, log2)

df7 <- df_log2

df_adjusted <- df7
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df7 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  expression(Delta * "nfxB" * Delta * "mexA"),
  expression(Delta * "nfxB" * Delta * "mexB"),
  expression(Delta * "nfxB" * Delta * "mexAB"),
  expression(Delta * "nfxB" * Delta * "mexXY"),
  expression(Delta * "nfxB" * Delta * "mexAB" * Delta * "mexXY")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df7[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df7, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(15, "mm"), 
        height = nrow(df1)*unit(15, "mm"),
        row_split = rep(c("A","B"),each=1), 
        column_split = rep(c("A","B","C",'D','a','b','c'),each=1),
        heatmap_legend_param = list(
          at = c(-4,0),
          labels = c("Susceptible",
                     "Intermediate"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of delta_nfxB_mexAB in PAO1

#################################################################################








rm(list = ls())

getwd()


########################################################

ATLAS<-readxl::read_xlsx("Atlas_Reuse_Data.xlsx")

Antibiotic_fam<-read.delim('antibiotic_groups.txt')

ATLAS<-ATLAS[5:dim(ATLAS)[1],]

colnames(ATLAS)<-ATLAS[1,]

ATLAS<-ATLAS[2:dim(ATLAS)[1],]

ATLAS1<-ATLAS %>% 
  pivot_longer(names_to = "Variable", values_to = "Value", -c(`Isolate Id`:Phenotype)) %>% 
  filter(!is.na(Value)) %>% 
  mutate(Species=case_when( Species== "Enterobacter aerogenes" ~ "Klebsiella aerogenes", 
                            Species == "Klebsiella (Enterobacter) aerogenes" ~ "Klebsiella aerogenes",
                            TRUE ~ Species))

Key_MICs<-read.table("key_MIC.txt")
colnames(Key_MICs)<- c("Value", "MICs_transformed", "Log2_MIC") 
ATLAS2<-ATLAS1 %>% left_join(Key_MICs)

ATLAS2<-ATLAS2[grepl("_I", ATLAS2$Variable)==F,]


####################################################################


ATLAS3<- ATLAS2 %>%
  left_join(Antibiotic_fam, by=c("Variable"="Antibiotic")) %>%
  mutate(Variable=Short) %>%
  rename(Isolate=`Isolate Id`) %>%
  filter(Species=="Pseudomonas aeruginosa") %>% 
  group_by(Isolate, Variable) %>%
  summarise(Log2_MIC=max(Log2_MIC, na.rm=T)) %>% 
  pivot_wider(names_from = Variable, values_from = Log2_MIC, values_fill = NA) %>%
  column_to_rownames("Isolate")  

#################################
levofloxacin_R_amikacin_S




ATLAS_Levofloxacin_R<- ATLAS3[grepl("^[3-9]*$", ATLAS3$LEV)==T,]

#ATLAS_Amikacin_S<- ATLAS3[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS3$AMK)==T,]

ATLAS_Levofloxacin_R_Amikacin_S<- ATLAS_Levofloxacin_R[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS_Levofloxacin_R$AMK)==T,]
ATLAS_Levofloxacin_R_Amikacin_S$LEV<- ATLAS_Levofloxacin_R_Amikacin_S$LEV - 2  
ATLAS_Levofloxacin_R_Amikacin_S$AMK<- ATLAS_Levofloxacin_R_Amikacin_S$AMK - 4 
ATLAS_Levofloxacin_R_Amikacin_S$AZM<- ATLAS_Levofloxacin_R_Amikacin_S$AZM - 3
ATLAS_Levofloxacin_R_Amikacin_S$CEF<- ATLAS_Levofloxacin_R_Amikacin_S$CEF - 3
ATLAS_Levofloxacin_R_Amikacin_S$CAZ<- ATLAS_Levofloxacin_R_Amikacin_S$CAZ - 3
ATLAS_Levofloxacin_R_Amikacin_S$`CAZ+AVI`<- ATLAS_Levofloxacin_R_Amikacin_S$`CAZ+AVI` - 3
ATLAS_Levofloxacin_R_Amikacin_S$`CTO+TAZ`<- ATLAS_Levofloxacin_R_Amikacin_S$`CTO+TAZ` - 2
ATLAS_Levofloxacin_R_Amikacin_S$IMI<- ATLAS_Levofloxacin_R_Amikacin_S$IMI - 1
ATLAS_Levofloxacin_R_Amikacin_S$MER<- ATLAS_Levofloxacin_R_Amikacin_S$MER - 1
ATLAS_Levofloxacin_R_Amikacin_S$DOR<- ATLAS_Levofloxacin_R_Amikacin_S$DOR - 1
ATLAS_Levofloxacin_R_Amikacin_S$`PIP+TAZ`<- ATLAS_Levofloxacin_R_Amikacin_S$`PIP+TAZ` - 4
ATLAS_Levofloxacin_R_Amikacin_S$COL<- ATLAS_Levofloxacin_R_Amikacin_S$COL - 1
ATLAS_Levofloxacin_R_Amikacin_S$MIN<- ATLAS_Levofloxacin_R_Amikacin_S$MIN - 3
ATLAS_Levofloxacin_R_Amikacin_S$TIG<- ATLAS_Levofloxacin_R_Amikacin_S$TIG - 3
ATLAS_Levofloxacin_R_Amikacin_S$`COLP80 `<- ATLAS_Levofloxacin_R_Amikacin_S$`COLP80 ` - 1


new_ATLAS_Levofloxacin_R_Amikacin_S <- select(ATLAS_Levofloxacin_R_Amikacin_S,-c(AMC,AMP,CTX,`AZM+AVI`,CAR,`CAR+AVI`,ERT,SUL,AZI))




if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(grid)
library(ComplexHeatmap)

n_new_ATLAS_Levofloxacin_R_Amikacin_S <- mutate_all(new_ATLAS_Levofloxacin_R_Amikacin_S, ~replace(., is.na(.), 0))
mat_new_ATLAS_Levofloxacin_R_Amikacin_S <- as.matrix(n_new_ATLAS_Levofloxacin_R_Amikacin_S) 
colnames(mat_new_ATLAS_Levofloxacin_R_Amikacin_S) <- c("Amikacin","Ceftazidime","Cefepime","Levofloxacin","Meropenem","Minocycline","Piperacillin tazobactam","Tigecycline","Aztreonam","Ceftazidime avibactam","ColistinP80","Doripenem","Imipenem","Colistin","Ceftolozane tazobactam") 
mat_new_ATLAS_Levofloxacin_R_Amikacin_S_reordered <- mat_new_ATLAS_Levofloxacin_R_Amikacin_S[,c("Amikacin","Aztreonam","Ceftazidime","Cefepime","Piperacillin tazobactam","Ceftazidime avibactam","Ceftolozane tazobactam","Meropenem","Doripenem","Imipenem","Colistin","ColistinP80","Minocycline","Tigecycline","Levofloxacin")] 


data <- mat_new_ATLAS_Levofloxacin_R_Amikacin_S_reordered


primary_groups <- c(rep("Aminoglucosides", 1), rep("Beta-lactams", 9),rep("Polymixin", 2),rep("Tetracycline", 2),rep("Quinolones", 1))

secondary_groups <- c(rep("Aminoglucosides", 1), rep("Monobactams", 1), rep("Cephalosporins", 2),rep("Beta-lactam+inhibitor", 3),rep("Carbapenems", 3),rep("Polymixin", 2),rep("Tetracycline", 2),rep("Quinolones", 1))



primary_annotation <- HeatmapAnnotation(
  PrimaryGroup = primary_groups,
  col = list(PrimaryGroup = c("Aminoglucosides" = "#0000CD","Beta-lactams" = "#32CD32", "Polymixin" = "#20B2AA", "Tetracycline" = "#C71585", "Quinolones" = "#DB7093")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

secondary_annotation <- HeatmapAnnotation(
  SecondaryGroup = secondary_groups,
  col = list(SecondaryGroup = c("Aminoglucosides" = "#00008B","Monobactams" = "#ADD8E6", "Cephalosporins" = "#87CEFA", "Beta-lactam+inhibitor" = "#E1FFFF", "Carbapenems" = "#008B8B", "Polymixin" = "#00FFFF", "Tetracycline" = "#FF1493", "Quinolones" = "#DC143C")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# 合并注释
#combined_annotation <- primary_annotation + secondary_annotation








Heatmap(mat_new_ATLAS_Levofloxacin_R_Amikacin_S_reordered, name = "Collateral sensitivity mode", 
        show_row_names = FALSE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        row_title = "Pseudomonas aeruginosa\nisolates",
        row_title_side = "left",
        row_title_gp = gpar(fontface = "italic"),
        column_title = "Antibiotics",
        column_title_side = "bottom",
        column_names_rot = 60,
        #column_names_gp = grid::gpar(fontsize = 12),
        #row_names_gp = grid::gpar(fontsize = 10),
        column_split = list(primary_groups, secondary_groups),  
        top_annotation = c(primary_annotation, secondary_annotation),
        heatmap_legend_param = list(
          at = c(-10,0,10),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))


#Collateral sensitivity mode in Pseudomonas aeruginosa isolates when levofloxacin_R_amikacin_S



#############################


tigecycline_R_amikacin_S


ATLAS_Tigecycline_R<- ATLAS3[grepl("^[4-9]*$", ATLAS3$TIG)==T,]

#ATLAS_Amikacin_S<- ATLAS3[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS3$AMK)==T,]

ATLAS_Tigecycline_R_Amikacin_S<- ATLAS_Tigecycline_R[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS_Tigecycline_R$AMK)==T,]
ATLAS_Tigecycline_R_Amikacin_S$LEV<- ATLAS_Tigecycline_R_Amikacin_S$LEV - 2  
ATLAS_Tigecycline_R_Amikacin_S$AMK<- ATLAS_Tigecycline_R_Amikacin_S$AMK - 4 
ATLAS_Tigecycline_R_Amikacin_S$AZM<- ATLAS_Tigecycline_R_Amikacin_S$AZM - 3
ATLAS_Tigecycline_R_Amikacin_S$CEF<- ATLAS_Tigecycline_R_Amikacin_S$CEF - 3
ATLAS_Tigecycline_R_Amikacin_S$CAZ<- ATLAS_Tigecycline_R_Amikacin_S$CAZ - 3
ATLAS_Tigecycline_R_Amikacin_S$`CAZ+AVI`<- ATLAS_Tigecycline_R_Amikacin_S$`CAZ+AVI` - 3
ATLAS_Tigecycline_R_Amikacin_S$`CTO+TAZ`<- ATLAS_Tigecycline_R_Amikacin_S$`CTO+TAZ` - 2
ATLAS_Tigecycline_R_Amikacin_S$IMI<- ATLAS_Tigecycline_R_Amikacin_S$IMI - 1
ATLAS_Tigecycline_R_Amikacin_S$MER<- ATLAS_Tigecycline_R_Amikacin_S$MER - 1
ATLAS_Tigecycline_R_Amikacin_S$DOR<- ATLAS_Tigecycline_R_Amikacin_S$DOR - 1
ATLAS_Tigecycline_R_Amikacin_S$`PIP+TAZ`<- ATLAS_Tigecycline_R_Amikacin_S$`PIP+TAZ` - 4
ATLAS_Tigecycline_R_Amikacin_S$COL<- ATLAS_Tigecycline_R_Amikacin_S$COL - 1
ATLAS_Tigecycline_R_Amikacin_S$MIN<- ATLAS_Tigecycline_R_Amikacin_S$MIN - 3
ATLAS_Tigecycline_R_Amikacin_S$TIG<- ATLAS_Tigecycline_R_Amikacin_S$TIG - 3
ATLAS_Tigecycline_R_Amikacin_S$`COLP80 `<- ATLAS_Tigecycline_R_Amikacin_S$`COLP80 ` - 1


new_ATLAS_Tigecycline_R_Amikacin_S <- select(ATLAS_Tigecycline_R_Amikacin_S,-c(AMC,AMP,CTX,`AZM+AVI`,CAR,`CAR+AVI`,ERT,SUL,AZI))




if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(grid)
library(ComplexHeatmap)

n_new_ATLAS_Tigecycline_R_Amikacin_S <- mutate_all(new_ATLAS_Tigecycline_R_Amikacin_S, ~replace(., is.na(.), 0))
mat_new_ATLAS_Tigecycline_R_Amikacin_S <- as.matrix(n_new_ATLAS_Tigecycline_R_Amikacin_S) 
colnames(mat_new_ATLAS_Tigecycline_R_Amikacin_S) <- c("Amikacin","Ceftazidime","Cefepime","Levofloxacin","Meropenem","Minocycline","Piperacillin tazobactam","Tigecycline","Aztreonam","Ceftazidime avibactam","ColistinP80","Doripenem","Imipenem","Colistin","Ceftolozane tazobactam") 
mat_new_ATLAS_Tigecycline_R_Amikacin_S_reordered <- mat_new_ATLAS_Tigecycline_R_Amikacin_S[,c("Amikacin","Aztreonam","Ceftazidime","Cefepime","Piperacillin tazobactam","Ceftazidime avibactam","Ceftolozane tazobactam","Meropenem","Doripenem","Imipenem","Colistin","ColistinP80","Levofloxacin","Minocycline","Tigecycline")] 


data <- mat_new_ATLAS_Tigecycline_R_Amikacin_S_reordered


primary_groups <- c(rep("Aminoglucosides", 1), rep("Beta-lactams", 9),rep("Polymixin", 2),rep("Quinolones", 1),rep("Tetracycline", 2))

secondary_groups <- c(rep("Aminoglucosides", 1), rep("Monobactams", 1), rep("Cephalosporins", 2),rep("Beta-lactam+inhibitor", 3),rep("Carbapenems", 3),rep("Polymixin", 2),rep("Quinolones", 1),rep("Tetracycline", 2))



primary_annotation <- HeatmapAnnotation(
  PrimaryGroup = primary_groups,
  col = list(PrimaryGroup = c("Aminoglucosides" = "#0000CD","Beta-lactams" = "#32CD32", "Polymixin" = "#20B2AA", "Quinolones" = "#DB7093", "Tetracycline" = "#C71585")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

secondary_annotation <- HeatmapAnnotation(
  SecondaryGroup = secondary_groups,
  col = list(SecondaryGroup = c("Aminoglucosides" = "#00008B","Monobactams" = "#ADD8E6", "Cephalosporins" = "#87CEFA", "Beta-lactam+inhibitor" = "#E1FFFF", "Carbapenems" = "#008B8B", "Polymixin" = "#00FFFF", "Quinolones" = "#DC143C", "Tetracycline" = "#FF1493")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)



Heatmap(mat_new_ATLAS_Tigecycline_R_Amikacin_S_reordered, name = "Collateral sensitivity mode", 
        show_row_names = FALSE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        row_title = "Pseudomonas aeruginosa\nisolates",
        row_title_gp = gpar(fontface = "italic"),
        row_title_side = "left",
        column_title = "Antibiotics",
        column_title_side = "bottom",
        column_names_rot = 60,
        #column_names_gp = grid::gpar(fontsize = 12),
        #row_names_gp = grid::gpar(fontsize = 10),
        column_split = list(primary_groups, secondary_groups),  
        top_annotation = c(primary_annotation, secondary_annotation),
        heatmap_legend_param = list(
          at = c(-10,0,10),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))


#Collateral sensitivity mode in Pseudomonas aeruginosa isolates when tigecycline_R_amikacin_S




#########################

Minocycline_R_amikacin_S



ATLAS_Minocycline_R<- ATLAS3[grepl("^[4-9]*$", ATLAS3$MIN)==T,]

#ATLAS_Amikacin_S<- ATLAS3[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS3$AMK)==T,]

ATLAS_Minocycline_R_Amikacin_S<- ATLAS_Minocycline_R[grepl("^-[0-9]*[1-9][0-9]*|[0-3]$", ATLAS_Minocycline_R$AMK)==T,]
ATLAS_Minocycline_R_Amikacin_S$LEV<- ATLAS_Minocycline_R_Amikacin_S$LEV - 2  
ATLAS_Minocycline_R_Amikacin_S$AMK<- ATLAS_Minocycline_R_Amikacin_S$AMK - 4 
ATLAS_Minocycline_R_Amikacin_S$AZM<- ATLAS_Minocycline_R_Amikacin_S$AZM - 3
ATLAS_Minocycline_R_Amikacin_S$CEF<- ATLAS_Minocycline_R_Amikacin_S$CEF - 3
ATLAS_Minocycline_R_Amikacin_S$CAZ<- ATLAS_Minocycline_R_Amikacin_S$CAZ - 3
ATLAS_Minocycline_R_Amikacin_S$`CAZ+AVI`<- ATLAS_Minocycline_R_Amikacin_S$`CAZ+AVI` - 3
ATLAS_Minocycline_R_Amikacin_S$`CTO+TAZ`<- ATLAS_Minocycline_R_Amikacin_S$`CTO+TAZ` - 2
ATLAS_Minocycline_R_Amikacin_S$IMI<- ATLAS_Minocycline_R_Amikacin_S$IMI - 1
ATLAS_Minocycline_R_Amikacin_S$MER<- ATLAS_Minocycline_R_Amikacin_S$MER - 1
ATLAS_Minocycline_R_Amikacin_S$DOR<- ATLAS_Minocycline_R_Amikacin_S$DOR - 1
ATLAS_Minocycline_R_Amikacin_S$`PIP+TAZ`<- ATLAS_Minocycline_R_Amikacin_S$`PIP+TAZ` - 4
ATLAS_Minocycline_R_Amikacin_S$COL<- ATLAS_Minocycline_R_Amikacin_S$COL - 1
ATLAS_Minocycline_R_Amikacin_S$MIN<- ATLAS_Minocycline_R_Amikacin_S$MIN - 3
ATLAS_Minocycline_R_Amikacin_S$TIG<- ATLAS_Minocycline_R_Amikacin_S$TIG - 3
ATLAS_Minocycline_R_Amikacin_S$`COLP80 `<- ATLAS_Minocycline_R_Amikacin_S$`COLP80 ` - 1


new_ATLAS_Minocycline_R_Amikacin_S <- select(ATLAS_Minocycline_R_Amikacin_S,-c(AMC,AMP,CTX,`AZM+AVI`,CAR,`CAR+AVI`,ERT,SUL,AZI))




if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(grid)
library(ComplexHeatmap)

n_new_ATLAS_Minocycline_R_Amikacin_S <- mutate_all(new_ATLAS_Minocycline_R_Amikacin_S, ~replace(., is.na(.), 0))
mat_new_ATLAS_Minocycline_R_Amikacin_S <- as.matrix(n_new_ATLAS_Minocycline_R_Amikacin_S) 
colnames(mat_new_ATLAS_Minocycline_R_Amikacin_S) <- c("Amikacin","Ceftazidime","Cefepime","Levofloxacin","Meropenem","Minocycline","Piperacillin tazobactam","Tigecycline","Aztreonam","Ceftazidime avibactam","ColistinP80","Doripenem","Imipenem","Colistin","Ceftolozane tazobactam") 
mat_new_ATLAS_Minocycline_R_Amikacin_S_reordered <- mat_new_ATLAS_Minocycline_R_Amikacin_S[,c("Amikacin","Aztreonam","Ceftazidime","Cefepime","Piperacillin tazobactam","Ceftazidime avibactam","Ceftolozane tazobactam","Meropenem","Doripenem","Imipenem","Colistin","ColistinP80","Levofloxacin","Tigecycline","Minocycline")] 


data <- mat_new_ATLAS_Minocycline_R_Amikacin_S_reordered


primary_groups <- c(rep("Aminoglucosides", 1), rep("Beta-lactams", 9),rep("Polymixin", 2),rep("Quinolones", 1),rep("Tetracycline", 2))

secondary_groups <- c(rep("Aminoglucosides", 1), rep("Monobactams", 1), rep("Cephalosporins", 2),rep("Beta-lactam+inhibitor", 3),rep("Carbapenems", 3),rep("Polymixin", 2),rep("Quinolones", 1),rep("Tetracycline", 2))



primary_annotation <- HeatmapAnnotation(
  PrimaryGroup = primary_groups,
  col = list(PrimaryGroup = c("Aminoglucosides" = "#0000CD","Beta-lactams" = "#32CD32", "Polymixin" = "#20B2AA", "Quinolones" = "#DB7093", "Tetracycline" = "#C71585")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

secondary_annotation <- HeatmapAnnotation(
  SecondaryGroup = secondary_groups,
  col = list(SecondaryGroup = c("Aminoglucosides" = "#00008B","Monobactams" = "#ADD8E6", "Cephalosporins" = "#87CEFA", "Beta-lactam+inhibitor" = "#E1FFFF", "Carbapenems" = "#008B8B", "Polymixin" = "#00FFFF", "Quinolones" = "#DC143C", "Tetracycline" = "#FF1493")),
  show_annotation_name = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)



Heatmap(mat_new_ATLAS_Minocycline_R_Amikacin_S_reordered, name = "Collateral sensitivity mode", 
        show_row_names = FALSE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        row_title = "Pseudomonas aeruginosa\nisolates",
        row_title_gp = gpar(fontface = "italic"),
        row_title_side = "left",
        column_title = "Antibiotics",
        column_title_side = "bottom",
        column_names_rot = 60,
        #column_names_gp = grid::gpar(fontsize = 12),
        #row_names_gp = grid::gpar(fontsize = 10),
        column_split = list(primary_groups, secondary_groups),  
        top_annotation = c(primary_annotation, secondary_annotation),
        heatmap_legend_param = list(
          at = c(-10,0,10),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))


#Collateral sensitivity mode in Pseudomonas aeruginosa isolates when Minocycline_R_amikacin_S

#Collateral sensitivity mode in Pseudomonas aeruginosa isolates when MIN_or_TIG_or_LEV_R_amikacin_S




###



################################################################################


df8 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin","Tobramycin"),
  delta_n = c(8,8,1),
  empty_vector = c(4,4,0.5),
  mexAB_oprM = c(8,4,1),
  mexXY_oprM = c(16,8,2)
)

rownames(df8) <- df8$Antibiotics

df8 <- df8[ , -1]

df_log2 <- df8

df_log2[,] <- lapply(df_log2, log2)

df8 <- df_log2

df_adjusted <- df8
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df8 <- as.matrix(df_adjusted)



col_labels <- c(
  expression(Delta * "nfxB"),
  "empty_vector",
  "mexAB_oprM", 
  "mexXY_oprM")



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df8[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df8, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df8)*unit(15, "mm"), 
        height = nrow(df8)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C",'D'),each=1),
        heatmap_legend_param = list(
          at = c(-1,0,1),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of delta_nfxB over expression mexAB_oprM or mexXY_oprM in Gen or Amk or Tob 
###################################################################

df9 <- data.frame(
  Antibiotics = c("Amikacin","Tobramycin","Carbenicillin"),
  PAO1 = c(16,2,128),
  empty_vector = c(4,1,64),
  mexA = c(4,1,64),
  mexB = c(4,1,64),
  mexX = c(4,1,64),
  mexY = c(4,1,64)
)

rownames(df9) <- df9$Antibiotics

df9 <- df9[ , -1]

df_log2 <- df9

df_log2[,] <- lapply(df_log2, log2)

df9 <- df_log2

df_adjusted <- df9
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df9 <- as.matrix(df_adjusted)



col_labels <- c(
  "PAO1",
  "empty_vector",
  "mexA", 
  "mexB",
  "mexX",
  "mexY")



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df9[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df9, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df9)*unit(15, "mm"), 
        height = nrow(df9)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C","D","E","F"),each=1),
        heatmap_legend_param = list(
          at = c(-2,0),
          labels = c("Susceptible",
                     "Intermediate"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of PAO1 over expression mexA or mexB or mexX or mexY in Car or Amk or Tob 

################################################################################

df10 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin","Tobramycin","Streptomycin","Kanamycin","Ciprofloxacin"),
  PAO1 = c(16,8,2,256,256,0.25),
  delta_m = c(1,2,0.5,8,128,0.0625),
  delta_n = c(4,4,1,64,128,4),
  delta_nm = c(2,2,0.5,16,64,4)
)

rownames(df10) <- df10$Antibiotics

df10 <- df10[ , -1]

df_log2 <- df10

df_log2[,] <- lapply(df_log2, log2)

df10 <- df_log2

df_adjusted <- df10
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df10 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "oprM"),
  expression(Delta * "nfxB"),
  expression(Delta * "nfxB" * Delta * "oprM")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df10[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df10, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df10)*unit(15, "mm"), 
        height = nrow(df10)*unit(15, "mm"),
        row_split = rep(c("a","b","c","d","e","f"),each=1), 
        column_split = rep(c("A","B","C",'D'),each=1),
        heatmap_legend_param = list(
          at = c(-4,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of delta_nfxB and delta_oprM and delta_nfxB_oprM in PAO1

#####################################################################

df10 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin","Tobramycin","Streptomycin","Kanamycin","Ciprofloxacin"),
  PAO1 = c(16,8,2,256,256,0.25),
  delta_r = c(16,8,0.5,8,128,0.0625),
  delta_z = c(32,32,1,64,128,4),
  delta_n = c(2,2,0.5,16,64,4),
  delta_nr = c(16,8,0.5,16,64,4),
  delta_nz = c(16,8,0.5,16,64,4)
)

rownames(df10) <- df10$Antibiotics

df10 <- df10[ , -1]

df_log2 <- df10

df_log2[,] <- lapply(df_log2, log2)

df10 <- df_log2

df_adjusted <- df10
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df10 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "oprM"),
  expression(Delta * "nfxB"),
  expression(Delta * "nfxB" * Delta * "oprM")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df10[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df10, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df10)*unit(15, "mm"), 
        height = nrow(df10)*unit(15, "mm"),
        row_split = rep(c("a","b","c","d","e","f"),each=1), 
        column_split = rep(c("A","B","C",'D'),each=1),
        heatmap_legend_param = list(
          at = c(-4,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of delta_nfxB and delta_mexR and delta_mexZ and delta_nfxB_mexR and delta_nfxB_mexZ in PAO1

####################################################################################################

df11 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(8,2,0.25),
  delta_n = c(2,0.5,4),
  empty_vector = c(2,1,4),
  oprM = c(2,1,4)
)

rownames(df11) <- df11$Antibiotics

df11 <- df11[ , -1]

df_log2 <- df11

df_log2[,] <- lapply(df_log2, log2)

df11 <- df_log2

df_adjusted <- df11
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df11 <- as.matrix(df_adjusted)



col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  expression("empty_vector" * (Delta * "nfxB")),
  expression("oprM" * (Delta * "nfxB")))



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df11[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df11, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df11)*unit(15, "mm"), 
        height = nrow(df11)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C","D"),each=1),
        heatmap_legend_param = list(
          at = c(-2,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))
#MIC of delta_nfxB over expression oprM in Gen or Tob or Cip

############################################################################

df12 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  PAO1 = c(0,0,0),
  delta_n = c(-2,-1,4),
  empty_vector = c(-2,-1,4),
  mexAB = c(-4,-3,4),
  mexXY = c(-2,-2,4)
)

rownames(df12) <- df12$Antibiotics

df12 <- df12[ , -1]

df12 <- as.matrix(df12)



col_labels <- c(
  "PAO1",
  expression(Delta * "nfxB"),
  expression("empty_vector" * (Delta * "nfxB")),
  expression("mexAB" * (Delta * "nfxB")),
  expression("mexXY" * (Delta * "nfxB"))
  )



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df12[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df12, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        column_names_gp = gpar(fontface = "italic"),
        width = ncol(df11)*unit(15, "mm"), 
        height = nrow(df11)*unit(15, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B","C","D","E"),each=1),
        heatmap_legend_param = list(
          at = c(-4,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))
#MIC of delta_nfxB over expression mexXY or mexAB in PAO1

################################################################################



df1 <- data.frame(
  Antibiotics = c("Gentamicin","Tobramycin","Ciprofloxacin"),
  wild_type = c(8,4,0.25),
  delta_n = c(2,1,4)
)

rownames(df1) <- df1$Antibiotics

df1 <- df1[ , -1]

df_log2 <- df1

df_log2[,] <- lapply(df_log2, log2)

df1 <- df_log2

df_adjusted <- df1
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df1 <- as.matrix(df_adjusted)

col_labels <- c(
  "wild type",
  expression(italic(Delta) * italic("nfxB"))
)


cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df1[i, j]), x, y, gp = gpar(fontsize = 10))
}

tiff("MIC of delta_nfxB in Tob or Gen or Cip.tiff", width = 5.7, height = 3, units = "in", res = 300)  # 设置宽高为英寸，分辨率为300 DPI
pdf("MIC of delta_nfxB in Tob or Gen or Cip.pdf", width = 5.7, height = 3)
svg("MIC of delta_nfxB in Tob or Gen or Cip.svg", width = 5.7, height = 3)


Heatmap(df1, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE,
        row_names_side = "left",
        show_column_names = TRUE, 
        column_names_side = "top",
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        #row_title = "Antibiotics",
        #row_title_side = "left",
        #column_title = "Mutants",
        #column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df1)*unit(10, "mm"), 
        height = nrow(df1)*unit(10, "mm"),
        row_split = rep(c("A","B","C"),each=1), 
        column_split = rep(c("A","B"),each=1),
        row_title = NULL,                  # 移除行分割标签
        column_title = NULL,               # 移除列分割标签
        heatmap_legend_param = list(
          at = c(-2,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

dev.off()







###################################################################################

df10 <- data.frame(
  Antibiotics = c("Gentamicin","Amikacin","Tobramycin","Streptomycin","Kanamycin","Ciprofloxacin"),
  PAO1 = c(16,8,2,256,256,0.25),
  delta_r = c(16,8,0.5,8,128,0.0625),
  delta_z = c(32,32,1,64,128,4),
  delta_n = c(2,2,0.5,16,64,4),
  delta_nr = c(16,8,0.5,16,64,4),
  delta_nz = c(16,8,0.5,16,64,4)
)

rownames(df10) <- df10$Antibiotics

df10 <- df10[ , -1]

df_log2 <- df10

df_log2[,] <- lapply(df_log2, log2)

df10 <- df_log2

df_adjusted <- df10
df_adjusted[,] <- df_adjusted[,] - df_adjusted[, 1]

df10 <- as.matrix(df_adjusted)

col_labels <- c(
  "PAO1",
  expression(Delta * "oprM"),
  expression(Delta * "nfxB"),
  expression(Delta * "nfxB" * Delta * "oprM")
)



cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%d", df10[i, j]), x, y, gp = gpar(fontsize = 10))
}


Heatmap(df10, name = "Collateral sensitivity mode",
        cell_fun = cell_fun,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,  
        show_column_dend = FALSE,  
        row_title = "Antibiotics",
        row_title_side = "left",
        column_title = "Mutants",
        column_title_side = "bottom",
        column_names_rot = 60,
        column_labels = col_labels,
        #column_names_gp = gpar(fontface = "italic"),
        width = ncol(df10)*unit(15, "mm"), 
        height = nrow(df10)*unit(15, "mm"),
        row_split = rep(c("a","b","c","d","e","f"),each=1), 
        column_split = rep(c("A","B","C",'D'),each=1),
        heatmap_legend_param = list(
          at = c(-4,0,4),
          labels = c("Susceptible",
                     "Intermediate",
                     "Resistant"),
          title = "Antibiotic sensitivity", 
          legend_gp = gpar(fill = 1:2)))

#MIC of delta_nfxB and delta_mexR and delta_mexZ and delta_nfxB_mexR and delta_nfxB_mexZ in PAO1
