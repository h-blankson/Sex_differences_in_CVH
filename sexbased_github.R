#Script for differential expression analysis based on sex
#source:https://ucdavis-bioinformatics-training.github.io/2022-April-GGI-DE-in-R/data_analysis/enrichment_with_quizzes_fixed 
#edited by Harriet Blankson 09-29-2026

#load libraries
library(smooth)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(matrixStats)
library(Rfast)
library(ggplot2)
library(devtools)
library(ggfortify)
library(edgeR)
library(limma)
library(MatchIt)
library(SummarizedExperiment)
library(RColorBrewer)
library(limma)
library(edgeR)
library(EnhancedVolcano)
library(dplyr)
library(readxl)
library(pheatmap)
library("statmod")
library(ComplexHeatmap)
library(grid)

rm(list = ls())
gc() #free up memory and report memory usage

#create a new directory 
dir.create("sexbased")
#set working directory
 setwd("/pathtoworkingdirectory")
################################################################################
## STEP1. Read in dataset
Counts<-read.csv("/path/to/countmatrix/Counts.csv",
                 check.names = FALSE, row.names = 1)
Counts[1:10,1:10]
dim(Counts) #  60779   443

#read in phenodata
pheno_data <- read.csv("/path/to/phenodata.csv")
pheno_data[1:10,1:10]
dim(pheno_data) #443  22

str(pheno_data) #chk structure
pheno_data$age <- as.numeric (pheno_data$age) #make sure age is numeric
pheno_data$ls7total <- as.numeric(as.character(pheno_data$ls7total))

#create bilaterals
pheno_data$ls7 <- ifelse (pheno_data$ls7total < 10, "Low", "High") # combining low and intermediate groups
table(pheno_data$ls7) # high 90 low 353

pheno_data$age <- as.numeric(as.character(pheno_data$age))
pheno_data$Batch <- factor(pheno_data$Batch)

#set groups
pheno_data$group <- as.factor(paste(pheno_data$gender, 
                                    pheno_data$ls7, sep="."))

ls7 <- as.factor(paste( pheno_data$ls7))

#check
table(pheno_data$group) 
#female.High 47  female.Low 227  male.High  43  male.Low 126

##matchit for propensity scoring to correct for age and sex for the 2 group analyses
targetsb <- pheno_data
targetsb$ls7 <- factor(targetsb$ls7, levels = c("Low", "High"))  # define reference

#remove missing data
targetsb <- targetsb[!is.na(targetsb$ls7) & !is.na(targetsb$age), ]

#Perform matching on age
match_out <- matchit(ls7 ~ age  , data = targetsb, 
                     method = "full",
                     estimand = "ATT")

# Extract matched data
matched_data <- match.data(match_out) # matched data for age for sex analysis
pheno_data$weights <-matched_data$weights[match(pheno_data$record_id, matched_data$record_id)] # phenodata with weightfor  the sex based dataset analyses
head(pheno_data)

###############################
## STEP.2 Create design and contrast
###############################
## now create contrast matrix and design matrix
design <- model.matrix(~ 0 + group +  monocyte + nkc + dendritic + plasmablast 
                        + lymphocyte + thymocyte + platelet + Batch, 
                        data = pheno_data) 

head(design)

## Now create a contrasts matrix 
contr.matrix <- makeContrasts(femaleLow_vs_High = groupfemale.Low - groupfemale.High,
                               maleLow_vs_High= groupmale.Low - groupmale.High,
                               levels = colnames(design))

contr.matrix

head(pheno_data)
head(Counts)
#check to ensure pheno and count data are in the same order
identical(as.numeric(names(Counts)), as.numeric(pheno_data$record_id) ) 

###############################
## STEP.3 Let's LIMMA baby
###############################
#  Create DGEList
dge <- DGEList(counts = Counts)

#calculate CPM 
cpm_vals <- edgeR::cpm(dge)

#filter low expression genes
min_group_size <- min(table(pheno_data$group)) #find the min group size
keep <- rowSums(cpm_vals > 1) >= min_group_size
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge_filtered) #17158   373

#normalize with TMM
dge_filtered <- calcNormFactors(dge_filtered)

logcpm <- edgeR::cpm (dge_filtered, log = TRUE, prior.count = 1) #for  gene heatmap

#save log normalized values
write.csv(logcpm, "logcpm_sex.csv")
#############################################
# # 4. Proceed with voom for sex analysis
v <- voom(dge_filtered, design = design, plot = TRUE, normalize = "quantile")

#get expression data from voom plot
expr_data <- v$E

#save expression data
write.csv(expr_data, "expr_data_sex.csv")

# # then fit to the linear model using the design  bilateral
fit <- lmFit(v, design, weights = pheno_data$weights)
# fit to contrast matix, to identify contrast variables bilateral
cfit <- contrasts.fit(fit, contrasts=contr.matrix)
# finally apply Bayesian correction
efit <- eBayes(cfit)
# Mason paper used 1.5 fold changes
tfit <- treat(efit, lfc=(log2(1)))
dt <- decideTests(tfit)
summary(dt)
# femaleLow_vs_High maleLow_vs_High
# Down                 370             112
# NotSig             16637           17347
# Up                   502              50

#get DEGs
top_groupfemale_High_groupfemale_Low <- topTable(efit, coef = "femaleLow_vs_High", n = Inf, adjust.method = "BH")
top_groupmale_High_groupmale_Low <- topTable(efit, coef = "maleLow_vs_High", n = Inf, adjust.method = "BH")
#save toptable
write.csv(top_groupfemale_High_groupfemale_Low, file = "top_groupfemale_High_groupfemale_Low.csv")
write.csv(top_groupmale_High_groupmale_Low, file = "top_groupmale_High_groupmale_Low.csv")

#save degs list
filtered_femaleHigh_vs_femaleLow <- top_groupfemale_High_groupfemale_Low[top_groupfemale_High_groupfemale_Low$adj.P.Val < 0.05 & abs(top_groupfemale_High_groupfemale_Low$logFC) > log2(1.2), ]
output_file <- "degs_femaleLow_vsHigh.csv"
write.csv(filtered_femaleHigh_vs_femaleLow, file = output_file, row.names = TRUE)

filtered_maleHigh_vs_maleLow <- top_groupmale_High_groupmale_Low[top_groupmale_High_groupmale_Low$adj.P.Val < 0.05 & abs(top_groupmale_High_groupmale_Low$logFC) > log2(1.2), ]
output_file <- "deg_maleLow_vsHigh.csv"
write.csv(filtered_maleHigh_vs_maleLow, file = output_file, row.names = TRUE)

#######################################################
#######################################################
# 5. Create a volcano plot
#######################################################
# enhanced volcano male
EnhancedVolcano(top_groupmale_High_groupmale_Low,
                lab = rownames(top_groupmale_High_groupmale_Low), # Or a column name with gene labels
                x = 'logFC',
                y = 'adj.P.Val', 
                ylim = c(0,4),
                xlim = c(-3,3),
                pCutoff = 0.05, 
                FCcutoff = log2(1.2))+ 
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

ggsave("volcano_lvhenhanced_malelvh.tiff", width = 6, height = 6, units = "in")
ggsave ("volcano_lvhenhanced_malelvh_lvh.svg")

# enhanced volcano female
EnhancedVolcano(top_groupfemale_High_groupfemale_Low,
                lab = rownames(top_groupfemale_High_groupfemale_Low), # Or a column name with gene labels
                x = 'logFC',
                y = 'adj.P.Val', 
                ylim = c(0,4),
                xlim = c(-3,3),
                pCutoff = 0.05, 
                FCcutoff = log2(1.2))+ theme(
                  plot.title = element_blank(),
                  plot.subtitle = element_blank(),
                  plot.caption = element_blank()
                )

ggsave("volcano_lvhenhanced_femalelvh.tiff", width = 6, height = 6, units = "in")
ggsave ("volcano_lvhenhanced_femalelvh_lvh.svg")



