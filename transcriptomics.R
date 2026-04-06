### MECA Project R script for data manipulation and DEG analysis
# Edit by Harriet Blankson
#https://ucdavis-bioinformatics-training.github.io/2022-April-GGI-DE-in-R/data_analysis/enrichment_with_quizzes_fixed 

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

#set working directory
#dir.create("newdata")
#setwd("~/newdata")
################################################################################
## STEP1. Read gene count matrix from visit one samples only
counts<-read.csv("~/gene_count_matrix.csv")
counts[1:10,1:10]
dim(counts) 

#separate gene id and assign stringtie name
exp_counts <- separate(counts, col=gene_id, into= c("stringtie_id", "gene_id"),
                      extra = "merge", sep ="\\|", fill = "left")
exp_counts[1:10,1:10]
exp_counts$stringtie_id <- NULL
dim(exp_counts) 

#remove rRNAs
smallnuclist <- c("U2|RNU5A-1|RNU5A-2|RNU5A-3|RNU5A-4|RNU5A-5|RNR_5_8S|RNR_18S_|RNR_28S|MT-RNR1|MT-RNR2|RNU5A-1P|RNU5A-2P|RNU18SP|RNU28SP|RNA5S|RPS|RNVU") 
Counts <- exp_counts[!grepl(smallnuclist, exp_counts$gene_id),]
dim(Counts) 

# Extract record_id/sample_id to simple number
record_id <- sub(".*_(\\d+)$", "\\1", names(Counts))

#insert the record numbers back into the count matrix
names(Counts)<-record_id
names(Counts)# 

rownames(Counts) <- make.unique(Counts$gene_id)
Counts$gene_id <- NULL
dim(Counts)


#read in annotation file
anno <- read.delim("/ensembl_hm_112.txt")
head(anno)
#check and remove duplicates
any(duplicated(anno$Gene.stable.ID))
dup <-duplicated(anno$Gene.stable.ID)
df_unique <- anno[!duplicated(anno$Gene.stable.ID), ]
print(df_unique)
# save for later

################################################################################
## STEP 2. Read in phenotypic data
clinical_data <- read.csv("~/clin_data.csv")
clinical_data[1:10,1:10]
dim(clinical_data) 

#remove  unqualified participant
clinical_data <- clinical_data %>% filter(record_id != 60186) #probably an error age (minor)
dim(clinical_data) 

LS7_data<-cbind(record_id=clinical_data$record_id,
                age=clinical_data$age_atenrollment,
                gender=clinical_data$male,
                bmi=clinical_data$ls7_bmi,
                bp=clinical_data$ls7_bpsubcomp,
                glucose=clinical_data$ls7_glucosesubcomp,
                ch=clinical_data$ls7_chsubcomp,
                exercise=clinical_data$ls7_exercise2,
                diet=clinical_data$ls7_diet2,
                smoke=clinical_data$smoke_enrol,
                ls7total=clinical_data$ls7_total, 
                glulev=clinical_data$glucose_level_categories,
                ls7_tertile = clinical_data$ls7_bi, 
                ls7_tertile_2 = clinical_data$ls7_tertile_2, 
                ls7_tertile_3 = clinical_data$ls7_tertile_3)

LS7_data[LS7_data<0]<-NA #replace negative values with NA's
LS7_data<-na.omit(LS7_data) #omit NA's which were actully negative
LS7_data <- data.frame(LS7_data)
dim(LS7_data) # 393 x 11

## assign gender and or omit females
LS7_data$gender <-ifelse(LS7_data$gender  == "1", "male", "female")
table(LS7_data$gender)

# Need to match samples so IDs are in correct order as Count matrix
rows = match(names(Counts),LS7_data$record_id)
LS7_data_subset <- LS7_data[rows,]
dim(Counts) 
dim(LS7_data_subset) 

# now remove NA variable in record_id
LS7_data_subset <- LS7_data_subset[!is.na(LS7_data_subset$record_id),]
cols <- match(LS7_data_subset$record_id,names(Counts))
Counts_subset <- Counts[,cols]
dim(Counts_subset); dim(LS7_data_subset)

head(Counts_subset)


#check 
identical(as.numeric(names(Counts_subset)), as.numeric(LS7_data_subset$record_id) ) 
        
## Change names back to Pheno_data and counts
pheno_data <- LS7_data_subset
Counts <- Counts_subset

#reassign names
str(pheno_data)
pheno_data$age <- as.numeric (pheno_data$age)
pheno_data$ls7total <- as.numeric(as.character(pheno_data$ls7total))
pheno_data$ls7 <- ifelse (pheno_data$ls7total < 10, "Low", "High") # combining low and intermediate groups

table(pheno_data$ls7) #

#set groups
pheno_data$group <- as.factor(paste(pheno_data$gender, pheno_data$ls7, sep="."))
ls7 <- as.factor(paste( pheno_data$ls7))

table(pheno_data$group)
table(ls7)

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

## now create contrast matrix and design matrix
design5 <- model.matrix(~ 0 + group , data = pheno_data)
head(design5)

## Now create a contrasts matrix and save it
# contrast for by sex comparison
contr.matrix2 <- makeContrasts(femaleLow_vs_High = groupfemale.Low - groupfemale.High,
                               maleLow_vs_High= groupmale.Low - groupmale.High,
                               levels = colnames(design5))

contr.matrix2

#save pheno and count datasets that have been filtered
head(pheno_data)
head(Counts)
write.csv(pheno_data, "pheno_data.csv")
write.csv(Counts, "Counts.csv")
################################################################################
## STEP.3 Let's LIMMA baby
###############################
#  Create DGEList
dge <- DGEList(counts = Counts)

#remove outliers
max_cpm_values <- apply(cpm(dge), 1, max)
outlier_threshold <- quantile(max_cpm_values, 0.99)
non_outliers <- which(max_cpm_values <= outlier_threshold)
d1 <- dge[non_outliers, ]
dim(d1) # 60181   373

#calculate CPM 
cpm_vals <- cpm(d1)

#filter low expression genes
min_group_size <- min(table(pheno_data$group)) #find the min group size
keep <- rowSums(cpm_vals > 1) >= min_group_size
dge_filtered <- d1[keep, , keep.lib.sizes = FALSE]
dim(dge_filtered) #16550   373

#normalize with TMM
dge_filtered <- calcNormFactors(dge_filtered)

logcpm <- cpm (dge_filtered, log = TRUE, prior.count = 1) #can be used for  gene heatmap

write.csv(logcpm, "logcpm.csv")

#############################################
# STEP 4. Proceed with voom for sex analysis
v <- voom(dge_filtered, design = design5, plot = TRUE, normalize = "quantile")
#get expression data from voom plot
expr_data <- v$E

write.csv(expr_data, "expr_data_sex.csv")

# then fit to the linear model using the design  
fit <- lmFit(v, design5, weights = pheno_data$weights)

# fit to contrast matix, to identify contrast variables
cfit <- contrasts.fit(fit, contrasts=contr.matrix2)
# finally apply Bayesian correction
efit <- eBayes(cfit)
# Mason paper used 1.5 fold changes 
tfit <- treat(efit, lfc=(log2(1)))
dt <- decideTests(tfit)
summary(dt)

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

# Create a volcano plot
volvano_top_groupmale_High_groupmale_Low <- top_groupmale_High_groupmale_Low %>%
  mutate(Color = case_when(
    adj.P.Val >= 0.05 ~ "grey",  
    TRUE ~ "red",     # Non-significant genes                        # Other category (if any)
  ))

ggplot(volcano_top_groupmale_High_groupmale_Low, aes(x = logFC, y = -log10(adj.P.Val), color = Color)) +
  geom_point() +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "orange" = "orange", "black" = "black", "forestgreen" = "forestgreen")) +
  theme_minimal() +
  labs(
    x = expression(log[2]~"FC"),  # Corrected x-axis label
    y = expression(-log[10]~"Adj P-Value"),  # Corrected y-axis label
    title = "Male Low vs High"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "green") +
  geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "green") +  # Add green vertical lines
  scale_y_continuous(limits = c(0, 4)) +
  theme(
    axis.line = element_line(),          # ensure axis lines show
    plot.title = element_text(size = 22),                # Title size
    axis.title.x = element_text(size = 22),              # X-axis title size
    axis.title.y = element_text(size = 22),              # Y-axis title size
    axis.text.x = element_text(size = 22),               # X-axis text size
    axis.text.y = element_text(size = 22),               # Y-axis text size
    plot.background = element_rect(fill = "white"),      # Set the background to white
    legend.position = "none"                             # No legend
  )
ggsave("volcano_male_Low_vs_High.tiff", width = 6.67, height = 6.67, units = "in")
ggsave ("volcano_male_lvh.svg")

# enhanced volcano male
EnhancedVolcano(top_groupmale_High_groupmale_Low,
                lab = rownames(top_groupmale_High_groupmale_Low), # Or a column name with gene labels
                x = 'logFC',
                y = 'adj.P.Val', 
                ylim = c(0,5),
                pCutoff = 0.05, 
                FCcutoff = log2(1.2))+ 
                theme(
                  plot.title = element_blank(),
                  plot.subtitle = element_blank(),
                  plot.caption = element_blank()
                )

ggsave("volcano_lvhenhanced_malelvh.tiff", width = 6, height = 6, units = "in")
ggsave ("volcano_lvhenhanced_malelvh_lvh.svg")


# Create a volcano plot
volcano_top_groupfemale_High_groupfemale_Low  <- top_groupfemale_High_groupfemale_Low  %>%
  mutate(Color = case_when(
    adj.P.Val >= 0.05 ~ "grey",  
    TRUE ~ "red",     # Non-significant genes                        # Other category (if any)
  ))

ggplot(volcano_top_groupfemale_High_groupfemale_Low , aes(x = logFC, y = -log10(adj.P.Val), color = Color)) +
  geom_point() +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "orange" = "orange", "black" = "black", "forestgreen" = "forestgreen")) +
  theme_minimal() +
  labs(
    x = expression(log[2]~"FC"),  # Corrected x-axis label
    y = expression(-log[10]~"Adj P-Value"),  # Corrected y-axis label
    title = "Female Low vs High"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "green") +
    geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "green") +  # Add green vertical lines
  scale_y_continuous(limits = c(0, 4)) +
  theme(
    axis.line = element_line(),          # ensure axis lines show
    plot.title = element_text(size = 22),                # Title size
    axis.title.x = element_text(size = 22),              # X-axis title size
    axis.title.y = element_text(size = 22),              # Y-axis title size
    axis.text.x = element_text(size = 22),               # X-axis text size
    axis.text.y = element_text(size = 22),               # Y-axis text size
    plot.background = element_rect(fill = "white"),      # Set the background to white
    legend.position = "none"                             # No legend
  )

ggsave("volcano_female_lvh.tiff", width = 6.67, height = 6.67, units = "in")
ggsave ("volcano_female_lvh.svg")

# enhanced volcano female
EnhancedVolcano(top_groupfemale_High_groupfemale_Low,
                lab = rownames(top_groupfemale_High_groupfemale_Low), # Or a column name with gene labels
                x = 'logFC',
                y = 'adj.P.Val', 
                ylim = c(0,5),
                pCutoff = 0.05, 
                FCcutoff = log2(1.2))+ theme(
                plot.title = element_blank(),
                plot.subtitle = element_blank(),
                plot.caption = element_blank()
                )

ggsave("volcano_lvhenhanced_femalelvh.tiff", width = 6, height = 6, units = "in")
ggsave ("volcano_lvhenhanced_femalelvh_lvh.svg")


