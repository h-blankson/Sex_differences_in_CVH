#Script to clean up data
#edited by Harriet Blankson

#Load libraries
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
 
getwd()
#create a new directory 
dir.create("newdirectory")

#set working directory
setwd("/path/to/working/directory/")
################################################################################
## STEP1. Read gene count matrix from visit one samples only
counts<-read.csv("/path/to/countmatrix.csv")
counts[1:10,1:10]
dim(counts) 

#remove rRNAs
smallnuclist <- c("U2|RNU5A-1|RNU5A-2|RNU5A-3|RNU5A-4|RNU5A-5|RNR_5_8S|RNR_18S_|RNR_28S|MT-RNR1|MT-RNR2|RNU5A-1P|RNU5A-2P|RNU18SP|RNU28SP|RNA5S|RPS|RNVU") 
Counts <- counts[!grepl(smallnuclist, counts$gene_id),]
dim(Counts) 

#check for duplications because we will have to make gene_id
#rownames and we cannot have duplicated rownames
sum(duplicated(Counts$gene_id)) #0

# Extract record_id/sample_id to simple number
sample_cols <- colnames(Counts)[1:ncol(Counts)]
record_ids <- sub(".*_([0-9]+)\\.Visit\\.1_.*", "\\1",
  sample_cols)

head(record_ids)

#check if extraction worked
head(data.frame(original_name = sample_cols,
                record_id = record_ids), 10)

#rename name Counts
counts_v1 <- Counts

#rename sample names to the record_id
colnames(counts_v1)[1:ncol(counts_v1)] <- record_ids
colnames(counts_v1)[1:10]
anyDuplicated(record_ids) 

#keep annotations
gene_annot <- counts_v1[, c("gene_id")]
counts <- counts_v1[, 1:ncol(counts_v1)]
rownames(counts) <- counts_v1$gene_id
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

#now the meta data
################################################################################
## STEP 2. Read in phenotypic data
clinical_data <- read.csv("/path/to/metadata.csv")
clinical_data[1:10,1:10]
dim(clinical_data) 

#remove  unqualified participant
meca_visit1 <- clinical_data %>% filter(redcap_event_name %in%
                                          c("baseline_arm_2", "baseline_arm_3", "baseline_arm_4", "baseline_arm_5"))
dim(meca_visit1)  

#read batch data
batchdata <- read.csv("/path/to/batchdata.csv")
batchdata[1:10,1:10]
dim(batchdata) # 
 
#read in cell proportions from MUSIC
cells <- read.csv("/path/to/celltypedata.csv")
head(cells)
dim(cells)

#merge cell props data to pheno data
new_df <- left_join(meca_visit1,batchdata %>% 
                     dplyr:: select(record_id, Repository, Batch),by = "record_id")

head(new_df)

#non.classical.monocyte was omited, not used as reference
new_df1 <- left_join(new_df,cells %>%
                       dplyr::select(record_id, classical.monocyte, natural.killer.cell,
                                     plasmacytoid.dendritic.cell..human, plasmablast, lymphocyte,
                                     double.negative.thymocyte, platelet ), by = "record_id")

head(new_df1)
dim(new_df1) 

#add the older datasheet to get the age info
old_data <- read.csv("/path/to/olddataset.csv")
head(old_data)

old_data$age
new_df1$age <- NULL #delete

#merge date of birth
meta_visit1 <- left_join(new_df1, old_data %>%
                           dplyr::select(record_id,age), 
                         by = "record_id")

#check
meta_visit1$age

# #calculate age
# library(lubridate)
# 
# meta_visit1 <- meta_visit1 %>%
#   mutate(
#     date_birth = as.Date(date_birth),
#     date_enrollment = as.Date(enrol_date),
#     age = floor(time_length(
#       interval(date_birth, date_enrollment),
#       "years"
#     ))
#   )

# create a new dataframe of selected variables
LS7_data <- data.frame(
  record_id = as.character(meta_visit1$record_id),
  age = as.numeric(meta_visit1$age),
 # age = as.numeric(meta_visit1$age_atenrollment),
  gender = meta_visit1$male,
  bmi = as.numeric(meta_visit1$ls7_bmi),
  bp = as.numeric(meta_visit1$ls7_bpsubcomp),
  glucose = as.numeric(meta_visit1$ls7_glucosesubcomp),
  ch = as.numeric(meta_visit1$ls7_chsubcomp),
  exercise = as.numeric(meta_visit1$ls7_exercise2),
  diet = as.numeric(meta_visit1$ls7_diet2),
  smoke = as.numeric(meta_visit1$smoke_enrol),
  ls7total = as.numeric(meta_visit1$ls7_total),
 bmi_raw = as.numeric(meta_visit1$bmi ),
 sysbp = as.numeric(meta_visit1$sysbp_enrol),
 diabp = as.numeric(meta_visit1$diabp_enrol),
 glucose_raw = as.numeric(meta_visit1$glucose_enrol),
 yrs_smoked = as.numeric(meta_visit1$chol_enrol),
 ch_raw = as.numeric(meta_visit1$chol_enrol),
 yrs_smoked = as.numeric(meta_visit1$yrssmoked),
 hdl = as.numeric(meta_visit1$hdl),
 ldl = as.numeric(meta_visit1$ldl),
 trigl = as.numeric(meta_visit1$trigly),
  consent = meta_visit1$informedconsent,
  Repository = factor(meta_visit1$Repository),
  Batch = factor(meta_visit1$Batch),
  monocyte = as.numeric(meta_visit1$classical.monocyte),
  nkc = as.numeric(meta_visit1$natural.killer.cell),
  dendritic = as.numeric(meta_visit1$plasmacytoid.dendritic.cell..human),
  plasmablast = as.numeric(meta_visit1$plasmablast),
  lymphocyte = as.numeric(meta_visit1$lymphocyte),
  thymocyte = as.numeric(meta_visit1$double.negative.thymocyte),
  platelet = as.numeric(meta_visit1$platelet),
  stringsAsFactors = FALSE
)
dim(LS7_data) #571  31

LS7_data <- data.frame(LS7_data)
dim(LS7_data) # 571  31

## assign gender 
LS7_data$gender <-ifelse(LS7_data$gender  == "1", "male", "female")
table(LS7_data$gender)
#female   male 
#363    208

#check informed consent
table(LS7_data$consent) #571

# Need to match samples so IDs are in correct order as Count matrix
rows = match(colnames(counts),as.character(LS7_data$record_id))
LS7_data_subset <- LS7_data[rows,]
dim(counts) #60779   471
dim(LS7_data_subset) #471  31


# now remove NA variable in record_id and other areas
LS7_data_subset <- LS7_data_subset %>%
  filter(!is.na(record_id), !is.na(age) , !is.na(ls7total), !is.na(monocyte)
         )

dim(LS7_data_subset) #443  31

# #take out participants with descrepant LS7
# different_ls7 <- read.csv("decrepant_ids.csv")
# dim(different_ls7) #25  5
# 
# LS7_data_subset <- LS7_data_subset %>% 
#   dplyr::filter(!record_id %in% different_ls7$record_id)

#how many are left
nrow(LS7_data_subset) ; dim(LS7_data_subset)

#match
cols <- match(LS7_data_subset$record_id,colnames(counts))
Counts_subset <- counts[,cols]
dim(Counts_subset); dim(LS7_data_subset) 
# 60779   443;  443  31
head(Counts_subset)

#check to ensure pheno and count data are in the same order
identical(as.numeric(colnames(Counts_subset)), as.numeric(LS7_data_subset$record_id) ) 

#save data
write.csv(LS7_data_subset, "pheno_data.csv")
write.csv(Counts_subset, "Counts.csv")

#############################################################################
#############################################################################
#############################################################################
#do the sheets have the same ids?
setdiff(as.character(old_data$record_id),
        as.character(new_df1$record_id))

setdiff(as.character(old_data$record_id),
        as.character(new_df1$record_id))


#Comparing  LS7 for the diff datasheets
ls7_compare <- old_data %>%
  dplyr::select(record_id, ls7_total) %>%
  dplyr::rename(LS7TOTAL_old = ls7_total) %>%
  inner_join(
    new_df1 %>%
      dplyr::select(record_id, ls7_total) %>%
      dplyr::rename(LS7TOTAL_new = ls7_total),
    by = "record_id"
  ) %>%
  mutate(
    difference = LS7TOTAL_old - LS7TOTAL_new,
    same = LS7TOTAL_old == LS7TOTAL_new
  )

#compare 
table(ls7_compare$same, useNA = "ifany")

 #look at the ones that differ
ls7_compare %>% filter(same == FALSE)

#summerize differences
table(ls7_compare$difference, useNA = "ifany")

#recalculate LS7
old_data <- old_data %>%
  mutate(
    ls7_recalc = ls7_exercise +
      ls7_diet2 +
      smoke_enrol +
      ls7_bpsubcomp +
      ls7_glucosesubcomp +
      ls7_chsubcomp +
      ls7_bmi
  )

new_df1 <- new_df1 %>%
  mutate(
    ls7_recalc = ls7_exercise +
      ls7_diet2 +
      smoke_enrol +
      ls7_bpsubcomp +
      ls7_glucosesubcomp +
      ls7_chsubcomp +
      ls7_bmi
  )

#compare recalculated values
ls7_compare_recalc <- old_data %>%
  dplyr::select(record_id, ls7_recalc) %>%
  dplyr::rename(ls7_old = ls7_recalc) %>%
  inner_join(
    new_df1 %>%
      dplyr::select(record_id, ls7_recalc) %>%
      dplyr::rename(ls7_new = ls7_recalc),
    by = "record_id"
  ) %>%
  mutate(
    difference = ls7_old - ls7_new,
    same = ls7_old == ls7_new
  )

#check agreements
table(ls7_compare_recalc$same, useNA = "ifany")

#check differences
ls7_compare_recalc %>% filter(same == FALSE)

#compare recalculated to existing
table(old_data$ls7_total == old_data$ls7_recalc, useNA = "ifany")
table(new_df1$ls7_total == new_df1$ls7_recalc, useNA = "ifany")

#get descrpant ids
#LS7 descrepancy and size
different_ls7 <- new_df1 %>%
  filter(
    !is.na(ls7_total),
    !is.na(ls7_recalc),
    ls7_total != ls7_recalc
  ) %>%
  mutate(difference = ls7_total - ls7_recalc) %>%
  dplyr::select(
    record_id,
    ls7_total,
    ls7_recalc,
    difference
  )

different_ls7

#count
nrow(different_ls7)

write.csv(different_ls7, "decrepant_ids.csv")




