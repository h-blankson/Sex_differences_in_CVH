#script for correlatio matrix
#edited by Harriet Blankson 09/29/2026

#load libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(Hmisc)
library(GGally)
library(corrgram)
library(psych)
library(gridExtra)
library(ggplot2)
library(corrplot)
library(corrplot)
library(grid)
library(gridGraphics)
library(RColorBrewer)
library(ggpubr)
library(readr)
library(Cairo)

#load pheno data
MECA_data <-read.csv("pheno_data.csv")

#separate male and female data
male_data <- filter(MECA_data, gender == "male")
female_data <- filter(MECA_data, gender == "female")
dim(male_data)
dim(female_data)

#create dataframe for females
female_ls7 <- cbind(
  BP = female_data$bp,
  BMI = female_data$bmi,
  Exercise = female_data$exercise,
  Smoking = female_data$smoke,
  Glucose = female_data$glucose,
  Cholesterol = female_data$ch,
  Diet = female_data$diet,
  'LS7 Total' = female_data$ls7total
)

#calculate correlation for females
female_ls7_corr <- cor(
  female_ls7,
  method = "spearman",
  use = "pairwise.complete.obs"
)

sort(female_ls7_corr, decreasing = TRUE)

#plot and save matrix for demales
svg("female_cor_plot.svg",  width = 5, height = 5)
tiff("female_cor_plot.tiff", units ="in", res = 300,  width = 5, height = 5)

#plot the matrix
 corrplot(
  female_ls7_corr,
  type = "lower",
  method = "color",
  addCoef.col = "black",
  number.cex = 0.8,
  col = brewer.pal(n = 8, name = "RdYlBu"),
  col.lim = c(-1, 1),
  diag = FALSE
)
dev.off()

#create dataframe for MALES
male_ls7 <- cbind(
  BP = male_data$bp,
  BMI = male_data$bmi,
  Exercise = male_data$exercise,
  Smoking = male_data$smoke,
  Glucose = male_data$glucose,
  Cholesterol = male_data$ch,
  Diet = male_data$diet,
  'LS7 Total' = male_data$ls7total
)

#calculate correlation for males
male_ls7_corr <- cor(
  male_ls7,
  method = "spearman",
  use = "pairwise.complete.obs"
)

sort(male_ls7_corr, decreasing = TRUE)

#plot and save matrix for males
svg("male_cor_plot.svg",  width = 5, height = 5)
tiff("male_cor_plot.tiff", units ="in", res = 300,  width = 5, height = 5)
#plot the matrix
corrplot(
  male_ls7_corr,
  type = "lower",
  method = "color",
  addCoef.col = "black",
  number.cex = 0.8,
  col = brewer.pal(n = 8, name = "RdYlBu"),
  col.lim = c(-1, 1),
  diag = FALSE
)

dev.off()


