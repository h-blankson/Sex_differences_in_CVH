#Script to calculate the Interaction of DEGs
#edited by Harriet Blankson 09-29-2026

#load libraries
library(dplyr)
library(tidydr)

rm(list = ls())
gc() #free up memory and report memory usage


#read in degs/toptable genes
interaction_all <-read.csv("/path/to/toptable_interactions.csv")
head(interaction_all)
male <-read.csv("path/t/omale/top_table.csv")
female <- read.csv("path/to/female/toptable.csv")

head(male)
head(female)

interaction_sig <- interaction_all[interaction_all$adj.P.Val < 0.05 &
    abs(interaction_all$logFC) > log2(1.2),]

head(interaction_sig)

#compare the interaction genes
interaction_compare <- interaction_sig %>%
  dplyr::select(gene_id, 
                Interaction_logFC = logFC,
                Interaction_FDR = adj.P.Val) %>%
  left_join(female %>%
              dplyr::select(gene_id, 
                  Female_logFC = logFC,
                  Female_FDR = adj.P.Val) %>%
      left_join(male %>%
                  dplyr::select(gene_id,
                                Male_logFC = logFC,
                                Male_FDR = adj.P.Val),
                by = "gene_id"))


#arrange the columns 
interaction_compare <- interaction_compare %>%
  select(Gene = gene_id, Female_logFC, Female_FDR,  Male_logFC,
    Male_FDR, Interaction_logFC, Interaction_FDR)

head(interaction_compare)

#check interaction algebra
all.equal(interaction_compare$Interaction_logFC,
  interaction_compare$Female_logFC - interaction_compare$Male_logFC,
  tolerance = 1e-8)
#answer should be TRUE


#classify interaction genes
interaction_compare$Pattern <- with(interaction_compare, ifelse(
    Female_logFC > 0 & Male_logFC < 0,"Opposite directions",
    ifelse(Female_logFC < 0 & Male_logFC > 0, "Opposite directions",
      ifelse(abs(Female_logFC) > abs(Male_logFC),
        "Stronger in females",
        "Stronger in males" ))))

table(interaction_compare$Pattern) 
# Opposite directions Stronger in females   Stronger in males 
# 258                   3                   1 

#random chek of genes
tail(interaction_compare[, c("Gene", "Female_logFC", "Male_logFC",
  "Interaction_logFC")], 20)


#scatter plot to show interaction
library(ggplot2)
library(ggrepel)

ggplot(interaction_compare,
       aes(x = Female_logFC, y = Male_logFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,  linetype = "dashed", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0,  colour = "red", linetype = "dotted") +
  geom_point(size = 2.5, colour = "#2C7FB8", alpha = 0.8) +
  theme_classic(base_size = 14) + labs(
    x = "Female log2 Fold Change (Low vs High)",
    y = "Male log2 Fold Change (Low vs High)",
    title = "Sex interaction genes"
  )


#if you want to label top interaction genes
topGenes <- interaction_compare[
  order(interaction_compare$Interaction_FDR), ][1:20, ]

ggplot(interaction_compare, aes(Female_logFC, Male_logFC)) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  geom_abline(slope = 1,  intercept = 0,  colour = "grey50", linetype = "dotted") +
  geom_point(size = 2.5, colour = "#2C7FB8") +
  geom_text_repel(data = topGenes, aes(label = Gene), size = 3,
    max.overlaps = Inf ) +  theme_classic(base_size = 14)

#calculate Spearman correlation
rho_int<- cor.test(interaction_compare$Female_logFC,
  interaction_compare$Male_logFC,method = "spearman")

#color coded interaction scatter plot
ggplot(interaction_compare, aes(Female_logFC, Male_logFC, colour = Pattern)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0,   linetype = "dashed") +
  geom_abline(slope = 1,  intercept = 0,  colour = "grey40",  linetype = "dotted") +
  geom_point(size = 2, alpha = 0.85) +
  theme_classic(base_size = 14) +
  annotate(  "text",
    x = min(interaction_compare$Female_logFC) + 0.1,
    y = max(interaction_compare$Male_logFC) - 0.1,
    hjust = 0,
    size = 5,
    label = paste0(
      "n = ", nrow(interaction_compare),
      " genes\n",
      "Spearman's \u03C1 = ",
      round(rho_int$estimate, 3),
      "\nP < 2.2 \u00D7 10\u207B\u00B9\u2076" ) ) +
  labs( x = "Female log2 Fold Change", y = "Male log2 Fold Change",
    colour = "Interaction Pattern")

ggsave("scatter_sexinteraction.tiff", width = 8, height = 6, dpi = 600)

#add more direction in sexes
interaction_compare <- interaction_compare %>%
  mutate(
    Pattern = case_when(
      Female_logFC < 0 & Male_logFC > 0 ~ "Female down / Male up",
      Female_logFC > 0 & Male_logFC < 0 ~ "Female up / Male down",
      Female_logFC > 0 & Male_logFC > 0 ~ "Both positive",
      Female_logFC < 0 & Male_logFC < 0 ~ "Both negative",
      TRUE ~ "Other"
    )
  )

table(interaction_compare$Pattern)
#Female down / Male up  Female up / Male down 
#51                    58

#save
write.csv (interaction_compare, "opposite effect genes.csv")


#######################
#Comparing all genes 
all_gene_compare <- female %>%
  dplyr::select(  gene_id,  Female_logFC = logFC,  Female_FDR = adj.P.Val ) %>%
  dplyr::inner_join(  male %>%  dplyr::select(gene_id, Male_logFC = logFC,
   Male_FDR = adj.P.Val ), by = "gene_id" )

#correlation of all genes
cor.test(all_gene_compare$Female_logFC,
  all_gene_compare$Male_logFC,
  method = "spearman")


#plot
rho <- cor.test(
  all_gene_compare$Female_logFC,
  all_gene_compare$Male_logFC,
  method = "spearman"
)$estimate

ggplot(all_gene_compare,
       aes(Female_logFC, Male_logFC)) +geom_point(alpha = 0.12,size = 0.5,
    colour = "black") +
   geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed",
    linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  annotate( "text",
    x = min(all_gene_compare$Female_logFC) + 0.2,
    y = max(all_gene_compare$Male_logFC) - 0.2,
    hjust = 0,
    size = 5,
    label = paste0(
      "Spearman's \u03C1 = ",
      round(rho,3),
      "\nP < 2.2 \u00D7 10\u207B\u00B9\u2076"
    )
  ) +
  labs(
    x = expression("Female LS7-associated log"[2]*" fold change"),
    y = expression("Male LS7-associated log"[2]*" fold change")
  ) +
  theme_classic(base_size = 16)


label = paste0(
  "n = ", nrow(all_gene_compare),
  " genes\n",
  "Spearman's \u03C1 = ",
  round(rho,3),
  "\nP < 2.2 \u00D7 10\u207B\u00B9\u2076"
)

ggsave("scatter_all_sexinteraction.tiff", width = 6, height = 5, dpi = 600)




