#Functional gene set enrichment analysis
#Edited by Harriet Blankson

rm(list = ls())
gc() #free up memory and report memory usage

getwd()
setwd("/newdata")

#source : https://biostatsquid.com/fgsea-tutorial-gsea/ 

#BiocManager::install("fgsea")
#install.packages("tidyverse")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("enrichplot")
#BiocManager::install("clusterProfiler")

library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(clusterProfiler)
library(readr)
library(ComplexHeatmap)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)

gene_list <- c("ADGRA3", "AKR1C3", "B4GAT1", "DNAJC6", "KANK2", "KIR2DL4",
              "KLRF1", "MSTRG.22508", "SH2D1B", "SPATC1L", "SPTB", 
              "WAPL-DT", "ZNF595")

gene_ids <- bitr(gene_list,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

gene_ids

ego <- enrichGO(gene = gene_ids$ENTREZID,   # <- ONLY this column!
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)


# Network of genes ↔ enriched terms
svg("networkcommon.svg")
cnetplot(ego, showCategory = 10, circular = FALSE, colorEdge = TRUE)
dev.off()


# Set relevant paths
list.files()
in_path <- "/newdata/"
#dir.create ("~/Results/")
out_path <- "/newdata/"
bg_path <- "~/"

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}


## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){

  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}



### ANALYSIS
###################################################################
##  Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'top_groupmale_High_groupmale_Low.csv'))

#go

##  Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$Counts.gene_id
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

#analysis by positional gene sets (human chromosome) C1
bg_genes <- prepare_gmt(gmt_files[3], my_genes, savefile = TRUE)

rankings <- df$t # we will use the t statistics as ranking

names(rankings) <- df$Counts.gene_id # genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#look at the rankings
max(rankings)
min(rankings)

ggplot(data.frame(Counts.gene_id = names(rankings)[1:50], ranks = rankings[1:50]), aes(Counts.gene_id, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####remove duplicates and assign IDs with not gene symbols unique identifiers
# Identify indices of empty names
empty_name_indices <- which(names(rankings) == "")

# Generate unique identifiers for empty names
new_names <- paste0("Unknown_", seq_along(empty_name_indices))

# Assign these new names to the empty slots
names(rankings)[empty_name_indices] <- new_names

# Verify no empty names
sum(names(rankings) == "")

# Create a dataframe from the rankings
rankings_df <- data.frame(gene = names(rankings), rank = rankings)

# Identify and keep only the highest ranked entry for each gene
rankings_df <- rankings_df[order(-rankings_df$rank), ]
rankings_df <- rankings_df[!duplicated(rankings_df$gene), ]

# Create a unique rankings vector again
rankings <- rankings_df$rank
names(rankings) <- rankings_df$gene

# Verify no duplicates
any(duplicated(names(rankings)))  # Should return FALSE
## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

#there will be a warning message, this is ok, should not have a big results

head(GSEAres)
str(GSEAres)
gsea_results <- as.data.frame(GSEAres)

##Visualize
## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by  adj.p-val, only use p-value if you do not have for adjust p-value)
head(GSEAres[order(padj), ])
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.01])

#now that you know the pathways that are sign set for the number to plot
number_of_top_pathways_up <-25
number_of_top_pathways_down <- 25

#topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]

#topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top35pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

## 5. Save the results -----------------------------------------------
name_of_comparison <- 'male_lvh_kegglegacy_pathwayS'
background_genes <- 'all'
filename <- paste0(out_path, '', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Sort by NES for main pathways
GSEAres_ordered  <- GSEAres[order(-GSEAres$NES), ]


# Assuming GSEAres is your dataframe with GSEA results containing NES and padj columns
# Filter for significant pathways 
significant_pathways <- GSEAres %>% filter(padj < 0.05)


# Top 20 positive NES
top_positive_NES <- significant_pathways %>%
  arrange(desc(NES)) %>%
  head(10)

# Top 10 negative NES
top_negative_NES <- significant_pathways %>%
  arrange(NES) %>%
  head(10)

# Combine both datasets for plotting
top_pathways <- bind_rows(top_positive_NES, top_negative_NES)

top_pathways$pathway <-  gsub('GOBP_|KEGG_|REACTOME_|HP_|GOMF_|GOCC_', '', top_pathways$pathway)

#dotplot
male_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(male_go_plot)

ggsave( "enrichdot_male_lvh_kegglegacy.svg", width =12, height= 8, unit = "in")
ggsave( "enrichdot_male_lvh_kegglegacy.tiff", width =12, height= 3, unit = "in")


#dotplot male kegg no label
male_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 0),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(male_go_plot)

ggsave( "enrichdot_male_lvh_kegglegacy_nolab.svg", width =5, height= 3, unit = "in")
ggsave( "enrichdot_male_lvh_kegglegacy_nolab.tiff", width =5, height= 3, unit = "in")

##reactome pathway
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'top_groupmale_High_groupmale_Low.csv'))

#go

## . Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$Counts.gene_id
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

#analysis by positional gene sets (human chromosome) C1
#bg_genes <- prepare_gmt(gmt_files[2], my_genes, savefile = TRUE)

# Get Reactome pathways for humans
reactome_sets <- msigdbr(species = "Homo sapiens", 
                         category = "C2", 
                         subcategory = "CP:REACTOME")

# Optional: get only unique gene symbols
my_bg_genes <- unique(reactome_sets$gene_symbol)

# If you want it in GMT format for prepare_gmt()
gmt_ready <- reactome_sets %>%
  group_by(gs_name) %>%
  summarise(genes = paste(gene_symbol, collapse = "\t")) %>%
  mutate(description = "", .before = genes)

write.table(gmt_ready, "reactome_only.gmt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Now run your prepare_gmt
bg_genes <- prepare_gmt("reactome_only.gmt", my_genes, savefile = TRUE)


#rankings <- sign(df$logFC)*(-log10(df$P.Val)) # we will use the signed p values from spatial DGE as ranking
rankings <- df$t # we will use the t statistics as ranking

names(rankings) <- df$Counts.gene_id # genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#look at the rankings
max(rankings)
min(rankings)

ggplot(data.frame(Counts.gene_id = names(rankings)[1:50], ranks = rankings[1:50]), aes(Counts.gene_id, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####remove duplicates and assign IDs with not gene symbols unique identifiers
# Identify indices of empty names
empty_name_indices <- which(names(rankings) == "")

# Generate unique identifiers for empty names
new_names <- paste0("Unknown_", seq_along(empty_name_indices))

# Assign these new names to the empty slots
names(rankings)[empty_name_indices] <- new_names

# Verify no empty names
sum(names(rankings) == "")

# Create a dataframe from the rankings
rankings_df <- data.frame(gene = names(rankings), rank = rankings)

# Identify and keep only the highest ranked entry for each gene
rankings_df <- rankings_df[order(-rankings_df$rank), ]
rankings_df <- rankings_df[!duplicated(rankings_df$gene), ]

# Create a unique rankings vector again
rankings <- rankings_df$rank
names(rankings) <- rankings_df$gene

# Verify no duplicates
any(duplicated(names(rankings)))  # Should return FALSE
## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

#there will be a warning message, this is ok, should not have a big results


head(GSEAres)
str(GSEAres)
gsea_results <- as.data.frame(GSEAres)

##Visualize
## Check results ------------------------------------------------------
# Top  enriched pathways (ordered by  adj.p-val, only use p-value if you do not have for adjust p-value)
head(GSEAres[order(padj), ])
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.01])

#now that you know the pathways that are sign set for the number to plot
number_of_top_pathways_up <-25
number_of_top_pathways_down <- 25

topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]

topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top35pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

## 5. Save the results -----------------------------------------------
name_of_comparison <- 'male_lvh_reactome_pathwayS'
background_genes <- 'all'
filename <- paste0(out_path, '', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Sort by NES for main pathways
GSEAres_ordered  <- GSEAres[order(-GSEAres$NES), ]


# Assuming GSEAres is your dataframe with GSEA results containing NES and padj columns
# Filter for significant pathways 
significant_pathways <- GSEAres %>% filter(padj < 0.05)


# Top 20 positive NES
top_positive_NES <- significant_pathways %>%
  arrange(desc(NES)) %>%
  head(10)

# Top 10 negative NES
top_negative_NES <- significant_pathways %>%
  arrange(NES) %>%
  head(10)

# Combine both datasets for plotting
top_pathways <- bind_rows(top_positive_NES, top_negative_NES)

top_pathways$pathway <-  gsub('GOBP_|KEGG_|REACTOME_|HP_|GOMF_|GOCC_', '', top_pathways$pathway)

#dotplot
male_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(male_go_plot)

ggsave( "enrichdot_male_lvh_reactome.svg", width =14, height= 8, unit = "in")
ggsave( "enrichdot_male_lvh_reactome.tiff", width =14, height= 8, unit = "in")


#########################################################################################
##Female analysis
#########################################################################################
# Analysis ====================================================

## Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'top_groupfemale_High_groupfemale_Low.csv'))

#go

## Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$Counts.gene_id
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

#analysis by positional gene sets (human chromosome) C1
bg_genes <- prepare_gmt(gmt_files[7], my_genes, savefile = TRUE)


#rankings <- sign(df$logFC)*(-log10(df$P.Val)) # we will use the signed p values from spatial DGE as ranking
rankings <- df$t # we will use the t statistics as ranking

names(rankings) <- df$Counts.gene_id # genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#look at the rankings
max(rankings)
min(rankings)

ggplot(data.frame(Counts.gene_id = names(rankings)[1:50], ranks = rankings[1:50]), aes(Counts.gene_id, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####remove duplicates and assign IDs with not gene symbols unique idenitfiers
# Identify indices of empty names
empty_name_indices <- which(names(rankings) == "")

# Generate unique identifiers for empty names
new_names <- paste0("Unknown_", seq_along(empty_name_indices))

# Assign these new names to the empty slots
names(rankings)[empty_name_indices] <- new_names

# Verify no empty names
sum(names(rankings) == "")

# Create a dataframe from the rankings
rankings_df <- data.frame(gene = names(rankings), rank = rankings)

# Identify and keep only the highest ranked entry for each gene
rankings_df <- rankings_df[order(-rankings_df$rank), ]
rankings_df <- rankings_df[!duplicated(rankings_df$gene), ]

# Create a unique rankings vector again
rankings <- rankings_df$rank
names(rankings) <- rankings_df$gene

# Verify no duplicates
any(duplicated(names(rankings)))  # Should return FALSE
## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

#there will be a warning message, this is ok, should not have a big results

head(GSEAres)
str(GSEAres)
gsea_results <- as.data.frame(GSEAres)

##Visualize
## Check results ------------------------------------------------------
# Top enriched pathways (ordered by  adj.p-val, only use p-value if you do not have for adjust p-value)
head(GSEAres[order(padj), ])
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.01])

#now that you know the pathways that are sign set for the number to plot
number_of_top_pathways_up <-25
number_of_top_pathways_down <- 25

#topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]

#topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top35pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

## Save the results -----------------------------------------------
name_of_comparison <- 'female_lvh_go_pathwayS'
background_genes <- 'all'
filename <- paste0(out_path, '', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Sort by NES for main pathways
GSEAres_ordered  <- GSEAres[order(-GSEAres$NES), ]

# Assuming GSEAres is your dataframe with GSEA results containing NES and padj columns
# Filter for significant pathways 
significant_pathways <- GSEAres %>% filter(padj < 0.05)

# Top positive NES
top_positive_NES <- significant_pathways %>%
  arrange(desc(NES)) %>%
  head(10)

# Top 10 negative NES
top_negative_NES <- significant_pathways %>%
  arrange(NES) %>%
  head(10)

# Combine both datasets for plotting
top_pathways <- bind_rows(top_positive_NES, top_negative_NES)

top_pathways$pathway <-  gsub('GOBP_|KEGG_|REACTOME_|HP_|GOMF_|GOCC_', '', top_pathways$pathway)

#dotplot
female_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(female_go_plot)

ggsave( "enrichdot_female_lvh_go.svg", width =12, height= 8, unit = "in")
ggsave( "enrichdot_female_lvh_go.tiff", width =12, height= 8, unit = "in")


###KEGG PATHWAY
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'top_groupfemale_High_groupfemale_Low.csv'))

#go

## Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$Counts.gene_id
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

#analysis by positional gene sets (human chromosome) C1
bg_genes <- prepare_gmt(gmt_files[3], my_genes, savefile = TRUE)

rankings <- df$t # we will use the t statistics as ranking

names(rankings) <- df$Counts.gene_id # genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#look at the rankings
max(rankings)
min(rankings)

ggplot(data.frame(Counts.gene_id = names(rankings)[1:50], ranks = rankings[1:50]), aes(Counts.gene_id, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####remove duplicates and assign IDs with not gene symbols unique idenitfiers
# Identify indices of empty names
empty_name_indices <- which(names(rankings) == "")

# Generate unique identifiers for empty names
new_names <- paste0("Unknown_", seq_along(empty_name_indices))

# Assign these new names to the empty slots
names(rankings)[empty_name_indices] <- new_names

# Verify no empty names
sum(names(rankings) == "")

# Create a dataframe from the rankings
rankings_df <- data.frame(gene = names(rankings), rank = rankings)

# Identify and keep only the highest ranked entry for each gene
rankings_df <- rankings_df[order(-rankings_df$rank), ]
rankings_df <- rankings_df[!duplicated(rankings_df$gene), ]

# Create a unique rankings vector again
rankings <- rankings_df$rank
names(rankings) <- rankings_df$gene

# Verify no duplicates
any(duplicated(names(rankings)))  # Should return FALSE
## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

#there will be a warning message, this is ok, should not have a big results


head(GSEAres)
str(GSEAres)
gsea_results <- as.data.frame(GSEAres)

##Visualize
## Check results ------------------------------------------------------
# Top enriched pathways (ordered by  adj.p-val, only use p-value if you do not have for adjust p-value)
head(GSEAres[order(padj), ])
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.01])

#now that you know the pathways that are sign set for the number to plot
number_of_top_pathways_up <-25
number_of_top_pathways_down <- 25

#topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]

#topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top35pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

## Save the results -----------------------------------------------
name_of_comparison <- 'female_lvh_kegglegacy_pathwayS'
background_genes <- 'all'
filename <- paste0(out_path, '', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Sort by NES for main pathways
GSEAres_ordered  <- GSEAres[order(-GSEAres$NES), ]

# Assuming GSEAres is your dataframe with GSEA results containing NES and padj columns
# Filter for significant pathways 
significant_pathways <- GSEAres %>% filter(padj < 0.05)

# Top positive NES
top_positive_NES <- significant_pathways %>%
  arrange(desc(NES)) %>%
  head(10)

# Top 10 negative NES
top_negative_NES <- significant_pathways %>%
  arrange(NES) %>%
  head(10)

# Combine both datasets for plotting
top_pathways <- bind_rows(top_positive_NES, top_negative_NES)

top_pathways$pathway <-  gsub('GOBP_|KEGG_|REACTOME_|HP_|GOMF_|GOCC_', '', top_pathways$pathway)

#dotplot
female_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(female_go_plot)

ggsave( "enrichdot_female_lvh_kegglegacy.svg", width =12, height= 8, unit = "in")
ggsave( "enrichdot_female_lvh_kegglegacy.tiff", width =12, height= 8, unit = "in")


##REACTOME PATHWAY
list.files(in_path)
df <- read.csv(paste0(in_path, 'top_groupfemale_High_groupfemale_Low.csv'))

## Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$Counts.gene_id
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

# Get Reactome pathways for humans
reactome_sets <- msigdbr(species = "Homo sapiens", 
                         category = "C2", 
                         subcategory = "CP:REACTOME")

# Optional: get only unique gene symbols
my_bg_genes <- unique(reactome_sets$gene_symbol)

# If you want it in GMT format for prepare_gmt()
gmt_ready <- reactome_sets %>%
  group_by(gs_name) %>%
  summarise(genes = paste(gene_symbol, collapse = "\t")) %>%
  mutate(description = "", .before = genes)

write.table(gmt_ready, "reactome_only.gmt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Now run your prepare_gmt
bg_genes <- prepare_gmt("reactome_only.gmt", my_genes, savefile = TRUE)

rankings <- df$t # we will use the t statistics as ranking

names(rankings) <- df$Counts.gene_id # genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#look at the rankings
max(rankings)
min(rankings)

ggplot(data.frame(Counts.gene_id = names(rankings)[1:50], ranks = rankings[1:50]), aes(Counts.gene_id, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####remove duplicates and assign IDs with not gene symbols unique idenitfiers
# Identify indices of empty names
empty_name_indices <- which(names(rankings) == "")

# Generate unique identifiers for empty names
new_names <- paste0("Unknown_", seq_along(empty_name_indices))

# Assign these new names to the empty slots
names(rankings)[empty_name_indices] <- new_names

# Verify no empty names
sum(names(rankings) == "")

# Create a dataframe from the rankings
rankings_df <- data.frame(gene = names(rankings), rank = rankings)

# Identify and keep only the highest ranked entry for each gene
rankings_df <- rankings_df[order(-rankings_df$rank), ]
rankings_df <- rankings_df[!duplicated(rankings_df$gene), ]

# Create a unique rankings vector again
rankings <- rankings_df$rank
names(rankings) <- rankings_df$gene

# Verify no duplicates
any(duplicated(names(rankings)))  # Should return FALSE
## Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

#there will be a warning message, this is ok, should not have a big results


head(GSEAres)
str(GSEAres)
gsea_results <- as.data.frame(GSEAres)

##Visualize
## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by  adj.p-val, only use p-value if you do not have for adjust p-value)
head(GSEAres[order(padj), ])
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.01])

#now that you know the pathways that are sign set for the number to plot
number_of_top_pathways_up <-25
number_of_top_pathways_down <- 25

topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]

topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top35pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

## Save the results -----------------------------------------------
name_of_comparison <- 'female_lvh_reactome_pathwayS'
background_genes <- 'all'
filename <- paste0(out_path, '', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Sort by NES for main pathways
GSEAres_ordered  <- GSEAres[order(-GSEAres$NES), ]


# Assuming GSEAres is your dataframe with GSEA results containing NES and padj columns
# Filter for significant pathways 
significant_pathways <- GSEAres %>% filter(padj < 0.05)


# Top positive NES
top_positive_NES <- significant_pathways %>%
  arrange(desc(NES)) %>%
  head(10)

# Top 10 negative NES
top_negative_NES <- significant_pathways %>%
  arrange(NES) %>%
  head(10)

# Combine both datasets for plotting
top_pathways <- bind_rows(top_positive_NES, top_negative_NES)

top_pathways$pathway <-  gsub('GOBP_|KEGG_|REACTOME_|HP_|GOMF_|GOCC_', '', top_pathways$pathway)

#dotplot
female_go_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), fill = -log10(padj), color = -log10(padj)),
             shape = 21, stroke = 1.5) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Enrichment of GO Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Functional Category",
    size = expression(-log[10](padj)),
    fill = expression(-log[10](padj)),
    color = expression(-log[10](padj))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    #  panel.grid.major.y = element_blank(),
    #  panel.grid.minor.y = element_blank()
  )

print(female_go_plot)

ggsave( "enrichdot_female_lvh_reactome.svg", width =12, height= 8, unit = "in")
ggsave( "enrichdot_female_lvh_reactome.tiff", width =12, height= 8, unit = "in")


