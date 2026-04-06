#Edited by Harriet Blankson

library(dplyr); library(metafor)

#back calculate Standard Error (SE)
df <- read.csv("common_male_females.csv")
head(df)

df <- df %>%
  mutate(SE_m = abs(logFC_m / qnorm(pval_m / 2, lower.tail = FALSE)))

df <- df %>%
  mutate(SE_f = abs(logFC_f / qnorm(pval_f / 2, lower.tail = FALSE)))
write.csv(df, "common_male_females_all.csv")

library(purrr)
# df_m / df_f: data frames with columns: gene, logFC, SE (or CI to back-calc SE)
meta_results <- pmap_dfr(df, function(gene, logFC_m, SE_m, logFC_f, SE_f, ...) {
  res <- rma(
    yi  = c(logFC_m, logFC_f),
    sei = c(SE_m,    SE_f),
    method = "FE"
  )
  tibble(
    gene = gene,
    meta_logFC = res$b[1],
    meta_se = res$se,
    meta_p = res$pval,
    dir_consistent = sign(logFC_m) == sign(logFC_f)
  )
})

meta_results <- meta_results %>% mutate(meta_FDR = p.adjust(meta_p, "fdr"))

#sanity check
head(meta_results)
summary(meta_results$meta_p)
write.csv(meta_results, "meta_results.csv")

#forest plot
#for one gene
g <- "ADGRA3"
d <- df %>% filter(gene == g)
res <- rma(yi = c(d$logFC_m, d$logFC_f),
           sei = c(d$SE_m, d$SE_f),
           slab = c("Male", "Female"),
           method = "FE")

forest(res, xlab = "log2 Fold Change", mlab = paste("Fixed-effects model:", g))

#forest plot for all the genes
# Output folder (create if it doesn't exist)
out_dir <- "forestplots_svg"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Loop through unique genes
for (g in unique(df$gene)) {
  
  d <- df %>% filter(gene == g)
  
  # Perform meta-analysis for each gene
  res <- rma(
    yi  = c(d$logFC_m, d$logFC_f),
    sei = c(d$SE_m, d$SE_f),
    slab = c("Male", "Female"),
    method = "FE"
  )
  
  # Set output file name
  file_path <- file.path(out_dir, paste0(g, "_forest.svg"))
  
  # Open TIFF device
  svg(filename = file_path, width = 6, height = 5
      #, units = "in", res = 300
      )
  par(mar = c(5, 5, 4, 2))
  
  # Plot forest
  forest(
    res,
    xlab = "log2 Fold Change",
    mlab = paste("Fixed-effects model:", g),
    main = g,
    refline = 0,
    xlim = c(-2.5, 2.5),
    alim = c(-2.5, 2.5),
    cex = 1.1
  )

  
  dev.off()  # close the svg device
}

cat("✅ All forest plots saved in:", out_dir, "\n")
