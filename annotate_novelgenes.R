#script to annonate novel genes

#load libraries
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------
# EDIT THESE PATHS
# ------------------------------------------------------------


novel_id_file <- "novel_genes.txt"
stringtie_gtf <- "/path/to/merged/stringtiematrix.gtf"
reference_gtf <- "/pathtoreferenceuseforannotation.gtf"
output_prefix <- "novel_gene_annotation"
getwd()


#set working directory
#setwd("~/Library/CloudStorage/OneDrive-morehouseschoolofmedicine/Documents/MECA/meca_new/visit1_v1")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
get_metadata_column <- function(gr, possible_names, default = NA_character_) {
  available <- colnames(mcols(gr))
  selected <- possible_names[possible_names %in% available]
  
  if (length(selected) == 0) {
    return(rep(default, length(gr)))
  }
  
  as.character(mcols(gr)[[selected[1]]])
}

collapse_values <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(sort(x), collapse = "; ")
}

# ------------------------------------------------------------
# Read requested MSTRG IDs
# ------------------------------------------------------------

novel_ids <- read_lines(novel_id_file) |>
  trimws() |>
  unique()

novel_ids <- novel_ids[nzchar(novel_ids)]

if (length(novel_ids) == 0) {
  stop("No IDs were found in: ", novel_id_file)
}

message("Requested novel IDs: ", length(novel_ids))

# ------------------------------------------------------------
# Import GTF files
# ------------------------------------------------------------

message("Importing StringTie GTF...")
stringtie <- import(stringtie_gtf)

message("Importing reference GTF...")
reference <- import(reference_gtf)

# ------------------------------------------------------------
# Extract StringTie transcript coordinates
# ------------------------------------------------------------

stringtie_type <- get_metadata_column(stringtie, "type")
stringtie_gene_id <- get_metadata_column(
  stringtie,
  c("gene_id", "geneID", "gene")
)
stringtie_transcript_id <- get_metadata_column(
  stringtie,
  c("transcript_id", "transcriptID")
)

is_transcript <- stringtie_type == "transcript"
is_requested <- stringtie_gene_id %in% novel_ids

novel_transcripts <- stringtie[is_transcript & is_requested]

if (length(novel_transcripts) == 0) {
  stop(
    "None of the requested IDs were found in transcript rows of ",
    stringtie_gtf
  )
}

novel_transcripts$gene_id <- stringtie_gene_id[
  is_transcript & is_requested
]

novel_transcripts$transcript_id <- stringtie_transcript_id[
  is_transcript & is_requested
]

# Collapse transcript isoforms into one gene-level interval
novel_df <- as.data.frame(novel_transcripts) |>
  transmute(
    Novel_ID = gene_id,
    Transcript_ID = transcript_id,
    Chromosome = as.character(seqnames),
    Start = start,
    End = end,
    Strand = as.character(strand)
  )

novel_gene_df <- novel_df |>
  group_by(Novel_ID, Chromosome, Strand) |>
  summarise(
    Start = min(Start),
    End = max(End),
    Width_bp = End - Start + 1L,
    Transcript_count = n_distinct(Transcript_ID),
    Transcript_IDs = collapse_values(Transcript_ID),
    .groups = "drop"
  )

novel_gr <- makeGRangesFromDataFrame(
  novel_gene_df,
  seqnames.field = "Chromosome",
  start.field = "Start",
  end.field = "End",
  strand.field = "Strand",
  keep.extra.columns = TRUE
)

# ------------------------------------------------------------
# Extract reference genes
# ------------------------------------------------------------

reference_type <- get_metadata_column(reference, "type")
reference_gene_id <- get_metadata_column(
  reference,
  c("gene_id", "geneID", "ID")
)
reference_gene_symbol <- get_metadata_column(
  reference,
  c("gene_name", "gene_symbol", "Name")
)
reference_gene_type <- get_metadata_column(
  reference,
  c("gene_biotype", "gene_type", "biotype")
)

reference_genes <- reference[reference_type == "gene"]

reference_genes$Reference_gene_ID <- reference_gene_id[
  reference_type == "gene"
]

reference_genes$Reference_gene_symbol <- reference_gene_symbol[
  reference_type == "gene"
]

reference_genes$Reference_gene_type <- reference_gene_type[
  reference_type == "gene"
]

missing_symbol <- is.na(reference_genes$Reference_gene_symbol) |
  reference_genes$Reference_gene_symbol == ""

reference_genes$Reference_gene_symbol[missing_symbol] <-
  reference_genes$Reference_gene_ID[missing_symbol]

# Keep standard chromosomes shared between the two files
common_seqlevels <- intersect(
  seqlevels(novel_gr),
  seqlevels(reference_genes)
)

if (length(common_seqlevels) == 0) {
  stop(
    "No chromosome names match between the GTF files. ",
    "One file may use 'chr1' while the other uses '1'."
  )
}

novel_gr <- keepSeqlevels(
  novel_gr,
  common_seqlevels,
  pruning.mode = "coarse"
)

reference_genes <- keepSeqlevels(
  reference_genes,
  common_seqlevels,
  pruning.mode = "coarse"
)

# ------------------------------------------------------------
# Find nearest annotated gene
# ------------------------------------------------------------

nearest_index <- nearest(
  novel_gr,
  reference_genes,
  ignore.strand = TRUE
)

nearest_distance <- distance(
  novel_gr,
  reference_genes[nearest_index],
  ignore.strand = TRUE
)

nearest_df <- tibble(
  Novel_ID = novel_gr$Novel_ID,
  Nearest_gene_ID =
    reference_genes$Reference_gene_ID[nearest_index],
  Nearest_gene_symbol =
    reference_genes$Reference_gene_symbol[nearest_index],
  Nearest_gene_type =
    reference_genes$Reference_gene_type[nearest_index],
  Nearest_gene_chromosome =
    as.character(seqnames(reference_genes[nearest_index])),
  Nearest_gene_start =
    start(reference_genes[nearest_index]),
  Nearest_gene_end =
    end(reference_genes[nearest_index]),
  Nearest_gene_strand =
    as.character(strand(reference_genes[nearest_index])),
  Distance_to_nearest_gene_bp = nearest_distance
)

# ------------------------------------------------------------
# Find same-strand overlaps
# ------------------------------------------------------------

same_hits <- findOverlaps(
  novel_gr,
  reference_genes,
  ignore.strand = FALSE
)

same_overlap_df <- tibble(
  Novel_ID = novel_gr$Novel_ID[queryHits(same_hits)],
  Gene_ID =
    reference_genes$Reference_gene_ID[subjectHits(same_hits)],
  Gene_symbol =
    reference_genes$Reference_gene_symbol[subjectHits(same_hits)],
  Gene_type =
    reference_genes$Reference_gene_type[subjectHits(same_hits)]
) |>
  group_by(Novel_ID) |>
  summarise(
    Same_strand_overlapping_gene_ID =
      collapse_values(Gene_ID),
    Same_strand_overlapping_gene_symbol =
      collapse_values(Gene_symbol),
    Same_strand_overlapping_gene_type =
      collapse_values(Gene_type),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Find opposite-strand overlaps
# ------------------------------------------------------------

any_hits <- findOverlaps(
  novel_gr,
  reference_genes,
  ignore.strand = TRUE
)

opposite_mask <-
  as.character(strand(novel_gr[queryHits(any_hits)])) !=
  as.character(strand(reference_genes[subjectHits(any_hits)]))

opposite_hits <- any_hits[opposite_mask]

opposite_overlap_df <- tibble(
  Novel_ID = novel_gr$Novel_ID[queryHits(opposite_hits)],
  Gene_ID =
    reference_genes$Reference_gene_ID[subjectHits(opposite_hits)],
  Gene_symbol =
    reference_genes$Reference_gene_symbol[subjectHits(opposite_hits)],
  Gene_type =
    reference_genes$Reference_gene_type[subjectHits(opposite_hits)]
) |>
  group_by(Novel_ID) |>
  summarise(
    Opposite_strand_overlapping_gene_ID =
      collapse_values(Gene_ID),
    Opposite_strand_overlapping_gene_symbol =
      collapse_values(Gene_symbol),
    Opposite_strand_overlapping_gene_type =
      collapse_values(Gene_type),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Combine results and classify transcripts
# ------------------------------------------------------------

annotation <- novel_gene_df |>
  left_join(nearest_df, by = "Novel_ID") |>
  left_join(same_overlap_df, by = "Novel_ID") |>
  left_join(opposite_overlap_df, by = "Novel_ID") |>
  mutate(
    Classification = case_when(
      !is.na(Same_strand_overlapping_gene_symbol) ~
        "Overlaps annotated gene on same strand",
      
      is.na(Same_strand_overlapping_gene_symbol) &
        !is.na(Opposite_strand_overlapping_gene_symbol) ~
        "Antisense overlap with annotated gene",
      
      Distance_to_nearest_gene_bp == 0 ~
        "Overlaps annotated gene",
      
      TRUE ~
        "Intergenic"
    )
  ) |>
  arrange(Chromosome, Start, Novel_ID)

# ------------------------------------------------------------
# Check for requested IDs not found
# ------------------------------------------------------------

found_ids <- unique(annotation$Novel_ID)
missing_ids <- setdiff(novel_ids, found_ids)

# ------------------------------------------------------------
# Export files
# ------------------------------------------------------------

write_tsv(
  annotation,
  paste0(output_prefix, ".tsv"),
  na = "NA"
)

write_csv(
  annotation,
  paste0(output_prefix, ".csv"),
  na = "NA"
)

write_lines(
  missing_ids,
  paste0(output_prefix, "_missing_ids.txt")
)

message("")
message("Annotation complete.")
message("Annotated IDs: ", length(found_ids))
message("Missing IDs: ", length(missing_ids))
message("Output: ", output_prefix, ".tsv")
message("Output: ", output_prefix, ".csv")
message("Missing-ID file: ", output_prefix, "_missing_ids.txt")

