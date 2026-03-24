suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggVennDiagram)
  library(patchwork)
})

# =========================================================
# 1. SETTINGS
# =========================================================

BASE_DIR <- "/home/Azu/Pulpit/data_NGS"

SPECIES_COL <- "group"      # PRD / PTSD
COMPARTMENT_COL <- "type"   # ENV / SILK / EGGS
SAMPLE_ID_COL <- "sample-id"

LEVEL_PRD  <- "PRD"
LEVEL_PTSD <- "PTSD"

LEVEL_ENV  <- "ENV"
LEVEL_SILK <- "SILK"
LEVEL_EGGS <- "EGGS"

TAX_LEVEL <- "family"   # "genus" or "family"

COUNT_THRESHOLD <- 10
PREVALENCE_THRESHOLD <- 0.66

PNG_DPI <- 600
WIDTH_3SET <- 7
HEIGHT_3SET <- 6
WIDTH_2SET <- 6
HEIGHT_2SET <- 5

BASE_TEXT_SIZE <- 12
TITLE_SIZE <- 13

COLOR_PRD_HIGH  <- "#B02A2A"
COLOR_PTSD_HIGH <- "#3137BD"
COLOR_ENV_HIGH  <- "#237A3C"
COLOR_SILK_HIGH <- "#3E455C"
COLOR_EGGS_HIGH <- "#F58E27"

COLOR_LOW <- "#f2f2f2"

# =========================================================
# 2. PATHS
# =========================================================

METADATA_FILE <- file.path(BASE_DIR, "01_Metadata", "metadata_plik.tsv")
FEATURE_TABLE_FILE <- file.path(BASE_DIR, "06_Exports", "deseq2_input", "feature-table.tsv")
TAXONOMY_FILE <- file.path(BASE_DIR, "06_Exports", "deseq2_input", "taxonomy_export", "taxonomy.tsv")

OUT_DIR <- file.path(BASE_DIR, "Core Microbiome", "Venn diagram")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PLOTS_DIR <- file.path(OUT_DIR, "plots")
TABLES_DIR <- file.path(OUT_DIR, "supplementary_tables")
RAW_SETS_DIR <- file.path(OUT_DIR, "raw_sets")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_SETS_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 3. HELPERS
# =========================================================

extract_family <- function(tax) {
  fam <- str_match(tax, "f__([^;]+)")[, 2]
  fam[is.na(fam) | fam == ""] <- "Unclassified_family"
  fam
}

extract_genus <- function(tax) {
  gen <- str_match(tax, "g__([^;]+)")[, 2]
  gen[is.na(gen) | gen == ""] <- "Unclassified_genus"
  gen
}

filter_taxa_for_venn <- function(taxa_vec) {
  taxa_vec <- as.character(taxa_vec)
  taxa_vec <- trimws(taxa_vec)
  taxa_vec <- taxa_vec[!is.na(taxa_vec)]
  taxa_vec <- taxa_vec[taxa_vec != ""]
  taxa_vec <- taxa_vec[!taxa_vec %in% c("Unassigned", "unassigned", "Unknown", "unknown", "NA")]
  taxa_vec <- taxa_vec[!grepl("^d__Bacteria(\\.[0-9]+)?$", taxa_vec, ignore.case = TRUE)]
  taxa_vec <- taxa_vec[!grepl("uncultured", taxa_vec, ignore.case = TRUE)]
  taxa_vec <- taxa_vec[!grepl("^Unknown($|[_[:punct:]])", taxa_vec, ignore.case = TRUE)]
  taxa_vec <- taxa_vec[!grepl("^Subgroup[_-]?[A-Za-z0-9]+$", taxa_vec, ignore.case = TRUE)]
  sort(unique(taxa_vec))
}

save_vector_csv <- function(x, file) {
  out <- data.frame(taxon = sort(unique(x)), stringsAsFactors = FALSE)
  write.csv(out, file, row.names = FALSE)
}

save_two_set_table <- function(set_a, set_b, name_a, name_b, out_prefix) {
  shared <- intersect(set_a, set_b)
  only_a <- setdiff(set_a, set_b)
  only_b <- setdiff(set_b, set_a)
  
  save_vector_csv(shared, paste0(out_prefix, "_shared.csv"))
  save_vector_csv(only_a, paste0(out_prefix, "_only_", name_a, ".csv"))
  save_vector_csv(only_b, paste0(out_prefix, "_only_", name_b, ".csv"))
  
  summary_df <- data.frame(
    category = c("shared", paste0("only_", name_a), paste0("only_", name_b)),
    n = c(length(shared), length(only_a), length(only_b))
  )
  write.csv(summary_df, paste0(out_prefix, "_summary_counts.csv"), row.names = FALSE)
}

save_three_set_table <- function(set_a, set_b, set_c, name_a, name_b, name_c, out_prefix) {
  abc <- Reduce(intersect, list(set_a, set_b, set_c))
  ab_only <- setdiff(intersect(set_a, set_b), set_c)
  ac_only <- setdiff(intersect(set_a, set_c), set_b)
  bc_only <- setdiff(intersect(set_b, set_c), set_a)
  
  a_only <- setdiff(set_a, union(set_b, set_c))
  b_only <- setdiff(set_b, union(set_a, set_c))
  c_only <- setdiff(set_c, union(set_a, set_b))
  
  save_vector_csv(abc,     paste0(out_prefix, "_shared_all_three.csv"))
  save_vector_csv(ab_only, paste0(out_prefix, "_shared_", name_a, "_", name_b, "_only.csv"))
  save_vector_csv(ac_only, paste0(out_prefix, "_shared_", name_a, "_", name_c, "_only.csv"))
  save_vector_csv(bc_only, paste0(out_prefix, "_shared_", name_b, "_", name_c, "_only.csv"))
  save_vector_csv(a_only,  paste0(out_prefix, "_only_", name_a, ".csv"))
  save_vector_csv(b_only,  paste0(out_prefix, "_only_", name_b, ".csv"))
  save_vector_csv(c_only,  paste0(out_prefix, "_only_", name_c, ".csv"))
  
  summary_df <- data.frame(
    category = c(
      "shared_all_three",
      paste0("shared_", name_a, "_", name_b, "_only"),
      paste0("shared_", name_a, "_", name_c, "_only"),
      paste0("shared_", name_b, "_", name_c, "_only"),
      paste0("only_", name_a),
      paste0("only_", name_b),
      paste0("only_", name_c)
    ),
    n = c(
      length(abc),
      length(ab_only),
      length(ac_only),
      length(bc_only),
      length(a_only),
      length(b_only),
      length(c_only)
    )
  )
  
  write.csv(summary_df, paste0(out_prefix, "_summary_counts.csv"), row.names = FALSE)
}

make_venn_plot <- function(sets_list, title_text = NULL,
                           low_color = "#f2f2f2",
                           high_color = "steelblue") {
  p <- ggVennDiagram(
    sets_list,
    label = "count",
    label_alpha = 0
  ) +
    scale_fill_gradient(
      low = low_color,
      high = high_color
    ) +
    theme(
      text = element_text(size = BASE_TEXT_SIZE),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  if (!is.null(title_text)) {
    p <- p + labs(title = title_text)
  }
  
  p
}

save_plot_all_formats <- function(plot_obj, out_prefix, width, height, dpi = 600) {
  ggsave(paste0(out_prefix, ".png"), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(paste0(out_prefix, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(out_prefix, ".svg"), plot_obj, width = width, height = height)
}

get_core_taxa <- function(abund_df, meta_df, species_value, compartment_value,
                          count_threshold = 10, prevalence_threshold = 0.70) {
  
  meta_sub <- meta_df %>%
    filter(.data[[SPECIES_COL]] == species_value,
           .data[[COMPARTMENT_COL]] == compartment_value)
  
  sample_ids <- meta_sub[[SAMPLE_ID_COL]]
  
  if (length(sample_ids) == 0) {
    return(character(0))
  }
  
  df_sub <- abund_df %>%
    select(taxon, all_of(sample_ids))
  
  n_samples <- length(sample_ids)
  min_required <- ceiling(prevalence_threshold * n_samples)
  
  counts_mat <- as.matrix(df_sub[, sample_ids, drop = FALSE])
  
  n_present <- rowSums(counts_mat > count_threshold, na.rm = TRUE)
  keep <- n_present >= min_required
  
  taxa <- df_sub$taxon[keep]
  taxa <- filter_taxa_for_venn(taxa)
  
  taxa
}

# =========================================================
# 4. LOAD DATA
# =========================================================

meta <- read_tsv(METADATA_FILE, show_col_types = FALSE) %>%
  mutate(across(all_of(c(SAMPLE_ID_COL, SPECIES_COL, COMPARTMENT_COL)), as.character))

feature_table <- read_tsv(
  FEATURE_TABLE_FILE,
  skip = 1,
  show_col_types = FALSE
)
colnames(feature_table)[1] <- "FeatureID"

taxonomy <- read_tsv(TAXONOMY_FILE, show_col_types = FALSE)
colnames(taxonomy)[1] <- "FeatureID"
tax_col <- colnames(taxonomy)[2]

taxonomy <- taxonomy %>%
  mutate(
    taxonomy_string = .data[[tax_col]],
    family = extract_family(taxonomy_string),
    genus = extract_genus(taxonomy_string)
  ) %>%
  select(FeatureID, family, genus)

common_samples <- intersect(colnames(feature_table)[-1], meta[[SAMPLE_ID_COL]])

feature_table <- feature_table %>%
  select(FeatureID, all_of(common_samples))

meta <- meta %>%
  filter(.data[[SAMPLE_ID_COL]] %in% common_samples)

feature_long <- feature_table %>%
  pivot_longer(
    cols = -FeatureID,
    names_to = SAMPLE_ID_COL,
    values_to = "count"
  ) %>%
  mutate(count = as.numeric(count)) %>%
  left_join(meta, by = SAMPLE_ID_COL) %>%
  left_join(taxonomy, by = "FeatureID")

if (TAX_LEVEL == "genus") {
  abund_long <- feature_long %>%
    group_by(.data[[SAMPLE_ID_COL]], .data[[SPECIES_COL]], .data[[COMPARTMENT_COL]], genus) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    rename(taxon = genus)
  TAXON_LABEL <- "Genus"
} else if (TAX_LEVEL == "family") {
  abund_long <- feature_long %>%
    group_by(.data[[SAMPLE_ID_COL]], .data[[SPECIES_COL]], .data[[COMPARTMENT_COL]], family) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    rename(taxon = family)
  TAXON_LABEL <- "Family"
} else {
  stop("TAX_LEVEL must be 'genus' or 'family'")
}

abund_wide <- abund_long %>%
  select(all_of(SAMPLE_ID_COL), taxon, count) %>%
  pivot_wider(names_from = SAMPLE_ID_COL, values_from = count, values_fill = 0)

# =========================================================
# 5. BUILD 6 CORE SETS
# =========================================================

PRD_ENV <- get_core_taxa(abund_wide, meta, LEVEL_PRD, LEVEL_ENV,
                         COUNT_THRESHOLD, PREVALENCE_THRESHOLD)
PRD_SILK <- get_core_taxa(abund_wide, meta, LEVEL_PRD, LEVEL_SILK,
                          COUNT_THRESHOLD, PREVALENCE_THRESHOLD)
PRD_EGGS <- get_core_taxa(abund_wide, meta, LEVEL_PRD, LEVEL_EGGS,
                          COUNT_THRESHOLD, PREVALENCE_THRESHOLD)

PTSD_ENV <- get_core_taxa(abund_wide, meta, LEVEL_PTSD, LEVEL_ENV,
                          COUNT_THRESHOLD, PREVALENCE_THRESHOLD)
PTSD_SILK <- get_core_taxa(abund_wide, meta, LEVEL_PTSD, LEVEL_SILK,
                           COUNT_THRESHOLD, PREVALENCE_THRESHOLD)
PTSD_EGGS <- get_core_taxa(abund_wide, meta, LEVEL_PTSD, LEVEL_EGGS,
                           COUNT_THRESHOLD, PREVALENCE_THRESHOLD)

# =========================================================
# 6. SAVE RAW SETS
# =========================================================

save_vector_csv(PRD_ENV,  file.path(RAW_SETS_DIR,  paste0("core_PRD_ENV_", TAX_LEVEL, ".csv")))
save_vector_csv(PRD_SILK, file.path(RAW_SETS_DIR,  paste0("core_PRD_SILK_", TAX_LEVEL, ".csv")))
save_vector_csv(PRD_EGGS, file.path(RAW_SETS_DIR,  paste0("core_PRD_EGGS_", TAX_LEVEL, ".csv")))

save_vector_csv(PTSD_ENV,  file.path(RAW_SETS_DIR, paste0("core_PTSD_ENV_", TAX_LEVEL, ".csv")))
save_vector_csv(PTSD_SILK, file.path(RAW_SETS_DIR, paste0("core_PTSD_SILK_", TAX_LEVEL, ".csv")))
save_vector_csv(PTSD_EGGS, file.path(RAW_SETS_DIR, paste0("core_PTSD_EGGS_", TAX_LEVEL, ".csv")))

# =========================================================
# 7. SAVE TABLES FOR INTERSECTIONS
# =========================================================

save_three_set_table(
  PRD_ENV, PRD_SILK, PRD_EGGS,
  "ENV", "SILK", "EGGS",
  file.path(TABLES_DIR, paste0("core_PRD_ENV_SILK_EGGS_", TAX_LEVEL))
)

save_three_set_table(
  PTSD_ENV, PTSD_SILK, PTSD_EGGS,
  "ENV", "SILK", "EGGS",
  file.path(TABLES_DIR, paste0("core_PTSD_ENV_SILK_EGGS_", TAX_LEVEL))
)

save_two_set_table(
  PRD_ENV, PTSD_ENV,
  "PRD_ENV", "PTSD_ENV",
  file.path(TABLES_DIR, paste0("core_PRD_vs_PTSD_ENV_", TAX_LEVEL))
)

save_two_set_table(
  PRD_SILK, PTSD_SILK,
  "PRD_SILK", "PTSD_SILK",
  file.path(TABLES_DIR, paste0("core_PRD_vs_PTSD_SILK_", TAX_LEVEL))
)

save_two_set_table(
  PRD_EGGS, PTSD_EGGS,
  "PRD_EGGS", "PTSD_EGGS",
  file.path(TABLES_DIR, paste0("core_PRD_vs_PTSD_EGGS_", TAX_LEVEL))
)

# =========================================================
# 8. MAKE PLOTS
# =========================================================

p_prd <- make_venn_plot(
  sets_list = list(ENV = PRD_ENV, SILK = PRD_SILK, EGGS = PRD_EGGS),
  title_text = "PRD",
  low_color = COLOR_LOW,
  high_color = COLOR_PRD_HIGH
)

save_plot_all_formats(
  p_prd,
  file.path(PLOTS_DIR, paste0("Core_Venn_PRD_ENV_SILK_EGGS_", TAX_LEVEL)),
  WIDTH_3SET, HEIGHT_3SET, PNG_DPI
)

p_ptsd <- make_venn_plot(
  sets_list = list(ENV = PTSD_ENV, SILK = PTSD_SILK, EGGS = PTSD_EGGS),
  title_text = "PTSD",
  low_color = COLOR_LOW,
  high_color = COLOR_PTSD_HIGH
)

save_plot_all_formats(
  p_ptsd,
  file.path(PLOTS_DIR, paste0("Core_Venn_PTSD_ENV_SILK_EGGS_", TAX_LEVEL)),
  WIDTH_3SET, HEIGHT_3SET, PNG_DPI
)

p_env <- make_venn_plot(
  sets_list = list(PRD = PRD_ENV, PTSD = PTSD_ENV),
  title_text = "ENV",
  low_color = COLOR_LOW,
  high_color = COLOR_ENV_HIGH
)

save_plot_all_formats(
  p_env,
  file.path(PLOTS_DIR, paste0("Core_Venn_PRD_vs_PTSD_ENV_", TAX_LEVEL)),
  WIDTH_2SET, HEIGHT_2SET, PNG_DPI
)

p_silk <- make_venn_plot(
  sets_list = list(PRD = PRD_SILK, PTSD = PTSD_SILK),
  title_text = "SILK",
  low_color = COLOR_LOW,
  high_color = COLOR_SILK_HIGH
)

save_plot_all_formats(
  p_silk,
  file.path(PLOTS_DIR, paste0("Core_Venn_PRD_vs_PTSD_SILK_", TAX_LEVEL)),
  WIDTH_2SET, HEIGHT_2SET, PNG_DPI
)

p_eggs <- make_venn_plot(
  sets_list = list(PRD = PRD_EGGS, PTSD = PTSD_EGGS),
  title_text = "EGGS",
  low_color = COLOR_LOW,
  high_color = COLOR_EGGS_HIGH
)

save_plot_all_formats(
  p_eggs,
  file.path(PLOTS_DIR, paste0("Core_Venn_PRD_vs_PTSD_EGGS_", TAX_LEVEL)),
  WIDTH_2SET, HEIGHT_2SET, PNG_DPI
)

# =========================================================
# 9. SAVE SIMPLE SUMMARY
# =========================================================

summary_counts <- data.frame(
  set = c("PRD_ENV", "PRD_SILK", "PRD_EGGS", "PTSD_ENV", "PTSD_SILK", "PTSD_EGGS"),
  n_taxa = c(
    length(PRD_ENV), length(PRD_SILK), length(PRD_EGGS),
    length(PTSD_ENV), length(PTSD_SILK), length(PTSD_EGGS)
  )
)

write.csv(
  summary_counts,
  file.path(OUT_DIR, paste0("core_set_sizes_summary_", TAX_LEVEL, ".csv")),
  row.names = FALSE
)

# =========================================================
# 10. SAVE FILTER INFO
# =========================================================

filter_info <- data.frame(
  parameter = c("COUNT_THRESHOLD", "PREVALENCE_THRESHOLD", "rule"),
  value = c(
    COUNT_THRESHOLD,
    PREVALENCE_THRESHOLD,
    "core if count > 10 in a sample and present in >= ceiling(0.70 * n_samples)"
  )
)

write.csv(
  filter_info,
  file.path(OUT_DIR, paste0("core_filter_settings_used_", TAX_LEVEL, ".csv")),
  row.names = FALSE
)

# =========================================================
# 11. COMBINED FIGURES
# =========================================================

combined_3set <- p_prd + p_ptsd +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = paste("Overlap of core bacterial", tolower(TAXON_LABEL), "across compartments")
  ) &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

save_plot_all_formats(
  combined_3set,
  file.path(PLOTS_DIR, paste0("Combined_Core_Venn_3set_PRD_PTSD_", TAX_LEVEL)),
  width = 12, height = 6.5, dpi = PNG_DPI
)

combined_2set <- p_env + p_silk + p_eggs +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = paste("Shared and unique core bacterial", tolower(TAXON_LABEL), "between PRD and PTSD")
  ) &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

save_plot_all_formats(
  combined_2set,
  file.path(PLOTS_DIR, paste0("Combined_Core_Venn_2set_ENV_SILK_EGGS_", TAX_LEVEL)),
  width = 15, height = 5.5, dpi = PNG_DPI
)

combined_all <- (p_prd + p_ptsd) / (p_env + p_silk + p_eggs) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = paste("Overlap of core bacterial", tolower(TAXON_LABEL), "across compartments and species"),
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 14, face = "bold")
  )

save_plot_all_formats(
  combined_all,
  file.path(PLOTS_DIR, paste0("Combined_Core_Venn_All_5plots_", TAX_LEVEL)),
  width = 16, height = 10, dpi = PNG_DPI
)

cat("Done. Core microbiome Venn results saved to:\n", OUT_DIR, "\n")

