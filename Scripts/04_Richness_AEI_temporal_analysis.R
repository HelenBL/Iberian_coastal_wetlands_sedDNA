# ==========================================================
# RICHNESS (18S + COI) AND ANTHROPOGENIC ENRICHMENT INDEX
# Final corrected analysis
# ==========================================================

# Primary temporal summaries use the median within each bin.
# Means and maxima are exported only for sensitivity checks.

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "writexl",
  "purrr", "broom", "zoo", "nlme", "cowplot",
  "sf", "rnaturalearth", "rnaturalearthdata", "tibble",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(purrr)
library(broom)
library(zoo)
library(nlme)
library(cowplot)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# ==========================================================
# 0) PATHS AND PARAMETERS
# ==========================================================

base_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/Data"
base_meta <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/Metadatas"

metadata_file <- file.path(base_meta, "metadata_samples.xlsx")
file_18S <- file.path(base_dir, "All_Peninsula_18S_AbRel.xlsx")
file_COI <- file.path(base_dir, "All_Peninsula_COI_AbRel.xlsx")
xrf_results_file <- file.path(
  base_dir,
  "AEI_general_site_specific_comparison_1600.xlsx"
)

output_dir <- file.path(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/timelines",
  "richness_18S_COI_xrf_AEI_OK_final_NEW"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tax_cols_n <- 12
coi_multiplier <- 10
rolling_width <- 4
minimum_bins_lm <- 5
minimum_bins_adjusted <- 6
minimum_bins_CAR1 <- 8

site_order <- c(
  "DEE", "STO", "CCO", "COL",
  "CSD", "CMB", "HIT", "HIO"
)

short_bin_sites <- c("CCO", "COL", "HIT", "HIO")
bin_label <- "site_specific_bins_20yr_30yr_median"

site_colors <- c(
  "CMB" = "#458B74",
  "CSD" = "#76EEC6",
  "DEE" = "#CDCD00",
  "HIT" = "#CD96CD",
  "HIO" = "#8B0A50",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "COL" = "burlywood1"
)

marker_colors <- c(
  "18S" = "#0072B2",
  "COI" = "#D55E00"
)

event_colors <- c(
  "Richness decreased; AEI increased" = "#D55E00",
  "Richness increased; AEI decreased" = "#0072B2"
)

# ==========================================================
# 1) HELPERS
# ==========================================================

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- gsub(",", ".", as.character(x), fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

rename_COL <- function(x) {
  dplyr::recode(as.character(x), "CLR" = "COL")
}

assign_bin_width <- function(site) {
  dplyr::if_else(site %in% short_bin_sites, 20, 30)
}

safe_z <- function(x) {
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(x_sd) || x_sd == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

safe_spearman <- function(data, x, y, minimum_n = 4) {
  dat <- data %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  
  if (
    nrow(dat) < minimum_n ||
    n_distinct(dat[[x]]) < 2 ||
    n_distinct(dat[[y]]) < 2
  ) {
    return(tibble(
      n = nrow(dat),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  tt <- suppressWarnings(
    cor.test(dat[[x]], dat[[y]], method = "spearman", exact = FALSE)
  )
  
  tibble(
    n = nrow(dat),
    rho = unname(tt$estimate),
    p_value = tt$p.value
  )
}

safe_roll_cor <- function(x) {
  x <- as.data.frame(x)
  ok <- complete.cases(x[, 1], x[, 2])
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok, 1], x[ok, 2], method = "spearman"))
}

save_plot_both <- function(plot, filename_base, width, height, dpi = 600) {
  ggsave(
    file.path(output_dir, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    file.path(output_dir, paste0(filename_base, ".svg")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# ==========================================================
# 2) READ AND STANDARDISE INPUT DATA
# ==========================================================

metadata <- read_excel(metadata_file) %>%
  mutate(
    Site = rename_COL(Site),
    Year = to_num(Year)
  )

data_18S <- read_excel(file_18S)
data_COI <- read_excel(file_COI)

available_AEI_sheets <- readxl::excel_sheets(xrf_results_file)

if (!"AEI_general" %in% available_AEI_sheets) {
  stop(
    "Sheet 'AEI_general' was not found in the AEI workbook. Available sheets: ",
    paste(available_AEI_sheets, collapse = ", ")
  )
}

df_aei <- read_excel(xrf_results_file, sheet = "AEI_general") %>%
  mutate(
    Site = rename_COL(Site),
    Year = to_num(Year)
  )

required_AEI_columns <- c("Site", "Year", "AEI_general")

if (!all(required_AEI_columns %in% names(df_aei))) {
  stop(
    "The AEI_general sheet must contain: ",
    paste(required_AEI_columns, collapse = ", ")
  )
}

# ==========================================================
# 3) CALCULATE EUKARYOTIC ASV RICHNESS
# ==========================================================

calculate_richness <- function(df, marker_name) {
  if (!"domain" %in% names(df)) {
    stop("Column 'domain' is absent from the ", marker_name, " table.")
  }
  
  all_euk <- df %>%
    filter(domain == "Eukaryota")
  
  if (ncol(all_euk) <= tax_cols_n) {
    stop("No sample columns were found in the ", marker_name, " table.")
  }
  
  sample_cols <- names(all_euk)[seq.int(tax_cols_n + 1, ncol(all_euk))]
  
  abundance_matrix <- all_euk %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), to_num)) %>%
    as.matrix()
  
  abundance_matrix[!is.finite(abundance_matrix)] <- 0
  
  tibble(
    sample = sample_cols,
    richness_q0 = colSums(abundance_matrix > 0),
    Marker = marker_name
  )
}

richness_all <- bind_rows(
  calculate_richness(data_18S, "18S"),
  calculate_richness(data_COI, "COI")
)

plot_df_rich <- richness_all %>%
  left_join(metadata, by = "sample") %>%
  mutate(
    Site = rename_COL(Site),
    Year = to_num(Year),
    richness_plot = if_else(
      Marker == "COI",
      richness_q0 * coi_multiplier,
      as.numeric(richness_q0)
    ),
    Marker_plot = if_else(
      Marker == "COI",
      paste0("COI ×", coi_multiplier),
      "18S"
    ),
    Site = factor(Site, levels = site_order)
  ) %>%
  filter(
    !is.na(Site),
    is.finite(Year),
    is.finite(richness_q0)
  ) %>%
  arrange(Site, Marker, Year)

df_aei_plot <- df_aei %>%
  mutate(
    Site = factor(rename_COL(Site), levels = site_order),
    Year = to_num(Year),
    AEI_general = to_num(AEI_general)
  ) %>%
  filter(
    !is.na(Site),
    is.finite(Year),
    is.finite(AEI_general)
  ) %>%
  arrange(Site, Year)

# ==========================================================
# 4) SITE-SPECIFIC TEMPORAL BINNING
# ==========================================================

# The median is the primary bin summary. The mean and maximum
# are retained only for sensitivity analyses and data checking.

rich_binned <- plot_df_rich %>%
  mutate(
    Site = as.character(Site),
    bin_width = assign_bin_width(Site),
    time_bin = floor(Year / bin_width) * bin_width,
    time_bin_end = time_bin + bin_width,
    time_bin_mid = time_bin + bin_width / 2
  ) %>%
  group_by(
    Site, Marker, Marker_plot, bin_width,
    time_bin, time_bin_end, time_bin_mid
  ) %>%
  summarise(
    richness = median(richness_q0, na.rm = TRUE),
    richness_mean = mean(richness_q0, na.rm = TRUE),
    richness_max = max(richness_q0, na.rm = TRUE),
    n_rich = sum(is.finite(richness_q0)),
    .groups = "drop"
  ) %>%
  mutate(
    richness_plot = if_else(
      Marker == "COI",
      richness * coi_multiplier,
      richness
    )
  )

aei_binned <- df_aei_plot %>%
  mutate(
    Site = as.character(Site),
    bin_width = assign_bin_width(Site),
    time_bin = floor(Year / bin_width) * bin_width,
    time_bin_end = time_bin + bin_width,
    time_bin_mid = time_bin + bin_width / 2
  ) %>%
  group_by(
    Site, bin_width,
    time_bin, time_bin_end, time_bin_mid
  ) %>%
  summarise(
    AEI = median(AEI_general, na.rm = TRUE),
    AEI_mean = mean(AEI_general, na.rm = TRUE),
    AEI_max = max(AEI_general, na.rm = TRUE),
    n_AEI = sum(is.finite(AEI_general)),
    .groups = "drop"
  )

df_binned <- rich_binned %>%
  inner_join(
    aei_binned,
    by = c(
      "Site", "bin_width", "time_bin",
      "time_bin_end", "time_bin_mid"
    )
  ) %>%
  filter(is.finite(richness), is.finite(AEI)) %>%
  arrange(Site, Marker, time_bin_mid)

bin_width_check <- df_binned %>%
  count(Site, Marker, bin_width, name = "n_bins") %>%
  arrange(factor(Site, levels = site_order), Marker)

print(bin_width_check)

# ==========================================================
# 5) STANDARDISATION AND CONSECUTIVE CHANGES
# ==========================================================

df_analysis <- df_binned %>%
  group_by(Site, Marker) %>%
  arrange(time_bin_mid, .by_group = TRUE) %>%
  mutate(
    richness_z = safe_z(richness),
    AEI_z = safe_z(AEI),
    Year_z = safe_z(time_bin_mid),
    Year_century =
      (time_bin_mid - mean(time_bin_mid, na.rm = TRUE)) / 100,
    bin_index = row_number(),
    d_richness_z = richness_z - lag(richness_z),
    d_AEI_z = AEI_z - lag(AEI_z),
    interaction_type = case_when(
      d_richness_z < 0 & d_AEI_z > 0 ~
        "Richness decreased; AEI increased",
      d_richness_z > 0 & d_AEI_z < 0 ~
        "Richness increased; AEI decreased",
      d_richness_z > 0 & d_AEI_z > 0 ~ "Both increased",
      d_richness_z < 0 & d_AEI_z < 0 ~ "Both decreased",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

df_opposite_events <- df_analysis %>%
  filter(
    interaction_type %in% c(
      "Richness decreased; AEI increased",
      "Richness increased; AEI decreased"
    )
  )

# ==========================================================
# 6) PRIMARY ASSOCIATIONS: SPEARMAN
# ==========================================================

spearman_levels <- df_analysis %>%
  group_by(Site, Marker) %>%
  group_modify(~ safe_spearman(.x, "richness", "AEI")) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

spearman_changes <- df_analysis %>%
  group_by(Site, Marker) %>%
  group_modify(~ safe_spearman(.x, "d_richness_z", "d_AEI_z")) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

# ==========================================================
# 7) STANDARDISED RICHNESS–AEI MODELS
# ==========================================================

fit_standard_lm <- function(data) {
  dat <- data %>%
    filter(is.finite(richness_z), is.finite(AEI_z))
  
  if (
    nrow(dat) < minimum_bins_lm ||
    n_distinct(dat$richness_z) < 2 ||
    n_distinct(dat$AEI_z) < 2
  ) return(NULL)
  
  lm(richness_z ~ AEI_z, data = dat)
}

fit_time_adjusted_lm <- function(data) {
  dat <- data %>%
    filter(
      is.finite(richness_z),
      is.finite(AEI_z),
      is.finite(Year_z)
    )
  
  if (
    nrow(dat) < minimum_bins_adjusted ||
    n_distinct(dat$richness_z) < 2 ||
    n_distinct(dat$AEI_z) < 2 ||
    n_distinct(dat$Year_z) < 3
  ) return(NULL)
  
  lm(richness_z ~ AEI_z + Year_z, data = dat)
}

fit_CAR1 <- function(data) {
  dat <- data %>%
    filter(
      is.finite(richness_z),
      is.finite(AEI_z),
      is.finite(Year_z),
      is.finite(time_bin_mid)
    ) %>%
    arrange(time_bin_mid)
  
  if (
    nrow(dat) < minimum_bins_CAR1 ||
    n_distinct(dat$richness_z) < 2 ||
    n_distinct(dat$AEI_z) < 2 ||
    n_distinct(dat$time_bin_mid) < minimum_bins_CAR1
  ) return(NULL)
  
  tryCatch(
    nlme::gls(
      richness_z ~ AEI_z + Year_z,
      correlation = nlme::corCAR1(form = ~ time_bin_mid),
      data = dat,
      method = "REML",
      na.action = na.omit,
      control = nlme::glsControl(returnObject = TRUE)
    ),
    error = function(e) NULL
  )
}

tidy_or_empty <- function(model) {
  if (is.null(model)) return(tibble())
  broom::tidy(model, conf.int = TRUE)
}

glance_or_empty <- function(model) {
  
  if (is.null(model)) {
    return(tibble())
  }
  
  broom::glance(model)
}


# Extract results from nlme::gls models

tidy_CAR1_or_empty <- function(
    model,
    conf_level = 0.95
) {
  
  if (
    is.null(model) ||
    !inherits(model, "gls")
  ) {
    return(tibble())
  }
  
  coefficient_table <- tryCatch(
    summary(model)$tTable,
    error = function(e) NULL
  )
  
  if (
    is.null(coefficient_table) ||
    nrow(coefficient_table) == 0
  ) {
    return(tibble())
  }
  
  residual_df <- tryCatch(
    model$dims$N - model$dims$p,
    error = function(e) NA_real_
  )
  
  if (
    length(residual_df) != 1 ||
    !is.finite(residual_df) ||
    residual_df <= 0
  ) {
    residual_df <- Inf
  }
  
  critical_value <- if (
    is.finite(residual_df)
  ) {
    
    qt(
      1 - (1 - conf_level) / 2,
      df = residual_df
    )
    
  } else {
    
    qnorm(
      1 - (1 - conf_level) / 2
    )
  }
  
  tibble(
    term = rownames(coefficient_table),
    
    estimate = unname(
      coefficient_table[, "Value"]
    ),
    
    std.error = unname(
      coefficient_table[, "Std.Error"]
    ),
    
    statistic = unname(
      coefficient_table[, "t-value"]
    ),
    
    p.value = unname(
      coefficient_table[, "p-value"]
    ),
    
    conf.low =
      estimate -
      critical_value * std.error,
    
    conf.high =
      estimate +
      critical_value * std.error,
    
    df = residual_df
  )
}

models_site_marker <- df_analysis %>%
  group_by(Site, Marker) %>%
  nest() %>%
  mutate(
    model_lm = map(data, fit_standard_lm),
    model_adjusted = map(data, fit_time_adjusted_lm),
    model_CAR1 = map(data, fit_CAR1)
  )

lm_AEI_effect <- models_site_marker %>%
  transmute(Site, Marker, result = map(model_lm, tidy_or_empty)) %>%
  unnest(result) %>%
  filter(term == "AEI_z") %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))

lm_model_fit <- models_site_marker %>%
  transmute(Site, Marker, result = map(model_lm, glance_or_empty)) %>%
  unnest(result)

time_adjusted_AEI_effect <- models_site_marker %>%
  transmute(Site, Marker, result = map(model_adjusted, tidy_or_empty)) %>%
  unnest(result) %>%
  filter(term == "AEI_z") %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))

time_adjusted_model_fit <- models_site_marker %>%
  transmute(Site, Marker, result = map(model_adjusted, glance_or_empty)) %>%
  unnest(result)

CAR1_AEI_effect <- models_site_marker %>%
  transmute(
    Site,
    Marker,
    
    result = map(
      model_CAR1,
      tidy_CAR1_or_empty
    )
  ) %>%
  unnest(result) %>%
  filter(
    term == "AEI_z"
  ) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(
      p.value,
      method = "BH"
    )
  )

# ==========================================================
# 8) TEMPORAL TRENDS (SEPARATE FROM RICHNESS–AEI COUPLING)
# ==========================================================

fit_richness_time <- function(data) {
  dat <- data %>%
    filter(is.finite(richness_z), is.finite(Year_century))
  
  if (
    nrow(dat) < minimum_bins_lm ||
    n_distinct(dat$richness_z) < 2 ||
    n_distinct(dat$Year_century) < 3
  ) return(NULL)
  
  lm(richness_z ~ Year_century, data = dat)
}

richness_temporal_trends <- df_analysis %>%
  group_by(Site, Marker) %>%
  nest() %>%
  mutate(
    model = map(data, fit_richness_time),
    result = map(model, tidy_or_empty)
  ) %>%
  select(Site, Marker, result) %>%
  unnest(result) %>%
  filter(term == "Year_century") %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))

AEI_temporal_data <- df_analysis %>%
  distinct(Site, time_bin_mid, AEI, AEI_z) %>%
  group_by(Site) %>%
  arrange(time_bin_mid, .by_group = TRUE) %>%
  mutate(
    Year_century =
      (time_bin_mid - mean(time_bin_mid, na.rm = TRUE)) / 100
  ) %>%
  ungroup()

fit_AEI_time <- function(data) {
  dat <- data %>%
    filter(is.finite(AEI_z), is.finite(Year_century))
  
  if (
    nrow(dat) < minimum_bins_lm ||
    n_distinct(dat$AEI_z) < 2 ||
    n_distinct(dat$Year_century) < 3
  ) return(NULL)
  
  lm(AEI_z ~ Year_century, data = dat)
}

AEI_temporal_trends <- AEI_temporal_data %>%
  group_by(Site) %>%
  nest() %>%
  mutate(
    model = map(data, fit_AEI_time),
    result = map(model, tidy_or_empty)
  ) %>%
  select(Site, result) %>%
  unnest(result) %>%
  filter(term == "Year_century") %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))

# ==========================================================
# 9) DIRECTIONAL TRANSITIONS (DESCRIPTIVE ONLY)
# ==========================================================

transition_summary <- df_analysis %>%
  filter(!is.na(interaction_type)) %>%
  count(Site, Marker, interaction_type, name = "n_transitions") %>%
  group_by(Site, Marker) %>%
  mutate(
    n_total = sum(n_transitions),
    proportion = n_transitions / n_total,
    percentage = 100 * proportion
  ) %>%
  ungroup()

# ==========================================================
# 10) ROLLING CORRELATION (EXPLORATORY ONLY)
# ==========================================================

df_rolling <- df_analysis %>%
  arrange(Site, Marker, time_bin_mid) %>%
  group_by(Site, Marker) %>%
  mutate(
    rho_roll = zoo::rollapply(
      cbind(richness, AEI),
      width = rolling_width,
      FUN = safe_roll_cor,
      by.column = FALSE,
      fill = NA_real_,
      align = "center"
    )
  ) %>%
  ungroup()

p_rolling <- ggplot(
  df_rolling,
  aes(
    x = rho_roll,
    y = time_bin_mid,
    color = Marker,
    group = Marker
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_path(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Site, nrow = 1, scales = "free_y") +
  scale_color_manual(values = marker_colors) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Exploratory rolling richness–AEI association",
    subtitle = paste0(
      "Four-bin rolling Spearman correlations; ",
      "20-year bins for CCO, COL, HIT and HIO and 30-year bins elsewhere"
    ),
    x = "Rolling Spearman rho",
    y = "Bin midpoint (CE)",
    color = "Marker"
  )

save_plot_both(
  p_rolling,
  paste0("Fig.SX_rolling_richness_AEI_", bin_label),
  width = 12,
  height = 16
)

# ==========================================================
# 11) MANUSCRIPT STATISTICAL TABLE
# ==========================================================

table_richness_AEI <- spearman_levels %>%
  rename(
    n_bins = n,
    Spearman_rho = rho,
    Spearman_p = p_value,
    Spearman_q = p_adj_BH
  ) %>%
  left_join(
    lm_AEI_effect %>%
      transmute(
        Site, Marker,
        standardized_slope = estimate,
        slope_CI_low = conf.low,
        slope_CI_high = conf.high,
        slope_p = p.value,
        slope_q = p_adj_BH
      ),
    by = c("Site", "Marker")
  ) %>%
  left_join(
    time_adjusted_AEI_effect %>%
      transmute(
        Site, Marker,
        adjusted_AEI_slope = estimate,
        adjusted_CI_low = conf.low,
        adjusted_CI_high = conf.high,
        adjusted_p = p.value,
        adjusted_q = p_adj_BH
      ),
    by = c("Site", "Marker")
  ) %>%
  left_join(
    CAR1_AEI_effect %>%
      transmute(
        Site, Marker,
        CAR1_AEI_slope = estimate,
        CAR1_CI_low = conf.low,
        CAR1_CI_high = conf.high,
        CAR1_p = p.value,
        CAR1_q = p_adj_BH
      ),
    by = c("Site", "Marker")
  ) %>%
  left_join(
    spearman_changes %>%
      transmute(
        Site, Marker,
        change_rho = rho,
        change_p = p_value,
        change_q = p_adj_BH
      ),
    by = c("Site", "Marker")
  ) %>%
  left_join(
    bin_width_check,
    by = c("Site", "Marker", "n_bins")
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4)),
    Site = factor(Site, levels = site_order)
  ) %>%
  arrange(Site, Marker)

# ==========================================================
# 12) FOREST PLOT OF STANDARDISED ASSOCIATIONS
# ==========================================================

p_AEI_effects <- lm_AEI_effect %>%
  mutate(
    Site = factor(Site, levels = rev(site_order)),
    BH_result = if_else(
      is.finite(p_adj_BH) & p_adj_BH < 0.05,
      "q < 0.05",
      "q ≥ 0.05"
    )
  ) %>%
  ggplot(
    aes(
      x = estimate,
      y = Site,
      color = Marker,
      shape = BH_result
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.18,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = marker_colors) +
  scale_shape_manual(values = c("q ≥ 0.05" = 16, "q < 0.05" = 17)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Standardised richness response per 1-SD increase in AEI",
    y = NULL,
    color = "Marker",
    shape = "BH-adjusted result"
  )

save_plot_both(
  p_AEI_effects,
  "Fig.SX_site_specific_richness_AEI_effects",
  width = 8,
  height = 6
)

# ==========================================================
# 13) BINNED TIMELINE WITH OPPOSITE-DIRECTION EVENTS
# ==========================================================

site_scale <- df_analysis %>%
  group_by(Site) %>%
  summarise(
    scale_factor = max(richness_plot, na.rm = TRUE) /
      max(AEI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scale_factor = if_else(
      is.finite(scale_factor) & scale_factor > 0,
      scale_factor,
      1
    )
  )

df_timeline <- df_analysis %>%
  left_join(site_scale, by = "Site")

df_AEI_line <- df_timeline %>%
  distinct(Site, time_bin_mid, AEI, scale_factor)

p_trend_events <- ggplot() +
  geom_line(
    data = df_timeline,
    aes(
      x = time_bin_mid,
      y = richness_plot,
      color = Site,
      linetype = Marker_plot,
      group = interaction(Site, Marker)
    ),
    linewidth = 1.2
  ) +
  geom_line(
    data = df_AEI_line,
    aes(
      x = time_bin_mid,
      y = AEI * scale_factor,
      group = Site
    ),
    color = "black",
    linewidth = 1,
    alpha = 0.5
  ) +
  geom_point(
    data = df_timeline %>%
      filter(interaction_type %in% names(event_colors)),
    aes(
      x = time_bin_mid,
      y = richness_plot,
      fill = interaction_type
    ),
    shape = 21,
    color = "black",
    size = 3.5,
    stroke = 0.8
  ) +
  facet_wrap(~ factor(Site, levels = site_order), scales = "free") +
  scale_color_manual(values = site_colors, drop = FALSE) +
  scale_linetype_manual(
    values = c(
      "18S" = "solid",
      "COI ×10" = "dashed"
    )
  ) +
  scale_fill_manual(values = event_colors) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Richness and Anthropogenic Enrichment Index",
    subtitle = paste0(
      "Bin medians: 20 years for CCO, COL, HIT and HIO; ",
      "30 years for remaining sites. COI ×10 for display only."
    ),
    x = "Calendar year (bin midpoint)",
    y = "ASV richness",
    color = "Site",
    linetype = "Marker",
    fill = "Opposite-direction change"
  )

save_plot_both(
  p_trend_events,
  paste0("richness_18S_COI_AEI_opposite_trends_", bin_label),
  width = 23,
  height = 20
)

# ==========================================================
# 14) FIGURE 1: MAP AND RAW SITE TIMELINES
# ==========================================================

site_coords <- tibble(
  Site = site_order,
  Site_full = c(
    "Ebro Delta Encanyissada",
    "Santoña",
    "Cicero",
    "Colindres",
    "Cadiz Salina Dolores",
    "Cadiz Militar Base",
    "Huelva Isla Saltes Tinto River",
    "Huelva Isla Saltes Ocean"
  ),
  Latitude = c(
    40.653111, 43.4463765, 43.4227010, 43.4002387,
    36.4648962, 36.499135, 37.181788, 37.18080
  ),
  Longitude = c(
    0.685250, -3.4662690, -3.4831056, -3.4552713,
    -6.2434993, -6.193053, -6.937831, -6.938724
  )
)

site_coords_sf <- st_as_sf(
  site_coords,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

iberian_peninsula <- ne_countries(
  scale = "medium",
  country = c("Spain", "Portugal"),
  returnclass = "sf"
)

make_site_panel <- function(site_id) {
  rich_site <- plot_df_rich %>%
    filter(as.character(Site) == site_id)
  
  aei_site <- df_aei_plot %>%
    filter(as.character(Site) == site_id)
 
  event_site <- df_analysis %>%
    filter(
      Site == site_id,
      interaction_type %in% names(event_colors)
    ) %>%
    select(
      Site,
      Marker,
      time_bin,
      time_bin_mid,
      interaction_type
    )
  
  
  # Interpolate richness at the midpoint of each event bin
  
  interpolate_richness <- function(
    marker_value,
    event_year
  ) {
    
    marker_data <- rich_site %>%
      filter(
        Marker == marker_value,
        is.finite(Year),
        is.finite(richness_plot)
      ) %>%
      arrange(Year) %>%
      distinct(
        Year,
        .keep_all = TRUE
      )
    
    if (
      nrow(marker_data) < 2 ||
      event_year < min(marker_data$Year) ||
      event_year > max(marker_data$Year)
    ) {
      return(NA_real_)
    }
    
    stats::approx(
      x = marker_data$Year,
      y = marker_data$richness_plot,
      xout = event_year,
      method = "linear",
      rule = 1,
      ties = mean
    )$y
  }
  
  
  event_site <- event_site %>%
    rowwise() %>%
    mutate(
      y_event = interpolate_richness(
        marker_value = Marker,
        event_year = time_bin_mid
      )
    ) %>%
    ungroup() %>%
    filter(
      is.finite(time_bin_mid),
      is.finite(y_event)
    )
  
  max_rich <- max(rich_site$richness_plot, na.rm = TRUE)
  max_AEI <- max(aei_site$AEI_general, na.rm = TRUE)
  sf_site <- max_rich / max_AEI
  
  if (!is.finite(sf_site) || sf_site <= 0) sf_site <- 1
  
  site_col <- unname(site_colors[[site_id]])
  site_name <- site_coords$Site_full[site_coords$Site == site_id]
  
  ggplot() +
    geom_line(
      data = aei_site,
      aes(x = Year, y = AEI_general * sf_site),
      color = "black",
      linewidth = 1,
      alpha = 0.5
    ) +
    geom_line(
      data = rich_site,
      aes(
        x = Year,
        y = richness_plot,
        linetype = Marker_plot,
        group = Marker
      ),
      color = site_col,
      linewidth = 1.1
    ) +
    geom_point(
      data = event_site,
      aes(
        x = time_bin_mid,
        y = y_event,
        fill = interaction_type
      ),
      shape = 21,
      color = "black",
      size = 3.2,
      stroke = 0.8,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    scale_linetype_manual(
      values = c("18S" = "solid", "COI ×10" = "dashed"),
      drop = FALSE
    ) +
    scale_fill_manual(values = event_colors, drop = FALSE) +
    scale_y_continuous(
      name = "ASV richness",
      sec.axis = sec_axis(~ . / sf_site, name = "AEI")
    ) +
    labs(title = paste0(site_id, " - ", site_name), x = NULL) +
    theme_classic(base_size = 11, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.background = element_rect(
        fill = "white",
        color = site_col,
        linewidth = 2.2
      ),
      axis.text = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 17, face = "bold"),
      axis.title.y.right = element_text(size = 17, face = "bold"),
      plot.title = element_text(
        size = 22,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

site_panels <- setNames(
  lapply(site_order, make_site_panel),
  site_order
)

p_map <- ggplot() +
  geom_sf(
    data = iberian_peninsula,
    fill = "grey96",
    color = "grey45",
    linewidth = 0.35
  ) +
  geom_sf(
    data = site_coords_sf,
    aes(color = Site),
    size = 3
  ) +
  scale_color_manual(values = site_colors) +
  coord_sf(
    xlim = c(-10.0, 3.5),
    ylim = c(35.5, 44.2),
    expand = FALSE
  ) +
  theme_void() +
  theme(legend.position = "none")

panel_pos <- tibble(
  Site = c("STO", "CCO", "COL", "CSD", "CMB", "HIT", "HIO", "DEE"),
  x = c(0.05, 0.35, 0.69, 0.02, 0.34, 0.02, 0.02, 0.72),
  y = c(0.72, 0.72, 0.72, 0.08, 0.08, 0.48, 0.30, 0.38),
  w = c(0.27, 0.27, 0.27, 0.25, 0.27, 0.25, 0.25, 0.26),
  h = c(0.19, 0.19, 0.19, 0.17, 0.17, 0.17, 0.17, 0.18)
)

map_box <- list(x = 0.30, y = 0.24, w = 0.40, h = 0.45)

lonlat_to_canvas <- function(lon, lat) {
  xlim <- c(-10.0, 3.5)
  ylim <- c(35.5, 44.2)
  
  tibble(
    x = map_box$x + ((lon - xlim[1]) / diff(xlim)) * map_box$w,
    y = map_box$y + ((lat - ylim[1]) / diff(ylim)) * map_box$h
  )
}

site_canvas <- site_coords %>%
  rowwise() %>%
  mutate(canvas = list(lonlat_to_canvas(Longitude, Latitude))) %>%
  unnest(canvas) %>%
  ungroup()

arrow_df <- panel_pos %>%
  left_join(site_canvas, by = "Site", suffix = c("_panel", "_site")) %>%
  rowwise() %>%
  mutate(
    panel_center_x = x_panel + w / 2,
    panel_center_y = y_panel + h / 2,
    dx = panel_center_x - x_site,
    dy = panel_center_y - y_site,
    use_vertical_edge = abs(dx / w) > abs(dy / h),
    xend = if_else(
      use_vertical_edge,
      if_else(dx > 0, x_panel, x_panel + w),
      panel_center_x
    ),
    yend = if_else(
      use_vertical_edge,
      panel_center_y,
      if_else(dy > 0, y_panel, y_panel + h)
    )
  ) %>%
  ungroup()

legend_plot <- ggplot() +
  geom_line(
    data = tibble(
      x = rep(c(1, 2), 3),
      y = rep(c(1, 1), 3),
      type = rep(
        c("18S", "COI ×10", "Anthropogenic Enrichment Index (AEI)"),
        each = 2
      )
    ),
    aes(x = x, y = y, linetype = type, group = type),
    linewidth = 1.5
  ) +
  geom_point(
    data = tibble(
      x = c(1, 2),
      y = c(1, 1),
      interaction_type = names(event_colors)
    ),
    aes(x = x, y = y, fill = interaction_type),
    shape = 21,
    size = 4,
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      "18S" = "solid",
      "COI ×10" = "dashed",
      "Anthropogenic Enrichment Index (AEI)" = "dotted"
    ),
    name = NULL
  ) +
  scale_fill_manual(values = event_colors, name = "Opposite-direction change") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 13, family = "Arial"),
    legend.title = element_text(size = 13, face = "bold", family = "Arial"),
    legend.key.width = grid::unit(2.2, "cm")
  )

legend_grob <- cowplot::get_legend(legend_plot)

final_map_figure <- ggdraw() +
  draw_label(
    "18S and COI richness vs Anthropogenic Enrichment Index",
    x = 0.5,
    y = 0.985,
    size = 22,
    fontface = "bold"
  ) +
  draw_label(
    paste0(
      "COI richness is multiplied by ", coi_multiplier,
      " for visualisation only; analyses use untransformed richness."
    ),
    x = 0.5,
    y = 0.955,
    size = 17
  ) +
  draw_plot(
    p_map,
    x = map_box$x,
    y = map_box$y,
    width = map_box$w,
    height = map_box$h
  )

for (i in seq_len(nrow(panel_pos))) {
  ss <- panel_pos$Site[i]
  final_map_figure <- final_map_figure +
    draw_plot(
      site_panels[[ss]],
      x = panel_pos$x[i],
      y = panel_pos$y[i],
      width = panel_pos$w[i],
      height = panel_pos$h[i]
    )
}

for (i in seq_len(nrow(arrow_df))) {
  ss <- arrow_df$Site[i]
  final_map_figure <- final_map_figure +
    draw_line(
      x = c(arrow_df$x_site[i], arrow_df$xend[i]),
      y = c(arrow_df$y_site[i], arrow_df$yend[i]),
      color = unname(site_colors[[ss]]),
      linewidth = 0.8
    )
}

final_map_figure <- final_map_figure +
  draw_grob(
    legend_grob,
    x = 0.67,
    y = 0.005,
    width = 0.31,
    height = 0.13
  )

save_plot_both(
  final_map_figure,
  "Fig.1_Spain_map_richness_18S_COI_AEI_TRENDS",
  width = 28,
  height = 21
)

# ==========================================================
# 15) EXPORT TABLES
# ==========================================================

writexl::write_xlsx(
  list(
    bin_width_check = bin_width_check,
    richness_raw = plot_df_rich,
    AEI_raw = df_aei_plot,
    richness_binned = rich_binned,
    AEI_binned = aei_binned,
    analysis_data = df_analysis,
    Spearman_levels = spearman_levels,
    Spearman_changes = spearman_changes,
    LM_AEI_effect = lm_AEI_effect,
    LM_model_fit = lm_model_fit,
    Time_adjusted_AEI = time_adjusted_AEI_effect,
    Time_adjusted_fit = time_adjusted_model_fit,
    CAR1_AEI_effect = CAR1_AEI_effect,
    Richness_time_trends = richness_temporal_trends,
    AEI_time_trends = AEI_temporal_trends,
    Directional_transitions = transition_summary,
    Rolling_exploratory = df_rolling,
    Manuscript_table = table_richness_AEI
  ),
  file.path(
    output_dir,
    paste0("Richness_AEI_final_statistics_", bin_label, ".xlsx")
  )
)

cat(
  "DONE\n",
  "Primary bin summary: median.\n",
  "CCO, COL, HIT and HIO: 20-year bins.\n",
  "DEE, STO, CSD and CMB: 30-year bins.\n",
  "Results and figures saved in: ", output_dir, "\n",
  sep = ""
)
