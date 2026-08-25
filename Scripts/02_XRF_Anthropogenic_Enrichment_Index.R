# ============================================================
# XRF CLR-based pollution analysis
# Reads dataset_clean from ALL_PENINSULA_XRF_DATASETS.xlsx
# Uses clr-transformed cleaned elemental data + Euclidean distance
# ============================================================

# =========================
# 0) Packages
# =========================

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "vegan",
  "pheatmap", "writexl", "purrr", "tibble", "viridis"
)

to_install <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(pheatmap)
library(writexl)
library(purrr)
library(tibble)
library(viridis)

# =========================
# 1) Paths
# =========================

input_datasets_file <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/Data/ALL_PENINSULA_XRF_DATASETS.xlsx"

output_stats_file <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/XRF/XRF_CLR_pollution_importance_tables.xlsx"

plot_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/XRF/plots_xrf_clr_NEW"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2) Read clean dataset
# =========================

dataset_clean <- readxl::read_excel(input_datasets_file, sheet = "dataset_clean") |>
  mutate(Site = as.character(Site))

# =========================
# 3) Elements
# =========================

elements_pollution <- c("Pb", "Zn", "Cu", "Sb", "Cd", "Hg", "Sn", "Cr", "Ni", "As", "P")

clean_cols <- paste0(elements_pollution, "_clean")
clean_cols <- clean_cols[clean_cols %in% names(dataset_clean)]

elements_available <- gsub("_clean$", "", clean_cols)

message("Elements used: ", paste(elements_available, collapse = ", "))

# =========================
# 4) Helper functions
# =========================

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

replace_missing_nonpositive <- function(x) {
  pos <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(pos) == 0) return(rep(NA_real_, length(x)))
  
  repl <- min(pos) / 2
  
  x2 <- x
  x2[is.na(x2) | !is.finite(x2) | x2 <= 0] <- repl
  x2
}

clr_matrix <- function(df_elements) {
  X <- as.data.frame(lapply(df_elements, replace_missing_nonpositive))
  M <- as.matrix(X)
  
  clr <- t(apply(M, 1, function(row) {
    gm <- exp(mean(log(row)))
    log(row / gm)
  }))
  
  clr <- as.data.frame(clr)
  names(clr) <- gsub("_clean$", "_clr", names(df_elements))
  clr
}

add_period_bin <- function(df) {
  df |>
    group_by(Site) |>
    mutate(
      site_has_year = sum(!is.na(Year)) > 0,
      Century = ifelse(
        site_has_year & !is.na(Year),
        paste0(floor(Year / 100) * 100, "s"),
        NA_character_
      ),
      Depth_bin = ifelse(
        !site_has_year & !is.na(Depth),
        paste0(
          floor(Depth / 20) * 20,
          "-",
          floor(Depth / 20) * 20 + 20,
          " cm"
        ),
        NA_character_
      ),
      Period_bin = case_when(
        site_has_year ~ Century,
        !site_has_year ~ Depth_bin,
        TRUE ~ NA_character_
      )
    ) |>
    ungroup()
}

# =========================
# 5) CLR transformation
# =========================

id_cols <- c("Site", "Depth", "Year", "X_plot", "X_source")
id_cols <- id_cols[id_cols %in% names(dataset_clean)]

clr_df <- clr_matrix(dataset_clean[, clean_cols, drop = FALSE])

dataset_clr <- bind_cols(
  dataset_clean |> select(all_of(id_cols)),
  clr_df
) |>
  add_period_bin() |>
  filter(!is.na(Period_bin))

clr_cols <- names(clr_df)

# =========================
# 6) Long format for importance
# =========================

pollution_long_clr <- dataset_clr |>
  select(Site, Depth, Year, Period_bin, all_of(clr_cols)) |>
  pivot_longer(
    cols = all_of(clr_cols),
    names_to = "Element",
    values_to = "Value"
  ) |>
  mutate(Element = gsub("_clr$", "", Element))

# =========================
# 7) Robust baseline
# =========================

pollution_long_clr <- pollution_long_clr |>
  group_by(Site) |>
  mutate(
    has_year = sum(!is.na(Year)) > 0,
    depth_cutoff = quantile(Depth, probs = 0.80, na.rm = TRUE),
    baseline_flag_initial = case_when(
      has_year & !is.na(Year) ~ Year < 1600,
      !has_year & !is.na(Depth) ~ Depth >= depth_cutoff,
      TRUE ~ FALSE
    )
  ) |>
  ungroup()

baseline_stats_clr <- pollution_long_clr |>
  group_by(Site, Element) |>
  summarise(
    n_total = sum(!is.na(Value)),
    n_baseline_initial = sum(baseline_flag_initial & !is.na(Value)),
    baseline_median_initial = median(Value[baseline_flag_initial], na.rm = TRUE),
    baseline_mad_initial = mad(Value[baseline_flag_initial], na.rm = TRUE, constant = 1),
    all_median = median(Value, na.rm = TRUE),
    all_mad = mad(Value, na.rm = TRUE, constant = 1),
    .groups = "drop"
  ) |>
  mutate(
    baseline_median = case_when(
      n_baseline_initial >= 3 & is.finite(baseline_median_initial) ~ baseline_median_initial,
      TRUE ~ all_median
    ),
    baseline_mad = case_when(
      n_baseline_initial >= 3 & is.finite(baseline_mad_initial) & baseline_mad_initial > 0 ~ baseline_mad_initial,
      is.finite(all_mad) & all_mad > 0 ~ all_mad,
      TRUE ~ NA_real_
    )
  ) |>
  select(Site, Element, n_total, n_baseline_initial, baseline_median, baseline_mad)

pollution_scores_clr <- pollution_long_clr |>
  left_join(baseline_stats_clr, by = c("Site", "Element")) |>
  mutate(
    robust_z = if_else(
      !is.na(Value) & !is.na(baseline_median) & !is.na(baseline_mad) & baseline_mad > 0,
      (Value - baseline_median) / (1.4826 * baseline_mad),
      NA_real_
    ),
    robust_z_pos = pmax(robust_z, 0, na.rm = FALSE)
  )


# ============================================================
# Transparent description of site-specific tracer enrichment
# Run after `pollution_scores_clr` has been created.
#
# This block does NOT modify PI_general. It replaces the composite
# importance score with two directly interpretable quantities:
#   1) magnitude: mean positive robust z-score;
#   2) recurrence: percentage of observations with robust z > 3.
# ============================================================

required_objects <- c("pollution_scores_clr", "plot_dir")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("Run the main XRF script first. Missing object(s): ",
       paste(missing_objects, collapse = ", "))
}

required_pkgs_transparent <- c("dplyr", "tidyr", "ggplot2", "writexl", "tidytext")
missing_pkgs_transparent <- required_pkgs_transparent[
  !vapply(required_pkgs_transparent, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs_transparent) > 0) {
  stop("Install the following package(s): ",
       paste(missing_pkgs_transparent, collapse = ", "))
}

# Operational definition used only for the complementary site-specific
# analysis. Change this value for sensitivity analyses (e.g. 1, 2, or 3).
strong_z_threshold <- 3
min_strong_events <- 2L

# -------------------------
# 1) Site x element summary
# -------------------------
tracer_enrichment_site <- pollution_scores_clr |>
  dplyr::group_by(Site, Element) |>
  dplyr::summarise(
    n_observations = sum(is.finite(robust_z)),
    n_positive = sum(robust_z > 0, na.rm = TRUE),
    mean_positive_z = ifelse(
      n_positive > 0,
      mean(robust_z[robust_z > 0], na.rm = TRUE),
      0
    ),
    median_positive_z = ifelse(
      n_positive > 0,
      median(robust_z[robust_z > 0], na.rm = TRUE),
      0
    ),
    n_strong_events = sum(robust_z > strong_z_threshold, na.rm = TRUE),
    frequency_strong_pct = ifelse(
      n_observations > 0,
      100 * n_strong_events / n_observations,
      NA_real_
    ),
    max_z = ifelse(
      n_observations > 0,
      max(robust_z, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    selected_site_specific = n_strong_events >= min_strong_events,
    selection_rule = paste0(
      "robust z > ", strong_z_threshold,
      " in at least ", min_strong_events, " samples"
    )
  ) |>
  dplyr::arrange(Site, dplyr::desc(mean_positive_z),
                 dplyr::desc(frequency_strong_pct), Element)

# Order elements within each site using the magnitude metric. This is done
# before pivot_longer(), exactly as in the original working version, so both
# metric panels inherit the same element order.
tracer_enrichment_plot_data <- tracer_enrichment_site |>
  dplyr::group_by(Site) |>
  dplyr::arrange(mean_positive_z, Element, .by_group = TRUE) |>
  dplyr::mutate(
    # A unique integer rank avoids unstable ordering when two or more elements
    # have the same magnitude (often zero). Higher ranks appear at the top
    # after coord_flip().
    Magnitude_rank = dplyr::row_number(),
    Element_order = tidytext::reorder_within(Element, Magnitude_rank, Site)
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    Site, Element, Element_order, Magnitude_rank, selected_site_specific,
    mean_positive_z, frequency_strong_pct
  ) |>
  tidyr::pivot_longer(
    cols = c(mean_positive_z, frequency_strong_pct),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("mean_positive_z", "frequency_strong_pct"),
      labels = c(
        "Magnitude: mean positive robust z-score",
        "Recurrence: observations with robust z > 3 (%)"
      )
    ),
    Site = factor(
      Site,
      levels = c("DEE", "STO", "CCO", "CLR","CSD", "CMB", "HIT", "HIO")
    )
  )

# Colour-blind-friendly qualitative palette. The first eight colours are from
# the Okabe-Ito palette; the additional dark colours extend it for 10-11
# tracers while retaining strong contrast on a white background.
element_colours_cb <- c(
  Pb = "#0072B2", # blue
  Zn = "#009E73", # bluish green
  Cu = "#D55E00", # vermillion
  Sb = "#CC79A7", # reddish purple
  Cd = "#56B4E9", # sky blue
  Hg = "#E69F00", # orange
  Sn = "#F0E442", # yellow
  Cr = "#000000", # black
  Ni = "#6F4E7C", # dark purple
  As = "#8C6D31", # dark ochre
  P  = "#2F6B3C"  # dark green
)

# Add an invisible value of 10 only to the magnitude facets. This fixes the
# displayed magnitude range at approximately 0-10 without changing any bar or
# any underlying calculation.
magnitude_width_spacer <- tracer_enrichment_plot_data |>
  dplyr::filter(Metric == "Magnitude: mean positive robust z-score") |>
  dplyr::group_by(Site) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::mutate(Value = 10)

p_tracer_enrichment <- ggplot2::ggplot(
  tracer_enrichment_plot_data,
  ggplot2::aes(
    x = Element_order,
    y = Value,
    fill = Element,
    alpha = selected_site_specific
  )
) +
  ggplot2::geom_blank(
    data = magnitude_width_spacer,
    ggplot2::aes(x = Element_order, y = Value),
    inherit.aes = FALSE
  ) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::coord_flip() +
  # Panel widths are assigned explicitly after building the plot (1:2), so
  # they are independent of the numerical ranges of the two metrics.
  ggplot2::facet_grid(Site ~ Metric, scales = "free", space = "fixed") +
  tidytext::scale_x_reordered() +
  ggplot2::scale_fill_manual(
    values = element_colours_cb,
    guide = "none"
  ) +
  ggplot2::scale_alpha_manual(
    values = c(`FALSE` = 0.35, `TRUE` = 1),
    breaks = c(FALSE, TRUE),
    labels = c("Not retained", "Retained"),
    name = "Complementary\nsite-specific rule"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Magnitude and recurrence of candidate tracer enrichment by site",
    subtitle = paste0(
      "Full colour: robust z > ", strong_z_threshold,
      " in at least ", min_strong_events,
      " samples; faded colour: criterion not met"
    )
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 13),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 17, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 11, face = "bold")
  )

# Convert to a graphical table and assign one third of the available panel
# width to magnitude and two thirds to recurrence. This widens the first facet
# without extending its x-axis beyond 10.
p_tracer_enrichment_grob <- ggplot2::ggplotGrob(p_tracer_enrichment)
panel_columns <- sort(unique(
  p_tracer_enrichment_grob$layout$l[
    grepl("^panel", p_tracer_enrichment_grob$layout$name)
  ]
))

if (length(panel_columns) != 2) {
  stop("Expected two metric panel columns, but found ", length(panel_columns))
}

p_tracer_enrichment_grob$widths[panel_columns] <-
  grid::unit(c(1, 2), "null")

grid::grid.newpage()
grid::grid.draw(p_tracer_enrichment_grob)

ggplot2::ggsave(
  filename = file.path(plot_dir, "Fig_S5X_tracer_enrichment_magnitude_recurrence.png"),
  plot = p_tracer_enrichment_grob,
  width = 14,
  height = 15,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(plot_dir, "Fig_S5X_tracer_enrichment_magnitude_recurrence.svg"),
  plot = p_tracer_enrichment_grob,
  width = 14,
  height = 15,
  bg = "white"
)

# -------------------------
# 2) Explicit selection table
# -------------------------

selected_tracers_summary_transparent <- tracer_enrichment_site |>
  dplyr::filter(selected_site_specific) |>
  dplyr::select(
    Site, Element, n_observations, n_positive, mean_positive_z,
    median_positive_z, n_strong_events, frequency_strong_pct,
    max_z, selection_rule
  )

# Sensitivity table: shows how the retained set changes if the minimum number
# of strong events is changed. This is preferable to choosing a threshold after
# seeing only the preferred result.
selection_sensitivity <- tidyr::crossing(
  tracer_enrichment_site,
  min_events_tested = 1:3
) |>
  dplyr::mutate(retained = n_strong_events >= min_events_tested) |>
  dplyr::select(
    Site, Element, min_events_tested, retained,
    n_observations, n_strong_events, frequency_strong_pct
  )

writexl::write_xlsx(
  list(
    tracer_metrics_all = tracer_enrichment_site,
    selected_tracers = selected_tracers_summary_transparent,
    threshold_sensitivity = selection_sensitivity
  ),
  file.path(plot_dir, "Site_specific_tracer_enrichment_transparent_1600.xlsx")
)

# -------------------------
# 3) Optional complementary site-specific PI
# -------------------------
# PI_general is intentionally not recalculated here. In the main analysis it
# remains the mean positive robust z-score across all 11 candidate tracers.

pollution_index_site_specific_transparent <- pollution_scores_clr |>
  dplyr::inner_join(
    selected_tracers_summary_transparent |>
      dplyr::select(Site, Element),
    by = c("Site", "Element")
  ) |>
  dplyr::mutate(z_pos = pmax(robust_z, 0, na.rm = FALSE)) |>
  dplyr::group_by(Site, Depth, Year, Period_bin) |>
  dplyr::summarise(
    PI_site_specific = ifelse(
      sum(!is.na(z_pos)) > 0,
      mean(z_pos, na.rm = TRUE),
      NA_real_
    ),
    n_elements_site_specific = sum(!is.na(z_pos)),
    elements_used = paste(sort(unique(Element[!is.na(z_pos)])), collapse = ", "),
    .groups = "drop"
  )

message(
  "Transparent tracer summary completed. PI_general was not modified."
)


# ============================================================
# General versus site-specific Anthropogenic Enrichment Index (AEI)
#
# Run after:
#   1) the main XRF script (creates pollution_scores_clr), and
#   2) tracer_enrichment_transparent.R (creates
#      selected_tracers_summary_transparent).
#
# The analysis distinguishes association (similar temporal ranking) from
# agreement (numerically interchangeable values).
# ============================================================

required_objects <- c(
  "pollution_scores_clr",
  "selected_tracers_summary_transparent",
  "plot_dir"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]
if (length(missing_objects) > 0) {
  stop(
    "Run the main XRF script and tracer_enrichment_transparent.R first. ",
    "Missing object(s): ", paste(missing_objects, collapse = ", ")
  )
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "writexl", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install package(s): ", paste(missing_packages, collapse = ", "))
}

site_order <- c("DEE", "STO", "CCO", "CLR", "CSD", "CMB", "HIT", "HIO")

site_colors <- c(
  "CMB" = "#458B74",
  "CSD" = "#76EEC6",
  "DEE" = "#CDCD00",
  "HIT" = "#CD96CD",
  "HIO" = "#8B0A50",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "CLR" = "burlywood1"
)

# ============================================================
# 1) Calculate both AEIs from the elemental robust z-scores
# ============================================================
#
# For sample j and element e:
#   enrichment_je = max(0, robust_z_je)
#
# AEI_general,j = arithmetic mean of enrichment_je across every available
# a-priori candidate tracer.
#
# AEI_site_specific,j = the same arithmetic mean, restricted to tracers that
# show robust_z > 3 in at least two samples at that site.
#
# Both indices therefore have the same units and aggregation rule. They differ
# only in the set of included tracers.

sample_keys <- c("Site", "Depth", "Year", "Period_bin")

aei_scores <- pollution_scores_clr |>
  dplyr::mutate(element_enrichment = pmax(robust_z, 0, na.rm = FALSE))

elements_expected <- aei_scores |>
  dplyr::distinct(Element) |>
  dplyr::arrange(Element) |>
  dplyr::pull(Element)

aei_general <- aei_scores |>
  dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
  dplyr::summarise(
    AEI_general = ifelse(
      all(is.finite(element_enrichment)),
      mean(element_enrichment),
      NA_real_
    ),
    n_elements_general = sum(is.finite(element_enrichment)),
    elements_general = paste(
      sort(Element[is.finite(element_enrichment)]),
      collapse = ", "
    ),
    .groups = "drop"
  ) |>
  # A fixed denominator is essential for direct comparison among samples.
  dplyr::mutate(
    complete_general = n_elements_general == length(elements_expected),
    AEI_general = ifelse(complete_general, AEI_general, NA_real_)
  )

aei_site_specific <- aei_scores |>
  dplyr::inner_join(
    selected_tracers_summary_transparent |>
      dplyr::select(Site, Element),
    by = c("Site", "Element")
  ) |>
  dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
  dplyr::summarise(
    AEI_site_specific = ifelse(
      all(is.finite(element_enrichment)),
      mean(element_enrichment),
      NA_real_
    ),
    n_elements_site_specific = sum(is.finite(element_enrichment)),
    elements_site_specific = paste(
      sort(Element[is.finite(element_enrichment)]),
      collapse = ", "
    ),
    .groups = "drop"
  )

index_comparison <- aei_general |>
  dplyr::inner_join(
    aei_site_specific,
    by = sample_keys
  ) |>
  dplyr::filter(is.finite(AEI_general), is.finite(AEI_site_specific)) |>
  dplyr::filter(is.finite(Year)) |>
  dplyr::mutate(
    Site = factor(Site, levels = site_order),
    X_value = as.numeric(Year),
    difference = AEI_site_specific - AEI_general,
    pair_mean = (AEI_site_specific + AEI_general) / 2
  )

# -------------------------
# Statistical functions
# -------------------------

lin_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  
  covariance_xy <- mean((x - mean(x)) * (y - mean(y)))
  2 * covariance_xy /
    (stats::var(x) * (length(x) - 1) / length(x) +
       stats::var(y) * (length(y) - 1) / length(y) +
       (mean(x) - mean(y))^2)
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

summarise_comparison <- function(data) {
  x <- data$AEI_general
  y <- data$AEI_site_specific
  difference <- y - x
  
  tibble::tibble(
    n_pairs = sum(is.finite(x) & is.finite(y)),
    mean_general = mean(x, na.rm = TRUE),
    mean_site_specific = mean(y, na.rm = TRUE),
    mean_difference = mean(difference, na.rm = TRUE),
    relative_mean_difference_pct =
      100 * mean(difference, na.rm = TRUE) / mean(x, na.rm = TRUE),
    lower_limit_agreement = mean(difference, na.rm = TRUE) -
      1.96 * stats::sd(difference, na.rm = TRUE),
    upper_limit_agreement = mean(difference, na.rm = TRUE) +
      1.96 * stats::sd(difference, na.rm = TRUE),
    pearson_r = safe_cor(x, y, "pearson"),
    spearman_rho = safe_cor(x, y, "spearman"),
    lin_ccc = lin_ccc(x, y),
    mean_absolute_difference = mean(abs(difference), na.rm = TRUE),
    root_mean_squared_difference = sqrt(mean(difference^2, na.rm = TRUE))
  )
}

comparison_overall <- summarise_comparison(index_comparison) |>
  dplyr::mutate(Site = "All sites", .before = 1)

comparison_by_site <- index_comparison |>
  dplyr::group_by(Site) |>
  dplyr::group_modify(~ summarise_comparison(.x)) |>
  dplyr::ungroup()

# Compare temporal shapes after removing site-specific location and scale.
# This is a sensitivity analysis of pattern similarity, not raw agreement.
index_comparison_standardised <- index_comparison |>
  dplyr::group_by(Site) |>
  dplyr::mutate(
    AEI_general_standardised = as.numeric(scale(AEI_general)),
    AEI_site_specific_standardised = as.numeric(scale(AEI_site_specific))
  ) |>
  dplyr::ungroup()

shape_comparison <- index_comparison_standardised |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    pearson_r = safe_cor(
      AEI_general_standardised,
      AEI_site_specific_standardised,
      "pearson"
    ),
    spearman_rho = safe_cor(
      AEI_general_standardised,
      AEI_site_specific_standardised,
      "spearman"
    ),
    lin_ccc = lin_ccc(
      AEI_general_standardised,
      AEI_site_specific_standardised
    )
  )

# -------------------------
# Figure A: temporal profiles
# -------------------------

timeline_data <- index_comparison |>
  dplyr::select(Site, Year, AEI_general, AEI_site_specific) |>
  tidyr::pivot_longer(
    cols = c(AEI_general, AEI_site_specific),
    names_to = "Index",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    Index = factor(
      Index,
      levels = c("AEI_general", "AEI_site_specific"),
      labels = c("General", "Site-specific")
    )
  ) |>
  dplyr::arrange(Site, Index, Year)

p_index_timeline <- ggplot2::ggplot(
  timeline_data,
  ggplot2::aes(x = Year, y = Value, colour = Index, group = Index)
) +
  ggplot2::geom_line(linewidth = 0.75, alpha = 0.90) +
  ggplot2::facet_wrap(~ Site, scales = "free", ncol = 2) +
  ggplot2::scale_colour_manual(
    values = c("General" = "#0072B2", "Site-specific" = "#D55E00")
  ) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Anthropogenic Enrichment Index (AEI)",
    colour = "Index",
    title = "General and site-specific enrichment histories"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

# -------------------------
# Figure B: agreement
# -------------------------

p_index_agreement <- ggplot2::ggplot(
  index_comparison,
  ggplot2::aes(x = AEI_general, y = AEI_site_specific, colour = Site)
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", colour = "grey35", linewidth = 0.7
  ) +
  ggplot2::geom_point(alpha = 0.65, size = 1.8) +
  # Fixed scales are required with coord_equal() and also make distances from
  # the 1:1 identity line directly comparable among sites.
  ggplot2::facet_wrap(~ Site, scales = "fixed") +
  ggplot2::coord_equal() +
  ggplot2::guides(colour = "none") +
  ggplot2::labs(
    x = "General AEI",
    y = "Site-specific AEI",
    title = "Agreement between general and site-specific AEI",
    subtitle = "Dashed line indicates numerical identity"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

# -------------------------
# Figure C: Bland-Altman differences
# -------------------------

overall_bias <- comparison_overall$mean_difference
overall_lower <- comparison_overall$lower_limit_agreement
overall_upper <- comparison_overall$upper_limit_agreement

p_index_difference <- ggplot2::ggplot(
  index_comparison,
  ggplot2::aes(x = pair_mean, y = difference, colour = Site)
) +
  ggplot2::geom_hline(yintercept = overall_bias, linewidth = 0.8) +
  ggplot2::geom_hline(
    yintercept = c(overall_lower, overall_upper),
    linetype = "dashed", linewidth = 0.7
  ) +
  ggplot2::geom_point(alpha = 0.65, size = 1.8) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
  ggplot2::labs(
    x = "Mean of paired indices",
    y = "Site-specific AEI minus general AEI",
    colour = "Site",
    title = "Difference between the two composite indices"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

# -------------------------
# Figure D: paired timeline + agreement mosaic
# -------------------------
# Each site contributes two adjacent panels. Colour identifies the site;
# line type and point shape distinguish the two AEI formulations.

make_site_pair <- function(site_code) {
  site_timeline <- timeline_data |>
    dplyr::filter(as.character(Site) == site_code)
  
  site_pairs <- index_comparison |>
    dplyr::filter(as.character(Site) == site_code)
  
  site_limit <- max(
    c(site_pairs$AEI_general, site_pairs$AEI_site_specific),
    na.rm = TRUE
  )
  site_limit <- site_limit * 1.05
  
  timeline_panel <- ggplot2::ggplot(
    site_timeline,
    ggplot2::aes(
      x = Year,
      y = Value,
      linetype = Index,
      shape = Index,
      group = Index
    )
  ) +
    ggplot2::geom_line(
      colour = unname(site_colors[site_code]),
      linewidth = 0.75,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      colour = unname(site_colors[site_code]),
      size = 1.25,
      alpha = 0.80
    ) +
    ggplot2::scale_linetype_manual(
      values = c("General" = "solid", "Site-specific" = "dashed")
    ) +
    ggplot2::scale_shape_manual(
      values = c("General" = 16, "Site-specific" = 1)
    ) +
    ggplot2::labs(
      x = "Year (CE)",
      y = "AEI",
      title = site_code,
      subtitle = "Temporal profiles",
      linetype = "Index",
      shape = "Index"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold", size = 13,
        colour = unname(site_colors[site_code])
      ),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  agreement_panel <- ggplot2::ggplot(
    site_pairs,
    ggplot2::aes(x = AEI_general, y = AEI_site_specific)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      colour = unname(site_colors[site_code]),
      size = 1.7,
      alpha = 0.70
    ) +
    ggplot2::scale_x_continuous(limits = c(0, site_limit)) +
    ggplot2::scale_y_continuous(limits = c(0, site_limit)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "General AEI",
      y = "Site-specific AEI",
      title = site_code,
      subtitle = "Numerical agreement"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold", size = 13,
        colour = unname(site_colors[site_code])
      ),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 10)
    )
  
  list(timeline_panel, agreement_panel)
}

# Two sites per row; within each site, timeline is immediately followed by
# agreement. The requested site order is preserved.
paired_plot_list <- unlist(
  lapply(site_order, make_site_pair),
  recursive = FALSE
)

p_AEI_mosaic <- patchwork::wrap_plots(
  paired_plot_list,
  ncol = 4,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "General and site-specific Anthropogenic Enrichment Indices",
    subtitle = paste0(
      "For each site, temporal profiles are shown alongside agreement with ",
      "the 1:1 line"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

print(p_index_timeline)
print(p_index_agreement)
print(p_index_difference)
print(p_AEI_mosaic)

ggplot2::ggsave(
  file.path(plot_dir, "Fig_SX_AEI_general_site_specific_timeline.svg"),
  p_index_timeline, width = 11, height = 13, bg = "white"
)
ggplot2::ggsave(
  file.path(plot_dir, "Fig_SX_AEI_general_site_specific_agreement.svg"),
  p_index_agreement, width = 10, height = 8, bg = "white"
)
ggplot2::ggsave(
  file.path(plot_dir, "Fig_SX_AEI_general_site_specific_difference.svg"),
  p_index_difference, width = 9, height = 7, bg = "white"
)
ggplot2::ggsave(
  file.path(plot_dir, "Fig_SX_OK_AEI_timeline_agreement_mosaic.svg"),
  p_AEI_mosaic, width = 20, height = 17, bg = "white"
)
ggplot2::ggsave(
  file.path(plot_dir, "Fig_SX_ok_AEI_timeline_agreement_mosaic.png"),
  p_AEI_mosaic, width = 20, height = 17, dpi = 600, bg = "white"
)

writexl::write_xlsx(
  list(
    AEI_general = aei_general,
    AEI_site_specific = aei_site_specific,
    paired_indices = index_comparison,
    overall_comparison = comparison_overall,
    site_comparison = comparison_by_site,
    temporal_shape_comparison = shape_comparison
  ),
  file.path(plot_dir, "AEI_general_site_specific_comparison_1600.xlsx")
)

print(comparison_overall)
print(comparison_by_site)
print(shape_comparison)
