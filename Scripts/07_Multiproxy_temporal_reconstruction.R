# ==========================================================
# MULTIPROXY TIMELINES:
# Combined 18S + COI richness, stable isotopes and general AEI
# ==========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(purrr)
library(janitor)

# ==========================================================
# 0) PATHS
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

iso_file <- file.path(base_dir,"EB_Stable_Isotopes_Spanish_Cores_FInal.xlsx")

output_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/timelines/richness_18S_COI_xrf_AEI_OK_final"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tax_cols_n <- 12
year_col <- "Year"

site_order <- c(
  "DEE", "STO", "CCO", "CLR",
  "CSD", "CMB", "HIT", "HIO"
)

# ==========================================================
# 1) HELPERS
# ==========================================================

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

save_plot_both <- function(plot, filename_base, width, height, dpi = 600) {
  ggsave(
    filename = file.path(output_dir, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(filename_base, ".svg")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

add_year_axis <- function(df) {
  df %>%
    mutate(Year = to_num(Year)) %>%
    filter(!is.na(Year)) %>%
    mutate(
      X_type = "Year",
      Y_plot = Year
    )
}

y_breaks_fun <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) return(NULL)
  
  if (max(abs(rng), na.rm = TRUE) < 500) {
    return(pretty(rng, n = 6))
  } else {
    return(seq(
      floor(rng[1] / 100) * 100,
      ceiling(rng[2] / 100) * 100,
      by = 100
    ))
  }
}

# ==========================================================
# 2) READ METADATA
# ==========================================================

metadata <- read_excel(metadata_file) %>%
  mutate(
    Year = to_num(.data[[year_col]])
  )

site_region <- metadata %>%
  distinct(Site, Region)

# ==========================================================
# 3) RICHNESS 18S + COI
# ==========================================================

calculate_richness <- function(df, marker_name) {
  
  all_euk <- df %>%
    filter(domain == "Eukaryota")
  
  sample_cols <- colnames(all_euk)[(tax_cols_n + 1):ncol(all_euk)]
  
  abundance_df <- all_euk[, sample_cols, drop = FALSE]
  
  for (cc in sample_cols) {
    abundance_df[[cc]] <- as.numeric(abundance_df[[cc]])
    abundance_df[[cc]][is.na(abundance_df[[cc]])] <- 0
  }
  
  data.frame(
    sample = sample_cols,
    richness_q0 = colSums(abundance_df > 0, na.rm = TRUE),
    Marker = marker_name,
    stringsAsFactors = FALSE
  )
}

data_18S <- read_excel(file_18S)
data_COI <- read_excel(file_COI)

richness_all <- bind_rows(
  calculate_richness(data_18S, "18S"),
  calculate_richness(data_COI, "COI")
)

richness_long <- richness_all %>%
  left_join(metadata, by = "sample") %>%
  mutate(
    Year = to_num(Year),
    Site = as.character(Site)
  ) %>%
  filter(!is.na(Site), !is.na(Region)) %>%
  add_year_axis() %>%
  filter(!is.na(Y_plot)) %>%
  mutate(
    Marker_clean = toupper(
      trimws(
        as.character(Marker)
      )
    ),
    
    richness_q0_plot = if_else(
      Marker_clean == "COI",
      as.numeric(richness_q0) * 10,
      as.numeric(richness_q0)
    )
  ) %>%
  transmute(
    Site,
    Y_plot,
    X_type,
    Proxy = "ASV Richness",
    Series = Marker_clean,
    Value_raw = as.numeric(richness_q0),
    Multiplier = if_else(Marker_clean == "COI", 10, 1),
    Value = richness_q0_plot,
    Data_type = "ASV Richness"
  )

coi_check <- richness_long %>%
  filter(Series == "COI")

if (nrow(coi_check) == 0) {
  stop("COI was not found in richness_long.")
}

if (!all(
  abs(coi_check$Value - coi_check$Value_raw * 10) < 1e-10,
  na.rm = TRUE
)) {
  stop("COI was not multiplied by 10 correctly.")
}

message(
  "COI ×10 verified. Raw range: ",
  paste(range(coi_check$Value_raw, na.rm = TRUE), collapse = "–"),
  "; plotted range: ",
  paste(range(coi_check$Value, na.rm = TRUE), collapse = "–")
)

# ==========================================================
# 4) STABLE ISOTOPES
# ==========================================================

sheets_iso <- excel_sheets(iso_file)

df_iso <- map_dfr(sheets_iso, function(sh) {
  read_excel(iso_file, sheet = sh) %>%
    clean_names() %>%
    mutate(Site = toupper(sh))
})

df_iso <- df_iso %>%
  rename(
    Year  = year,
    Sample = name,
    N_pct = n_percent,
    C_pct = c_percent,
    d15N  = d15n_air,
    d13C  = d13c_vpdb,
    Run   = run,
    CN    = c_n
  ) %>%
  mutate(
    Site  = as.character(Site),
    Year  = to_num(Year),
    d15N  = to_num(d15N),
    d13C  = to_num(d13C),
    CN    = to_num(CN)
  ) %>%
  left_join(site_region, by = "Site") %>%
  add_year_axis() %>%
  filter(!is.na(Y_plot))

iso_long <- df_iso %>%
  select(Site, Y_plot, X_type, d13C, d15N, CN) %>%
  pivot_longer(
    cols = c(d13C, d15N, CN),
    names_to = "Proxy",
    values_to = "Value"
  ) %>%
  mutate(
    Proxy = recode(
      Proxy,
      d13C = "δ13C",
      d15N = "δ15N",
      CN = "C:N"
    ),
    Series = Proxy,
    Data_type = "Stable isotopes"
  ) %>%
  filter(!is.na(Value))

# ==========================================================
# 5) ANTHROPOGENIC ENRICHMENT INDEX (GENERAL AEI)
# ==========================================================

aei_general <- read_excel(xrf_results_file, sheet = "AEI_general") %>%
  mutate(
    Site = as.character(Site),
    Year = to_num(Year),
    AEI_general = to_num(AEI_general)
  )

aei_long <- aei_general %>%
  add_year_axis() %>%
  filter(!is.na(Y_plot)) %>%
  select(Site, Y_plot, X_type, AEI_general) %>%
  pivot_longer(
    cols = AEI_general,
    names_to = "Series",
    values_to = "Value"
  ) %>%
  mutate(
    Proxy = "AEI",
    Series = recode(
      Series,
      AEI_general = "General AEI"
    ),
    Data_type = "AEI"
  ) %>%
  filter(!is.na(Value))

# ==========================================================
# 6) COMBINE DATASETS
# ==========================================================

timeline_aei <- bind_rows(
  richness_long,
  iso_long,
  aei_long
)

proxy_order_aei <- c(
  "ASV Richness",
  "δ13C",
  "δ15N",
  "C:N",
  "AEI"
)

timeline_aei <- timeline_aei %>%
  mutate(
    Proxy = factor(Proxy, levels = proxy_order_aei),
    Site = factor(Site, levels = site_order)
  ) %>%
  arrange(Site, Proxy, Series, Y_plot)

# ----------------------------------------------------------
# SYNCHRONIZED TEMPORAL TRENDS
# ----------------------------------------------------------
# Each series is standardized independently within each site.
# Therefore, multiplying COI richness by 10 changes its graphical
# representation but does not affect the synchronized-trend analysis.

# Window length used for each sediment record.
# CCO and CLR use 50 years because their records span less than 100 years.
trend_window_by_site <- c(
  "DEE" = 100,
  "STO" = 100,
  "CCO" = 25,
  "CLR" = 25,
  "CSD" = 100,
  "CMB" = 100,
  "HIT" = 25,
  "HIO" = 25
)

# Minimum change required during the complete temporal window,
# expressed in within-series standard deviations.
min_standardized_change <- 0.25

# Minimum number of proxies changing in the same direction.
min_proxies_synchronous <- 3


# ==========================================================
# STANDARDISE EACH TEMPORAL SERIES
# ==========================================================

timeline_trend_standardized <- timeline_aei %>%
  group_by(
    Site,
    Proxy,
    Series
  ) %>%
  mutate(
    Value_mean = mean(
      Value,
      na.rm = TRUE
    ),
    
    Value_sd = sd(
      Value,
      na.rm = TRUE
    ),
    
    Value_z = case_when(
      is.finite(Value_sd) &
        Value_sd > 0 ~
        (Value - Value_mean) / Value_sd,
      
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(
    is.finite(Value_z),
    is.finite(Y_plot)
  ) %>%
  select(
    -Value_mean,
    -Value_sd
  )


# ==========================================================
# FUNCTION TO CALCULATE LOCAL SLOPES
# ==========================================================

calculate_local_proxy_slopes <- function(site_df) {
  
  # Identify the site.
  site_here <- as.character(
    site_df$Site[[1]]
  )
  
  # Obtain the site-specific window length.
  trend_window_here <- unname(
    trend_window_by_site[site_here]
  )
  
  if (
    length(trend_window_here) != 1 ||
    is.na(trend_window_here)
  ) {
    stop(
      "No temporal trend window was defined for site: ",
      site_here
    )
  }
  
  # Half of the complete window is used before and after
  # each central year.
  trend_half_window_here <-
    trend_window_here / 2
  
  # Evaluate the trend for every year covered by the site.
  years_here <- seq(
    floor(
      min(
        site_df$Y_plot,
        na.rm = TRUE
      )
    ),
    
    ceiling(
      max(
        site_df$Y_plot,
        na.rm = TRUE
      )
    ),
    
    by = 1
  )
  
  local_slopes <- tidyr::crossing(
    Year_center = years_here,
    
    Proxy = unique(
      as.character(
        site_df$Proxy
      )
    )
  ) %>%
    mutate(
      fit = purrr::map2(
        Year_center,
        Proxy,
        
        function(
    year_here,
    proxy_here
        ) {
          
          # Select observations within the temporal window.
          window_df <- site_df %>%
            filter(
              as.character(Proxy) ==
                proxy_here,
              
              Y_plot >=
                year_here -
                trend_half_window_here,
              
              Y_plot <=
                year_here +
                trend_half_window_here
            )
          
          # A regression requires at least:
          # - three observations;
          # - three different years;
          # - two different standardized values.
          if (
            nrow(window_df) < 3 ||
            n_distinct(
              window_df$Y_plot
            ) < 3 ||
            n_distinct(
              window_df$Value_z
            ) < 2
          ) {
            return(
              tibble(
                slope = NA_real_,
                p_value = NA_real_,
                n_observations = nrow(
                  window_df
                ),
                n_years = n_distinct(
                  window_df$Y_plot
                )
              )
            )
          }
          
          # Fit the local linear regression.
          model_here <- lm(
            Value_z ~ Y_plot,
            data = window_df
          )
          
          coefficient_table <-
            summary(
              model_here
            )$coefficients
          
          tibble(
            slope = unname(
              coefficient_table[
                "Y_plot",
                "Estimate"
              ]
            ),
            
            p_value = unname(
              coefficient_table[
                "Y_plot",
                "Pr(>|t|)"
              ]
            ),
            
            n_observations = nrow(
              window_df
            ),
            
            n_years = n_distinct(
              window_df$Y_plot
            )
          )
        }
      ),
    
    # Store the window used for each result.
    trend_window_years =
      trend_window_here
    ) %>%
    unnest(
      cols = fit
    )
  
  return(
    local_slopes
  )
}


# ==========================================================
# CALCULATE LOCAL SLOPES FOR EVERY SITE
# ==========================================================

local_proxy_slopes <-
  timeline_trend_standardized %>%
  group_by(Site) %>%
  group_split() %>%
  map_dfr(
    function(site_df) {
      
      site_here <- as.character(
        site_df$Site[[1]]
      )
      
      calculate_local_proxy_slopes(
        site_df
      ) %>%
        mutate(
          Site = site_here,
          .before = 1
        )
    }
  )


# ==========================================================
# IDENTIFY SYNCHRONIZED YEARS
# ==========================================================

synchronized_years <-
  local_proxy_slopes %>%
  mutate(
    # Estimated standardized change across the complete
    # site-specific temporal window.
    standardized_change =
      slope *
      trend_window_years,
    
    positive_trend =
      is.finite(
        standardized_change
      ) &
      standardized_change >=
      min_standardized_change,
    
    negative_trend =
      is.finite(
        standardized_change
      ) &
      standardized_change <=
      -min_standardized_change
  ) %>%
  group_by(
    Site,
    Year_center
  ) %>%
  summarise(
    n_positive = sum(
      positive_trend,
      na.rm = TRUE
    ),
    
    n_negative = sum(
      negative_trend,
      na.rm = TRUE
    ),
    
    n_proxies_available = sum(
      is.finite(
        standardized_change
      )
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Direction = case_when(
      n_positive >=
        min_proxies_synchronous ~
        "Synchronous increase",
      
      n_negative >=
        min_proxies_synchronous ~
        "Synchronous decrease",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Direction)
  ) %>%
  arrange(
    Site,
    Year_center
  )


# ==========================================================
# GROUP CONSECUTIVE SYNCHRONIZED YEARS INTO INTERVALS
# ==========================================================

synchronized_years <-
  synchronized_years %>%
  group_by(Site) %>%
  mutate(
    new_interval =
      row_number() == 1 |
      
      Direction != lag(
        Direction,
        default = first(Direction)
      ) |
      
      Year_center -
      lag(
        Year_center,
        default = first(Year_center)
      ) > 1,
    
    interval_id =
      cumsum(
        new_interval
      )
  ) %>%
  ungroup()


synchronized_intervals <-
  synchronized_years %>%
  group_by(
    Site,
    Direction,
    interval_id
  ) %>%
  summarise(
    ymin =
      min(
        Year_center
      ) - 0.5,
    
    ymax =
      max(
        Year_center
      ) + 0.5,
    
    n_years_interval =
      n(),
    
    maximum_agreement =
      max(
        pmax(
          n_positive,
          n_negative
        )
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = site_order
    )
  )


# ==========================================================
# DEFINE A COMMON Y-AXIS WITHIN EACH SITE
# ==========================================================
# The facet grid uses one temporal y-axis per site row. These anchor points
# ensure that the complete temporal range of each site is retained.

site_y_ranges <-
  timeline_aei %>%
  group_by(Site) %>%
  summarise(
    site_ymin = min(
      Y_plot,
      na.rm = TRUE
    ),
    
    site_ymax = max(
      Y_plot,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


panel_axis_anchors <-
  timeline_aei %>%
  group_by(
    Site,
    Proxy
  ) %>%
  summarise(
    x_anchor = mean(
      range(
        Value,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    site_y_ranges,
    by = "Site"
  ) %>%
  pivot_longer(
    cols = c(
      site_ymin,
      site_ymax
    ),
    
    names_to = "axis_limit",
    values_to = "Y_plot"
  )

# ==========================================================
# 7) COLORS
# ==========================================================

series_colors <- c(
  "18S" = "#0072B2",
  "COI" = "#D55E00",
  "δ13C" = "#009E73",
  "δ15N" = "#CC79A7",
  "C:N" = "#F0E442",
  "General AEI" = "#777777"
)

# ==========================================================
# 8) COMMON THEME
# ==========================================================

theme_multiproxy <- function() {
  theme_classic(base_size = 15, base_family = "Arial") +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = 17,
        face = "bold",
        color = "black",
        margin = margin(b = 8)
      ),
      strip.text.y = element_text(
        size = 17,
        face = "bold",
        color = "black",
        angle = 0
      ),
      
      axis.text.x = element_text(
        size = 12,
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.title.x = element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      axis.title.y = element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
      axis.line.x.top = element_line(color = "black", linewidth = 0.5),
      axis.line.y.left = element_line(color = "black", linewidth = 0.5),
      axis.line.y.right = element_blank(),
      
      axis.ticks = element_line(color = "black"),
      axis.ticks.y.right = element_blank(),
      
      panel.grid = element_blank(),
      panel.spacing.x = unit(2.0, "lines"),
      panel.spacing.y = unit(1.2, "lines"),
      
      legend.position = "bottom",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15, face = "bold"),
      
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 16),
      plot.margin = margin(15, 20, 15, 20)
    )
}

# ==========================================================
# PUBLICATION THEME
# ==========================================================

x_labels_clean <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-8,
    as.character(round(x)),
    sprintf("%.2f", x)
  )
}

theme_multiproxy_pub <- function() {
  theme_classic(base_size = 22, base_family = "Arial") +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      
      # nombres de proxies arriba
      strip.text.x = element_text(
        size = 28,
        face = "bold",
        color = "black",
        margin = margin(b = 14)
      ),
      
      # nombres de sites a la izquierda
      strip.text.y.left = element_text(
        size = 30,
        face = "bold",
        color = "black",
        angle = 0,
        margin = margin(r = 16)
      ),
      
      axis.text.x = element_text(
        size = 20,
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 22,
        color = "black"
      ),
      
      axis.title.x = element_text(
        size = 28,
        face = "bold",
        color = "black",
        margin = margin(t = 14)
      ),
      axis.title.y = element_text(
        size = 30,
        face = "bold",
        color = "black",
        margin = margin(r = 14)
      ),
      
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.7),
      axis.line.x.top = element_line(color = "black", linewidth = 0.7),
      axis.line.y.left = element_line(color = "black", linewidth = 0.7),
      axis.line.y.right = element_blank(),
      
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.18, "cm"),
      
      panel.grid = element_blank(),
      panel.spacing.x = unit(2.4, "lines"),
      panel.spacing.y = unit(2.0, "lines"),
      
      legend.position = "bottom",
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 24, face = "bold"),
      legend.key.width = unit(1.5, "cm"),
      
      plot.title = element_text(size = 34, face = "bold"),
      plot.subtitle = element_text(size = 24),
      plot.margin = margin(20, 25, 20, 25)
    )
}

# ==========================================================
# 9) MULTIPROXY FIGURE: RICHNESS + ISOTOPES + GENERAL AEI
# ==========================================================

p_aei_multi <- ggplot(
  timeline_aei,
  aes(
    x = Value,
    y = Y_plot,
    color = Series,
    group = interaction(Site, Proxy, Series)
  )
) +
  geom_rect(
    data = synchronized_intervals,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    colour = NA,
    alpha = 0.22
  ) +
  geom_blank(
    data = panel_axis_anchors,
    aes(x = x_anchor, y = Y_plot),
    inherit.aes = FALSE
  ) +
  geom_path(linewidth = 1.15, alpha = 0.9, na.rm = TRUE) +
  geom_point(size = 2.4, alpha = 0.9, na.rm = TRUE) +
  ggplot2::facet_grid(
    rows = vars(Site),
    cols = vars(Proxy),
    scales = "free",
    switch = "y"
  ) +
  # Default temporal scale: ticks every 100 years.
  scale_y_continuous(
    name = "Year",
    labels = x_labels_clean,
    breaks = scales::breaks_width(100)
  ) +
  
  # Site-specific temporal scales.
  ggh4x::facetted_pos_scales(
    y = list(
      Site %in% c(
        "CCO",
        "CLR",
        "HIT",
        "HIO"
      ) ~ scale_y_continuous(
        name = "Year",
        labels = x_labels_clean,
        breaks = scales::breaks_width(20)
      )
    )
  ) +
  scale_x_continuous(
    labels = x_labels_clean,
    breaks = scales::breaks_pretty(n = 3)
  ) +
  scale_color_manual(
    values = series_colors,
    labels = function(x) {
      ifelse(
        x == "COI",
        "COI (×10)",
        x
      )
    },
    name = "Variable"
  ) +
  labs(
    x = NULL,
    y = "Year"
  ) +
  theme_multiproxy_pub() +
  theme(
    # Keep columns close together, but leave enough space for the five
    # independent x-axes in the final row to remain visually separated.
    panel.spacing.x = unit(0.70, "lines"),
    panel.spacing.y = unit(1.4, "lines"),
    
    # Keep a single y-axis on the left side of every site row.
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    
    # Reduce margins between adjacent panels.
    plot.margin = margin(18, 16, 18, 16)
  ) +
  ggh4x::force_panelsizes(
    cols = grid::unit(
      c(
        0.85,  # Richness
        0.85,  # delta 13C
        0.85,  # delta 15N
        0.85,  # C:N
        0.85   # AEI
      ),
      "null"
    )
  )

print(p_aei_multi)

save_plot_both(
  plot = p_aei_multi,
  filename_base = "Fig_multiproxy_richness_isotopes_AEI",
  width = 21,
  height = 32
)

# ==========================================================
# 10) EXPORT DATA USED FOR PLOT
# ==========================================================

write_xlsx(
  list(
    timeline_AEI = timeline_aei,
    standardized_trend_input = timeline_trend_standardized,
    local_proxy_slopes = local_proxy_slopes,
    synchronized_years = synchronized_years,
    synchronized_intervals = synchronized_intervals,
    richness_long = richness_long,
    isotopes_long = iso_long,
    AEI_long = aei_long
  ),
  file.path(output_dir, "Multiproxy_timeline_plot_data.xlsx")
)

cat("DONE: multiproxy figures and data exported to:\n", output_dir, "\n")




# ==========================================================
# 7) CONCURRENT MULTIVARIATE CHANGE
# ==========================================================
# A concurrent-change interval is identified when at least
# three of the five proxies show a sufficiently large local
# change, irrespective of whether they increase or decrease.
#
# The direction of every proxy is retained for subsequent
# ecological interpretation.

concurrent_change_details <-
  local_proxy_slopes %>%
  mutate(
    Proxy = factor(
      Proxy,
      levels = proxy_order_aei
    ),
    
    # Estimated standardized change over the complete
    # site-specific temporal window.
    standardized_change =
      slope * trend_window_years,
    
    # A proxy is considered active when its absolute change
    # is at least 0.25 within-series standard deviations.
    active_change =
      is.finite(standardized_change) &
      abs(standardized_change) >=
      min_standardized_change,
    
    Trend_direction = case_when(
      !active_change ~
        "No marked change",
      
      standardized_change > 0 ~
        "Increase",
      
      standardized_change < 0 ~
        "Decrease",
      
      TRUE ~
        "No marked change"
    ),
    
    Direction_symbol = case_when(
      Trend_direction == "Increase" ~ "↑",
      Trend_direction == "Decrease" ~ "↓",
      TRUE ~ "–"
    )
  ) %>%
  arrange(
    Site,
    Year_center,
    Proxy
  )


# ==========================================================
# COUNT ACTIVE PROXIES FOR EVERY SITE AND YEAR
# ==========================================================

concurrent_change_years <-
  concurrent_change_details %>%
  group_by(
    Site,
    Year_center
  ) %>%
  summarise(
    n_proxies_available = sum(
      is.finite(standardized_change)
    ),
    
    n_active_proxies = sum(
      active_change,
      na.rm = TRUE
    ),
    
    n_increasing = sum(
      Trend_direction == "Increase",
      na.rm = TRUE
    ),
    
    n_decreasing = sum(
      Trend_direction == "Decrease",
      na.rm = TRUE
    ),
    
    # Complete directional signature for ecological
    # interpretation.
    Change_signature = paste(
      as.character(Proxy),
      Direction_symbol,
      sep = " ",
      collapse = " | "
    ),
    
    .groups = "drop"
  ) %>%
  filter(
    n_active_proxies >=
      min_proxies_synchronous
  ) %>%
  arrange(
    Site,
    Year_center
  )


# ==========================================================
# GROUP CONSECUTIVE YEARS INTO INTERVALS
# ==========================================================

concurrent_change_years <-
  concurrent_change_years %>%
  group_by(Site) %>%
  mutate(
    new_interval =
      row_number() == 1 |
      
      Year_center -
      lag(
        Year_center,
        default = first(Year_center)
      ) > 1,
    
    interval_id =
      cumsum(new_interval)
  ) %>%
  ungroup()


concurrent_intervals <-
  concurrent_change_years %>%
  group_by(
    Site,
    interval_id
  ) %>%
  summarise(
    ymin =
      min(Year_center) - 0.5,
    
    ymax =
      max(Year_center) + 0.5,
    
    n_years_interval =
      n(),
    
    maximum_active_proxies =
      max(
        n_active_proxies,
        na.rm = TRUE
      ),
    
    maximum_increasing =
      max(
        n_increasing,
        na.rm = TRUE
      ),
    
    maximum_decreasing =
      max(
        n_decreasing,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = site_order
    )
  )


# ==========================================================
# TEMPORAL LIMITS SHARED WITHIN EACH SITE ROW
# ==========================================================

site_y_ranges <-
  timeline_aei %>%
  group_by(Site) %>%
  summarise(
    site_ymin = min(
      Y_plot,
      na.rm = TRUE
    ),
    
    site_ymax = max(
      Y_plot,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


panel_axis_anchors <-
  timeline_aei %>%
  group_by(
    Site,
    Proxy
  ) %>%
  summarise(
    x_anchor = mean(
      range(
        Value,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    site_y_ranges,
    by = "Site"
  ) %>%
  pivot_longer(
    cols = c(
      site_ymin,
      site_ymax
    ),
    
    names_to = "axis_limit",
    values_to = "Y_plot"
  )


# ==========================================================
# 8) COLOURS
# ==========================================================

series_colors <- c(
  "18S" = "#0072B2",
  "COI" = "#D55E00",
  "δ13C" = "#009E73",
  "δ15N" = "#CC79A7",
  "C:N" = "#E6AB02",
  "General AEI" = "#777777"
)


# ==========================================================
# 9) AXIS LABEL FUNCTION
# ==========================================================

x_labels_clean <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-8,
    as.character(round(x)),
    sprintf("%.2f", x)
  )
}


# ==========================================================
# 10) PUBLICATION THEME
# ==========================================================

theme_multiproxy_pub <- function() {
  
  theme_classic(
    base_size = 22,
    base_family = "Arial"
  ) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      
      # Proxy names above the columns.
      strip.text.x = element_text(
        size = 28,
        face = "bold",
        colour = "black",
        margin = margin(
          b = 10
        )
      ),
      
      # Site names on the left.
      strip.text.y.left = element_text(
        size = 30,
        face = "bold",
        colour = "black",
        angle = 0,
        margin = margin(
          r = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 20,
        colour = "black",
        angle = 45,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 22,
        colour = "black"
      ),
      
      axis.title.x = element_blank(),
      
      axis.title.y = element_text(
        size = 30,
        face = "bold",
        colour = "black",
        margin = margin(
          r = 12
        )
      ),
      
      axis.line.x.bottom = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.x.top = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.y.left = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.y.right =
        element_blank(),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks.length =
        unit(
          0.16,
          "cm"
        ),
      
      panel.grid =
        element_blank(),
      
      # Columns remain close, but their x-axes do not touch.
      panel.spacing.x =
        unit(
          0.70,
          "lines"
        ),
      
      panel.spacing.y =
        unit(
          1.4,
          "lines"
        ),
      
      legend.position =
        "bottom",
      
      legend.text = element_text(
        size = 22
      ),
      
      legend.title = element_text(
        size = 24,
        face = "bold"
      ),
      
      legend.key.width =
        unit(
          1.5,
          "cm"
        ),
      
      plot.margin = margin(
        18,
        16,
        18,
        16
      )
    )
}


# ==========================================================
# 11) MULTIPROXY FIGURE
# ==========================================================

p_aei_multi <- ggplot(
  timeline_aei,
  aes(
    x = Value,
    y = Y_plot,
    colour = Series,
    group = interaction(
      Site,
      Proxy,
      Series
    )
  )
) +
  
  # Grey intervals represent concurrent multivariate change.
  # They do not require every proxy to change in the same
  # direction.
  geom_rect(
    data = concurrent_intervals,
    
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax
    ),
    
    inherit.aes = FALSE,
    fill = "grey70",
    colour = NA,
    alpha = 0.22
  ) +
  
  # Invisible points preserve a common temporal range for
  # all five panels belonging to the same site.
  geom_blank(
    data = panel_axis_anchors,
    
    aes(
      x = x_anchor,
      y = Y_plot
    ),
    
    inherit.aes = FALSE
  ) +
  
  geom_path(
    linewidth = 1.15,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  
  geom_point(
    size = 2.4,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  
  facet_grid(
    rows = vars(Site),
    cols = vars(Proxy),
    scales = "free",
    switch = "y"
  ) +
  
  # Default temporal scale: ticks every 100 years.
  scale_y_continuous(
    name = "Year",
    labels = x_labels_clean,
    breaks = scales::breaks_width(100)
  ) +
  
  # Site-specific temporal scales.
  ggh4x::facetted_pos_scales(
    y = list(
      Site %in% c(
        "CCO",
        "CLR",
        "HIT",
        "HIO"
      ) ~ scale_y_continuous(
        name = "Year",
        labels = x_labels_clean,
        breaks = scales::breaks_width(20)
      )
    )
  ) +
  
  scale_x_continuous(
    labels = x_labels_clean,
    breaks = scales::breaks_pretty(
      n = 3
    ),
    
    expand = expansion(
      mult = c(
        0.04,
        0.04
      )
    )
  ) +
  
  scale_colour_manual(
    values = series_colors,
    
    labels = function(x) {
      ifelse(
        x == "COI",
        "COI (×10)",
        x
      )
    },
    
    name = "Variable"
  ) +
  
  labs(
    x = NULL,
    y = "Year"
  ) +
  
  theme_multiproxy_pub() +
  
  # Richness remains wider; the other four columns are
  # narrower and equal in width.
  ggh4x::force_panelsizes(
    cols = grid::unit(
      c(
        0.85,  # Richness
        0.85,  # delta 13C
        0.85,  # delta 15N
        0.85,  # C:N
        0.85   # AEI
      ),
      "null"
    )
  )


print(p_aei_multi)


# ==========================================================
# 12) SAVE FIGURE
# ==========================================================

save_plot_both(
  plot = p_aei_multi,
  
  filename_base =
    "Fig_multiproxy_richness_isotopes_AEI_concurrent_change",
  
  width = 21,
  height = 32
)


# ==========================================================
# 13) EXPORT DATA USED FOR THE FIGURE
# ==========================================================

write_xlsx(
  list(
    timeline_AEI =
      timeline_aei,
    
    standardized_trend_input =
      timeline_trend_standardized,
    
    local_proxy_slopes =
      local_proxy_slopes,
    
    concurrent_change_details =
      concurrent_change_details,
    
    concurrent_change_years =
      concurrent_change_years,
    
    concurrent_intervals =
      concurrent_intervals,
    
    richness_long =
      richness_long,
    
    isotopes_long =
      iso_long,
    
    AEI_long =
      aei_long
  ),
  
  file.path(
    output_dir,
    "Multiproxy_concurrent_change_plot_data.xlsx"
  )
)


cat(
  "DONE: concurrent multivariate-change figure and data exported to:\n",
  output_dir,
  "\n"
)




# ==========================================================
# 7) CONCURRENT MULTIVARIATE CHANGE
# ==========================================================

concurrent_change_details <-
  local_proxy_slopes %>%
  mutate(
    Proxy = factor(
      Proxy,
      levels = proxy_order_aei
    ),
    
    # Estimated change across the complete temporal window,
    # expressed in within-series standard deviations.
    standardized_change =
      slope * trend_window_years,
    
    # A proxy is considered active when the absolute change
    # is at least 0.25 standard deviations.
    active_change =
      is.finite(standardized_change) &
      abs(standardized_change) >=
      min_standardized_change,
    
    Trend_direction = case_when(
      !active_change ~
        "No marked change",
      
      standardized_change > 0 ~
        "Increase",
      
      standardized_change < 0 ~
        "Decrease",
      
      TRUE ~
        "No marked change"
    ),
    
    Direction_symbol = case_when(
      Trend_direction == "Increase" ~ "↑",
      Trend_direction == "Decrease" ~ "↓",
      TRUE ~ "–"
    )
  ) %>%
  arrange(
    Site,
    Year_center,
    Proxy
  )


# ==========================================================
# 8) COUNT ACTIVE PROXIES BY SITE AND YEAR
# ==========================================================

concurrent_change_years <-
  concurrent_change_details %>%
  group_by(
    Site,
    Year_center
  ) %>%
  summarise(
    n_proxies_available = sum(
      is.finite(standardized_change)
    ),
    
    n_active_proxies = sum(
      active_change,
      na.rm = TRUE
    ),
    
    n_increasing = sum(
      Trend_direction == "Increase",
      na.rm = TRUE
    ),
    
    n_decreasing = sum(
      Trend_direction == "Decrease",
      na.rm = TRUE
    ),
    
    Change_signature = paste(
      as.character(Proxy),
      Direction_symbol,
      sep = " ",
      collapse = " | "
    ),
    
    .groups = "drop"
  ) %>%
  filter(
    n_active_proxies >=
      min_proxies_synchronous
  ) %>%
  arrange(
    Site,
    Year_center
  )


# ==========================================================
# 9) DIRECTIONAL ECOLOGICAL CLASSIFICATION
# ==========================================================

directional_change_years <-
  concurrent_change_details %>%
  filter(
    as.character(Proxy) %in%
      c(
        "ASV Richness",
        "AEI"
      )
  ) %>%
  mutate(
    # Use stable internal names for pivot_wider(). This avoids
    # spaces in the generated column names.
    Proxy_key = case_when(
      as.character(Proxy) == "ASV Richness" ~ "Richness",
      as.character(Proxy) == "AEI" ~ "AEI",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Site,
    Year_center,
    Proxy_key,
    standardized_change
  ) %>%
  pivot_wider(
    names_from = Proxy_key,
    values_from = standardized_change,
    names_prefix = "change_"
  ) %>%
  left_join(
    concurrent_change_years %>%
      select(
        Site,
        Year_center,
        n_proxies_available,
        n_active_proxies,
        n_increasing,
        n_decreasing,
        Change_signature
      ),
    
    by = c(
      "Site",
      "Year_center"
    )
  ) %>%
  filter(
    n_active_proxies >=
      min_proxies_synchronous
  ) %>%
  mutate(
    Change_class = case_when(
      # Potential intensification:
      # AEI increases and richness decreases.
      is.finite(change_AEI) &
        is.finite(change_Richness) &
        change_AEI >=
        min_standardized_change &
        change_Richness <=
        -min_standardized_change ~
        "Potential intensification",
      
      # Potential relaxation:
      # AEI decreases and richness increases.
      is.finite(change_AEI) &
        is.finite(change_Richness) &
        change_AEI <=
        -min_standardized_change &
        change_Richness >=
        min_standardized_change ~
        "Potential relaxation",
      
      # At least three proxies change, but the AEI–richness
      # combination does not match either previous pattern.
      TRUE ~
        "Other concurrent change"
    )
  ) %>%
  arrange(
    Site,
    Year_center
  )


# ==========================================================
# 10) GROUP CONSECUTIVE YEARS INTO COLOURED INTERVALS
# ==========================================================

directional_change_years <-
  directional_change_years %>%
  group_by(Site) %>%
  mutate(
    new_interval =
      row_number() == 1 |
      
      Change_class != lag(
        Change_class,
        default = first(Change_class)
      ) |
      
      Year_center -
      lag(
        Year_center,
        default = first(Year_center)
      ) > 1,
    
    interval_id =
      cumsum(new_interval)
  ) %>%
  ungroup()


concurrent_intervals <-
  directional_change_years %>%
  group_by(
    Site,
    Change_class,
    interval_id
  ) %>%
  summarise(
    ymin =
      min(Year_center) - 0.5,
    
    ymax =
      max(Year_center) + 0.5,
    
    n_years_interval =
      n(),
    
    maximum_active_proxies =
      max(
        n_active_proxies,
        na.rm = TRUE
      ),
    
    maximum_increasing =
      max(
        n_increasing,
        na.rm = TRUE
      ),
    
    maximum_decreasing =
      max(
        n_decreasing,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = site_order
    ),
    
    Change_class = factor(
      Change_class,
      levels = c(
        "Potential intensification",
        "Potential relaxation",
        "Other concurrent change"
      )
    )
  )


# ==========================================================
# 11) COMMON TEMPORAL RANGE WITHIN EACH SITE
# ==========================================================

site_y_ranges <-
  timeline_aei %>%
  group_by(Site) %>%
  summarise(
    site_ymin = min(
      Y_plot,
      na.rm = TRUE
    ),
    
    site_ymax = max(
      Y_plot,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


panel_axis_anchors <-
  timeline_aei %>%
  group_by(
    Site,
    Proxy
  ) %>%
  summarise(
    x_anchor = mean(
      range(
        Value,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    site_y_ranges,
    by = "Site"
  ) %>%
  pivot_longer(
    cols = c(
      site_ymin,
      site_ymax
    ),
    
    names_to = "axis_limit",
    values_to = "Y_plot"
  )


# ==========================================================
# 12) COLOURS
# ==========================================================

series_colors <- c(
  "18S" = "#0072B2",
  "COI" = "#D55E00",
  "δ13C" = "#009E73",
  "δ15N" = "#CC79A7",
  "C:N" = "#E6AB02",
  "General AEI" = "#777777"
)


# Colour-blind-friendly interval colours.
change_class_colors <- c(
  "Potential intensification" = "#D55E00",
  "Potential relaxation" = "#0072B2",
  "Other concurrent change" = "grey70"
)


# ==========================================================
# 13) AXIS-LABEL FUNCTION
# ==========================================================

x_labels_clean <- function(x) {
  
  ifelse(
    abs(x - round(x)) < 1e-8,
    
    as.character(
      round(x)
    ),
    
    sprintf(
      "%.2f",
      x
    )
  )
}


# ==========================================================
# 14) PUBLICATION THEME
# ==========================================================

theme_multiproxy_pub <- function() {
  
  theme_classic(
    base_size = 22,
    base_family = "Arial"
  ) +
    theme(
      strip.placement =
        "outside",
      
      strip.background =
        element_blank(),
      
      # Proxy names.
      strip.text.x = element_text(
        size = 28,
        face = "bold",
        colour = "black",
        margin = margin(
          b = 10
        )
      ),
      
      # Site names.
      strip.text.y.left = element_text(
        size = 30,
        face = "bold",
        colour = "black",
        angle = 0,
        margin = margin(
          r = 12
        )
      ),
      
      axis.text.x = element_text(
        size = 20,
        colour = "black",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      axis.text.y = element_text(
        size = 22,
        colour = "black"
      ),
      
      axis.title.x =
        element_blank(),
      
      axis.title.y = element_text(
        size = 30,
        face = "bold",
        colour = "black",
        margin = margin(
          r = 12
        )
      ),
      
      axis.line.x.bottom = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.x.top = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.y.left = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      
      axis.line.y.right =
        element_blank(),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks.length =
        unit(
          0.16,
          "cm"
        ),
      
      panel.grid =
        element_blank(),
      
      # Horizontal separation between proxy columns.
      panel.spacing.x =
        unit(
          0.65,
          "lines"
        ),
      
      panel.spacing.y =
        unit(
          1.4,
          "lines"
        ),
      
      legend.position =
        "bottom",
      
      legend.box =
        "vertical",
      
      legend.text = element_text(
        size = 20
      ),
      
      legend.title = element_text(
        size = 22,
        face = "bold"
      ),
      
      legend.key.width =
        unit(
          1.4,
          "cm"
        ),
      
      legend.spacing.y =
        unit(
          0.20,
          "cm"
        ),
      
      plot.margin = margin(
        18,
        16,
        18,
        16
      )
    )
}


# ==========================================================
# 15) MULTIPROXY FIGURE
# ==========================================================

p_aei_multi <- ggplot(
  timeline_aei,
  
  aes(
    x = Value,
    y = Y_plot,
    colour = Series,
    
    group = interaction(
      Site,
      Proxy,
      Series
    )
  )
) +
  
  # Coloured intervals of concurrent multivariate change.
  geom_rect(
    data = concurrent_intervals,
    
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax,
      fill = Change_class
    ),
    
    inherit.aes = FALSE,
    colour = NA,
    alpha = 0.24
  ) +
  
  # Invisible anchors maintain the same temporal extent
  # across all five panels belonging to the same site.
  geom_blank(
    data = panel_axis_anchors,
    
    aes(
      x = x_anchor,
      y = Y_plot
    ),
    
    inherit.aes = FALSE
  ) +
  
  geom_path(
    linewidth = 1.15,
    alpha = 0.90,
    na.rm = TRUE
  ) +
  
  geom_point(
    size = 2.4,
    alpha = 0.90,
    na.rm = TRUE
  ) +
  
  facet_grid(
    rows = vars(Site),
    cols = vars(Proxy),
    scales = "free",
    switch = "y"
  ) +
  
  # Default scale: ticks separated by 100 years.
  scale_y_continuous(
    name = "Year",
    labels = x_labels_clean,
    breaks = scales::breaks_width(100)
  ) +
  
  # CCO, CLR, HIT and HIO: ticks separated by 20 years.
  ggh4x::facetted_pos_scales(
    y = list(
      Site %in% c(
        "CCO",
        "CLR",
        "HIT",
        "HIO"
      ) ~ scale_y_continuous(
        name = "Year",
        labels = x_labels_clean,
        breaks = scales::breaks_width(20)
      )
    )
  ) +
  
  scale_x_continuous(
    labels = x_labels_clean,
    
    breaks = scales::breaks_pretty(
      n = 3
    ),
    
    expand = expansion(
      mult = c(
        0.04,
        0.04
      )
    )
  ) +
  
  scale_colour_manual(
    values = series_colors,
    
    labels = function(x) {
      ifelse(
        x == "COI",
        "COI (×10)",
        x
      )
    },
    
    name = "Variable"
  ) +
  
  scale_fill_manual(
    values = change_class_colors,
    
    breaks = c(
      "Potential intensification",
      "Potential relaxation",
      "Other concurrent change"
    ),
    
    name = "Concurrent change",
    
    drop = FALSE
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      title = "Variable"
    ),
    
    fill = guide_legend(
      order = 2,
      title = "Concurrent change",
      
      override.aes = list(
        alpha = 0.55
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Year"
  ) +
  
  theme_multiproxy_pub() +
  
  # Richness is wider; the remaining four proxies are
  # narrower and have the same width.
  ggh4x::force_panelsizes(
    cols = grid::unit(
      c(
        1.35,  # Richness
        0.85,  # δ13C
        0.85,  # δ15N
        0.85,  # C:N
        0.85   # AEI
      ),
      
      "null"
    )
  )


print(p_aei_multi)


# ==========================================================
# 16) SAVE FIGURE
# ==========================================================

save_plot_both(
  plot = p_aei_multi,
  
  filename_base =
    "Fig_multiproxy_AEI_directional_concurrent_change",
  
  width = 21,
  height = 32
)


# ==========================================================
# 17) EXPORT DATA
# ==========================================================

write_xlsx(
  list(
    timeline_AEI =
      timeline_aei,
    
    standardized_trend_input =
      timeline_trend_standardized,
    
    local_proxy_slopes =
      local_proxy_slopes,
    
    concurrent_change_details =
      concurrent_change_details,
    
    concurrent_change_years =
      concurrent_change_years,
    
    directional_change_years =
      directional_change_years,
    
    concurrent_intervals =
      concurrent_intervals,
    
    richness_long =
      richness_long,
    
    isotopes_long =
      iso_long,
    
    AEI_long =
      aei_long
  ),
  
  file.path(
    output_dir,
    "Multiproxy_directional_concurrent_change_data.xlsx"
  )
)


cat(
  "DONE: directional concurrent-change figure and data exported to:\n",
  output_dir,
  "\n"
)