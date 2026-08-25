# ============================================================
# TAXONOMIC RICHNESS THROUGH TIME
# Total eukaryotic richness + top 5 phyla + pie chart
# One figure per marker and sampling site
# ============================================================

# Install if necessary:
# install.packages(c(
#   "readxl", "dplyr", "tidyr", "ggplot2",
#   "stringr", "purrr", "patchwork",
#   "scales", "writexl"
# ))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(patchwork)
library(scales)
library(writexl)

# ============================================================
# 0) PATHS
# ============================================================

base_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/Data"
)

base_meta <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/Metadatas"
)

out_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/Composition"
)

metadata_file <- file.path(
  base_meta,
  "metadata_samples.xlsx"
)

file_18S <- file.path(
  base_dir,
  "All_Peninsula_18S_AbRel.xlsx"
)

file_COI <- file.path(
  base_dir,
  "All_Peninsula_COI_AbRel.xlsx"
)

output_dir <- file.path(
  out_dir,
  "Taxonomic_richness_top_phyla"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Sites to represent by depth instead of year.
# Remove or add sites here if necessary.
depth_sites <- c()

# ============================================================
# 1) HELPER FUNCTIONS
# ============================================================

# Convert decimal commas to numeric values.
to_num <- function(x) {
  
  if (is.numeric(x)) {
    return(x)
  }
  
  x <- as.character(x)
  x <- str_replace_all(x, ",", ".")
  
  suppressWarnings(as.numeric(x))
}


# Standardise sample names.
# CMB_36, cmb_36 and " CMB_36 " all become CMB_36.
standardise_sample <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_to_upper(x)
  
  x
}


safe_filename <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9_-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove_all("^_|_$")
}


save_plot_both <- function(
    plot,
    filename_base,
    width = 14,
    height = 7.5,
    dpi = 600
) {
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(filename_base, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(filename_base, ".svg")
    ),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# ============================================================
# 2) READ METADATA
# ============================================================

metadata <- read_excel(metadata_file)

# Check expected columns.
required_metadata_cols <- c(
  "sample",
  "Site",
  "Depth",
  "Year",
  "Region"
)

missing_metadata_cols <- setdiff(
  required_metadata_cols,
  colnames(metadata)
)

if (length(missing_metadata_cols) > 0) {
  
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

metadata <- metadata %>%
  transmute(
    sample = standardise_sample(sample),
    Site = standardise_sample(Site),
    Depth = to_num(Depth),
    Year = to_num(Year),
    Region = as.character(Region),
    Sea = if ("Sea" %in% names(metadata)) {
      as.character(Sea)
    } else {
      NA_character_
    },
    Era = if ("Era" %in% names(metadata)) {
      as.character(Era)
    } else {
      NA_character_
    }
  ) %>%
  distinct(sample, .keep_all = TRUE)

cat(
  "\nNumber of metadata samples:",
  nrow(metadata),
  "\n"
)

cat(
  "First metadata sample names:\n"
)

print(head(metadata$sample))

# ============================================================
# 3) READ ASV TABLES
# ============================================================

read_asv_table <- function(
    file,
    marker_name
) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  df <- read_excel(file)
  
  # Required taxonomy columns.
  required_tax_cols <- c(
    "ASV",
    "domain",
    "phylum",
    "species"
  )
  
  missing_tax_cols <- setdiff(
    required_tax_cols,
    colnames(df)
  )
  
  if (length(missing_tax_cols) > 0) {
    
    stop(
      marker_name,
      ": missing taxonomy columns: ",
      paste(missing_tax_cols, collapse = ", ")
    )
  }
  
  # Sample columns are everything after species.
  species_position <- match(
    "species",
    colnames(df)
  )
  
  if (is.na(species_position)) {
    stop(
      marker_name,
      ": column 'species' was not found."
    )
  }
  
  sample_cols_original <- colnames(df)[
    (species_position + 1):ncol(df)
  ]
  
  if (length(sample_cols_original) == 0) {
    
    stop(
      marker_name,
      ": no abundance columns found after species."
    )
  }
  
  # Standardised sample names corresponding to table columns.
  sample_cols_standardised <- standardise_sample(
    sample_cols_original
  )
  
  # Check for duplicated names after standardisation.
  if (anyDuplicated(sample_cols_standardised) > 0) {
    
    duplicated_samples <- unique(
      sample_cols_standardised[
        duplicated(sample_cols_standardised)
      ]
    )
    
    stop(
      marker_name,
      ": duplicated sample names after standardisation: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }
  
  # Rename abundance columns using standardised sample names.
  names(df)[
    match(
      sample_cols_original,
      names(df)
    )
  ] <- sample_cols_standardised
  
  # Keep only columns represented in metadata.
  sample_cols_matching <- intersect(
    sample_cols_standardised,
    metadata$sample
  )
  
  sample_cols_not_in_metadata <- setdiff(
    sample_cols_standardised,
    metadata$sample
  )
  
  metadata_samples_not_in_table <- setdiff(
    metadata$sample,
    sample_cols_standardised
  )
  
  cat(
    "\n-----------------------------------\n",
    marker_name,
    "\n-----------------------------------\n"
  )
  
  cat(
    "Sample columns in ASV table:",
    length(sample_cols_standardised),
    "\n"
  )
  
  cat(
    "Matching metadata samples:",
    length(sample_cols_matching),
    "\n"
  )
  
  if (length(sample_cols_not_in_metadata) > 0) {
    
    cat(
      "ASV-table samples absent from metadata:",
      length(sample_cols_not_in_metadata),
      "\n"
    )
    
    print(
      head(sample_cols_not_in_metadata, 20)
    )
  }
  
  if (length(metadata_samples_not_in_table) > 0) {
    
    cat(
      "Metadata samples absent from ASV table:",
      length(metadata_samples_not_in_table),
      "\n"
    )
  }
  
  if (length(sample_cols_matching) == 0) {
    
    stop(
      marker_name,
      ": no abundance columns matched metadata sample names."
    )
  }
  
  # Convert abundance values to numeric.
  # This also handles decimal commas.
  df <- df %>%
    mutate(
      across(
        all_of(sample_cols_matching),
        to_num
      )
    )
  
  # Replace missing abundance values with zero.
  df <- df %>%
    mutate(
      across(
        all_of(sample_cols_matching),
        ~ replace_na(.x, 0)
      )
    )
  
  # Retain only required taxonomy and matching sample columns.
  df %>%
    transmute(
      ASV = as.character(ASV),
      Marker = marker_name,
      domain = as.character(domain),
      kingdom = as.character(kingdom),
      phylum = as.character(phylum),
      across(
        all_of(sample_cols_matching)
      )
    )
}


data_18S <- read_asv_table(
  file = file_18S,
  marker_name = "18S"
)

data_COI <- read_asv_table(
  file = file_COI,
  marker_name = "COI"
)

# ============================================================
# 4) CONVERT ASV TABLES TO LONG FORMAT
# ============================================================

asv_to_long <- function(df) {
  
  taxonomy_cols <- c(
    "ASV",
    "Marker",
    "domain",
    "kingdom",
    "phylum"
  )
  
  sample_cols <- setdiff(
    names(df),
    taxonomy_cols
  )
  
  df %>%
    filter(
      str_to_lower(
        str_trim(domain)
      ) == "eukaryota"
    ) %>%
    mutate(
      phylum = case_when(
        is.na(phylum) ~ "Unassigned",
        str_trim(phylum) == "" ~ "Unassigned",
        str_to_lower(str_trim(phylum)) %in% c(
          "na",
          "n/a",
          "unknown",
          "unclassified",
          "unassigned"
        ) ~ "Unassigned",
        TRUE ~ str_trim(phylum)
      )
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    mutate(
      sample = standardise_sample(sample),
      abundance = to_num(abundance)
    ) %>%
    filter(
      !is.na(abundance),
      abundance > 0
    ) %>%
    left_join(
      metadata,
      by = "sample"
    ) %>%
    filter(
      !is.na(Site)
    ) %>%
    select(
      Marker,
      ASV,
      domain,
      kingdom,
      phylum,
      sample,
      abundance,
      Site,
      Depth,
      Year,
      Region,
      Sea,
      Era
    )
}


asv_long <- bind_rows(
  asv_to_long(data_18S),
  asv_to_long(data_COI)
)

cat(
  "\nNumber of positive ASV detections:",
  nrow(asv_long),
  "\n"
)

# ============================================================
# 5) DEFINE TEMPORAL AXIS
# ============================================================

asv_long <- asv_long %>%
  mutate(
    Axis_type = case_when(
      Site %in% depth_sites ~ "Depth",
      !is.na(Year) ~ "Year",
      TRUE ~ "Depth"
    ),
    Axis_value = case_when(
      Axis_type == "Year" ~ Year,
      Axis_type == "Depth" ~ Depth,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Axis_value)
  )

# ============================================================
# 6) PHYLA SELECTED A PRIORI
# ============================================================

selected_phyla <- c(
  "Streptophyta",
  "Chordata",
  "Cnidaria",
  "Mollusca",
  "Porifera",
  "Chlorophyta",
  "Rhodophyta"
)

# Check that names are present in the dataset
missing_selected_phyla <- setdiff(
  selected_phyla,
  unique(asv_long$phylum)
)

if (length(missing_selected_phyla) > 0) {
  warning(
    "Selected phyla absent from asv_long: ",
    paste(missing_selected_phyla, collapse = ", ")
  )
}

# ============================================================
# 7) TOTAL EUKARYOTIC RICHNESS PER SAMPLE
# ============================================================

total_richness <- asv_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    name = "Richness"
  ) %>%
  mutate(
    Series = "Total richness"
  )

# ============================================================
# 8) RICHNESS OF THE SELECTED PHYLA
# ============================================================

# Observed richness: number of distinct ASVs belonging to each
# selected phylum in each sample.
selected_richness_observed <- asv_long %>%
  filter(
    phylum %in% selected_phyla
  ) %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    name = "Richness"
  ) %>%
  rename(
    Series = phylum
  )

# Generate all sample × selected-phylum combinations.
# This makes absent phyla equal to zero for that sample.
selected_richness_grid <- total_richness %>%
  select(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  distinct() %>%
  tidyr::crossing(
    Series = selected_phyla
  )

# Add observed richness and replace true absences with zero.
selected_richness <- selected_richness_grid %>%
  left_join(
    selected_richness_observed,
    by = c(
      "Marker",
      "Site",
      "sample",
      "Axis_type",
      "Axis_value",
      "Series"
    )
  ) %>%
  mutate(
    Richness = replace_na(Richness, 0L)
  )

# Final dataset for the temporal plots.
richness_plot_data <- bind_rows(
  total_richness,
  selected_richness
) %>%
  mutate(
    Series = factor(
      Series,
      levels = c(
        "Total richness",
        selected_phyla
      )
    )
  ) %>%
  arrange(
    Marker,
    Site,
    Series,
    Axis_value
  )

# ============================================================
# 9) PIE-CHART DATA
# ============================================================

# Unique ASVs per assigned phylum and marker-site combination.
# Unassigned phyla are excluded from the denominator here.
phylum_summary <- asv_long %>%
  filter(
    !is.na(phylum),
    !tolower(trimws(phylum)) %in% c(
      "",
      "unassigned",
      "unknown",
      "unclassified",
      "na",
      "n/a"
    )
  ) %>%
  distinct(
    Marker,
    Site,
    phylum,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    phylum,
    name = "n_unique_ASVs"
  )

pie_data <- phylum_summary %>%
  mutate(
    Pie_group = if_else(
      phylum %in% selected_phyla,
      phylum,
      "Other"
    )
  ) %>%
  group_by(
    Marker,
    Site,
    Pie_group
  ) %>%
  summarise(
    n_unique_ASVs = sum(n_unique_ASVs),
    .groups = "drop"
  ) %>%
  group_by(
    Marker,
    Site
  ) %>%
  mutate(
    Percentage = 100 *
      n_unique_ASVs /
      sum(n_unique_ASVs),
    
    Percentage_label = if_else(
      Percentage >= 3,
      paste0(
        round(Percentage, 1),
        "%"
      ),
      ""
    )
  ) %>%
  ungroup()

# ============================================================
# 10) COLOURS
# ============================================================

all_top_phyla <- sort(
  unique(top_phyla$phylum)
)

phylum_colours <- c(
  "Streptophyta" = "#009E73",
  "Chordata"     = "#0072B2",
  "Cnidaria"     = "#CC79A7",
  "Mollusca"     = "#D55E00",
  "Porifera"     = "#56B4E9",
  "Chlorophyta"  = "#E69F00",
  "Rhodophyta"   = "#A6761D"
)

line_colours <- c(
  "Total richness" = "grey40",
  phylum_colours
)

pie_colours <- c(
  phylum_colours,
  "Other" = "grey80"
)

# ============================================================
# 11) PLOT THEMES
# ============================================================

theme_richness <- function() {
  
  theme_classic(
    base_size = 15,
    base_family = "Arial"
  ) +
    theme(
      axis.title = element_text(
        size = 25,
        face = "bold",
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 22,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      legend.position = "bottom",
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = 20
      ),
      
      legend.key.width = unit(
        1.3,
        "cm"
      ),
      
      plot.title = element_text(
        size = 28,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 13
      ),
      
      plot.margin = margin(
        12,
        12,
        12,
        12
      )
    )
}


theme_pie <- function() {
  
  theme_void(
    base_family = "Arial"
  ) +
    theme(
      legend.position = "none",
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}

# ============================================================
# 12) FUNCTION TO CREATE ONE FIGURE
# ============================================================

make_site_marker_plot <- function(
    marker_name,
    site_name
) {
  
  line_df <- richness_plot_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  pie_df <- pie_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  if (nrow(line_df) == 0) {
    
    warning(
      "No richness data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  if (nrow(pie_df) == 0) {
    
    warning(
      "No pie-chart data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  axis_type_here <- unique(
    line_df$Axis_type
  )
  
  if (length(axis_type_here) > 1) {
    
    stop(
      marker_name,
      " - ",
      site_name,
      ": more than one axis type was found."
    )
  }
  
  top_phyla_here <- selected_phyla
  
  series_order <- c(
    "Total richness",
    top_phyla_here
  )
  
  line_df <- line_df %>%
    mutate(
      Series = factor(
        Series,
        levels = series_order
      )
    ) %>%
    arrange(
      Series,
      Axis_value
    )
  
  # Temporal richness plot.
  p_line <- ggplot(
    line_df,
    aes(
      x = Axis_value,
      y = Richness,
      colour = Series,
      group = Series
    )
  ) +
    geom_line(
      aes(
        linewidth = Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      aes(
        size = Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = phylum_colours,
      breaks = selected_phyla,
      limits = selected_phyla,
      drop = FALSE,
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = c(
        `FALSE` = 1,
        `TRUE` = 1.7
      ),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(
        `FALSE` = 2.2,
        `TRUE` = 3
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(
        base = 10,
        sigma = 1
      ),
      breaks = c(
        0, 1, 2, 5, 10, 20, 50,
        100, 200, 500, 1000, 2000
      ),
      labels = scales::label_number(),
      expand = expansion(
        mult = c(0.02, 0.08)
      )
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 7)
    ) +
    labs(
      title = paste0(
        marker_name,
        " — ",
        site_name
      ),
      
      subtitle = paste0(
        "Total eukaryotic richness and richness ",
        "of the five dominant phyla"
      ),
      
      x = ifelse(
        axis_type_here == "Year",
        "Year",
        "Depth (cm)"
      ),
      
      y = "ASV richness (q = 0)"
    ) +
    theme_richness()
  
  # Pie chart.
  pie_order <- c(
    top_phyla_here,
    "Other"
  )
  
  pie_df <- pie_df %>%
    mutate(
      Pie_group = factor(
        Pie_group,
        levels = pie_order
      )
    ) %>%
    arrange(Pie_group)
  
  p_pie <- ggplot(
    pie_df,
    aes(
      x = "",
      y = Percentage,
      fill = Pie_group
    )
  ) +
    geom_col(
      width = 1,
      colour = "white",
      linewidth = 0.5
    ) +
    coord_polar(
      theta = "y"
    ) +
    geom_text(
      aes(
        label = scales::number(
          Percentage,
          accuracy = 1
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 6,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = pie_colours,
      breaks = pie_order,
      limits = pie_order,
      drop = FALSE,
      guide = "none"
    ) +
    labs(
      title = "Unique ASV composition"
    ) +
    theme_pie()
  
  # Combine both panels.
  combined_plot <- p_line + p_pie +
    plot_layout(
      widths = c(
        3.5,
        1.3
      ),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom"
    )
  
  combined_plot
}

make_general_legend <- function() {
  
  legend_df <- tibble(
    Series = factor(
      c("Total richness", selected_phyla, "Other"),
      levels = c("Total richness", selected_phyla, "Other")
    ),
    x = seq_along(c("Total richness", selected_phyla, "Other")),
    y = 1
  )
  
  legend_colours <- c(
    "Total richness" = "grey40",
    phylum_colours,
    "Other" = "grey80"
  )
  
  ggplot(
    legend_df,
    aes(
      x = x,
      y = y,
      colour = Series
    )
  ) +
    geom_point(size = 4) +
    scale_colour_manual(
      values = legend_colours,
      breaks = c(
        "Total richness",
        selected_phyla,
        "Other"
      ),
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          shape = 16,
          size = 4,
          linewidth = 1.3
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(
        size = 15,
        family = "Arial"
      ),
      legend.key.width = unit(1.1, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      plot.margin = margin(0, 5, 0, 5)
    )
}

# ============================================================
# 13) CREATE ALL FIGURES
# ============================================================

marker_site_combinations <- richness_plot_data %>%
  distinct(
    Marker,
    Site
  ) %>%
  arrange(
    Marker,
    Site
  )


all_plots <- pmap(
  marker_site_combinations,
  function(Marker, Site) {
    
    message(
      "Creating plot: ",
      Marker,
      " - ",
      Site
    )
    
    p <- make_site_marker_plot(
      marker_name = Marker,
      site_name = Site
    )
    
    if (!is.null(p)) {
      
      filename_base <- paste0(
        "Richness_top5_phyla_",
        safe_filename(Marker),
        "_",
        safe_filename(Site)
      )
      
      save_plot_both(
        plot = p,
        filename_base = filename_base,
        width = 14,
        height = 7.5,
        dpi = 600
      )
    }
    
    p
  }
)


names(all_plots) <- paste(
  marker_site_combinations$Marker,
  marker_site_combinations$Site,
  sep = "_"
)

# Display the first plot in RStudio.
valid_plot_positions <- which(
  !vapply(
    all_plots,
    is.null,
    logical(1)
  )
)

if (length(valid_plot_positions) > 0) {
  
  print(
    all_plots[[valid_plot_positions[1]]]
  )
}

# ============================================================
# 14) OPTIONAL COMBINED FIGURE FOR EACH MARKER
# ============================================================

for (
  marker_name in unique(
    marker_site_combinations$Marker
  )
) {
  
  selected_names <- names(all_plots)[
    str_starts(
      names(all_plots),
      paste0(marker_name, "_")
    )
  ]
  
  marker_plots <- all_plots[
    selected_names
  ]
  
  marker_plots <- marker_plots[
    !vapply(
      marker_plots,
      is.null,
      logical(1)
    )
  ]
  
  if (length(marker_plots) > 0) {
    
    # Two columns to avoid an excessively tall figure
    n_columns <- 2
    
    n_rows <- ceiling(
      length(marker_plots) / n_columns
    )
    
    general_legend <- make_general_legend()
    
    for (
      marker_name in unique(
        marker_site_combinations$Marker
      )
    ) {
      
      selected_names <- names(all_plots)[
        str_starts(
          names(all_plots),
          paste0(marker_name, "_")
        )
      ]
      
      marker_plots <- all_plots[
        selected_names
      ]
      
      marker_plots <- marker_plots[
        !vapply(
          marker_plots,
          is.null,
          logical(1)
        )
      ]
      
      if (length(marker_plots) > 0) {
        
        n_columns <- 2
        n_rows <- ceiling(
          length(marker_plots) / n_columns
        )
        
        site_panel <- wrap_plots(
          marker_plots,
          ncol = n_columns
        ) +
          plot_annotation(
            title = paste0(
              marker_name,
              ": total eukaryotic richness ",
              "and selected phyla by site"
            ),
            theme = theme(
              plot.title = element_text(
                size = 30,
                face = "bold",
                family = "Arial",
                hjust = 0.5
              )
            )
          )
        
        combined_marker_plot <- site_panel / general_legend +
          plot_layout(
            heights = c(20, 1.2)
          )
        
        save_plot_both(
          plot = combined_marker_plot,
          filename_base = paste0(
            "Richness_selected_phyla_all_sites_",
            marker_name
          ),
          width = 24,
          height = 7.5 * n_rows + 2,
          dpi = 600
        )
      }
    }
  }}

# ============================================================
# 15) EXPORT SOURCE DATA
# ============================================================

write_xlsx(
  list(
    richness_time_series = richness_plot_data,
    
    top5_phyla_by_site = top5_phyla,
    
    unique_ASVs_by_phylum = phylum_summary,
    
    pie_chart_data = pie_data,
    
    metadata_used = metadata
  ),
  
  file.path(
    output_dir,
    "Taxonomic_richness_top5_phyla_plot_data.xlsx"
  )
)

cat(
  "\nDONE\n",
  "Figures and source data exported to:\n",
  output_dir,
  "\n"
)



# ============================================================
# INDEPENDENT ANALYSIS:
# 100% STACKED TAXONOMIC RICHNESS THROUGH TIME
#
# This block does not modify the previous richness figures.
# It uses asv_long and selected_phyla already created.
# ============================================================

stacked_output_dir <- file.path(
  out_dir,
  "Taxonomic_richness_100_percent"
)

dir.create(
  stacked_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 1) ASSIGN EACH ASV TO A PLOTTING CATEGORY
# ============================================================

asv_composition_long <- asv_long %>%
  mutate(
    phylum_clean = str_trim(as.character(phylum)),
    
    Composition_group = case_when(
      is.na(phylum_clean) ~ "Unassigned",
      
      str_to_lower(phylum_clean) %in% c(
        "",
        "unassigned",
        "unknown",
        "unclassified",
        "na",
        "n/a"
      ) ~ "Unassigned",
      
      phylum_clean %in% selected_phyla ~ phylum_clean,
      
      TRUE ~ "Other"
    )
  )

# ============================================================
# 2) CALCULATE UNIQUE-ASV RICHNESS PER SAMPLE AND GROUP
# ============================================================

stacked_richness_observed <- asv_composition_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    Composition_group,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    Composition_group,
    name = "Group_richness"
  )

# ============================================================
# 3) COMPLETE ABSENT GROUPS WITH ZERO
# ============================================================

composition_levels <- c(
  selected_phyla,
  "Other",
  "Unassigned"
)

stacked_sample_grid <- asv_composition_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  tidyr::crossing(
    Composition_group = composition_levels
  )

stacked_richness_complete <- stacked_sample_grid %>%
  left_join(
    stacked_richness_observed,
    by = c(
      "Marker",
      "Site",
      "sample",
      "Axis_type",
      "Axis_value",
      "Composition_group"
    )
  ) %>%
  mutate(
    Group_richness = replace_na(
      Group_richness,
      0L
    )
  )

# ============================================================
# 4) CONVERT RICHNESS TO PERCENTAGE
# ============================================================

stacked_percentage_data <- stacked_richness_complete %>%
  group_by(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  mutate(
    Total_richness_check = sum(
      Group_richness,
      na.rm = TRUE
    ),
    
    Percentage = if_else(
      Total_richness_check > 0,
      100 * Group_richness / Total_richness_check,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    Composition_group = factor(
      Composition_group,
      levels = composition_levels
    )
  )

# Check that percentages sum to 100.
percentage_check <- stacked_percentage_data %>%
  group_by(
    Marker,
    Site,
    sample
  ) %>%
  summarise(
    Percentage_sum = sum(
      Percentage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  summary(
    percentage_check$Percentage_sum
  )
)

# ============================================================
# 5) COLOURS
# ============================================================

stacked_colours <- c(
  phylum_colours,
  "Other" = "grey70",
  "Unassigned" = "grey90"
)

# ============================================================
# 6) FUNCTION FOR ONE MARKER × SITE
# ============================================================

make_stacked_percentage_plot <- function(
    marker_name,
    site_name
) {
  
  plot_df <- stacked_percentage_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    ) %>%
    arrange(
      Axis_value,
      Composition_group
    )
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  
  axis_type_here <- unique(
    plot_df$Axis_type
  )
  
  if (length(axis_type_here) != 1) {
    
    stop(
      marker_name,
      " - ",
      site_name,
      ": more than one axis type detected."
    )
  }
  
  ggplot(
    plot_df,
    aes(
      x = Axis_value,
      y = Percentage,
      fill = Composition_group
    )
  ) +
    geom_area(
      position = "stack",
      alpha = 0.9,
      colour = NA,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = stacked_colours,
      breaks = composition_levels,
      limits = composition_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(
        0,
        100,
        by = 20
      ),
      labels = function(x) {
        paste0(x, "%")
      },
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 7)
    ) +
    labs(
      title = paste0(
        marker_name,
        " — ",
        site_name
      ),
      subtitle = paste0(
        "Relative contribution of selected phyla ",
        "to total ASV richness"
      ),
      x = ifelse(
        axis_type_here == "Year",
        "Year",
        "Depth (cm)"
      ),
      y = "Proportion of total ASV richness"
    ) +
    theme_classic(
      base_size = 15,
      base_family = "Arial"
    ) +
    theme(
      axis.title = element_text(
        size = 16,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      legend.position = "bottom",
      legend.text = element_text(
        size = 12
      ),
      legend.key.width = unit(
        1.1,
        "cm"
      ),
      plot.title = element_text(
        size = 20,
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 13
      )
    )
}

# ============================================================
# 7) CREATE INDIVIDUAL FIGURES
# ============================================================

stacked_combinations <- stacked_percentage_data %>%
  distinct(
    Marker,
    Site
  ) %>%
  arrange(
    Marker,
    Site
  )

stacked_plots <- pmap(
  stacked_combinations,
  function(Marker, Site) {
    
    message(
      "Creating 100% stacked plot: ",
      Marker,
      " - ",
      Site
    )
    
    p <- make_stacked_percentage_plot(
      marker_name = Marker,
      site_name = Site
    )
    
    if (!is.null(p)) {
      
      ggsave(
        filename = file.path(
          stacked_output_dir,
          paste0(
            "Richness_composition_100pct_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".png"
          )
        ),
        plot = p,
        width = 12,
        height = 7,
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        filename = file.path(
          stacked_output_dir,
          paste0(
            "Richness_composition_100pct_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".svg"
          )
        ),
        plot = p,
        width = 12,
        height = 7,
        bg = "white"
      )
    }
    
    p
  }
)

names(stacked_plots) <- paste(
  stacked_combinations$Marker,
  stacked_combinations$Site,
  sep = "_"
)

# ============================================================
# 8) MULTIPANEL PER MARKER WITH ONE GENERAL LEGEND
# ============================================================

for (
  marker_name in unique(
    stacked_combinations$Marker
  )
) {
  
  selected_names <- names(stacked_plots)[
    str_starts(
      names(stacked_plots),
      paste0(marker_name, "_")
    )
  ]
  
  marker_stacked_plots <- stacked_plots[
    selected_names
  ]
  
  marker_stacked_plots <- marker_stacked_plots[
    !vapply(
      marker_stacked_plots,
      is.null,
      logical(1)
    )
  ]
  
  if (length(marker_stacked_plots) > 0) {
    
    # Remove legends from individual panels.
    marker_stacked_plots_no_legend <- lapply(
      marker_stacked_plots,
      function(p) {
        p + theme(
          legend.position = "none"
        )
      }
    )
    
    # General legend.
    stacked_legend_df <- tibble(
      Composition_group = factor(
        composition_levels,
        levels = composition_levels
      ),
      x = seq_along(composition_levels),
      y = 1
    )
    
    stacked_general_legend <- ggplot(
      stacked_legend_df,
      aes(
        x = x,
        y = y,
        fill = Composition_group
      )
    ) +
      geom_col() +
      scale_fill_manual(
        values = stacked_colours,
        breaks = composition_levels,
        drop = FALSE,
        name = NULL
      ) +
      guides(
        fill = guide_legend(
          nrow = 2,
          byrow = TRUE
        )
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(
          size = 15,
          family = "Arial"
        )
      )
    
    n_columns <- 2
    n_rows <- ceiling(
      length(marker_stacked_plots_no_legend) /
        n_columns
    )
    
    stacked_site_panel <- wrap_plots(
      marker_stacked_plots_no_legend,
      ncol = n_columns
    ) +
      plot_annotation(
        title = paste0(
          marker_name,
          ": relative taxonomic composition ",
          "of ASV richness through time"
        ),
        theme = theme(
          plot.title = element_text(
            size = 30,
            face = "bold",
            hjust = 0.5,
            family = "Arial"
          )
        )
      )
    
    combined_stacked_plot <-
      stacked_site_panel /
      stacked_general_legend +
      plot_layout(
        heights = c(20, 1.2)
      )
    
    ggsave(
      filename = file.path(
        stacked_output_dir,
        paste0(
          "Richness_composition_100pct_all_sites_",
          marker_name,
          ".png"
        )
      ),
      plot = combined_stacked_plot,
      width = 24,
      height = 7.2 * n_rows + 2,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(
        stacked_output_dir,
        paste0(
          "Richness_composition_100pct_all_sites_",
          marker_name,
          ".svg"
        )
      ),
      plot = combined_stacked_plot,
      width = 24,
      height = 7.2 * n_rows + 2,
      bg = "white",
      limitsize = FALSE
    )
  }
}

# ============================================================
# 9) EXPORT DATA
# ============================================================

write_xlsx(
  list(
    stacked_richness = stacked_richness_complete,
    stacked_percentages = stacked_percentage_data,
    percentage_check = percentage_check
  ),
  file.path(
    stacked_output_dir,
    "Richness_composition_100pct_data.xlsx"
  )
)

cat(
  "\nDONE: 100% stacked richness figures exported to:\n",
  stacked_output_dir,
  "\n"
)


# ============================================================
# INDEPENDENT ANALYSIS:
# TOTAL RICHNESS + FIVE MOST ASV-RICH PHYLA PER SITE
#
# This section is independent from the previous selected-phyla
# analysis. It requires the existing object:
#
#   asv_long
#
# Expected columns:
# Marker, Site, sample, Axis_type, Axis_value, ASV, phylum
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(patchwork)
library(scales)
library(writexl)

# ============================================================
# 0) OUTPUT DIRECTORY
# ============================================================

auto5_output_dir <- file.path(
  out_dir,
  "Taxonomic_richness_automatic_top5_phyla"
)

dir.create(
  auto5_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 1) CLEAN PHYLA
# ============================================================

auto5_asv_long <- asv_long %>%
  mutate(
    phylum_clean = str_trim(
      as.character(phylum)
    )
  )

# Categories excluded when selecting the five dominant phyla.
auto5_unassigned_terms <- c(
  "",
  "unassigned",
  "unknown",
  "unclassified",
  "unidentified",
  "incertae sedis",
  "na",
  "n/a"
)

# ============================================================
# 2) TOTAL EUKARYOTIC RICHNESS PER SAMPLE
# ============================================================

auto5_total_richness <- auto5_asv_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    name = "Richness"
  ) %>%
  mutate(
    Series = "Total richness"
  )

# ============================================================
# 3) NUMBER OF UNIQUE ASVs PER PHYLUM, MARKER AND SITE
# ============================================================
#
# A phylum is ranked using the number of different ASVs detected
# at least once in the complete sediment record of that site.
# Unassigned categories are excluded.
# ============================================================

auto5_phylum_summary <- auto5_asv_long %>%
  filter(
    !is.na(phylum_clean),
    !str_to_lower(phylum_clean) %in%
      auto5_unassigned_terms
  ) %>%
  distinct(
    Marker,
    Site,
    phylum_clean,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    phylum_clean,
    name = "n_unique_ASVs"
  ) %>%
  rename(
    phylum = phylum_clean
  )

# ============================================================
# 4) SELECT THE FIVE MOST ASV-RICH PHYLA PER MARKER AND SITE
# ============================================================

auto5_top_phyla <- auto5_phylum_summary %>%
  group_by(
    Marker,
    Site
  ) %>%
  arrange(
    desc(n_unique_ASVs),
    phylum,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 5
  ) %>%
  mutate(
    Rank = row_number()
  ) %>%
  ungroup()

cat(
  "\nFive most ASV-rich phyla by marker and site:\n"
)

print(
  auto5_top_phyla,
  n = Inf
)

# ============================================================
# 5) OBSERVED RICHNESS OF THE TOP-FIVE PHYLA PER SAMPLE
# ============================================================

auto5_richness_observed <- auto5_asv_long %>%
  mutate(
    phylum = phylum_clean
  ) %>%
  inner_join(
    auto5_top_phyla %>%
      select(
        Marker,
        Site,
        phylum,
        Rank
      ),
    by = c(
      "Marker",
      "Site",
      "phylum"
    )
  ) %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    Rank,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    Rank,
    name = "Richness"
  ) %>%
  rename(
    Series = phylum
  )

# ============================================================
# 6) COMPLETE ABSENT PHYLA WITH ZERO
# ============================================================
#
# Each sample is crossed only with the five phyla selected for
# its own marker-site combination.
# ============================================================

auto5_richness_grid <- auto5_total_richness %>%
  select(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  distinct() %>%
  inner_join(
    auto5_top_phyla %>%
      select(
        Marker,
        Site,
        phylum,
        Rank
      ),
    by = c(
      "Marker",
      "Site"
    )
  ) %>%
  rename(
    Series = phylum
  )

auto5_selected_richness <- auto5_richness_grid %>%
  left_join(
    auto5_richness_observed,
    by = c(
      "Marker",
      "Site",
      "sample",
      "Axis_type",
      "Axis_value",
      "Series",
      "Rank"
    )
  ) %>%
  mutate(
    Richness = replace_na(
      Richness,
      0L
    )
  )

auto5_richness_plot_data <- bind_rows(
  auto5_total_richness %>%
    mutate(
      Rank = 0L
    ),
  auto5_selected_richness
) %>%
  arrange(
    Marker,
    Site,
    Rank,
    Axis_value
  )

# ============================================================
# 7) PIE-CHART DATA
# ============================================================
#
# The pie chart contains:
# - the five dominant phyla of the corresponding site;
# - Other, containing every other assigned phylum.
#
# Unassigned ASVs are excluded from the denominator.
# ============================================================

auto5_pie_data <- auto5_phylum_summary %>%
  left_join(
    auto5_top_phyla %>%
      select(
        Marker,
        Site,
        phylum
      ) %>%
      mutate(
        Is_top5 = TRUE
      ),
    by = c(
      "Marker",
      "Site",
      "phylum"
    )
  ) %>%
  mutate(
    Pie_group = if_else(
      replace_na(
        Is_top5,
        FALSE
      ),
      phylum,
      "Other"
    )
  ) %>%
  group_by(
    Marker,
    Site,
    Pie_group
  ) %>%
  summarise(
    n_unique_ASVs = sum(
      n_unique_ASVs,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(
    Marker,
    Site
  ) %>%
  mutate(
    Percentage = 100 *
      n_unique_ASVs /
      sum(
        n_unique_ASVs,
        na.rm = TRUE
      ),
    
    Percentage_label = if_else(
      Percentage >= 3,
      paste0(
        round(
          Percentage,
          1
        ),
        "%"
      ),
      ""
    )
  ) %>%
  ungroup()

# ============================================================
# 8) GLOBAL COLOURS
# ============================================================
#
# A phylum receives the same colour in every site and marker.
# ============================================================

auto5_all_phyla <- sort(
  unique(
    auto5_top_phyla$phylum
  )
)

auto5_phylum_colours <- setNames(
  grDevices::hcl.colors(
    n = length(
      auto5_all_phyla
    ),
    palette = "Dark 3"
  ),
  auto5_all_phyla
)

auto5_line_colours <- c(
  "Total richness" = "grey40",
  auto5_phylum_colours
)

auto5_pie_colours <- c(
  auto5_phylum_colours,
  "Other" = "grey80"
)

# ============================================================
# 9) THEMES
# ============================================================

auto5_theme_richness <- function() {
  
  theme_classic(
    base_size = 15,
    base_family = "Arial"
  ) +
    theme(
      axis.title = element_text(
        size = 16,
        face = "bold",
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 22,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      legend.position = "none",
      
      plot.title = element_text(
        size = 25,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 13
      ),
      
      plot.margin = margin(
        12,
        12,
        12,
        12
      )
    )
}


auto5_theme_pie <- function() {
  
  theme_void(
    base_family = "Arial"
  ) +
    theme(
      legend.position = "none",
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}

# ============================================================
# 10) FUNCTION FOR ONE MARKER × SITE
# ============================================================

auto5_make_site_plot <- function(
    marker_name,
    site_name
) {
  
  line_df <- auto5_richness_plot_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  pie_df <- auto5_pie_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  top_phyla_here <- auto5_top_phyla %>%
    filter(
      Marker == marker_name,
      Site == site_name
    ) %>%
    arrange(
      Rank
    ) %>%
    pull(
      phylum
    )
  
  if (
    nrow(line_df) == 0 ||
    length(top_phyla_here) == 0
  ) {
    
    warning(
      "Insufficient data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  axis_type_here <- unique(
    line_df$Axis_type
  )
  
  if (length(axis_type_here) != 1) {
    
    stop(
      marker_name,
      " - ",
      site_name,
      ": more than one axis type was detected."
    )
  }
  
  series_order <- c(
    "Total richness",
    top_phyla_here
  )
  
  line_df <- line_df %>%
    mutate(
      Series = factor(
        Series,
        levels = series_order
      )
    ) %>%
    arrange(
      Series,
      Axis_value
    )
  
  colours_here <- auto5_line_colours[
    series_order
  ]
  
  # ----------------------------------------------------------
  # Temporal richness plot
  # ----------------------------------------------------------
  
  p_line <- ggplot(
    line_df,
    aes(
      x = Axis_value,
      y = Richness,
      colour = Series,
      group = Series
    )
  ) +
    geom_line(
      aes(
        linewidth =
          Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      aes(
        size =
          Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = colours_here,
      breaks = series_order,
      limits = series_order,
      drop = FALSE,
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = c(
        `FALSE` = 1,
        `TRUE` = 1.7
      ),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(
        `FALSE` = 2.2,
        `TRUE` = 3
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(
        base = 10,
        sigma = 1
      ),
      breaks = c(
        0,
        1,
        2,
        5,
        10,
        20,
        50,
        100,
        200,
        500,
        1000,
        2000
      ),
      labels = scales::label_number(),
      expand = expansion(
        mult = c(
          0.02,
          0.08
        )
      )
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(
        n = 7
      )
    ) +
    labs(
      title = paste0(
        marker_name,
        " — ",
        site_name
      ),
      
      subtitle = paste0(
        "Total eukaryotic richness and the five ",
        "most ASV-rich phyla"
      ),
      
      x = ifelse(
        axis_type_here == "Year",
        "Year",
        "Depth (cm)"
      ),
      
      y = "ASV richness (pseudo-log scale)"
    ) +
    auto5_theme_richness()
  
  # ----------------------------------------------------------
  # Pie chart
  # ----------------------------------------------------------
  
  pie_order <- c(
    top_phyla_here,
    "Other"
  )
  
  pie_df <- pie_df %>%
    mutate(
      Pie_group = factor(
        Pie_group,
        levels = pie_order
      )
    ) %>%
    arrange(
      Pie_group
    )
  
  pie_colours_here <- auto5_pie_colours[
    pie_order
  ]
  
  p_pie <- ggplot(
    pie_df,
    aes(
      x = "",
      y = Percentage,
      fill = Pie_group
    )
  ) +
    geom_col(
      width = 1,
      colour = "white",
      linewidth = 0.5
    ) +
    coord_polar(
      theta = "y"
    ) +
    geom_text(
      aes(
        label = scales::number(
          Percentage,
          accuracy = 1
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 6,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = pie_colours_here,
      breaks = pie_order,
      limits = pie_order,
      drop = FALSE,
      guide = "none"
    ) +
    labs(
      title = "Unique ASV composition"
    ) +
    auto5_theme_pie()
  
  p_line + p_pie +
    plot_layout(
      widths = c(
        3.5,
        1.3
      )
    )
}

# ============================================================
# 11) CREATE INDIVIDUAL FIGURES
# ============================================================

auto5_marker_site_combinations <-
  auto5_richness_plot_data %>%
  distinct(
    Marker,
    Site
  ) %>%
  arrange(
    Marker,
    Site
  )

auto5_all_plots <- pmap(
  auto5_marker_site_combinations,
  function(Marker, Site) {
    
    message(
      "Creating automatic top-five plot: ",
      Marker,
      " - ",
      Site
    )
    
    p <- auto5_make_site_plot(
      marker_name = Marker,
      site_name = Site
    )
    
    if (!is.null(p)) {
      
      ggsave(
        filename = file.path(
          auto5_output_dir,
          paste0(
            "Richness_automatic_top5_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".png"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        filename = file.path(
          auto5_output_dir,
          paste0(
            "Richness_automatic_top5_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".svg"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        bg = "white"
      )
    }
    
    p
  }
)

names(auto5_all_plots) <- paste(
  auto5_marker_site_combinations$Marker,
  auto5_marker_site_combinations$Site,
  sep = "_"
)

# Display first valid plot.
auto5_valid_positions <- which(
  !vapply(
    auto5_all_plots,
    is.null,
    logical(1)
  )
)



if (length(auto5_valid_positions) > 0) {
  
  print(
    auto5_all_plots[[auto5_valid_positions[1]]]
  )
}

# ============================================================
# 12) GENERAL LEGEND PER MARKER
# ============================================================
#
# Because the five phyla vary among sites, the general legend
# contains the union of all phyla selected for that marker.
# ============================================================

auto5_make_marker_legend <- function(
    marker_name
) {
  
  marker_phyla <- auto5_top_phyla %>%
    filter(
      Marker == marker_name
    ) %>%
    distinct(
      phylum
    ) %>%
    arrange(
      phylum
    ) %>%
    pull(
      phylum
    )
  
  legend_levels <- c(
    "Total richness",
    marker_phyla,
    "Other"
  )
  
  legend_colours <- c(
    "Total richness" = "grey40",
    auto5_phylum_colours[
      marker_phyla
    ],
    "Other" = "grey80"
  )
  
  legend_df <- tibble(
    Series = factor(
      legend_levels,
      levels = legend_levels
    ),
    x = seq_along(
      legend_levels
    ),
    y = 1
  )
  
  ggplot(
    legend_df,
    aes(
      x = x,
      y = y,
      colour = Series
    )
  ) +
    geom_point(
      size = 4
    ) +
    scale_colour_manual(
      values = legend_colours,
      breaks = legend_levels,
      limits = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          shape = 16,
          size = 4
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 14,
        family = "Arial"
      ),
      
      legend.key.width = unit(
        1.1,
        "cm"
      ),
      
      legend.spacing.x = unit(
        0.2,
        "cm"
      )
    )
}

# ============================================================
# 13) MULTIPANEL FIGURE PER MARKER
# ============================================================

for (
  marker_name in unique(
    auto5_marker_site_combinations$Marker
  )
) {
  
  selected_plot_names <- names(
    auto5_all_plots
  )[
    str_starts(
      names(auto5_all_plots),
      paste0(
        marker_name,
        "_"
      )
    )
  ]
  
  marker_plots <- auto5_all_plots[
    selected_plot_names
  ]
  
  marker_plots <- marker_plots[
    !vapply(
      marker_plots,
      is.null,
      logical(1)
    )
  ]
  
  if (length(marker_plots) > 0) {
    
    n_columns <- 2
    
    n_rows <- ceiling(
      length(marker_plots) /
        n_columns
    )
    
    marker_legend <-
      auto5_make_marker_legend(
        marker_name
      )
    
    site_panel <- wrap_plots(
      marker_plots,
      ncol = n_columns
    ) +
      plot_annotation(
        title = paste0(
          marker_name,
          ": total richness and five most ",
          "ASV-rich phyla per site"
        ),
        
        theme = theme(
          plot.title = element_text(
            size = 30,
            face = "bold",
            family = "Arial",
            hjust = 0.5
          )
        )
      )
    
    combined_marker_plot <-
      site_panel /
      marker_legend +
      plot_layout(
        heights = c(
          20,
          1.5
        )
      )
    
    ggsave(
      filename = file.path(
        auto5_output_dir,
        paste0(
          "Richness_automatic_top5_all_sites_",
          marker_name,
          ".png"
        )
      ),
      plot = combined_marker_plot,
      width = 24,
      height = 7.5 * n_rows + 2.5,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(
        auto5_output_dir,
        paste0(
          "Richness_automatic_top5_all_sites_",
          marker_name,
          ".svg"
        )
      ),
      plot = combined_marker_plot,
      width = 24,
      height = 7.5 * n_rows + 2.5,
      bg = "white",
      limitsize = FALSE
    )
  }
}

# ============================================================
# 14) EXPORT SOURCE DATA
# ============================================================

write_xlsx(
  list(
    top5_phyla_by_site =
      auto5_top_phyla,
    
    richness_time_series =
      auto5_richness_plot_data,
    
    unique_ASVs_by_phylum =
      auto5_phylum_summary,
    
    pie_chart_data =
      auto5_pie_data
  ),
  
  file.path(
    auto5_output_dir,
    "Automatic_top5_phyla_plot_data.xlsx"
  )
)

cat(
  "\nDONE\n",
  "Automatic top-five phylum figures exported to:\n",
  auto5_output_dir,
  "\n"
)



# ============================================================
# INDEPENDENT ANALYSIS:
# KINGDOM-LEVEL RICHNESS THROUGH TIME
#
# Requires an existing object:
# asv_long
#
# Expected columns:
# Marker, Site, sample, Axis_type, Axis_value,
# ASV, domain, kingdom
# ============================================================
# ============================================================
# KINGDOM-LEVEL RICHNESS THROUGH TIME
# Total richness + richness by kingdom + pie chart
#
# Includes "Unassigned"
# One individual figure per Marker × Site
# One combined figure per marker with a general legend
#
# Requires:
#   asv_long
#   out_dir
#   safe_filename()
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(patchwork)
library(scales)
library(writexl)

# ============================================================
# 0) OUTPUT DIRECTORY
# ============================================================

kingdom_output_dir <- file.path(
  out_dir,
  "Taxonomic_richness_all_kingdoms"
)

dir.create(
  kingdom_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 1) CHECK REQUIRED COLUMNS
# ============================================================

required_kingdom_columns <- c(
  "Marker",
  "Site",
  "sample",
  "Axis_type",
  "Axis_value",
  "ASV",
  "kingdom"
)

missing_kingdom_columns <- setdiff(
  required_kingdom_columns,
  names(asv_long)
)

if (length(missing_kingdom_columns) > 0) {
  
  stop(
    "Missing required columns in asv_long: ",
    paste(
      missing_kingdom_columns,
      collapse = ", "
    )
  )
}

# ============================================================
# 2) CLEAN KINGDOM NAMES
# ============================================================

unassigned_kingdom_terms <- c(
  "",
  "unassigned",
  "unknown",
  "unclassified",
  "unidentified",
  "na",
  "n/a"
)

kingdom_asv_long <- asv_long %>%
  mutate(
    kingdom_clean = str_trim(
      as.character(kingdom)
    ),
    
    kingdom_clean = case_when(
      is.na(kingdom_clean) ~ "Unassigned",
      
      str_to_lower(kingdom_clean) %in%
        unassigned_kingdom_terms ~ "Unassigned",
      
      TRUE ~ kingdom_clean
    )
  )

# Keep only Eukaryota when domain is available
if ("domain" %in% names(kingdom_asv_long)) {
  
  kingdom_asv_long <- kingdom_asv_long %>%
    filter(
      str_to_lower(
        str_trim(
          as.character(domain)
        )
      ) == "eukaryota"
    )
}

# ============================================================
# 3) LIST ALL KINGDOMS, INCLUDING UNASSIGNED
# ============================================================

all_kingdoms <- kingdom_asv_long %>%
  distinct(
    kingdom_clean
  ) %>%
  arrange(
    kingdom_clean == "Unassigned",
    kingdom_clean
  ) %>%
  pull(
    kingdom_clean
  )

cat(
  "\nKingdoms detected in the dataset:\n"
)

print(all_kingdoms)

# ============================================================
# 4) TOTAL EUKARYOTIC RICHNESS PER SAMPLE
# ============================================================

kingdom_total_richness <- kingdom_asv_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    name = "Richness"
  ) %>%
  mutate(
    Series = "Total richness"
  )

# ============================================================
# 5) OBSERVED RICHNESS PER KINGDOM AND SAMPLE
# ============================================================
#
# "Unassigned" is retained as a kingdom category.
# ============================================================

kingdom_richness_observed <- kingdom_asv_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    kingdom_clean,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    kingdom_clean,
    name = "Richness"
  ) %>%
  rename(
    Series = kingdom_clean
  )

# ============================================================
# 6) KINGDOMS PRESENT IN EACH MARKER × SITE
# ============================================================

kingdoms_by_site <- kingdom_asv_long %>%
  distinct(
    Marker,
    Site,
    kingdom_clean
  ) %>%
  rename(
    Series = kingdom_clean
  )

# ============================================================
# 7) COMPLETE ABSENT KINGDOMS WITH ZERO
# ============================================================
#
# A kingdom is included in all samples of a site if it occurs
# at least once in that Marker × Site combination.
# ============================================================

kingdom_richness_grid <- kingdom_total_richness %>%
  select(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  distinct() %>%
  inner_join(
    kingdoms_by_site,
    by = c(
      "Marker",
      "Site"
    )
  )

kingdom_richness_complete <- kingdom_richness_grid %>%
  left_join(
    kingdom_richness_observed,
    by = c(
      "Marker",
      "Site",
      "sample",
      "Axis_type",
      "Axis_value",
      "Series"
    )
  ) %>%
  mutate(
    Richness = replace_na(
      Richness,
      0L
    )
  )

kingdom_richness_plot_data <- bind_rows(
  kingdom_total_richness,
  kingdom_richness_complete
) %>%
  arrange(
    Marker,
    Site,
    Series,
    Axis_value
  )

# ============================================================
# 8) KINGDOM PIE-CHART DATA
# ============================================================
#
# Percentages are calculated from unique ASVs detected across
# the complete record of each Marker × Site combination.
#
# "Unassigned" is included in both the pie chart and denominator.
# ============================================================

kingdom_pie_data <- kingdom_asv_long %>%
  distinct(
    Marker,
    Site,
    kingdom_clean,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    kingdom_clean,
    name = "n_unique_ASVs"
  ) %>%
  group_by(
    Marker,
    Site
  ) %>%
  mutate(
    Percentage = 100 *
      n_unique_ASVs /
      sum(
        n_unique_ASVs,
        na.rm = TRUE
      ),
    
    Percentage_label = if_else(
      Percentage >= 3,
      paste0(
        round(
          Percentage,
          1
        ),
        "%"
      ),
      ""
    )
  ) %>%
  ungroup() %>%
  rename(
    Kingdom = kingdom_clean
  )

# Check that pie-chart percentages sum to 100
kingdom_pie_check <- kingdom_pie_data %>%
  group_by(
    Marker,
    Site
  ) %>%
  summarise(
    Percentage_sum = sum(
      Percentage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(kingdom_pie_check)

# ============================================================
# 9) COLOURS
# ============================================================
#
# Assigned kingdoms receive colours from Dark 3.
# Unassigned is always light grey.
# Total richness is always dark grey.
# ============================================================

assigned_kingdoms <- setdiff(
  all_kingdoms,
  "Unassigned"
)

assigned_kingdom_colours <- setNames(
  grDevices::hcl.colors(
    n = max(
      length(assigned_kingdoms),
      1
    ),
    palette = "Dark 3"
  )[seq_along(assigned_kingdoms)],
  assigned_kingdoms
)

kingdom_colours <- c(
  assigned_kingdom_colours,
  "Unassigned" = "grey75"
)

kingdom_line_colours <- c(
  "Total richness" = "grey30",
  kingdom_colours
)

# ============================================================
# 10) THEMES
# ============================================================

theme_kingdom_richness <- function() {
  
  theme_classic(
    base_size = 25,
    base_family = "Arial"
  ) +
    theme(
      axis.title = element_text(
        size = 32,
        face = "bold",
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 30,
        colour = "black"
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      # Legends are removed from individual panels.
      legend.position = "none",
      
      plot.title = element_text(
        size = 34,
        face = "bold"
      ),
      
      plot.subtitle = element_blank(),
      
      plot.margin = margin(
        18,
        34,
        16,
        16
      )
    )
}


theme_kingdom_pie <- function() {
  theme_void(base_family = "Arial") +
    theme(
      text = element_text(
        family = "Arial",
        face = "bold",
        size = 26
      ),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}

# ============================================================
# 11) FUNCTION:
# KINGDOM RICHNESS THROUGH TIME + PIE CHART
# ============================================================

make_kingdom_richness_plot <- function(
    marker_name,
    site_name
) {
  
  line_df <- kingdom_richness_plot_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  pie_df <- kingdom_pie_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  if (nrow(line_df) == 0) {
    
    warning(
      "No kingdom richness data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  if (nrow(pie_df) == 0) {
    
    warning(
      "No kingdom pie-chart data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  axis_type_here <- unique(
    line_df$Axis_type
  )
  
  if (length(axis_type_here) != 1) {
    
    stop(
      marker_name,
      " - ",
      site_name,
      ": more than one axis type detected."
    )
  }
  
  kingdoms_here <- kingdoms_by_site %>%
    filter(
      Marker == marker_name,
      Site == site_name
    ) %>%
    arrange(
      Series == "Unassigned",
      Series
    ) %>%
    pull(
      Series
    )
  
  series_order <- c(
    "Total richness",
    kingdoms_here
  )
  
  line_df <- line_df %>%
    mutate(
      Series = factor(
        Series,
        levels = series_order
      )
    ) %>%
    arrange(
      Series,
      Axis_value
    )
  
  colours_here <- kingdom_line_colours[
    series_order
  ]
  
  # ----------------------------------------------------------
  # Temporal richness plot
  # ----------------------------------------------------------
  
  p_line <- ggplot(
    line_df,
    aes(
      x = Axis_value,
      y = Richness,
      colour = Series,
      group = Series
    )
  ) +
    geom_line(
      aes(
        linewidth =
          Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      aes(
        size =
          Series == "Total richness"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = colours_here,
      breaks = series_order,
      limits = series_order,
      drop = FALSE,
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = c(
        `FALSE` = 1.6,
        `TRUE` = 2.5
      ),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(
        `FALSE` = 2.8,
        `TRUE` = 3.8
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(
        base = 10,
        sigma = 1
      ),
      breaks = c(
        0,
        1,
        2,
        5,
        10,
        20,
        50,
        100,
        200,
        500,
        1000,
        2000
      ),
      labels = scales::label_number(),
      expand = expansion(
        mult = c(
          0.02,
          0.08
        )
      )
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(
        n = 7
      )
    ) +
    labs(
      title = site_name,
      subtitle = NULL,
      
      x = ifelse(
        axis_type_here == "Year",
        "Year",
        "Depth (cm)"
      ),
      
      y = "ASV richness"
    ) +
    theme_kingdom_richness()
  
  # ----------------------------------------------------------
  # Pie chart
  # ----------------------------------------------------------
  
  pie_order <- kingdoms_here
  
  pie_df <- pie_df %>%
    mutate(
      Kingdom = factor(
        Kingdom,
        levels = pie_order
      )
    ) %>%
    arrange(
      Kingdom
    )
  
  pie_colours_here <- kingdom_colours[
    pie_order
  ]
  
  p_pie <- ggplot(
    pie_df,
    aes(
      x = "",
      y = Percentage,
      fill = Kingdom
    )
  ) +
    geom_col(
      width = 1,
      colour = "white",
      linewidth = 0.5
    ) +
    coord_polar(
      theta = "y"
    ) +
    geom_text(
      aes(
        label = scales::number(
          Percentage,
          accuracy = 1
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 6,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = pie_colours_here,
      breaks = pie_order,
      limits = pie_order,
      drop = FALSE,
      guide = "none"
    ) +
    labs(title = NULL) +
    theme_kingdom_pie()
  
  # ----------------------------------------------------------
  # Place a large pie chart at the upper-right corner of each timeline.
  # The inset deliberately crosses the panel boundary slightly.
  # ----------------------------------------------------------
  
  combined_plot <- p_line +
    patchwork::inset_element(
      p_pie,
      left = 0.69,
      bottom = 0.58,
      right = 1.07,
      top = 1.07,
      align_to = "panel",
      clip = FALSE,
      on_top = TRUE
    )
  
  combined_plot
}

# ============================================================
# 12) CREATE INDIVIDUAL FIGURES
# ============================================================

kingdom_combinations <- kingdom_richness_plot_data %>%
  distinct(
    Marker,
    Site
  ) %>%
  arrange(
    Marker,
    Site
  )

kingdom_richness_plots <- pmap(
  kingdom_combinations,
  function(Marker, Site) {
    
    message(
      "Creating kingdom richness plot: ",
      Marker,
      " - ",
      Site
    )
    
    p <- make_kingdom_richness_plot(
      marker_name = Marker,
      site_name = Site
    )
    
    if (!is.null(p)) {
      
      ggsave(
        filename = file.path(
          kingdom_output_dir,
          paste0(
            "Kingdom_richness_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".png"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        filename = file.path(
          kingdom_output_dir,
          paste0(
            "Kingdom_richness_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".svg"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        bg = "white"
      )
    }
    
    p
  }
)

names(kingdom_richness_plots) <- paste(
  kingdom_combinations$Marker,
  kingdom_combinations$Site,
  sep = "_"
)

# Display first valid plot
valid_kingdom_plot_positions <- which(
  !vapply(
    kingdom_richness_plots,
    is.null,
    logical(1)
  )
)


if (length(valid_kingdom_plot_positions) > 0) {
  
  print(
    kingdom_richness_plots[[valid_kingdom_plot_positions[1]]]
  )
}

# ============================================================
# 13) GENERAL LEGEND FOR EACH MARKER
# ============================================================

make_kingdom_marker_legend <- function(
    marker_name
) {
  
  kingdoms_marker <- kingdom_asv_long %>%
    filter(
      Marker == marker_name
    ) %>%
    distinct(
      kingdom_clean
    ) %>%
    arrange(
      kingdom_clean == "Unassigned",
      kingdom_clean
    ) %>%
    pull(
      kingdom_clean
    )
  
  legend_levels <- c(
    "Total richness",
    kingdoms_marker
  )
  
  legend_colours <- c(
    "Total richness" = "grey30",
    kingdom_colours[
      kingdoms_marker
    ]
  )
  
  legend_df <- tibble(
    Series = factor(
      legend_levels,
      levels = legend_levels
    ),
    
    x = seq_along(
      legend_levels
    ),
    
    y = 1
  )
  
  ggplot(
    legend_df,
    aes(
      x = x,
      y = y,
      colour = Series
    )
  ) +
    geom_point(
      size = 4
    ) +
    scale_colour_manual(
      values = legend_colours,
      breaks = legend_levels,
      limits = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        
        override.aes = list(
          shape = 16,
          size = 4,
          linewidth = 1.2
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 14,
        family = "Arial"
      ),
      
      legend.key.width = unit(
        1.1,
        "cm"
      ),
      
      legend.spacing.x = unit(
        0.2,
        "cm"
      ),
      
      plot.margin = margin(
        0,
        5,
        0,
        5
      )
    )
}

# ============================================================
# 14) COMBINED MULTIPANEL PER MARKER
# ONE GENERAL LEGEND
# ============================================================

for (
  marker_name in unique(
    kingdom_combinations$Marker
  )
) {
  
  selected_names <- names(
    kingdom_richness_plots
  )[
    str_starts(
      names(kingdom_richness_plots),
      paste0(
        marker_name,
        "_"
      )
    )
  ]
  
  marker_plots <- kingdom_richness_plots[
    selected_names
  ]
  
  marker_plots <- marker_plots[
    !vapply(
      marker_plots,
      is.null,
      logical(1)
    )
  ]
  
  if (length(marker_plots) > 0) {
    
    n_columns <- 2
    
    n_rows <- ceiling(
      length(marker_plots) /
        n_columns
    )
    
    # Explicitly remove legends from every individual plot
    marker_plots_no_legend <- lapply(
      marker_plots,
      function(p) {
        
        p & theme(
          legend.position = "none"
        )
      }
    )
    
    marker_legend <- make_kingdom_marker_legend(
      marker_name
    )
    
    site_panel <- wrap_plots(
      marker_plots_no_legend,
      ncol = n_columns
    ) +
      plot_annotation(
        title = paste0(
          marker_name,
          ": kingdom-level richness by site"
        ),
        
        theme = theme(
          plot.title = element_text(
            size = 30,
            face = "bold",
            family = "Arial",
            hjust = 0.5
          )
        )
      )
    
    combined_richness <-
      site_panel /
      marker_legend +
      plot_layout(
        heights = c(
          20,
          1.4
        )
      )
    
    ggsave(
      filename = file.path(
        kingdom_output_dir,
        paste0(
          "Kingdom_richness_all_sites_",
          marker_name,
          ".png"
        )
      ),
      plot = combined_richness,
      width = 24,
      height = 7.2 * n_rows + 2.5,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(
        kingdom_output_dir,
        paste0(
          "Kingdom_richness_all_sites_",
          marker_name,
          ".svg"
        )
      ),
      plot = combined_richness,
      width = 24,
      height = 7.2 * n_rows + 2.5,
      bg = "white",
      limitsize = FALSE
    )
  }
}

# ============================================================
# 15) EXPORT SOURCE DATA
# ============================================================

write_xlsx(
  list(
    kingdom_richness_time_series =
      kingdom_richness_plot_data,
    
    kingdoms_by_site =
      kingdoms_by_site,
    
    kingdom_pie_chart_data =
      kingdom_pie_data,
    
    kingdom_pie_percentage_check =
      kingdom_pie_check
  ),
  
  file.path(
    kingdom_output_dir,
    "Kingdom_level_richness_plot_data.xlsx"
  )
)

cat(
  "\nDONE\n",
  "Kingdom-level figures and source data exported to:\n",
  kingdom_output_dir,
  "\n"
)

# ============================================================
# 14) EXPORT SOURCE DATA
# ============================================================

write_xlsx(
  list(
    kingdom_richness_time_series =
      kingdom_richness_plot_data,
    
    kingdom_composition =
      kingdom_composition_data,
    
    kingdoms_by_site =
      kingdoms_by_site
  ),
  file.path(
    kingdom_output_dir,
    "Kingdom_level_richness_plot_data.xlsx"
  )
)

cat(
  "\nDONE\n",
  "Kingdom-level figures and source data exported to:\n",
  kingdom_output_dir,
  "\n"
)




# ============================================================
# METAZOAN RICHNESS + TOP 5 PHYLA PER SITE
#
# Outputs:
# 1. Individual timeline + pie chart per Marker × Site
# 2. Combined multipanel per marker with one common legend
#
# Unassigned phyla are excluded from:
# - total metazoan richness
# - top-five selection
# - timelines
# - pie-chart denominator
# - legends
#
# Requires:
#   asv_long
#   out_dir
#   safe_filename()
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(patchwork)
library(scales)
library(writexl)

# ============================================================
# 0) OUTPUT DIRECTORY
# ============================================================

metazoa_output_dir <- file.path(
  out_dir,
  "Metazoa_richness_top5_phyla"
)

dir.create(
  metazoa_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 1) CHECK REQUIRED COLUMNS
# ============================================================

required_metazoa_columns <- c(
  "Marker",
  "Site",
  "sample",
  "Axis_type",
  "Axis_value",
  "ASV",
  "kingdom",
  "phylum"
)

missing_metazoa_columns <- setdiff(
  required_metazoa_columns,
  names(asv_long)
)

if (length(missing_metazoa_columns) > 0) {
  
  stop(
    "Missing required columns in asv_long: ",
    paste(
      missing_metazoa_columns,
      collapse = ", "
    )
  )
}

# ============================================================
# 2) FILTER METAZOA AND CLEAN PHYLA
# ============================================================

unassigned_phylum_terms <- c(
  "",
  "unassigned",
  "unknown",
  "unclassified",
  "unidentified",
  "incertae sedis",
  "na",
  "n/a"
)

metazoa_asv_long <- asv_long %>%
  mutate(
    kingdom_clean = str_to_lower(
      str_trim(
        as.character(kingdom)
      )
    ),
    
    phylum_clean = str_trim(
      as.character(phylum)
    )
  ) %>%
  filter(
    kingdom_clean %in% c(
      "metazoa",
      "animalia"
    )
  ) %>%
  filter(
    !is.na(phylum_clean),
    !str_to_lower(phylum_clean) %in%
      unassigned_phylum_terms
  )

if (nrow(metazoa_asv_long) == 0) {
  
  stop(
    paste0(
      "No assigned metazoan ASVs were found. ",
      "Check the kingdom and phylum assignments."
    )
  )
}

cat(
  "\nNumber of positive, phylum-assigned metazoan ASV detections:",
  nrow(metazoa_asv_long),
  "\n"
)

cat(
  "\nMetazoan phyla retained:\n"
)

print(
  metazoa_asv_long %>%
    distinct(phylum_clean) %>%
    arrange(phylum_clean),
  n = Inf
)

# ============================================================
# 3) TOTAL ASSIGNED METAZOAN RICHNESS PER SAMPLE
# ============================================================
#
# Unassigned metazoan ASVs have already been removed.
# Therefore Total Metazoa is the total richness of ASVs
# assigned to a metazoan phylum.
# ============================================================

metazoa_total_richness <- metazoa_asv_long %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    name = "Richness"
  ) %>%
  mutate(
    Series = "Total Metazoa"
  )

# ============================================================
# 4) UNIQUE ASVs PER PHYLUM, MARKER AND SITE
# ============================================================

metazoa_phylum_summary <- metazoa_asv_long %>%
  distinct(
    Marker,
    Site,
    phylum_clean,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    phylum_clean,
    name = "n_unique_ASVs"
  ) %>%
  rename(
    phylum = phylum_clean
  )

# ============================================================
# 5) TOP FIVE PHYLA PER MARKER AND SITE
# ============================================================
#
# Ranking criterion:
# number of unique ASVs detected at least once across the
# complete record of each Marker × Site combination.
# ============================================================

metazoa_top5_phyla <- metazoa_phylum_summary %>%
  group_by(
    Marker,
    Site
  ) %>%
  arrange(
    desc(n_unique_ASVs),
    phylum,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 5
  ) %>%
  mutate(
    Rank = row_number()
  ) %>%
  ungroup()

cat(
  "\nTop five metazoan phyla per marker and site:\n"
)

print(
  metazoa_top5_phyla,
  n = Inf
)

# ============================================================
# 6) OBSERVED RICHNESS OF TOP-FIVE PHYLA PER SAMPLE
# ============================================================

metazoa_top5_observed <- metazoa_asv_long %>%
  mutate(
    phylum = phylum_clean
  ) %>%
  inner_join(
    metazoa_top5_phyla %>%
      select(
        Marker,
        Site,
        phylum,
        Rank
      ),
    by = c(
      "Marker",
      "Site",
      "phylum"
    )
  ) %>%
  distinct(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    Rank,
    ASV
  ) %>%
  count(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value,
    phylum,
    Rank,
    name = "Richness"
  ) %>%
  rename(
    Series = phylum
  )

# ============================================================
# 7) COMPLETE TRUE ABSENCES WITH ZERO
# ============================================================

metazoa_top5_grid <- metazoa_total_richness %>%
  select(
    Marker,
    Site,
    sample,
    Axis_type,
    Axis_value
  ) %>%
  distinct() %>%
  inner_join(
    metazoa_top5_phyla %>%
      select(
        Marker,
        Site,
        phylum,
        Rank
      ),
    by = c(
      "Marker",
      "Site"
    )
  ) %>%
  rename(
    Series = phylum
  )

metazoa_top5_richness <- metazoa_top5_grid %>%
  left_join(
    metazoa_top5_observed,
    by = c(
      "Marker",
      "Site",
      "sample",
      "Axis_type",
      "Axis_value",
      "Series",
      "Rank"
    )
  ) %>%
  mutate(
    Richness = replace_na(
      Richness,
      0L
    )
  )

metazoa_richness_plot_data <- bind_rows(
  metazoa_total_richness %>%
    mutate(
      Rank = 0L
    ),
  
  metazoa_top5_richness
) %>%
  arrange(
    Marker,
    Site,
    Rank,
    Axis_value
  )

# ============================================================
# 8) PIE-CHART DATA
# ============================================================
#
# The denominator includes all assigned metazoan ASVs.
#
# The five dominant phyla are represented individually.
# All remaining assigned phyla are grouped as "Other Metazoa".
# ============================================================

metazoa_pie_data <- metazoa_phylum_summary %>%
  left_join(
    metazoa_top5_phyla %>%
      select(
        Marker,
        Site,
        phylum
      ) %>%
      mutate(
        Is_top5 = TRUE
      ),
    by = c(
      "Marker",
      "Site",
      "phylum"
    )
  ) %>%
  mutate(
    Pie_group = if_else(
      replace_na(
        Is_top5,
        FALSE
      ),
      phylum,
      "Other Metazoa"
    )
  ) %>%
  group_by(
    Marker,
    Site,
    Pie_group
  ) %>%
  summarise(
    n_unique_ASVs = sum(
      n_unique_ASVs,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(
    Marker,
    Site
  ) %>%
  mutate(
    Percentage = 100 *
      n_unique_ASVs /
      sum(
        n_unique_ASVs,
        na.rm = TRUE
      ),
    
    Percentage_label = if_else(
      Percentage >= 3,
      paste0(
        round(
          Percentage,
          1
        ),
        "%"
      ),
      ""
    )
  ) %>%
  ungroup()

# Check that pie-chart percentages sum to 100
metazoa_pie_check <- metazoa_pie_data %>%
  group_by(
    Marker,
    Site
  ) %>%
  summarise(
    Percentage_sum = sum(
      Percentage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(metazoa_pie_check)

# ============================================================
# 9) GLOBAL COLOURS
# ============================================================
#
# The same phylum receives the same colour across all sites
# and both markers.
# ============================================================

metazoa_all_top_phyla <- sort(
  unique(
    metazoa_top5_phyla$phylum
  )
)

metazoa_phylum_colours <- setNames(
  grDevices::hcl.colors(
    n = max(
      length(metazoa_all_top_phyla),
      1
    ),
    palette = "Dark 3"
  )[seq_along(metazoa_all_top_phyla)],
  metazoa_all_top_phyla
)

metazoa_line_colours <- c(
  "Total Metazoa" = "grey35",
  metazoa_phylum_colours
)

metazoa_pie_colours <- c(
  metazoa_phylum_colours,
  "Other Metazoa" = "grey80"
)

# ============================================================
# 10) THEMES
# ============================================================

metazoa_theme_richness <- function() {
  
  theme_classic(
    base_size = 25,
    base_family = "Arial"
  ) +
    theme(
      axis.title = element_text(
        size = 32,
        face = "bold",
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 30,
        colour = "black"
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      legend.position = "none",
      
      plot.title = element_text(
        size = 34,
        face = "bold"
      ),
      
      plot.subtitle = element_blank(),
      
      plot.margin = margin(
        18,
        34,
        16,
        16
      )
    )
}

metazoa_theme_pie <- function() {
  theme_void(base_family = "Arial") +
    theme(
      text = element_text(
        family = "Arial",
        face = "bold",
        size = 18
      ),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}

# ============================================================
# 11) CREATE ONE TIMELINE + PIE CHART
# ============================================================

make_metazoa_site_plot <- function(
    marker_name,
    site_name
) {
  
  line_df <- metazoa_richness_plot_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  pie_df <- metazoa_pie_data %>%
    filter(
      Marker == marker_name,
      Site == site_name
    )
  
  top_phyla_here <- metazoa_top5_phyla %>%
    filter(
      Marker == marker_name,
      Site == site_name
    ) %>%
    arrange(
      Rank
    ) %>%
    pull(
      phylum
    )
  
  if (
    nrow(line_df) == 0 ||
    nrow(pie_df) == 0 ||
    length(top_phyla_here) == 0
  ) {
    
    warning(
      "Insufficient assigned metazoan data for ",
      marker_name,
      " - ",
      site_name
    )
    
    return(NULL)
  }
  
  axis_type_here <- unique(
    line_df$Axis_type
  )
  
  if (length(axis_type_here) != 1) {
    
    stop(
      marker_name,
      " - ",
      site_name,
      ": more than one axis type detected."
    )
  }
  
  series_order <- c(
    "Total Metazoa",
    top_phyla_here
  )
  
  line_df <- line_df %>%
    mutate(
      Series = factor(
        Series,
        levels = series_order
      )
    ) %>%
    arrange(
      Series,
      Axis_value
    )
  
  colours_here <- metazoa_line_colours[
    series_order
  ]
  
  # ----------------------------------------------------------
  # Timeline
  # ----------------------------------------------------------
  
  p_line <- ggplot(
    line_df,
    aes(
      x = Axis_value,
      y = Richness,
      colour = Series,
      group = Series
    )
  ) +
    geom_line(
      aes(
        linewidth =
          Series == "Total Metazoa"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      aes(
        size =
          Series == "Total Metazoa"
      ),
      alpha = 0.9,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = colours_here,
      breaks = series_order,
      limits = series_order,
      drop = FALSE,
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = c(
        `FALSE` = 1.6,
        `TRUE` = 2.5
      ),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(
        `FALSE` = 2.8,
        `TRUE` = 3.8
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(
        base = 10,
        sigma = 1
      ),
      breaks = c(
        0,
        1,
        2,
        5,
        10,
        20,
        50,
        100,
        200,
        500,
        1000,
        2000
      ),
      labels = scales::label_number(),
      expand = expansion(
        mult = c(
          0.02,
          0.08
        )
      )
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(
        n = 7
      )
    ) +
    labs(
      title = site_name,
      subtitle = NULL,
      
      x = ifelse(
        axis_type_here == "Year",
        "Year",
        "Depth (cm)"
      ),
      
      y = "ASV richness"
    ) +
    metazoa_theme_richness()
  
  # ----------------------------------------------------------
  # Pie chart
  # ----------------------------------------------------------
  
  pie_order <- c(
    top_phyla_here,
    "Other Metazoa"
  )
  
  pie_df <- pie_df %>%
    mutate(
      Pie_group = factor(
        Pie_group,
        levels = pie_order
      )
    ) %>%
    arrange(
      Pie_group
    )
  
  pie_colours_here <- metazoa_pie_colours[
    pie_order
  ]
  
  p_pie <- ggplot(
    pie_df,
    aes(
      x = "",
      y = Percentage,
      fill = Pie_group
    )
  ) +
    geom_col(
      width = 1,
      colour = "white",
      linewidth = 0.5
    ) +
    coord_polar(
      theta = "y"
    ) +
    geom_text(
      aes(
        label = scales::number(
          Percentage,
          accuracy = 1
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 6,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = pie_colours_here,
      breaks = pie_order,
      limits = pie_order,
      drop = FALSE,
      guide = "none"
    ) +
    labs(title = NULL) +
    metazoa_theme_pie()
  
  p_line +
    patchwork::inset_element(
      p_pie,
      left = 0.69,
      bottom = 0.58,
      right = 1.07,
      top = 1.07,
      align_to = "panel",
      clip = FALSE,
      on_top = TRUE
    )
}

# ============================================================
# 12) CREATE INDIVIDUAL FIGURES
# ============================================================

metazoa_combinations <- metazoa_richness_plot_data %>%
  distinct(
    Marker,
    Site
  ) %>%
  arrange(
    Marker,
    Site
  )

metazoa_all_plots <- pmap(
  metazoa_combinations,
  function(Marker, Site) {
    
    message(
      "Creating Metazoa timeline: ",
      Marker,
      " - ",
      Site
    )
    
    p <- make_metazoa_site_plot(
      marker_name = Marker,
      site_name = Site
    )
    
    if (!is.null(p)) {
      
      ggsave(
        filename = file.path(
          metazoa_output_dir,
          paste0(
            "Metazoa_richness_top5_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".png"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        filename = file.path(
          metazoa_output_dir,
          paste0(
            "Metazoa_richness_top5_",
            safe_filename(Marker),
            "_",
            safe_filename(Site),
            ".svg"
          )
        ),
        plot = p,
        width = 14,
        height = 7.5,
        bg = "white"
      )
    }
    
    p
  }
)

names(metazoa_all_plots) <- paste(
  metazoa_combinations$Marker,
  metazoa_combinations$Site,
  sep = "_"
)

# Display first valid plot
metazoa_valid_positions <- which(
  !vapply(
    metazoa_all_plots,
    is.null,
    logical(1)
  )
)


if (length(metazoa_valid_positions) > 0) {
  
  print(
    metazoa_all_plots[[metazoa_valid_positions[1]]]
  )
}

# ============================================================
# 13) COMMON LEGEND PER MARKER
# ============================================================
#
# The legend contains the union of all phyla selected among
# the sites of that marker.
# ============================================================

make_metazoa_marker_legend <- function(
    marker_name
) {
  
  marker_phyla <- metazoa_top5_phyla %>%
    filter(
      Marker == marker_name
    ) %>%
    distinct(
      phylum
    ) %>%
    arrange(
      phylum
    ) %>%
    pull(
      phylum
    )
  
  legend_levels <- c(
    "Total Metazoa",
    marker_phyla,
    "Other Metazoa"
  )
  
  legend_colours <- c(
    "Total Metazoa" = "grey35",
    
    metazoa_phylum_colours[
      marker_phyla
    ],
    
    "Other Metazoa" = "grey80"
  )
  
  legend_df <- tibble(
    Series = factor(
      legend_levels,
      levels = legend_levels
    ),
    
    x = seq_along(
      legend_levels
    ),
    
    y = 1
  )
  
  ggplot(
    legend_df,
    aes(
      x = x,
      y = y,
      colour = Series
    )
  ) +
    geom_point(
      size = 4
    ) +
    scale_colour_manual(
      values = legend_colours,
      breaks = legend_levels,
      limits = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        
        override.aes = list(
          shape = 16,
          size = 4,
          linewidth = 1.2
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 14,
        family = "Arial"
      ),
      
      legend.key.width = unit(
        1.1,
        "cm"
      ),
      
      legend.spacing.x = unit(
        0.2,
        "cm"
      ),
      
      plot.margin = margin(
        0,
        5,
        0,
        5
      )
    )
}

# ============================================================
# 14) COMBINED MULTIPANEL PER MARKER
# ============================================================

for (
  marker_name in unique(
    metazoa_combinations$Marker
  )
) {
  
  selected_plot_names <- names(
    metazoa_all_plots
  )[
    str_starts(
      names(metazoa_all_plots),
      paste0(
        marker_name,
        "_"
      )
    )
  ]
  
  marker_plots <- metazoa_all_plots[
    selected_plot_names
  ]
  
  marker_plots <- marker_plots[
    !vapply(
      marker_plots,
      is.null,
      logical(1)
    )
  ]
  
  if (length(marker_plots) > 0) {
    
    n_columns <- 2
    
    n_rows <- ceiling(
      length(marker_plots) /
        n_columns
    )
    
    # Remove any individual legends
    marker_plots_no_legend <- lapply(
      marker_plots,
      function(p) {
        
        p & theme(
          legend.position = "none"
        )
      }
    )
    
    marker_legend <- make_metazoa_marker_legend(
      marker_name
    )
    
    site_panel <- wrap_plots(
      marker_plots_no_legend,
      ncol = n_columns
    ) +
      plot_annotation(
        title = paste0(
          marker_name,
          ": assigned metazoan richness and ",
          "five dominant phyla per site"
        ),
        
        theme = theme(
          plot.title = element_text(
            size = 30,
            face = "bold",
            family = "Arial",
            hjust = 0.5
          )
        )
      )
    
    combined_marker_plot <-
      site_panel /
      marker_legend +
      plot_layout(
        heights = c(
          20,
          1.5
        )
      )
    
    ggsave(
      filename = file.path(
        metazoa_output_dir,
        paste0(
          "Metazoa_richness_top5_all_sites_",
          marker_name,
          ".png"
        )
      ),
      plot = combined_marker_plot,
      width = 24,
      height = 7.5 * n_rows + 2.5,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(
        metazoa_output_dir,
        paste0(
          "Metazoa_richness_top5_all_sites_",
          marker_name,
          ".svg"
        )
      ),
      plot = combined_marker_plot,
      width = 24,
      height = 7.5 * n_rows + 2.5,
      bg = "white",
      limitsize = FALSE
    )
  }
}

# ============================================================
# 15) EXPORT SOURCE DATA
# ============================================================

write_xlsx(
  list(
    top5_metazoan_phyla_by_site =
      metazoa_top5_phyla,
    
    metazoan_richness_time_series =
      metazoa_richness_plot_data,
    
    unique_ASVs_by_metazoan_phylum =
      metazoa_phylum_summary,
    
    metazoan_pie_chart_data =
      metazoa_pie_data,
    
    metazoan_pie_percentage_check =
      metazoa_pie_check
  ),
  
  file.path(
    metazoa_output_dir,
    "Metazoa_richness_top5_phyla_plot_data.xlsx"
  )
)

cat(
  "\nDONE\n",
  "Metazoan timelines, pie charts and source data exported to:\n",
  metazoa_output_dir,
  "\n"
)


# ============================================================
# FINAL "BOTH" MOSAICS: COI LEFT, 18S RIGHT
# ============================================================
# Each marker occupies a 2-column x 4-row block. Sites follow the
# manuscript order. The pie chart is already embedded in the upper-right
# corner of every site timeline by the plotting functions above.

final_site_order <- c(
  "DEE", "STO", "CCO", "CLR",
  "CSD", "CMB", "HIT", "HIO"
)

ordered_marker_plots <- function(plot_list, marker_name, site_order) {
  requested_names <- paste(marker_name, site_order, sep = "_")
  
  lapply(
    requested_names,
    function(plot_name) {
      candidate <- plot_list[[plot_name]]
      
      if (is.null(candidate)) {
        warning("Missing plot: ", plot_name)
        return(patchwork::plot_spacer())
      }
      
      candidate & ggplot2::theme(legend.position = "none")
    }
  )
}

make_marker_block <- function(
    plot_list,
    marker_name,
    site_order,
    legend_plot,
    block_title
) {
  site_plots <- ordered_marker_plots(
    plot_list = plot_list,
    marker_name = marker_name,
    site_order = site_order
  )
  
  site_grid <- patchwork::wrap_plots(
    site_plots,
    ncol = 2,
    nrow = 4,
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = block_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 40,
          face = "bold",
          family = "Arial",
          hjust = 0.5
        )
      )
    )
  
  site_grid /
    legend_plot +
    patchwork::plot_layout(
      heights = c(24, 2.2)
    )
}

# ------------------------------------------------------------
# Kingdom-level BOTH mosaic
# ------------------------------------------------------------

kingdom_COI_block <- make_marker_block(
  plot_list = kingdom_richness_plots,
  marker_name = "COI",
  site_order = final_site_order,
  legend_plot = make_kingdom_marker_legend("COI"),
  block_title = "COI"
)

kingdom_18S_block <- make_marker_block(
  plot_list = kingdom_richness_plots,
  marker_name = "18S",
  site_order = final_site_order,
  legend_plot = make_kingdom_marker_legend("18S"),
  block_title = "18S"
)

kingdom_Both <- kingdom_COI_block |
  kingdom_18S_block

print(kingdom_Both)

ggplot2::ggsave(
  filename = file.path(
    kingdom_output_dir,
    "Kingdom_richness_all_sites_Both.png"
  ),
  plot = kingdom_Both,
  width = 44,
  height = 32,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    kingdom_output_dir,
    "Kingdom_richness_all_sites_Both.svg"
  ),
  plot = kingdom_Both,
  width = 44,
  height = 32,
  bg = "white",
  limitsize = FALSE
)

# ------------------------------------------------------------
# Metazoa top-five-phylum BOTH mosaic
# ------------------------------------------------------------

metazoa_COI_block <- make_marker_block(
  plot_list = metazoa_all_plots,
  marker_name = "COI",
  site_order = final_site_order,
  legend_plot = make_metazoa_marker_legend("COI"),
  block_title = "COI"
)

metazoa_18S_block <- make_marker_block(
  plot_list = metazoa_all_plots,
  marker_name = "18S",
  site_order = final_site_order,
  legend_plot = make_metazoa_marker_legend("18S"),
  block_title = "18S"
)

metazoa_Both <- metazoa_COI_block |
  metazoa_18S_block

print(metazoa_Both)

ggplot2::ggsave(
  filename = file.path(
    metazoa_output_dir,
    "Metazoa_richness_top5_all_sites_Both.png"
  ),
  plot = metazoa_Both,
  width = 44,
  height = 32,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    metazoa_output_dir,
    "Metazoa_richness_top5_all_sites_Both.svg"
  ),
  plot = metazoa_Both,
  width = 44,
  height = 32,
  bg = "white",
  limitsize = FALSE
)