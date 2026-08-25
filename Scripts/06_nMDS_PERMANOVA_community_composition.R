# ==========================================================
# NMDS AND PERMANOVA
# ALL RETAINED ASVs, EUKARYOTA AND METAZOA
# COI AND 18S
# ==========================================================

# Install if required:
# install.packages(
#   c(
#     "readxl", "dplyr", "tibble", "ggplot2",
#     "ggrepel", "vegan", "patchwork", "writexl"
#   )
# )

library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(vegan)
library(patchwork)
library(writexl)


# ==========================================================
# 1. FILE PATHS
# ==========================================================

base_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/Data"
)

base_meta <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/Metadatas"
)

file_18S <- file.path(
  base_dir,
  "All_Peninsula_18S_AbRel.xlsx"
)

file_COI <- file.path(
  base_dir,
  "All_Peninsula_COI_AbRel.xlsx"
)

metadata_18S_file <- file.path(
  base_meta,
  "metadata_samples.xlsx"
)

metadata_COI_file <- file.path(
  base_meta,
  "metadata_samples.xlsx"
)

# Use the general metadata file if no COI-specific file exists.
if (!file.exists(metadata_COI_file)) {
  metadata_COI_file <- metadata_18S_file
}

base_output_dir <- paste0(
  "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/",
  "ALL_DATA/Analisis/nMDS_new/nMDS_PERMANOVA_final"
)

dir.create(
  base_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==========================================================
# 2. OUTPUT FOLDERS
# ==========================================================

output_folders <- c(
  "All retained ASVs" = "nMDS_all",
  "Eukaryota" = "nMDS_Euk",
  "Metazoa" = "nMDS_Metazoa"
)

for (folder_name in output_folders) {
  
  dir.create(
    file.path(
      base_output_dir,
      folder_name
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==========================================================
# 3. ANALYSIS PARAMETERS
# ==========================================================

random_seed <- 123

nmds_trymax_2D <- 100
nmds_trymax_3D <- 150

n_permutations <- 9999

# Within-site comparison between pre- and post-industrial eras.
run_within_site_era_permanova <- TRUE

# At least three samples are required in each Era level.
minimum_samples_per_era <- 2

# Pairwise PERMANOVA entails 28 comparisons per analysis.
# Leave FALSE initially to avoid a very long computation.
# Change to TRUE after checking the global PERMANOVA.
run_pairwise_permanova <- FALSE

# Minimum number of samples in which an ASV must occur.
# 1 retains every ASV present in at least one sample.
min_prevalence <- 1

# Ensure that Bray–Curtis uses relative composition.
standardize_to_relative <- TRUE

# FALSE reproduces the previous abundance treatment.
# TRUE downweights highly dominant ASVs.
square_root_abundance <- FALSE

# Do not exclude samples without a documented reason.
excluded_samples <- c(
  "CMB_40",
  "CSD_80"
)

# Example:
# excluded_samples <- c("HIT_52", "HIO_52")


site_order <- c(
  "DEE", "STO", "CCO", "COL",
  "CSD", "CMB", "HIT", "HIO"
)

site_colors <- c(
  "DEE" = "#CDCD00",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "COL" = "burlywood1",
  "CSD" = "#76EEC6",
  "CMB" = "#458B74",
  "HIT" = "#CD96CD",
  "HIO" = "#8B0A50"
)


# ==========================================================
# 4. GENERAL FUNCTIONS
# ==========================================================

find_column <- function(
    data,
    candidates,
    required = TRUE) {
  
  data_names_lower <- tolower(
    names(data)
  )
  
  candidate_positions <- match(
    tolower(candidates),
    data_names_lower
  )
  
  candidate_positions <- candidate_positions[
    !is.na(candidate_positions)
  ]
  
  if (length(candidate_positions) == 0) {
    
    if (required) {
      
      stop(
        paste0(
          "None of these columns was found: ",
          paste(
            candidates,
            collapse = ", "
          )
        )
      )
    }
    
    return(NULL)
  }
  
  names(data)[candidate_positions[1]]
}


safe_numeric <- function(x) {
  
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x <- trimws(
    as.character(x)
  )
  
  x[
    x %in% c(
      "", "NA", "NaN",
      "Inf", "-Inf"
    )
  ] <- NA_character_
  
  suppressWarnings(
    as.numeric(x)
  )
}


add_analysis_columns <- function(
    data,
    marker,
    community) {
  
  data %>%
    mutate(
      Marker = marker,
      Community = community,
      .before = 1
    )
}


# ==========================================================
# 5. READ AND CLEAN METADATA
# ==========================================================

read_and_clean_metadata <- function(path) {
  
  metadata <- read_excel(path) %>%
    as.data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
  sample_column <- find_column(
    metadata,
    c(
      "sample", "Sample",
      "sample_id", "Sample_ID"
    )
  )
  
  site_column <- find_column(
    metadata,
    c(
      "Site", "site",
      "Core", "core", "Core_ID"
    )
  )
  
  year_column <- find_column(
    metadata,
    c(
      "Year", "year",
      "Age", "age"
    ),
    required = FALSE
  )
  
  depth_column <- find_column(
    metadata,
    c(
      "Depth", "depth",
      "Depth_cm", "depth_cm"
    ),
    required = FALSE
  )
  
  era_column <- find_column(
    metadata,
    c(
      "Era", "era", "ERA"
    ),
    required = TRUE
  )
  
  names(metadata)[
    names(metadata) == sample_column
  ] <- "sample"
  
  names(metadata)[
    names(metadata) == site_column
  ] <- "Site"
  
  if (!is.null(year_column)) {
    
    names(metadata)[
      names(metadata) == year_column
    ] <- "Year"
  }
  
  if (!is.null(depth_column)) {
    
    names(metadata)[
      names(metadata) == depth_column
    ] <- "Depth"
  }
  
  names(metadata)[
    names(metadata) == era_column
  ] <- "Era"
  
  if (!"Year" %in% names(metadata)) {
    metadata$Year <- NA_real_
  }
  
  if (!"Depth" %in% names(metadata)) {
    metadata$Depth <- NA_real_
  }
  
  metadata <- metadata %>%
    mutate(
      sample = as.character(sample),
      
      Site = as.character(Site),
      
      Site = recode(
        Site,
        "CLR" = "COL"
      ),
      
      Site = factor(
        Site,
        levels = site_order
      ),
      
      Year = safe_numeric(Year),
      Depth = safe_numeric(Depth),
      
      Era = trimws(
        as.character(Era)
      ),
      
      Era = na_if(
        Era,
        ""
      )
    ) %>%
    filter(
      !is.na(sample),
      sample != "",
      !is.na(Site)
    ) %>%
    select(
      sample,
      Site,
      Year,
      Depth,
      Era,
      everything()
    )
  
  duplicated_samples <- metadata$sample[
    duplicated(metadata$sample)
  ]
  
  if (length(duplicated_samples) > 0) {
    
    stop(
      paste0(
        "Duplicated samples in metadata: ",
        paste(
          unique(duplicated_samples),
          collapse = ", "
        )
      )
    )
  }
  
  metadata
}


read_taxonomic_table <- function(path) {
  
  read_excel(path) %>%
    as.data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
}


# ==========================================================
# 6. CREATE COMMUNITY MATRIX
# ==========================================================

create_community_matrix <- function(
    taxonomic_data,
    metadata,
    taxonomic_group = c(
      "All retained ASVs",
      "Eukaryota",
      "Metazoa"
    ),
    min_prevalence = 1,
    excluded_samples = character(0),
    standardize_to_relative = TRUE,
    square_root_abundance = FALSE) {
  
  taxonomic_group <- match.arg(
    taxonomic_group
  )
  
  asv_column <- find_column(
    taxonomic_data,
    c(
      "ASV", "asv", "ASV_ID",
      "Feature.ID", "FeatureID"
    )
  )
  
  if (taxonomic_group == "All retained ASVs") {
    
    keep_rows <- rep(
      TRUE,
      nrow(taxonomic_data)
    )
    
  } else if (taxonomic_group == "Eukaryota") {
    
    domain_column <- find_column(
      taxonomic_data,
      c(
        "domain", "Domain",
        "superkingdom"
      )
    )
    
    keep_rows <- tolower(
      trimws(
        as.character(
          taxonomic_data[[domain_column]]
        )
      )
    ) == "eukaryota"
    
  } else {
    
    kingdom_column <- find_column(
      taxonomic_data,
      c(
        "Kingdom", "kingdom"
      )
    )
    
    keep_rows <- tolower(
      trimws(
        as.character(
          taxonomic_data[[kingdom_column]]
        )
      )
    ) == "metazoa"
  }
  
  keep_rows[is.na(keep_rows)] <- FALSE
  
  filtered_taxa <- taxonomic_data[
    keep_rows,
    ,
    drop = FALSE
  ]
  
  if (nrow(filtered_taxa) == 0) {
    
    stop(
      paste(
        "No ASVs were found for",
        taxonomic_group
      )
    )
  }
  
  # Detect sample columns from metadata instead of using
  # a fixed number of taxonomic columns.
  sample_columns <- metadata$sample[
    metadata$sample %in%
      names(filtered_taxa)
  ]
  
  sample_columns <- setdiff(
    sample_columns,
    excluded_samples
  )
  
  if (length(sample_columns) == 0) {
    
    stop(
      paste(
        "No sample columns matched the metadata for",
        taxonomic_group
      )
    )
  }
  
  message(
    taxonomic_group,
    ": ",
    nrow(filtered_taxa),
    " ASV rows; ",
    length(sample_columns),
    " matched sample columns."
  )
  
  # --------------------------------------------------------
  # Efficient numerical conversion
  # --------------------------------------------------------
  # Conversion occurs only once, before transposition.
  
  abundance_by_asv <- do.call(
    cbind,
    lapply(
      filtered_taxa[
        sample_columns
      ],
      safe_numeric
    )
  )
  
  abundance_by_asv <- as.matrix(
    abundance_by_asv
  )
  
  storage.mode(abundance_by_asv) <- "double"
  
  colnames(abundance_by_asv) <- sample_columns
  
  asv_names <- as.character(
    filtered_taxa[[asv_column]]
  )
  
  missing_asv_names <- (
    is.na(asv_names) |
      asv_names == ""
  )
  
  asv_names[missing_asv_names] <- paste0(
    "Unnamed_ASV_",
    which(missing_asv_names)
  )
  
  rownames(abundance_by_asv) <- make.unique(
    asv_names
  )
  
  abundance_by_asv[
    !is.finite(abundance_by_asv)
  ] <- 0
  
  abundance_by_asv[
    abundance_by_asv < 0
  ] <- 0
  
  # Samples become rows and ASVs become columns.
  community_matrix <- t(
    abundance_by_asv
  )
  
  rm(abundance_by_asv)
  gc(verbose = FALSE)
  
  # Remove empty samples.
  community_matrix <- community_matrix[
    rowSums(community_matrix) > 0,
    ,
    drop = FALSE
  ]
  
  # Remove empty ASVs.
  community_matrix <- community_matrix[
    ,
    colSums(community_matrix) > 0,
    drop = FALSE
  ]
  
  # Remove ASVs below the selected prevalence.
  prevalence <- colSums(
    community_matrix > 0
  )
  
  community_matrix <- community_matrix[
    ,
    prevalence >= min_prevalence,
    drop = FALSE
  ]
  
  if (nrow(community_matrix) < 3) {
    
    stop(
      paste(
        "Fewer than three non-empty samples remained for",
        taxonomic_group
      )
    )
  }
  
  if (ncol(community_matrix) < 2) {
    
    stop(
      paste(
        "Fewer than two ASVs remained for",
        taxonomic_group
      )
    )
  }
  
  analysis_metadata <- metadata[
    match(
      rownames(community_matrix),
      metadata$sample
    ),
    ,
    drop = FALSE
  ]
  
  if (
    !identical(
      rownames(community_matrix),
      analysis_metadata$sample
    )
  ) {
    
    stop(
      "The community matrix and metadata are not aligned."
    )
  }
  
  analysis_metadata$Site <- droplevels(
    analysis_metadata$Site
  )
  
  if (standardize_to_relative) {
    
    community_matrix <- decostand(
      community_matrix,
      method = "total",
      MARGIN = 1
    )
  }
  
  if (square_root_abundance) {
    
    community_matrix <- sqrt(
      community_matrix
    )
  }
  
  message(
    taxonomic_group,
    ": final matrix = ",
    nrow(community_matrix),
    " samples × ",
    ncol(community_matrix),
    " ASVs."
  )
  
  list(
    community = community_matrix,
    metadata = analysis_metadata
  )
}


# ==========================================================
# 7. NMDS FUNCTIONS
# ==========================================================

run_nmds <- function(
    community_matrix,
    dimensions,
    trymax,
    seed) {
  
  set.seed(seed)
  
  metaMDS(
    community_matrix,
    distance = "bray",
    k = dimensions,
    trymax = trymax,
    autotransform = FALSE,
    # Use extended dissimilarities when samples share no ASVs.
    # This prevents degenerate, line-like ordinations in sparse datasets
    # such as the COI Metazoa subset.
    noshare = TRUE,
    zerodist = "add",
    trace = 1
  )
}


extract_nmds_scores <- function(
    nmds_model,
    metadata,
    marker,
    taxonomic_group,
    dimensions) {
  
  if (!"Era" %in% names(metadata)) {
    
    stop(
      paste0(
        "The column 'Era' was not found in the metadata used for ",
        marker,
        " — ",
        taxonomic_group,
        "."
      )
    )
  }
  
  score_data <- scores(
    nmds_model,
    display = "sites",
    choices = seq_len(dimensions)
  ) %>%
    as.data.frame() %>%
    rownames_to_column(
      var = "sample"
    )
  
  metadata_minimal <- metadata %>%
    transmute(
      sample = as.character(sample),
      Site,
      Year,
      Depth,
      Era
    ) %>%
    distinct(
      sample,
      .keep_all = TRUE
    )
  
  score_data %>%
    left_join(
      metadata_minimal,
      by = "sample"
    ) %>%
    mutate(
      Marker = marker,
      Community = taxonomic_group,
      Dimensions = dimensions,
      Site = factor(
        Site,
        levels = site_order
      ),
      Era = factor(Era)
    ) %>%
    arrange(
      Site,
      Year,
      sample
    )
}


# ==========================================================
# 8. GLOBAL PERMANOVA AND DISPERSION
# ==========================================================

run_global_statistics <- function(
    community_matrix,
    metadata,
    permutations,
    seed) {
  
  bray_distance <- vegdist(
    community_matrix,
    method = "bray"
  )
  
  set.seed(seed)
  
  permanova_model <- adonis2(
    bray_distance ~ Site,
    data = metadata,
    permutations = permutations
  )
  
  permanova_table <- as.data.frame(
    permanova_model
  ) %>%
    rownames_to_column(
      var = "Term"
    )
  
  dispersion_model <- betadisper(
    bray_distance,
    group = metadata$Site,
    type = "median",
    bias.adjust = TRUE,
    add = "lingoes"
  )
  
  set.seed(seed)
  
  dispersion_test <- permutest(
    dispersion_model,
    permutations = permutations
  )
  
  dispersion_table <- as.data.frame(
    dispersion_test$tab
  ) %>%
    rownames_to_column(
      var = "Term"
    )
  
  dispersion_by_site <- tibble(
    sample = metadata$sample,
    Site = metadata$Site,
    distance_to_centroid =
      dispersion_model$distances
  ) %>%
    group_by(Site) %>%
    summarise(
      n = n(),
      
      mean_distance_to_centroid = mean(
        distance_to_centroid,
        na.rm = TRUE
      ),
      
      median_distance_to_centroid = median(
        distance_to_centroid,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  list(
    permanova = permanova_table,
    dispersion = dispersion_table,
    dispersion_by_site = dispersion_by_site
  )
}


# ==========================================================
# 9. WITHIN-SITE PERMANOVA: ERA
# ==========================================================

run_within_site_era_statistics <- function(
    community_matrix,
    metadata,
    permutations = 9999,
    seed = 123,
    minimum_samples_per_era = 3) {
  
  if (!"Era" %in% names(metadata)) {
    stop("The column 'Era' was not found in the analysis metadata.")
  }
  
  if (!identical(rownames(community_matrix), metadata$sample)) {
    stop(
      paste0(
        "Community matrix and metadata are not aligned in ",
        "run_within_site_era_statistics()."
      )
    )
  }
  
  metadata <- metadata %>%
    mutate(
      Era = trimws(as.character(Era)),
      Era = na_if(Era, ""),
      Site = droplevels(factor(Site, levels = site_order))
    )
  
  sites_present <- levels(droplevels(metadata$Site))
  
  site_results <- lapply(
    seq_along(sites_present),
    function(site_index) {
      
      current_site <- sites_present[site_index]
      
      keep_samples <- (
        metadata$Site == current_site &
          !is.na(metadata$Era)
      )
      
      site_metadata <- metadata[
        keep_samples,
        ,
        drop = FALSE
      ] %>%
        mutate(Era = droplevels(factor(Era)))
      
      site_community <- community_matrix[
        keep_samples,
        ,
        drop = FALSE
      ]
      
      if (nrow(site_community) > 0) {
        site_community <- site_community[
          ,
          colSums(site_community) > 0,
          drop = FALSE
        ]
      }
      
      era_counts <- table(site_metadata$Era)
      era_levels <- names(era_counts)
      
      era_counts_text <- paste(
        paste0(era_levels, " = ", as.integer(era_counts)),
        collapse = "; "
      )
      
      empty_result <- tibble(
        Site = current_site,
        n_samples = nrow(site_metadata),
        n_eras = length(era_levels),
        Era_counts = era_counts_text,
        minimum_n_per_era = if (length(era_counts) > 0) {
          min(as.integer(era_counts))
        } else {
          NA_integer_
        },
        PERMANOVA_F = NA_real_,
        PERMANOVA_R2 = NA_real_,
        PERMANOVA_p = NA_real_,
        Dispersion_F = NA_real_,
        Dispersion_p = NA_real_,
        Test_status = NA_character_
      )
      
      if (length(era_levels) != 2) {
        empty_result$Test_status <- paste0(
          "Not tested: expected two Era levels but found ",
          length(era_levels),
          "."
        )
        return(empty_result)
      }
      
      if (any(era_counts < minimum_samples_per_era)) {
        empty_result$Test_status <- paste0(
          "Not tested: fewer than ",
          minimum_samples_per_era,
          " samples in at least one Era."
        )
        return(empty_result)
      }
      
      if (nrow(site_community) < 2 * minimum_samples_per_era) {
        empty_result$Test_status <- paste0(
          "Not tested: fewer than ",
          2 * minimum_samples_per_era,
          " samples in total."
        )
        return(empty_result)
      }
      
      if (ncol(site_community) < 2) {
        empty_result$Test_status <- paste0(
          "Not tested: fewer than two ASVs were retained within the site."
        )
        return(empty_result)
      }
      
      site_distance <- tryCatch(
        vegdist(site_community, method = "bray"),
        error = function(e) NULL
      )
      
      if (is.null(site_distance)) {
        empty_result$Test_status <- paste0(
          "Not tested: Bray-Curtis distances could not be calculated."
        )
        return(empty_result)
      }
      
      distance_values <- as.numeric(site_distance)
      
      if (
        any(!is.finite(distance_values)) ||
        length(unique(distance_values)) < 2
      ) {
        empty_result$Test_status <- paste0(
          "Not tested: Bray-Curtis distances were non-finite or lacked variation."
        )
        return(empty_result)
      }
      
      set.seed(seed + site_index)
      
      permanova_model <- tryCatch(
        adonis2(
          site_distance ~ Era,
          data = site_metadata,
          permutations = permutations
        ),
        error = function(e) NULL
      )
      
      if (is.null(permanova_model)) {
        empty_result$Test_status <- "PERMANOVA failed for this site."
        return(empty_result)
      }
      
      dispersion_output <- tryCatch(
        {
          dispersion_model <- betadisper(
            site_distance,
            group = site_metadata$Era,
            type = "median",
            bias.adjust = TRUE,
            add = "lingoes"
          )
          
          set.seed(seed + 1000 + site_index)
          
          dispersion_test <- permutest(
            dispersion_model,
            permutations = permutations
          )
          
          list(
            F = dispersion_test$tab[1, "F"],
            p = dispersion_test$tab[1, "Pr(>F)"]
          )
        },
        error = function(e) {
          list(F = NA_real_, p = NA_real_)
        }
      )
      
      tibble(
        Site = current_site,
        n_samples = nrow(site_metadata),
        n_eras = length(era_levels),
        Era_counts = era_counts_text,
        minimum_n_per_era = min(as.integer(era_counts)),
        PERMANOVA_F = as.numeric(permanova_model["Era", "F"]),
        PERMANOVA_R2 = as.numeric(permanova_model["Era", "R2"]),
        PERMANOVA_p = as.numeric(permanova_model["Era", "Pr(>F)"]),
        Dispersion_F = as.numeric(dispersion_output$F),
        Dispersion_p = as.numeric(dispersion_output$p),
        Test_status = "Test completed"
      )
    }
  )
  
  bind_rows(site_results) %>%
    mutate(
      # These adjusted values are retained as diagnostic information.
      # Nominal p-values are also exported unchanged.
      PERMANOVA_p_BH = p.adjust(PERMANOVA_p, method = "BH"),
      Dispersion_p_BH = p.adjust(Dispersion_p, method = "BH"),
      PERMANOVA_interpretation = case_when(
        is.na(PERMANOVA_p) ~ "Not tested",
        PERMANOVA_p < 0.05 &
          !is.na(Dispersion_p) &
          Dispersion_p < 0.05 ~
          "Era difference detected, but dispersion also differed",
        PERMANOVA_p < 0.05 ~
          "Era difference detected; no significant dispersion difference",
        TRUE ~ "No significant Era difference detected"
      )
    ) %>%
    arrange(factor(Site, levels = site_order))
}


# ==========================================================
# 10. OPTIONAL PAIRWISE PERMANOVA
# ==========================================================

run_pairwise_statistics <- function(
    community_matrix,
    metadata,
    permutations,
    seed) {
  
  site_levels <- levels(
    droplevels(
      metadata$Site
    )
  )
  
  site_pairs <- combn(
    site_levels,
    2,
    simplify = FALSE
  )
  
  pairwise_results <- lapply(
    seq_along(site_pairs),
    function(i) {
      
      pair <- site_pairs[[i]]
      
      keep <- metadata$Site %in% pair
      
      pair_metadata <- droplevels(
        metadata[
          keep,
          ,
          drop = FALSE
        ]
      )
      
      pair_community <- community_matrix[
        keep,
        ,
        drop = FALSE
      ]
      
      pair_community <- pair_community[
        ,
        colSums(pair_community) > 0,
        drop = FALSE
      ]
      
      if (
        nrow(pair_community) < 4 ||
        ncol(pair_community) < 2
      ) {
        
        return(
          tibble(
            Site_1 = pair[1],
            Site_2 = pair[2],
            n_1 = sum(
              pair_metadata$Site == pair[1]
            ),
            n_2 = sum(
              pair_metadata$Site == pair[2]
            ),
            R2 = NA_real_,
            F = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      pair_distance <- vegdist(
        pair_community,
        method = "bray"
      )
      
      set.seed(seed + i)
      
      pair_model <- adonis2(
        pair_distance ~ Site,
        data = pair_metadata,
        permutations = permutations
      )
      
      tibble(
        Site_1 = pair[1],
        Site_2 = pair[2],
        
        n_1 = sum(
          pair_metadata$Site == pair[1]
        ),
        
        n_2 = sum(
          pair_metadata$Site == pair[2]
        ),
        
        R2 = pair_model[
          "Site",
          "R2"
        ],
        
        F = pair_model[
          "Site",
          "F"
        ],
        
        p_value = pair_model[
          "Site",
          "Pr(>F)"
        ]
      )
    }
  )
  
  bind_rows(
    pairwise_results
  ) %>%
    mutate(
      p_adjusted_BH = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    arrange(
      p_adjusted_BH,
      Site_1,
      Site_2
    )
}


# ==========================================================
# 10. ANALYSE ONE MARKER AND COMMUNITY
# ==========================================================

analyse_community <- function(
    taxonomic_data,
    metadata,
    marker,
    taxonomic_group) {
  
  message(
    "\n=========================================="
  )
  
  message(
    "Analysing ",
    marker,
    " — ",
    taxonomic_group
  )
  
  message(
    "=========================================="
  )
  
  prepared <- create_community_matrix(
    taxonomic_data = taxonomic_data,
    metadata = metadata,
    taxonomic_group = taxonomic_group,
    min_prevalence = min_prevalence,
    excluded_samples = excluded_samples,
    standardize_to_relative =
      standardize_to_relative,
    square_root_abundance =
      square_root_abundance
  )
  
  # --------------------------------------------------------
  # 2D nMDS
  # --------------------------------------------------------
  
  nmds_2D <- run_nmds(
    community_matrix = prepared$community,
    dimensions = 2,
    trymax = nmds_trymax_2D,
    seed = random_seed
  )
  
  scores_2D <- extract_nmds_scores(
    nmds_model = nmds_2D,
    metadata = prepared$metadata,
    marker = marker,
    taxonomic_group = taxonomic_group,
    dimensions = 2
  )
  
  # --------------------------------------------------------
  # 3D nMDS
  # --------------------------------------------------------
  
  nmds_3D <- run_nmds(
    community_matrix = prepared$community,
    dimensions = 3,
    trymax = nmds_trymax_3D,
    seed = random_seed + 1000
  )
  
  scores_3D <- extract_nmds_scores(
    nmds_model = nmds_3D,
    metadata = prepared$metadata,
    marker = marker,
    taxonomic_group = taxonomic_group,
    dimensions = 3
  )
  
  # --------------------------------------------------------
  # Statistics
  # --------------------------------------------------------
  
  global_statistics <- run_global_statistics(
    community_matrix = prepared$community,
    metadata = prepared$metadata,
    permutations = n_permutations,
    seed = random_seed
  )
  
  # --------------------------------------------------------
  # Within-site Era PERMANOVA and dispersion
  # --------------------------------------------------------
  
  if (run_within_site_era_permanova) {
    
    era_statistics <- run_within_site_era_statistics(
      community_matrix = prepared$community,
      metadata = prepared$metadata,
      permutations = n_permutations,
      seed = random_seed,
      minimum_samples_per_era = minimum_samples_per_era
    )
    
  } else {
    
    era_statistics <- tibble(
      Note = paste0(
        "Within-site Era PERMANOVA was not run. ",
        "Set run_within_site_era_permanova <- TRUE."
      )
    )
  }
  
  if (run_pairwise_permanova) {
    
    pairwise_statistics <- run_pairwise_statistics(
      community_matrix = prepared$community,
      metadata = prepared$metadata,
      permutations = n_permutations,
      seed = random_seed
    )
    
  } else {
    
    pairwise_statistics <- tibble(
      Note = paste0(
        "Pairwise PERMANOVA was not run. ",
        "Set run_pairwise_permanova <- TRUE ",
        "to calculate it."
      )
    )
  }
  
  sample_summary <- prepared$metadata %>%
    count(
      Site,
      name = "n_samples"
    )
  
  analysis_summary <- tibble(
    Marker = marker,
    Community = taxonomic_group,
    
    n_samples = nrow(
      prepared$community
    ),
    
    n_ASVs = ncol(
      prepared$community
    ),
    
    stress_2D = nmds_2D$stress,
    stress_3D = nmds_3D$stress,
    
    absolute_stress_reduction =
      nmds_2D$stress -
      nmds_3D$stress,
    
    relative_stress_reduction_percent =
      100 *
      (
        nmds_2D$stress -
          nmds_3D$stress
      ) /
      nmds_2D$stress
  )
  
  result <- list(
    marker = marker,
    community_name = taxonomic_group,
    
    scores_2D = scores_2D,
    scores_3D = scores_3D,
    
    stress_2D = nmds_2D$stress,
    stress_3D = nmds_3D$stress,
    
    analysis_summary = analysis_summary,
    sample_summary = sample_summary,
    
    permanova =
      global_statistics$permanova,
    
    dispersion =
      global_statistics$dispersion,
    
    dispersion_by_site =
      global_statistics$dispersion_by_site,
    
    era_permanova_by_site =
      era_statistics,
    
    pairwise_permanova =
      pairwise_statistics
  )
  
  # The large community matrix is not returned.
  rm(
    prepared,
    nmds_2D,
    nmds_3D
  )
  
  gc(verbose = FALSE)
  
  result
}


# ==========================================================
# 11. NMDS PLOT FUNCTION
# ==========================================================
create_nmds_plot <- function(
    analysis_results,
    model_dimensions = c("2D", "3D"),
    x_axis = "NMDS1",
    y_axis = "NMDS2",
    show_ellipses = TRUE,
    show_endpoint_labels = TRUE
) {
  
  model_dimensions <- match.arg(model_dimensions)
  
  # Select the scores and stress from the requested model
  if (model_dimensions == "2D") {
    
    plot_scores <- analysis_results$scores_2D
    model_stress <- analysis_results$stress_2D
    
  } else {
    
    plot_scores <- analysis_results$scores_3D
    model_stress <- analysis_results$stress_3D
    
  }
  
  # Check that the requested nMDS axes exist
  if (!all(c(x_axis, y_axis) %in% names(plot_scores))) {
    
    stop(
      paste0(
        "Requested axes were not found: ",
        x_axis,
        " and ",
        y_axis
      )
    )
  }
  
  # Check that Era is available in the ordination scores
  if (!"Era" %in% names(plot_scores)) {
    
    stop(
      paste0(
        "The variable 'Era' is not present in the nMDS scores. ",
        "Join the metadata column Era to scores_2D and scores_3D ",
        "before calling create_nmds_plot()."
      )
    )
  }
  
  # Prepare ordination scores
  plot_scores <- plot_scores %>%
    mutate(
      Era = factor(Era)
    ) %>%
    filter(
      is.finite(.data[[x_axis]]),
      is.finite(.data[[y_axis]])
    ) %>%
    arrange(
      Site,
      Year,
      sample
    )
  
  # Data used to calculate site ellipses
  ellipse_data <- plot_scores %>%
    group_by(Site) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  # Select four samples per site:
  # oldest, two approximately equally spaced intermediate samples,
  # and most recent.
  endpoint_labels <- plot_scores %>%
    filter(
      is.finite(Year)
    ) %>%
    arrange(
      Site,
      Year,
      sample
    ) %>%
    group_by(Site) %>%
    mutate(
      label_position = row_number(),
      n_site_samples = n()
    ) %>%
    filter(
      label_position == 1L |
        label_position == pmax(
          1L,
          round(1 + (n_site_samples - 1) / 3)
        ) |
        label_position == pmax(
          1L,
          round(1 + 2 * (n_site_samples - 1) / 3)
        ) |
        label_position == n_site_samples
    ) %>%
    ungroup() %>%
    mutate(
      Year_label = sprintf(
        "%.0f",
        Year
      )
    ) %>%
    select(
      -label_position,
      -n_site_samples
    )
  
  # Base plot
  p <- ggplot(
    plot_scores,
    aes(
      x = .data[[x_axis]],
      y = .data[[y_axis]],
      color = Site,
      shape = Era
    )
  )
  
  # Site ellipses
  if (
    show_ellipses &&
    nrow(ellipse_data) > 0
  ) {
    
    p <- p +
      stat_ellipse(
        data = ellipse_data,
        aes(
          group = Site
        ),
        type = "t",
        level = 0.95,
        linewidth = 1.25,
        alpha = 0.75,
        show.legend = FALSE,
        na.rm = TRUE
      )
  }
  
  # Sample points
  p <- p +
    geom_point(
      size = 3.2,
      stroke = 0.8,
      alpha = 0.90,
      na.rm = TRUE
    )
  
  # Year labels for four selected samples per site
  if (
    show_endpoint_labels &&
    nrow(endpoint_labels) > 0
  ) {
    
    p <- p +
      ggrepel::geom_text_repel(
        data = endpoint_labels,
        aes(
          label = Year_label,
          color = Site
        ),
        size = 4,
        family = "Arial",
        fontface = "bold",
        show.legend = FALSE,
        box.padding = 0.35,
        point.padding = 0.20,
        max.overlaps = Inf,
        
        # Suppress segments connecting labels and points
        segment.color = NA,
        
        seed = random_seed,
        na.rm = TRUE
      )
  }
  
  # Final plot
  p +
    scale_color_manual(
      values = site_colors,
      breaks = site_order,
      drop = FALSE
    ) +
    
    # Era is represented by point shape
    scale_shape_discrete(
      name = "Era",
      drop = FALSE
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    labs(
      title = analysis_results$marker,
      
      subtitle = paste0(
        model_dimensions,
        " Bray–Curtis nMDS; stress = ",
        format(
          round(model_stress, 3),
          nsmall = 3
        )
      ),
      
      x = x_axis,
      y = y_axis,
      color = "Site",
      shape = "Era"
    ) +
    
    theme_bw(
      base_size = 12,
      base_family = "Arial"
    ) +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 12,
        face = "bold"
      ),
      
      axis.text = element_text(
        size = 10,
        color = "black"
      ),
      
      panel.grid = element_blank(),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 10
      ),
      
      legend.box = "vertical",
      
      plot.margin = margin(
        15, 20, 15, 20
      )
    )
}


# ==========================================================
# 12. CREATE A COMBINED COI + 18S FIGURE
# ==========================================================

create_combined_marker_figure <- function(
    results_COI,
    results_18S,
    model_dimensions,
    x_axis,
    y_axis,
    figure_title) {
  
  p_COI <- create_nmds_plot(
    analysis_results = results_COI,
    model_dimensions = model_dimensions,
    x_axis = x_axis,
    y_axis = y_axis
  )
  
  p_18S <- create_nmds_plot(
    analysis_results = results_18S,
    model_dimensions = model_dimensions,
    x_axis = x_axis,
    y_axis = y_axis
  )
  
  (
    p_COI |
      p_18S
  ) +
    plot_layout(
      guides = "collect"
    ) +
    plot_annotation(
      title = figure_title
    ) &
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


# ==========================================================
# 13. SAVE RESULTS FOR ONE COMMUNITY GROUP
# ==========================================================

save_group_results <- function(
    results_COI,
    results_18S,
    folder_name,
    file_prefix,
    figure_title) {
  
  group_output_dir <- file.path(
    base_output_dir,
    folder_name
  )
  
  dir.create(
    group_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # --------------------------------------------------------
  # Create combined COI + 18S figures
  # --------------------------------------------------------
  
  figure_2D_12 <- create_combined_marker_figure(
    results_COI = results_COI,
    results_18S = results_18S,
    model_dimensions = "2D",
    x_axis = "NMDS1",
    y_axis = "NMDS2",
    figure_title = paste0(
      figure_title,
      ": two-dimensional nMDS"
    )
  )
  
  figure_3D_12 <- create_combined_marker_figure(
    results_COI = results_COI,
    results_18S = results_18S,
    model_dimensions = "3D",
    x_axis = "NMDS1",
    y_axis = "NMDS2",
    figure_title = paste0(
      figure_title,
      ": three-dimensional nMDS projected onto axes 1 and 2"
    )
  )
  
  figure_3D_13 <- create_combined_marker_figure(
    results_COI = results_COI,
    results_18S = results_18S,
    model_dimensions = "3D",
    x_axis = "NMDS1",
    y_axis = "NMDS3",
    figure_title = paste0(
      figure_title,
      ": three-dimensional nMDS projected onto axes 1 and 3"
    )
  )
  
  figure_3D_23 <- create_combined_marker_figure(
    results_COI = results_COI,
    results_18S = results_18S,
    model_dimensions = "3D",
    x_axis = "NMDS2",
    y_axis = "NMDS3",
    figure_title = paste0(
      figure_title,
      ": three-dimensional nMDS projected onto axes 2 and 3"
    )
  )
  
  figures <- list(
    NMDS_2D_axes_1_2 = figure_2D_12,
    NMDS_3D_axes_1_2 = figure_3D_12,
    NMDS_3D_axes_1_3 = figure_3D_13,
    NMDS_3D_axes_2_3 = figure_3D_23
  )
  
  for (
    figure_name in names(figures)
  ) {
    
    current_figure <- figures[[figure_name]]
    
    ggsave(
      filename = file.path(
        group_output_dir,
        paste0(
          file_prefix,
          "_",
          figure_name,
          ".svg"
        )
      ),
      plot = current_figure,
      width = 15,
      height = 8,
      bg = "white"
    )
    
    ggsave(
      filename = file.path(
        group_output_dir,
        paste0(
          file_prefix,
          "_",
          figure_name,
          ".png"
        )
      ),
      plot = current_figure,
      width = 15,
      height = 8,
      dpi = 400,
      bg = "white"
    )
  }
  
  # --------------------------------------------------------
  # Combine statistical outputs
  # --------------------------------------------------------
  
  analysis_summary <- bind_rows(
    results_COI$analysis_summary,
    results_18S$analysis_summary
  )
  
  sample_summary <- bind_rows(
    add_analysis_columns(
      results_COI$sample_summary,
      "COI",
      results_COI$community_name
    ),
    add_analysis_columns(
      results_18S$sample_summary,
      "18S",
      results_18S$community_name
    )
  )
  
  coordinates_2D <- bind_rows(
    results_COI$scores_2D,
    results_18S$scores_2D
  )
  
  coordinates_3D <- bind_rows(
    results_COI$scores_3D,
    results_18S$scores_3D
  )
  
  permanova_results <- bind_rows(
    add_analysis_columns(
      results_COI$permanova,
      "COI",
      results_COI$community_name
    ),
    add_analysis_columns(
      results_18S$permanova,
      "18S",
      results_18S$community_name
    )
  )
  
  dispersion_results <- bind_rows(
    add_analysis_columns(
      results_COI$dispersion,
      "COI",
      results_COI$community_name
    ),
    add_analysis_columns(
      results_18S$dispersion,
      "18S",
      results_18S$community_name
    )
  )
  
  dispersion_by_site <- bind_rows(
    add_analysis_columns(
      results_COI$dispersion_by_site,
      "COI",
      results_COI$community_name
    ),
    add_analysis_columns(
      results_18S$dispersion_by_site,
      "18S",
      results_18S$community_name
    )
  )
  
  # --------------------------------------------------------
  # Combine within-site Era results
  # --------------------------------------------------------
  
  if (run_within_site_era_permanova) {
    
    era_results <- bind_rows(
      add_analysis_columns(
        results_COI$era_permanova_by_site,
        "COI",
        results_COI$community_name
      ),
      add_analysis_columns(
        results_18S$era_permanova_by_site,
        "18S",
        results_18S$community_name
      )
    )
    
  } else {
    
    era_results <- tibble(
      Note = paste0(
        "Within-site Era PERMANOVA was not run. ",
        "Set run_within_site_era_permanova <- TRUE."
      )
    )
  }
  
  if (run_pairwise_permanova) {
    
    pairwise_results <- bind_rows(
      add_analysis_columns(
        results_COI$pairwise_permanova,
        "COI",
        results_COI$community_name
      ),
      add_analysis_columns(
        results_18S$pairwise_permanova,
        "18S",
        results_18S$community_name
      )
    )
    
  } else {
    
    pairwise_results <- tibble(
      Note = paste0(
        "Pairwise PERMANOVA was not run. ",
        "Set run_pairwise_permanova <- TRUE."
      )
    )
  }
  
  write_xlsx(
    list(
      NMDS_summary = analysis_summary,
      samples_by_site = sample_summary,
      coordinates_2D = coordinates_2D,
      coordinates_3D = coordinates_3D,
      PERMANOVA_global = permanova_results,
      dispersion_global = dispersion_results,
      dispersion_by_site = dispersion_by_site,
      PERMANOVA_Era_within_site = era_results,
      PERMANOVA_pairwise = pairwise_results
    ),
    path = file.path(
      group_output_dir,
      paste0(
        file_prefix,
        "_nMDS_PERMANOVA_results.xlsx"
      )
    )
  )
  
  cat(
    "\n\n==========================================\n"
  )
  
  cat(
    figure_title,
    "\n"
  )
  
  cat(
    "==========================================\n"
  )
  
  print(
    analysis_summary
  )
  
  cat(
    "\nGlobal PERMANOVA:\n"
  )
  
  print(
    permanova_results
  )
  
  cat(
    "\nDispersion test:\n"
  )
  
  print(
    dispersion_results
  )
  
  cat(
    "\nWithin-site Era PERMANOVA and dispersion:\n"
  )
  
  print(
    era_results
  )
  
  invisible(
    list(
      figures = figures,
      analysis_summary = analysis_summary,
      permanova = permanova_results,
      dispersion = dispersion_results,
      era_permanova_by_site = era_results
    )
  )
}


# ==========================================================
# 14. READ INPUT FILES
# ==========================================================

message(
  "Reading taxonomic tables..."
)

data_18S <- read_taxonomic_table(
  file_18S
)

data_COI <- read_taxonomic_table(
  file_COI
)

message(
  "Reading metadata..."
)

metadata_18S <- read_and_clean_metadata(
  metadata_18S_file
)

metadata_COI <- read_and_clean_metadata(
  metadata_COI_file
)


# ==========================================================
# 15. ALL RETAINED ASVs
# ==========================================================

results_COI_all <- analyse_community(
  taxonomic_data = data_COI,
  metadata = metadata_COI,
  marker = "COI",
  taxonomic_group = "All retained ASVs"
)

results_18S_all <- analyse_community(
  taxonomic_data = data_18S,
  metadata = metadata_18S,
  marker = "18S",
  taxonomic_group = "All retained ASVs"
)

saved_all <- save_group_results(
  results_COI = results_COI_all,
  results_18S = results_18S_all,
  folder_name = "nMDS_all",
  file_prefix = "All_ASVs_COI_18S",
  figure_title = "All retained ASVs"
)

rm(
  results_COI_all,
  results_18S_all
)

gc(verbose = FALSE)


# ==========================================================
# 16. EUKARYOTA
# ==========================================================

results_COI_euk <- analyse_community(
  taxonomic_data = data_COI,
  metadata = metadata_COI,
  marker = "COI",
  taxonomic_group = "Eukaryota"
)

results_18S_euk <- analyse_community(
  taxonomic_data = data_18S,
  metadata = metadata_18S,
  marker = "18S",
  taxonomic_group = "Eukaryota"
)

saved_euk <- save_group_results(
  results_COI = results_COI_euk,
  results_18S = results_18S_euk,
  folder_name = "nMDS_Euk",
  file_prefix = "Eukaryota_COI_18S",
  figure_title = "Eukaryotic community composition"
)

rm(
  results_COI_euk,
  results_18S_euk
)

gc(verbose = FALSE)


# ==========================================================
# 17. METAZOA
# ==========================================================

results_COI_metazoa <- analyse_community(
  taxonomic_data = data_COI,
  metadata = metadata_COI,
  marker = "COI",
  taxonomic_group = "Metazoa"
)

results_18S_metazoa <- analyse_community(
  taxonomic_data = data_18S,
  metadata = metadata_18S,
  marker = "18S",
  taxonomic_group = "Metazoa"
)

saved_metazoa <- save_group_results(
  results_COI = results_COI_metazoa,
  results_18S = results_18S_metazoa,
  folder_name = "nMDS_Metazoa",
  file_prefix = "Metazoa_COI_18S",
  figure_title = "Metazoan community composition"
)

rm(
  results_COI_metazoa,
  results_18S_metazoa
)

gc(verbose = FALSE)


# ==========================================================
# 18. COMBINED STRESS SUMMARY
# ==========================================================

combined_stress_summary <- bind_rows(
  saved_all$analysis_summary,
  saved_euk$analysis_summary,
  saved_metazoa$analysis_summary
) %>%
  arrange(
    Community,
    Marker
  )

write_xlsx(
  list(
    NMDS_stress_comparison =
      combined_stress_summary
  ),
  path = file.path(
    base_output_dir,
    "NMDS_2D_3D_stress_comparison_all_analyses.xlsx"
  )
)

# ==========================================================
# 19. COMBINED WITHIN-SITE ERA SUMMARY
# ==========================================================

if (run_within_site_era_permanova) {
  
  combined_era_summary <- bind_rows(
    saved_all$era_permanova_by_site,
    saved_euk$era_permanova_by_site,
    saved_metazoa$era_permanova_by_site
  ) %>%
    mutate(
      Community = factor(
        Community,
        levels = c(
          "All retained ASVs",
          "Eukaryota",
          "Metazoa"
        )
      ),
      Marker = factor(
        Marker,
        levels = c("COI", "18S")
      ),
      Site = factor(
        Site,
        levels = site_order
      )
    ) %>%
    arrange(
      Community,
      Marker,
      Site
    )
  
  write_xlsx(
    list(
      Era_PERMANOVA_all = combined_era_summary
    ),
    path = file.path(
      base_output_dir,
      paste0(
        "PERMANOVA_Era_within_site_",
        "all_markers_and_communities.xlsx"
      )
    )
  )
  
  cat(
    "\n\n==========================================\n"
  )
  
  cat(
    "WITHIN-SITE ERA PERMANOVA SUMMARY\n"
  )
  
  cat(
    "==========================================\n"
  )
  
  print(
    combined_era_summary
  )
}

cat(
  "\n\n==========================================\n"
)

cat(
  "COMPLETE NMDS STRESS COMPARISON\n"
)

cat(
  "==========================================\n"
)

print(
  combined_stress_summary
)

cat(
  "\nAnalysis completed successfully.\n"
)

cat(
  "Results saved in:\n",
  base_output_dir,
  "\n"
)