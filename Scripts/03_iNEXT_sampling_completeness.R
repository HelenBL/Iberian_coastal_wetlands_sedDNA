library(ggplot2)
library(openxlsx)
library(readr)
library(dplyr)
library(stringr)
library(svglite)

# ==================== CONFIG ====================

data_dir <- "C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/Analisis/Data"

outdir <- file.path(data_dir, "read_depth_diagnostics")
dir.create(outdir, showWarnings = FALSE)

coi_file <- file.path(data_dir, "all_data_COI_AbTot.csv")
s18_file <- file.path(data_dir, "All_Peninsula_18S_Reps_AbTot.xlsx")

thresholds <- c(10000, 4000000)

# ==================== FUNCIONES ====================

read_count_table <- function(file, marker) {
  
  if (grepl("\\.csv$", file)) {
    tab <- read.csv(file, check.names = FALSE)
  } else {
    tab <- read.xlsx(file)
  }
  
  rownames(tab) <- tab[[1]]
  tab[[1]] <- NULL
  
  tab[] <- lapply(tab, as.numeric)
  
  reads <- colSums(tab, na.rm = TRUE)
  
  df <- data.frame(
    PCR = names(reads),
    Reads = as.numeric(reads),
    Marker = marker,
    stringsAsFactors = FALSE
  )
  
  # Sample without PCR replicate suffix
  # Example: CCO.C1.45.1 -> CCO.C1.45
  df$Sample <- sub("(.*)\\.[0-9]+$", "\\1", df$PCR)
  
  df <- df %>%
    arrange(Reads) %>%
    mutate(Rank = row_number())
  
  return(df)
}

make_ordered_plot <- function(df, marker) {
  
  p <- ggplot(df, aes(x = Rank, y = Reads)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_hline(yintercept = thresholds, linetype = "dashed") +
    scale_y_log10() +
    labs(
      x = "PCR replicates ordered by read depth",
      y = "Total reads per PCR replicate (log10)",
      title = paste0(marker, ": sequencing depth across PCR replicates")
    ) +
    theme_bw()
  
  ggsave(
    file.path(outdir, paste0(marker, "_reads_ordered.png")),
    plot = p,
    width = 9, height = 6, dpi = 300
  )
  
  ggsave(
    file.path(outdir, paste0(marker, "_reads_ordered.svg")),
    plot = p,
    width = 9, height = 6,
    device = svglite::svglite
  )
  
  return(p)
}

make_sample_plot <- function(df, marker) {
  
  site_cols <- c(
    "CMB" = "#458B74",
    "CSD" = "#76EEC6",
    "DEE" = "#CDCD00",
    "HIT" = "#CD96CD",
    "HIO" = "#FFE1FF",
    "STO" = "#8B7355",
    "CCO" = "#CDAA7D",
    "CLR" = "burlywood1"
  )
  
  # ==================== LIMPIAR NOMBRES SIN RELEER TABLA ====================
  
  df <- df %>%
    mutate(
      PCR = as.character(PCR),
      
      # Quita punto final si existe: CCO_C3_3.2. -> CCO_C3_3.2
      PCR_clean = sub("\\.$", "", PCR),
      
      # Quita la réplica final: CCO_C3_3.2 -> CCO_C3_3
      Sample_clean = sub("\\.[0-9]+$", "", PCR_clean),
      
      # Reemplaza Sample por el nombre limpio
      Sample = Sample_clean,
      
      # Sitio = tres primeras letras
      Site = substr(Sample, 1, 3)
    ) %>%
    filter(Site %in% names(site_cols))
  
  # ==================== ORDENAR MUESTRAS ====================
  
  sample_order <- df %>%
    group_by(Sample) %>%
    summarise(
      median_reads = median(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_reads) %>%
    pull(Sample)
  
  df$Sample <- factor(df$Sample, levels = sample_order)
  
  # ==================== PLOT ====================
  
  p <- ggplot(df, aes(x = Sample, y = Reads)) +
    
    geom_jitter(
      aes(color = Site),
      width = 0.2,
      size = 1.5,
      alpha = 0.25
    ) +
    
    geom_boxplot(
      aes(fill = Site, group = Sample),
      outlier.shape = NA,
      alpha = 0.65,
      width = 0.6,
      color = "black"
    ) +
    
    geom_hline(yintercept = thresholds, linetype = "dashed") +
    scale_y_log10() +
    scale_color_manual(values = site_cols, drop = FALSE) +
    scale_fill_manual(values = site_cols, drop = FALSE) +
    labs(
      x = "Sediment sample",
      y = "Total reads per PCR replicate (log10)",
      color = "Site",
      fill = "Site",
      title = paste0(marker, ": read depth variability among PCR replicates within samples")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 7
      )
    )
  
  ggsave(
    file.path(outdir, paste0(marker, "_reads_by_sample.png")),
    plot = p,
    width = 13, height = 6, dpi = 300
  )
  
  ggsave(
    file.path(outdir, paste0(marker, "_reads_by_sample.svg")),
    plot = p,
    width = 13, height = 6,
    device = svglite::svglite
  )
  
  return(p)
}

make_summary <- function(df, marker) {
  
  summary_df <- df %>%
    summarise(
      Marker = marker,
      n_PCR = n(),
      min_reads = min(Reads, na.rm = TRUE),
      q1_reads = quantile(Reads, 0.25, na.rm = TRUE),
      median_reads = median(Reads, na.rm = TRUE),
      mean_reads = mean(Reads, na.rm = TRUE),
      q3_reads = quantile(Reads, 0.75, na.rm = TRUE),
      max_reads = max(Reads, na.rm = TRUE),
      n_below_1000 = sum(Reads < 1000, na.rm = TRUE),
      n_below_30000 = sum(Reads < 30000, na.rm = TRUE)
    )
  
  write.csv(
    summary_df,
    file.path(outdir, paste0(marker, "_read_depth_summary.csv")),
    row.names = FALSE
  )
  
  low_reads <- df %>%
    filter(Reads < max(thresholds)) %>%
    arrange(Reads)
  
  write.csv(
    low_reads,
    file.path(outdir, paste0(marker, "_PCRs_below_30000_reads.csv")),
    row.names = FALSE
  )
  
  return(summary_df)
}

# ==================== RUN ====================

coi_df <- read_count_table(coi_file, "COI")
coi_df <- coi_df[c(10:847),]
s18_df <- read_count_table(s18_file, "18S")
s18_df <- s18_df[c(10:840),]

coi_df <- coi_df %>%
  mutate(
    PCR_clean = sub("\\.$", "", as.character(PCR)),
    Sample = sub("\\.[0-9]+$", "", PCR_clean),
    Site = substr(Sample, 1, 3)
  )

s18_df <- s18_df %>%
  mutate(
    PCR_clean = sub("\\.$", "", as.character(PCR)),
    Sample = sub("\\.[0-9]+$", "", PCR_clean),
    Site = substr(Sample, 1, 3)
  )

write.csv(coi_df, file.path(outdir, "COI_read_depth_per_PCR.csv"), row.names = FALSE)
write.csv(s18_df, file.path(outdir, "18S_read_depth_per_PCR.csv"), row.names = FALSE)

make_ordered_plot(coi_df, "COI")
ggsave(
  filename = "coi_read_depth.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
make_sample_plot(coi_df, "COI")
ggsave(
  filename = "coi_read_depth_samples.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
coi_summary <- make_summary(coi_df, "COI")


make_ordered_plot(s18_df, "18S")
ggsave(
  filename = "18S_read_depth.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
make_sample_plot(s18_df, "18S")
ggsave(
  filename = "18s_read_depth_samples.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
s18_summary <- make_summary(s18_df, "18S")

summary_all <- bind_rows(coi_summary, s18_summary)

write.csv(
  summary_all,
  file.path(outdir, "read_depth_summary_COI_18S.csv"),
  row.names = FALSE
)

cat("Read-depth diagnostic plots completed.\n")
cat("Results saved in:", outdir, "\n")










#======== LUKE'S SCRIPT RAREFACTION =================
file_18 <- read.xlsx(s18_file)
row.names(file_18) <- file_18$ASV
file_18<-file_18[,c(13:841)]
file_18[is.na(file_18)] <- 0
file_18[] <- lapply(file_18, as.numeric)
### Alpha diversity 
#lets look at uneven read depth 
lib_sizes <- colSums(file_18) 
min(lib_sizes)                        
quantile(lib_sizes) 
#ignore samples below 60k
n_reads <- 20000
euk_filt_18 <- as.matrix(file_18[ , lib_sizes >= n_reads]) 


#some functions to help rarefy & calc
## Hill numbers for one count vector
hill_3 <- function(v){
  v <- v[v > 0]
  p <- v / sum(v)
  c(q0 = length(v),
    q1 = exp(-sum(p * log(p))),
    q2 = 1 / sum(p^2))
}

## Rarefy the whole matrix once and return a data-frame of Hill nums
one_draw <- function(mat, depth){
  rare <- t(rrarefy(t(mat), sample = depth))
  hills <- t(apply(rare, 2, hill_3))
  data.frame(sample = rownames(hills), hills, row.names = NULL)
}

euk_filt_18 <- round(euk_filt_18)
euk_filt_18 <- as.matrix(euk_filt_18)
storage.mode(euk_filt_18) <- "integer"

# Comprobaciones
stopifnot(!any(euk_filt_18 %% 1 != 0))
stopifnot(!any(euk_filt_18 < 0))

N <- 100
set.seed(12345)
cores <- 1   # Windows

out_list_18s <- pbmclapply(
  seq_len(N),
  function(i) one_draw(euk_filt_18, n_reads),
  mc.cores = cores,
  mc.preschedule = TRUE
)

clean_sample_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\.$", "", x)          # remove final dot if present
  x <- sub("\\.[0-9]+$", "", x)    # remove final PCR replicate number
  return(x)
}


alpha_all_18 <- bind_rows(out_list_18s, .id = "draw")

alpha_avg_18 <- alpha_all_18 %>% 
  group_by(sample) %>% 
  summarise(
    across(
      q0:q2,
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ==================== CLEAN SAMPLE NAMES ====================

alpha_avg_18 <- alpha_avg_18 %>%
  mutate(
    Sample_base = clean_sample_name(sample),
    Site = substr(Sample_base, 1, 3)
  )

alpha_avg_reps <- alpha_avg_18 %>% 
  group_by(Sample_base) %>% 
  summarise(
    across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
    across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))),
    .groups = "drop"
  )

# ==================== ORDER SAMPLES ====================

sample_order <- alpha_avg_18 %>%
  group_by(Sample_base) %>%
  summarise(
    median_richness = median(q0_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_richness) %>%
  pull(Sample_base)

alpha_avg_18$Sample_base <- factor(alpha_avg_18$Sample_base, levels = sample_order)

# ==================== COLOURS ====================

colors_site <- c(
  "CMB" = "#458B74",
  "CSD" = "#76EEC6",
  "DEE" = "#CDCD00",
  "HIT" = "#CD96CD",
  "HIO" = "#FFE1FF",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "CLR" = "burlywood1"
)
setdiff(unique(alpha_avg_18$Site), names(colors_site))
# ==================== PLOT ====================

p <- ggplot(alpha_avg_18, aes(x = Sample_base, y = q0_mean)) +
  geom_boxplot(
    aes(fill = Site, group = Sample_base),
    outlier.shape = NA,
    alpha = 0.45,
    colour = "black"
  ) +
  geom_jitter(
    aes(colour = Site),
    width = 0.12,
    size = 2,
    #alpha = 0.55
  ) +
  scale_fill_manual(values = colors_site, drop = FALSE) +
  scale_colour_manual(values = colors_site, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sediment sample",
    y = paste0("Rarefied ASV richness (", n_reads, " reads)"),
    fill = "Site",
    colour = "Site"
  )

p

ggsave(
  filename = "18s_rarefy_read_depth_ordered.svg",
  plot = p,
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# # ==================== SUMMARISE RAREFACTION ====================
# 
# alpha_all_18  <- bind_rows(out_list_18s, .id = "draw")
# 
# alpha_avg_18 <- alpha_all_18 %>%
#   group_by(sample) %>%
#   summarise(
#     across(
#       q0:q2,
#       list(
#         mean = \(x) mean(x, na.rm = TRUE),
#         sd   = \(x) sd(x, na.rm = TRUE)
#       ),
#       .names = "{.col}_{.fn}"
#     ),
#     .groups = "drop"
#   )

# ==================== OBSERVED RICHNESS ====================

observed <- data.frame(
  sample = colnames(euk_filt_18),
  Reads = colSums(euk_filt_18),
  Observed_q0 = colSums(euk_filt_18 > 0)
)

# ==================== MERGE OBSERVED + RAREFIED ====================

compare_18 <- merge(
  observed,
  alpha_avg_18[, c("sample", "q0_mean", "q0_sd")],
  by = "sample"
)

colnames(compare_18)[colnames(compare_18) == "q0_mean"] <- "Rarefied_q0"
colnames(compare_18)[colnames(compare_18) == "q0_sd"] <- "Rarefied_sd"

compare_18$Difference <- compare_18$Observed_q0 - compare_18$Rarefied_q0
compare_18$Percent_loss <- 100 * compare_18$Difference / compare_18$Observed_q0

# ==================== EXTRA METADATA ====================

compare_18$Sample_base <- clean_sample_name(compare_18$sample)
compare_18$Site <- substr(compare_18$Sample_base, 1, 3)

# Opcional: asegurar formato consistente
compare_18$Site <- factor(
  compare_18$Site,
  levels = names(colors_site)
)

# ==================== SAVE TABLE ====================

write.csv(
  compare_18,
  "18S_observed_vs_rarefied.csv",
  row.names = FALSE
)



# ==================== CORRELATIONS ====================

cor_obs_rare <- cor.test(
  compare_18$Observed_q0,
  compare_18$Rarefied_q0,
  method = "spearman"
)

print(cor_obs_rare)

cor_reads_rare <- cor.test(
  compare_18$Reads,
  compare_18$Rarefied_q0,
  method = "spearman"
)

print(cor_reads_rare)

lm_obs_rare <- lm(
  Rarefied_q0 ~ Observed_q0,
  data = compare_18
)

summary(lm_obs_rare)

lm_reads_rare <- lm(
  Rarefied_q0 ~ log10(Reads),
  data = compare_18
)

summary(lm_reads_rare)

# ==================== SUMMARY ====================

sink("18S_rarefaction_summary.txt")

cat("Read depth used for rarefaction:\n")
print(n_reads)

cat("\nObserved richness:\n")
print(summary(compare_18$Observed_q0))

cat("\nRarefied richness:\n")
print(summary(compare_18$Rarefied_q0))

cat("\nDifference (Observed - Rarefied):\n")
print(summary(compare_18$Difference))

cat("\nPercent richness loss:\n")
print(summary(compare_18$Percent_loss))

cat("\nMean richness loss (%):\n")
print(mean(compare_18$Percent_loss))

cat("\nMedian richness loss (%):\n")
print(median(compare_18$Percent_loss))

cat("Observed vs Rarefied richness\n\n")
print(cor_obs_rare)

cat("\n\n")
print(summary(lm_obs_rare))

cat("\n\nReads vs Rarefied richness\n\n")
print(cor_reads_rare)

cat("\n\n")
print(summary(lm_reads_rare))

sink()

colors_site <- c(
  "CMB" = "#458B74",
  "CSD" = "#76EEC6",
  "DEE" = "#CDCD00",
  "HIT" = "#CD96CD",
  "HIO" = "#FFE1FF",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "CLR" = "burlywood1"
)

# ==================== FIGURE 1 ====================
# Observed vs Rarefied

p1 <- ggplot(compare_18,
             aes(Observed_q0,
                 Rarefied_q0,
                 colour = Site)) +
  geom_point(size = 2.2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = TRUE
  ) +
  scale_colour_manual(values = colors_site) +
  theme_bw() +
  labs(
    x = "Observed richness",
    y = paste0("Rarefied richness (", n_reads, " reads)"),
    colour = "Site"
  )

ggsave(
  "18S_observed_vs_rarefied.svg",
  plot = p1,
  width = 7,
  height = 6
)

# ==================== FIGURE 2 ====================
# Richness loss

p2 <- ggplot(compare_18,
             aes(Difference)) +
  geom_histogram(
    bins = 30,
    fill = "grey75",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    x = "Observed − Rarefied richness",
    y = "Number of PCR replicates"
  )

ggsave(
  "18S_rarefaction_loss.svg",
  plot = p2,
  width = 7,
  height = 6
)

# ==================== FIGURE 3 ====================
# Rarefied richness vs read depth

p3 <- ggplot(compare_18,
             aes(Reads,
                 Rarefied_q0,
                 colour = Site)) +
  geom_point(size = 2.2) +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = TRUE
  ) +
  scale_x_log10() +
  scale_colour_manual(values = colors_site) +
  theme_bw() +
  labs(
    x = "Read depth",
    y = paste0("Rarefied richness (", n_reads, " reads)"),
    colour = "Site"
  )

ggsave(
  "18S_rarefied_vs_reads.svg",
  plot = p3,
  width = 7,
  height = 6
)


# ============================================================
# iNEXT RAREFACTION AND SAMPLE-COVERAGE CURVES
# ============================================================

library(iNEXT)
library(ggplot2)
library(dplyr)

# euk_filt_18 tiene:
# filas = ASVs
# columnas = PCR replicates

# Convertir cada PCR en un vector de abundancias > 0 para iNEXT
inext_input_18 <- lapply(
  seq_len(ncol(euk_filt_18)),
  function(i) {
    x <- euk_filt_18[, i]
    x[x > 0]
  }
)

names(inext_input_18) <- colnames(euk_filt_18)

# Eliminar cualquier PCR vacía, por seguridad
inext_input_18 <- inext_input_18[
  vapply(inext_input_18, sum, numeric(1)) > 0
]

# Ejecutar iNEXT solo para riqueza, q = 0
# endpoint = máximo observado; puedes aumentarlo si quieres extrapolación
curve_endpoint <- 100000

out_iNEXT_18 <- iNEXT(
  inext_input_18,
  q = 0,
  datatype = "abundance",
  endpoint = curve_endpoint,
  knots = 40,
  nboot = 0
)

saveRDS(
  out_iNEXT_18,
  file = "18S_iNEXT_curves_output.rds"
)

# Guardar las tablas generadas
write.csv(
  out_iNEXT_18$iNextEst$size_based,
  "18S_iNEXT_size_based_curves.csv",
  row.names = FALSE
)

write.csv(
  out_iNEXT_18$iNextEst$coverage_based,
  "18S_iNEXT_coverage_based_curves.csv",
  row.names = FALSE
)

write.csv(
  out_iNEXT_18$DataInfo,
  "18S_iNEXT_sample_coverage.csv",
  row.names = FALSE
)


# ==================== PREPARE SIZE-BASED DATA ====================

size_curve_18 <- out_iNEXT_18$iNextEst$size_based %>%
  mutate(
    Assemblage = as.character(Assemblage),
    
    # Nombre original de PCR
    PCR = Assemblage,
    
    # Nombre limpio de muestra
    Sample_base = clean_sample_name(PCR),
    
    # Las tres primeras letras
    Site = substr(Sample_base, 1, 3)
  ) %>%
  filter(Site %in% names(colors_site))

# Algunas versiones de iNEXT llaman m o t al número de reads
if ("m" %in% names(size_curve_18)) {
  size_curve_18$Reads_curve <- as.numeric(size_curve_18$m)
} else if ("t" %in% names(size_curve_18)) {
  size_curve_18$Reads_curve <- as.numeric(size_curve_18$t)
} else {
  stop(
    paste(
      "No se encontró la columna de tamaño de muestra.",
      "Columnas disponibles:",
      paste(names(size_curve_18), collapse = ", ")
    )
  )
}

p_rarefaction_18 <- ggplot(
  size_curve_18,
  aes(
    x = Reads_curve,
    y = qD,
    group = PCR,
    colour = Site
  )
) +
  geom_line(alpha = 0.28, linewidth = 0.45) +
  geom_vline(
    xintercept = n_reads,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = scales::label_number_si()
  ) +
  theme_bw() +
  labs(
    x = "Number of reads",
    y = "Estimated ASV richness",
    colour = "Site",
    title = "18S sample-size-based rarefaction curves"
  )

p_rarefaction_18

ggsave(
  filename = "18S_rarefaction_curves.svg",
  plot = p_rarefaction_18,
  width = 9,
  height = 7,
  bg = "white"
)

ggsave(
  filename = "18S_rarefaction_curves.png",
  plot = p_rarefaction_18,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ==================== PREPARE COVERAGE-BASED DATA ====================

coverage_curve_18 <- out_iNEXT_18$iNextEst$coverage_based %>%
  mutate(
    Assemblage = as.character(Assemblage),
    PCR = Assemblage,
    Sample_base = clean_sample_name(PCR),
    Site = substr(Sample_base, 1, 3),
    SC = as.numeric(SC),
    qD = as.numeric(qD)
  ) %>%
  filter(Site %in% names(colors_site))

p_coverage_18 <- ggplot(
  coverage_curve_18,
  aes(
    x = SC,
    y = qD,
    group = PCR,
    colour = Site
  )
) +
  geom_line(alpha = 0.28, linewidth = 0.45) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  theme_bw() +
  labs(
    x = "Sample coverage",
    y = "Estimated ASV richness",
    colour = "Site",
    title = "18S coverage-based rarefaction curves"
  )

p_coverage_18

ggsave(
  filename = "18S_coverage_based_curves.svg",
  plot = p_coverage_18,
  width = 9,
  height = 7,
  bg = "white"
)

ggsave(
  filename = "18S_coverage_based_curves.png",
  plot = p_coverage_18,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

coverage_points_18 <- out_iNEXT_18$DataInfo %>%
  mutate(
    Assemblage = as.character(Assemblage),
    PCR = Assemblage,
    Sample_base = clean_sample_name(PCR),
    Site = substr(Sample_base, 1, 3),
    SC = as.numeric(SC)
  ) %>%
  filter(Site %in% names(colors_site))

# Ordenar las PCR por cobertura
coverage_order_18 <- coverage_points_18 %>%
  arrange(SC) %>%
  pull(PCR)

coverage_points_18$PCR <- factor(
  coverage_points_18$PCR,
  levels = coverage_order_18
)

p_sample_coverage_18 <- ggplot(
  coverage_points_18,
  aes(
    x = PCR,
    y = SC,
    colour = Site
  )
) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      max(0, min(coverage_points_18$SC, na.rm = TRUE) - 0.01),
      1
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "PCR replicates ordered by sample coverage",
    y = "Observed sample coverage",
    colour = "Site",
    title = "18S sample coverage across PCR replicates"
  )

p_sample_coverage_18

ggsave(
  filename = "18S_sample_coverage.svg",
  plot = p_sample_coverage_18,
  width = 9,
  height = 6,
  bg = "white"
)

ggsave(
  filename = "18S_sample_coverage.png",
  plot = p_sample_coverage_18,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)




#####################    COI   ###################
#======== LUKE'S SCRIPT RAREFACTION =================
file_coi_rep <- read.csv(coi_file)
file_coi_rep<-file_coi_rep[,c(2:848)]

# Eliminar muestras (columnas) cuya suma sea < 1000
cols_keep <- c(
  1:11,  # Mantener siempre las primeras 11 columnas
  which(colSums(file_coi_rep[, 12:ncol(file_coi_rep)], na.rm = TRUE) >= 1000) + 11
)

file_coi_rep <- file_coi_rep[, cols_keep]

# Exportar a CSV
write.csv(file_coi_rep,
          "All_Peninsula_COI_Reps_AbTot.csv",
          row.names = TRUE)

# Convertir a abundancias relativas (solo las columnas de muestras)
file_coi_rep_rel <- file_coi_rep
file_coi_rep_rel[, 12:ncol(file_coi_rep_rel)] <-
  sweep(file_coi_rep_rel[, 12:ncol(file_coi_rep_rel)],
        2,
        colSums(file_coi_rep_rel[, 12:ncol(file_coi_rep_rel)]),
        "/")

# Guardar abundancias relativas
write.csv(file_coi_rep_rel,
          "All_Peninsula_COI_Reps_AbRel.csv",
          row.names = TRUE)

row.names(file_coi_rep) <- file_coi_rep$ASV
file_coi<-file_coi_rep[,c(12:843)]
file_coi[] <- lapply(file_coi, as.numeric)
### Alpha diversity 
#lets look at uneven read depth 
lib_sizes <- colSums(file_coi) 
min(lib_sizes)                        
quantile(lib_sizes) 
#ignore samples below 1k
n_reads <- 10000
euk_filt <- as.matrix(file_coi[ , lib_sizes >= n_reads]) 


#some functions to help rarefy & calc

## Rarefy the whole matrix once and return a data-frame of Hill nums
one_draw <- function(mat, depth){
  rare <- t(rrarefy(t(mat), sample = depth))
  hills <- t(apply(rare, 2, hill_3))
  data.frame(sample = rownames(hills), hills, row.names = NULL)
}

euk_filt <- round(euk_filt)
euk_filt <- as.matrix(euk_filt)
storage.mode(euk_filt) <- "integer"

# Comprobaciones
stopifnot(!any(euk_filt %% 1 != 0))
stopifnot(!any(euk_filt < 0))

N <- 100
set.seed(12345)
cores <- 1   # Windows

out_list <- pbmclapply(
  seq_len(N),
  function(i) one_draw(euk_filt, n_reads),
  mc.cores = cores,
  mc.preschedule = TRUE
)

alpha_all <- bind_rows(out_list, .id = "draw")

alpha_avg <- alpha_all %>% 
  group_by(sample) %>% 
  summarise(
    across(
      q0:q2,
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ==================== CLEAN SAMPLE NAMES ====================

alpha_avg <- alpha_avg %>%
  mutate(
    Sample_base = clean_sample_name(sample),
    Site = substr(Sample_base, 1, 3)
  )

alpha_avg_reps <- alpha_avg %>% 
  group_by(Sample_base) %>% 
  summarise(
    across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
    across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))),
    .groups = "drop"
  )

# ==================== ORDER SAMPLES ====================

sample_order <- alpha_avg %>%
  group_by(Sample_base) %>%
  summarise(
    median_richness = median(q0_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_richness) %>%
  pull(Sample_base)

alpha_avg$Sample_base <- factor(alpha_avg$Sample_base, levels = sample_order)

# ==================== COLOURS ====================

colors_site <- c(
  "CMB" = "#458B74",
  "CSD" = "#76EEC6",
  "DEE" = "#CDCD00",
  "HIT" = "#CD96CD",
  "HIO" = "#FFE1FF",
  "STO" = "#8B7355",
  "CCO" = "#CDAA7D",
  "CLR" = "burlywood1"
)

# ==================== PLOT ====================

p <- ggplot(alpha_avg, aes(x = Sample_base, y = q0_mean)) +
  geom_boxplot(
    aes(fill = Site, group = Sample_base),
    outlier.shape = NA,
    #alpha = 0.45,
    colour = "black"
  ) +
  geom_jitter(
    aes(colour = Site),
    width = 0.12,
    size = 2,
    alpha = 0.55
  ) +
  scale_fill_manual(values = colors_site, drop = FALSE) +
  scale_colour_manual(values = colors_site, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sediment sample",
    y = paste0("Rarefied ASV richness (", n_reads, " reads)"),
    fill = "Site",
    colour = "Site"
  )

p

ggsave(
  filename = "COI_rarefy_read_depth_ordered.svg",
  plot = p,
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
# ==================== SUMMARISE RAREFACTION ====================

alpha_all <- bind_rows(out_list, .id = "draw")

alpha_avg <- alpha_all %>%
  group_by(sample) %>%
  summarise(
    across(
      q0:q2,
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ==================== OBSERVED RICHNESS ====================

observed <- data.frame(
  sample = colnames(euk_filt),
  Reads = colSums(euk_filt),
  Observed_q0 = colSums(euk_filt > 0)
)

# ==================== MERGE OBSERVED + RAREFIED ====================


compare <- merge(
  observed,
  alpha_avg[, c("sample", "q0_mean", "q0_sd")],
  by = "sample"
)

colnames(compare)[colnames(compare) == "q0_mean"] <- "Rarefied_q0"
colnames(compare)[colnames(compare) == "q0_sd"] <- "Rarefied_sd"

compare$Difference <- compare$Observed_q0 - compare$Rarefied_q0
compare$Percent_loss <- 100 * compare$Difference / compare$Observed_q0

# ==================== EXTRA METADATA ====================

compare$Sample_base <- clean_sample_name(compare$sample)
compare$Site <- substr(compare$Sample_base, 1, 3)

# Opcional: asegurar formato consistente
compare$Site <- factor(
  compare$Site,
  levels = names(colors_site)
)

# ==================== SAVE TABLE ====================

write.csv(
  compare,
  "COI_observed_vs_rarefied.csv",
  row.names = FALSE
)



# ==================== CORRELATIONS ====================

cor_obs_rare <- cor.test(
  compare$Observed_q0,
  compare$Rarefied_q0,
  method = "spearman"
)

print(cor_obs_rare)

cor_reads_rare <- cor.test(
  compare$Reads,
  compare$Rarefied_q0,
  method = "spearman"
)

print(cor_reads_rare)

lm_obs_rare <- lm(
  Rarefied_q0 ~ Observed_q0,
  data = compare
)

summary(lm_obs_rare)

lm_reads_rare <- lm(
  Rarefied_q0 ~ log10(Reads),
  data = compare
)

summary(lm_reads_rare)


# ==================== SUMMARY ====================

sink("COI_rarefaction_summary.txt")

cat("Read depth used for rarefaction:\n")
print(n_reads)

cat("\nObserved richness:\n")
print(summary(compare$Observed_q0))

cat("\nRarefied richness:\n")
print(summary(compare$Rarefied_q0))

cat("\nDifference (Observed - Rarefied):\n")
print(summary(compare$Difference))

cat("\nPercent richness loss:\n")
print(summary(compare$Percent_loss))

cat("\nMean richness loss (%):\n")
print(mean(compare$Percent_loss))

cat("\nMedian richness loss (%):\n")
print(median(compare$Percent_loss))

cat("Observed vs Rarefied richness\n\n")
print(cor_obs_rare)

cat("\n\n")
print(summary(lm_obs_rare))

cat("\n\nReads vs Rarefied richness\n\n")
print(cor_reads_rare)

cat("\n\n")
print(summary(lm_reads_rare))

sink()

# ==================== FIGURE 1 ====================
# Observed vs Rarefied

p1 <- ggplot(compare,
             aes(Observed_q0,
                 Rarefied_q0,
                 colour = Site)) +
  geom_point(size = 2.2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = TRUE
  ) +
  scale_colour_manual(values = colors_site) +
  theme_bw() +
  labs(
    x = "Observed richness",
    y = paste0("Rarefied richness (", n_reads, " reads)"),
    colour = "Site"
  )

ggsave(
  "COI_observed_vs_rarefied.svg",
  plot = p1,
  width = 7,
  height = 6
)

# ==================== FIGURE 2 ====================
# Richness loss

p2 <- ggplot(compare,
             aes(Difference)) +
  geom_histogram(
    bins = 30,
    fill = "grey75",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    x = "Observed − Rarefied richness",
    y = "Number of PCR replicates"
  )

ggsave(
  "COI_rarefaction_loss.svg",
  plot = p2,
  width = 7,
  height = 6
)

# ==================== FIGURE 3 ====================
# Rarefied richness vs read depth

p3 <- ggplot(compare,
             aes(Reads,
                 Rarefied_q0,
                 colour = Site)) +
  geom_point(size = 2.2) +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = TRUE
  ) +
  scale_x_log10() +
  scale_colour_manual(values = colors_site) +
  theme_bw() +
  labs(
    x = "Read depth",
    y = paste0("Rarefied richness (", n_reads, " reads)"),
    colour = "Site"
  )

ggsave(
  "COI_rarefied_vs_reads.svg",
  plot = p3,
  width = 7,
  height = 6
)


# ============================================================
# iNEXT RAREFACTION AND SAMPLE-COVERAGE CURVES
# ============================================================

library(iNEXT)
library(ggplot2)
library(dplyr)

# euk_filt_coi tiene:
# filas = ASVs
# columnas = PCR replicates

# Convertir cada PCR en un vector de abundancias > 0 para iNEXT
inext_input_coi <- lapply(
  seq_len(ncol(euk_filt_coi)),
  function(i) {
    x <- euk_filt_coi[, i]
    x[x > 0]
  }
)

names(inext_input_coi) <- colnames(euk_filt_coi)

# Eliminar cualquier PCR vacía, por seguridad
inext_input_coi <- inext_input_coi[
  vapply(inext_input_coi, sum, numeric(1)) > 0
]

# Ejecutar iNEXT solo para riqueza, q = 0
# endpoint = máximo observado; puedes aumentarlo si quieres extrapolación
max_reads_coi <- max(colSums(euk_filt_coi))

out_iNEXT_coi <- iNEXT(
  inext_input_coi,
  q = 0,
  datatype = "abundance",
  endpoint = max_reads_coi,
  knots = 40,
  nboot = 0
)

saveRDS(
  out_iNEXT_coi,
  file = "COI_iNEXT_curves_output.rds"
)

# Guardar las tablas generadas
write.csv(
  out_iNEXT_coi$iNextEst$size_based,
  "COI_iNEXT_size_based_curves.csv",
  row.names = FALSE
)

write.csv(
  out_iNEXT_18$iNextEst$coverage_based,
  "COI_iNEXT_coverage_based_curves.csv",
  row.names = FALSE
)

write.csv(
  out_iNEXT_coi$DataInfo,
  "COI_iNEXT_sample_coverage.csv",
  row.names = FALSE
)


# ==================== PREPARE SIZE-BASED DATA ====================

size_curve_coi <- out_iNEXT_coi$iNextEst$size_based %>%
  mutate(
    Assemblage = as.character(Assemblage),
    
    # Nombre original de PCR
    PCR = Assemblage,
    
    # Nombre limpio de muestra
    Sample_base = clean_sample_name(PCR),
    
    # Las tres primeras letras
    Site = substr(Sample_base, 1, 3)
  ) %>%
  filter(Site %in% names(colors_site))

# Algunas versiones de iNEXT llaman m o t al número de reads
if ("m" %in% names(size_curve_coi)) {
  size_curve_coi$Reads_curve <- as.numeric(size_curve_coi$m)
} else if ("t" %in% names(size_curve_coi)) {
  size_curve_coi$Reads_curve <- as.numeric(size_curve_coi$t)
} else {
  stop(
    paste(
      "No se encontró la columna de tamaño de muestra.",
      "Columnas disponibles:",
      paste(names(size_curve_coi), collapse = ", ")
    )
  )
}

p_rarefaction_coi <- ggplot(
  size_curve_coi,
  aes(
    x = Reads_curve,
    y = qD,
    group = PCR,
    colour = Site
  )
) +
  geom_line(alpha = 0.28, linewidth = 0.45) +
  geom_vline(
    xintercept = n_reads,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = scales::label_number_si()
  ) +
  theme_bw() +
  labs(
    x = "Number of reads",
    y = "Estimated ASV richness",
    colour = "Site",
    title = "COI sample-size-based rarefaction curves"
  )

p_rarefaction_18

ggsave(
  filename = "COI_rarefaction_curves.svg",
  plot = p_rarefaction_coi,
  width = 9,
  height = 7,
  bg = "white"
)

ggsave(
  filename = "COI_rarefaction_curves.png",
  plot = p_rarefaction_coi,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ==================== PREPARE COVERAGE-BASED DATA ====================

coverage_curve_coi <- out_iNEXT_coi$iNextEst$coverage_based %>%
  mutate(
    Assemblage = as.character(Assemblage),
    PCR = Assemblage,
    Sample_base = clean_sample_name(PCR),
    Site = substr(Sample_base, 1, 3),
    SC = as.numeric(SC),
    qD = as.numeric(qD)
  ) %>%
  filter(Site %in% names(colors_site))

p_coverage_coi <- ggplot(
  coverage_curve_coi,
  aes(
    x = SC,
    y = qD,
    group = PCR,
    colour = Site
  )
) +
  geom_line(alpha = 0.28, linewidth = 0.45) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  theme_bw() +
  labs(
    x = "Sample coverage",
    y = "Estimated ASV richness",
    colour = "Site",
    title = "COI coverage-based rarefaction curves"
  )

p_coverage_coi

ggsave(
  filename = "COI_coverage_based_curves.svg",
  plot = p_coverage_coi,
  width = 9,
  height = 7,
  bg = "white"
)

ggsave(
  filename = "COI_coverage_based_curves.png",
  plot = p_coverage_coi,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

coverage_points_coi <- out_iNEXT_coi$DataInfo %>%
  mutate(
    Assemblage = as.character(Assemblage),
    PCR = Assemblage,
    Sample_base = clean_sample_name(PCR),
    Site = substr(Sample_base, 1, 3),
    SC = as.numeric(SC)
  ) %>%
  filter(Site %in% names(colors_site))

# Ordenar las PCR por cobertura
coverage_order_coi <- coverage_points_coi %>%
  arrange(SC) %>%
  pull(PCR)

coverage_points_coi$PCR <- factor(
  coverage_points_coi$PCR,
  levels = coverage_order_coi
)

p_sample_coverage_coi <- ggplot(
  coverage_points_coi,
  aes(
    x = PCR,
    y = SC,
    colour = Site
  )
) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  scale_colour_manual(
    values = colors_site,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      max(0, min(coverage_points_coi$SC, na.rm = TRUE) - 0.01),
      1
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "PCR replicates ordered by sample coverage",
    y = "Observed sample coverage",
    colour = "Site",
    title = "COI sample coverage across PCR replicates"
  )

p_sample_coverage_coi

ggsave(
  filename = "COI_sample_coverage.svg",
  plot = p_sample_coverage_coi,
  width = 9,
  height = 6,
  bg = "white"
)

ggsave(
  filename = "COI_sample_coverage.png",
  plot = p_sample_coverage_coi,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)
