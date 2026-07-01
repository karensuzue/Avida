# This script produces generation vs. fitness 
# and generation vs. mutation rate plots.
# PDF 1: one pair of pages per change_per_update condition
# PDF 2: 2 pages, one for mutation, one for fitness,
# faceted by change_per_update (all conditions on one page)

# .libPaths(c("~/R/library", .libPaths()))
 
library(data.table)
library(tidyverse)
library(patchwork)

# ---------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------- 
DATA_DIR <- "/mnt/scratch/suzuekar/ce-em-nn-data/" 
OUT_DIR <- "."
OUT_NAME <- "ce_em_nn_timeseries.pdf"
OUT_NAME_FACETED <- "ce_em_nn_timeseries_faceted.pdf"

# Reading every generation is heavy
# Set to 1 to read everything, or e.g. 1000 to keep every 1000th generation.
SUBSAMPLE_STEP <- 1000

# Filename pattern:
# ce_em_nn_{change_per_update}_{seed}.csv
REGEX_PATTERN <- "^ce_em_nn_([^_]+)_([0-9]+)\\.csv$"

ROMEO_COLS <- c("Update", "Fittest Organism ID", "Fittest Organism Error",
                "Fittest Organism Fitness", "Average Organism Error",
                "Average Organism Fitness", "Min Mutation Rate",
                "Average Mutation Rate", "Max Mutation Rate")

COLS_NEEDED <- c("Update", "Fittest Organism Fitness", "Average Organism Fitness",
                  "Min Mutation Rate", "Average Mutation Rate", "Max Mutation Rate")

LOG_FLOOR <- 1e-6

# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------
parse_filename <- function(file) {
    nm <- basename(file)
    m <- regexec(REGEX_PATTERN, nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 0) {
        print(paste("Error parsing filename:", nm))
        return(NULL)
    }
    parsed <- list(
        change_per_update = as.double(parts[2]),
        seed = as.integer(parts[3])
    )
    return(parsed)
}

read_file <- function(file) {
    meta <- parse_filename(file)

    # Skip header, take first 9 cols by position
    df <- fread(file, header = FALSE, skip = 1, 
                select = 1:length(ROMEO_COLS))
    setnames(df, ROMEO_COLS) # reattach header

    df <- df %>%
        filter(Update %% SUBSAMPLE_STEP == 0) %>%
        select(Update, `Fittest Organism Fitness`, `Average Organism Fitness`,
               `Min Mutation Rate`, `Average Mutation Rate`, `Max Mutation Rate`) %>%
        mutate(change_per_update = meta$change_per_update,
              seed = meta$seed)
    return(df)
}


make_fitness_plot <- function(df, cpu, ylim) {
    d <- df %>%
        filter(change_per_update == cpu) %>%
        pivot_longer(c(`Fittest Organism Fitness`, `Average Organism Fitness`),
                    names_to = "series", values_to = "value")
 
    # Spread across replicates (seeds) at each Update, per series
    d_summary <- d %>%
        group_by(Update, series) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  lo = min(value, na.rm = TRUE),
                  hi = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(mean = pmax(mean, LOG_FLOOR),
               lo = pmax(lo, LOG_FLOOR),
               hi = pmax(hi, LOG_FLOOR))
 
    ggplot(d_summary, aes(Update, mean, color = series, fill = series)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
        geom_line() +
        scale_y_log10(limits = ylim) +
        labs(title = paste("Fitness, change_per_update =", cpu),
             x = "Update", y = "Fitness", color = "Series", fill = "Series") +
        theme_bw()
}


make_mutrate_plot <- function(df, cpu, ylim) {
    d <- df %>%
        filter(change_per_update == cpu) %>%
        pivot_longer(c(`Min Mutation Rate`, `Average Mutation Rate`, `Max Mutation Rate`),
                     names_to = "series", values_to = "value")
 
    # Spread across replicates (seeds) at each Update, per series
    d_summary <- d %>%
        group_by(Update, series) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  lo = min(value, na.rm = TRUE),
                  hi = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(mean = pmax(mean, LOG_FLOOR),
               lo = pmax(lo, LOG_FLOOR),
               hi = pmax(hi, LOG_FLOOR))
 
    ggplot(d_summary, aes(Update, mean, color = series, fill = series)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
        geom_line() +
        scale_y_log10(limits = ylim) +
        labs(title = paste("Mutation rate, change_per_update =", cpu),
             x = "Update", y = "Mutation rate", color = "Series", fill = "Series") +
        theme_bw()
}

# Label facets as "change_per_update = X" instead of just the bare number
cpu_labeller <- function(value) paste("change_per_update =", value)

make_fitness_facet <- function(df, ylim) {
    d <- df %>%
        pivot_longer(c(`Fittest Organism Fitness`, `Average Organism Fitness`),
                    names_to = "series", values_to = "value")
 
    # Spread across replicates (seeds) at each Update, per series, per condition
    d_summary <- d %>%
        group_by(change_per_update, Update, series) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  lo = min(value, na.rm = TRUE),
                  hi = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(mean = pmax(mean, LOG_FLOOR),
               lo = pmax(lo, LOG_FLOOR),
               hi = pmax(hi, LOG_FLOOR))
 
    ggplot(d_summary, aes(Update, mean, color = series, fill = series)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
        geom_line() +
        scale_y_log10(limits = ylim) +
        facet_wrap(~ change_per_update, scales = "free_y",
                   labeller = labeller(change_per_update = cpu_labeller)) +
        labs(title = "Fitness across all change_per_update conditions",
             x = "Update", y = "Fitness", color = "Series", fill = "Series") +
        theme_bw()
}
 
make_mutrate_facet <- function(df, ylim) {
    d <- df %>%
        pivot_longer(c(`Min Mutation Rate`, `Average Mutation Rate`, `Max Mutation Rate`),
                     names_to = "series", values_to = "value")
 
    # Spread across replicates (seeds) at each Update, per series, per condition
    d_summary <- d %>%
        group_by(change_per_update, Update, series) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  lo = min(value, na.rm = TRUE),
                  hi = max(value, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(mean = pmax(mean, LOG_FLOOR),
               lo = pmax(lo, LOG_FLOOR),
               hi = pmax(hi, LOG_FLOOR))
 
    ggplot(d_summary, aes(Update, mean, color = series, fill = series)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
        geom_line() +
        scale_y_log10(limits = ylim) +
        facet_wrap(~ change_per_update, scales = "free_y",
                   labeller = labeller(change_per_update = cpu_labeller)) +
        labs(title = "Mutation rate across all change_per_update conditions",
             x = "Update", y = "Mutation rate", color = "Series", fill = "Series") +
        theme_bw()
}

# ---------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------
files <- list.files(DATA_DIR, pattern = REGEX_PATTERN, full.names = TRUE)
print(paste("Found", length(files), "files"))

all_data <- map_dfr(files, read_file)

cpu_values <- sort(unique(all_data$change_per_update))
print(paste("Conditions:", paste(cpu_values, collapse = ", ")))

# Global y-axis limits, shared across all plots
fitness_ylim <- c(
    pmax(min(all_data$`Average Organism Fitness`, na.rm = TRUE), LOG_FLOOR),
    max(all_data$`Fittest Organism Fitness`, na.rm = TRUE)
)
mutrate_ylim <- c(
    pmax(min(all_data$`Min Mutation Rate`, na.rm = TRUE), LOG_FLOOR),
    max(all_data$`Max Mutation Rate`, na.rm = TRUE)
)

out_path <- file.path(OUT_DIR, OUT_NAME)
pdf(out_path, width = 11, height = 7)
for (cpu in cpu_values) {
    print(make_fitness_plot(all_data, cpu))
    print(make_mutrate_plot(all_data, cpu))
}
dev.off()

print(paste("Wrote", length(cpu_values) * 2, "pages to", out_path))

out_path <- file.path(OUT_DIR, OUT_NAME_FACETED)
pdf(out_path, width = 14, height = 10)
print(make_fitness_facet(all_data))
print(make_mutrate_facet(all_data))
dev.off()

print(paste("Wrote 2 pages to", out_path))