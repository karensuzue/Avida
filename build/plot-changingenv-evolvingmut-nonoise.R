# This script produces a plot displaying both final fitness and
# per-gene mutation rate as a function of the rate of environment change.
# For each rate of environment change, mutation rate and fitness values
# were time-averaged over the last K updates of each of the last N complete
# environment turnovers, then averaged over 20 replicates.

# install.packages(
#     c("tidyverse", "data.table", "patchwork", "furrr"),
#     lib = "~/R/library",
#     repos = "https://cloud.r-project.org"
# )

# .libPaths(c("~/R/library", .libPaths()))


library(tidyverse)
library(data.table)
library(patchwork)
# library(furrr)

# ---------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------- 
DATA_DIR <- "/mnt/d/MINE/ce-em-nn-data/"
OUT_DIR <- "."

TURNOVER_COUNT <- 10

# Number of updates at the end of each turnover window to average over.
TAIL_WINDOW <- 10

# Filename pattern:
# ce_em_nn_{change_per_update}_{seed}.csv
REGEX_PATTERN <- "^ce_em_nn_([^_]+)_([0-9]+)\\.csv$"

# MAX_GEN <- 200000
# GENOME_LENGTH <- 100

# Use the cores Slurm actually allocated to this job, not the whole node's count
# n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
# plan(multisession, workers = max(1, n_cores))

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

# Grab the relevant rows of one file (last TURNOVER_COUNT turnovers),
# tagged with the rate of environment change they belong to.
#
# For each turnover, instead of taking a single row at the turnover
# boundary, takes the mean of the numeric columns over the last
# TAIL_WINDOW updates of that turnover (clamped to not reach into the
# previous turnover for very short windows).
#
# Uses data.table::fread instead of readr::read_csv for I/O speed
read_turnover_rows <- function(file) {
    meta <- parse_filename(file)
    if (is.null(meta)) {
        print(paste("Error reading turnover rows for:", file))
        return(NULL)
    }
    change_per_update <- meta$change_per_update

    # Number of turnovers = 2000 * change per update
    # Number of updates per turnover = 100 / change per update
    updates_per_turnover <- 100 / change_per_update
    num_turnovers <- 2000 * change_per_update

    # Cap at TURNOVER_COUNT, but don't request more turnovers than exist
    if (num_turnovers > TURNOVER_COUNT) num_turnovers <- TURNOVER_COUNT

    # Don't let the tail window exceed the turnover length itself
    # (relevant for very fast change_per_update, e.g. 1-2 updates/turnover)
    tail_window <- min(TAIL_WINDOW, floor(updates_per_turnover))
    
    # How many rows from the end of the file could we possibly need?
    rows_needed <- ceiling(num_turnovers * updates_per_turnover) + tail_window

    # Read only the tail of the file, plus the header separately
    df <- fread(cmd = paste("tail -n", rows_needed + 1, shQuote(file)), header = FALSE)
    setnames(df, names(fread(file, nrows = 0)))

    rows <- list()
    for (turnover in 0:(num_turnovers - 1)) {
        turnover_end_index <- nrow(df) - round(updates_per_turnover * turnover)
        turnover_start_index <- turnover_end_index - tail_window + 1

        if (turnover_start_index < 1) {
            print(paste("Not enough rows in", file, "for turnover", turnover))
            next
        }

        window_rows <- df[turnover_start_index:turnover_end_index, ]

        # Average the numeric columns over this turnover's tail window
        averaged_row <- window_rows[, lapply(.SD, mean, na.rm = TRUE), .SDcols = is.numeric]

        rows[[length(rows) + 1]] <- averaged_row
    }
    rows_df <- rbindlist(rows) # rows is a list of single-row data.tables, so we bind them

    # Tag rows with which rate of environment change and seed they came from
    rows_df$change_per_update <- change_per_update
    rows_df$seed <- meta$seed

    return(rows_df)
}

# Compile relevant rows across all replicates, then average per change_per_update.
# Resulting data frame:
# change_per_update | max_fitness | avg_fitness | min_mut_rate | avg_mut_rate | max_mut_rate | fittest_org_mut_rate
compile_average <- function(files) { # 'files' is a vector of filenames
    # rows <- future_map(files, read_turnover_rows)
    # rows <- rows[!sapply(rows, is.null)] # drop any files that failed to parse
    # all_rows <- rbindlist(rows)
    rows <- list()
    for (f in files) {
        turnover_rows <- read_turnover_rows(f)
        if (!is.null(turnover_rows)) {
            rows[[length(rows) + 1]] <- turnover_rows
        }
    }
    all_rows <- rbindlist(rows)

    summary_df <- all_rows %>%
        group_by(change_per_update) %>%
        summarise(
            max_fitness = mean(`Fittest Organism Fitness`, na.rm = TRUE),
            avg_fitness = mean(`Average Organism Fitness`, na.rm = TRUE),
            min_mut_rate = mean(`Min Mutation Rate`, na.rm = TRUE),
            avg_mut_rate = mean(`Average Mutation Rate`, na.rm = TRUE),
            max_mut_rate = mean(`Max Mutation Rate`, na.rm = TRUE),
            fittest_org_mut_rate = mean(`Fittest Organism Mutation Rate`, na.rm = TRUE),
            n_rows = n(),
            .groups = "drop"
        ) %>%
        arrange(change_per_update)

    return(summary_df)
}

# ---------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------
files <- list.files(DATA_DIR, pattern = REGEX_PATTERN, full.names = TRUE)
result <- compile_average(files)
result

write_csv(result, file.path(OUT_DIR, "turnover_summary.csv"))

# Reshape fitness columns into long format for a single legend
fitness_long <- result %>%
    select(change_per_update, max_fitness, avg_fitness) %>%
    pivot_longer(
        cols = c(max_fitness, avg_fitness),
        names_to = "series",
        values_to = "value"
    ) %>%
    mutate(series = recode(series,
        max_fitness = "Best organism",
        avg_fitness = "Average"
    ))

# Reshape mutation rate columns into long format
mutation_long <- result %>%
    select(change_per_update, min_mut_rate, avg_mut_rate, max_mut_rate, fittest_org_mut_rate) %>%
    pivot_longer(
        cols = c(min_mut_rate, avg_mut_rate, max_mut_rate, fittest_org_mut_rate),
        names_to = "series",
        values_to = "value"
    ) %>%
    mutate(series = recode(series,
        min_mut_rate = "Min",
        avg_mut_rate = "Average",
        max_mut_rate = "Max",
        fittest_org_mut_rate = "Best organism"
    ))

# ---------------------------------------------------------------------
# PLOT
# ---------------------------------------------------------------------
p_fitness <- ggplot(fitness_long, aes(x = change_per_update, y = value, color = series)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(
        x = NULL,
        y = "Fitness",
        color = NULL,
        title = "Fitness vs. rate of environment change"
    ) +
    theme_minimal()

p_mutation <- ggplot(mutation_long, aes(x = change_per_update, y = value, color = series)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(
        x = "Rate of environment change (genes per update)",
        y = "Mutation rate",
        color = NULL,
        title = "Mutation rate vs. rate of environment change"
    ) +
    theme_minimal()


# Stack into one page, save
combined_plot <- p_fitness / p_mutation

ggsave(
    file.path(OUT_DIR, "fitness_and_mutation_vs_change_rate.pdf"),
    plot = combined_plot, width = 8, height = 9
)

combined_plot