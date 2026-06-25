# Diagnostic plot: for a couple of representative files, show mutation rate
# over the last N turnover windows. Used to check whether mutation rate
# converges (flattens) within a turnover, which determines whether the
# main analysis should sample a single endpoint row per turnover or
# average across rows near the end of each window.

library(tidyverse)

# ---------------------------------------------------------------------
# CONFIG 
# ---------------------------------------------------------------------
DATA_DIR <- "/mnt/d/MINE/ce-em-nn-data"
OUT_DIR <- "."

# Files to inspect
FILES_TO_CHECK <- c(
    # file.path(DATA_DIR, "ce_em_nn_0.001_1.csv"),
    # file.path(DATA_DIR, "ce_em_nn_0.002_1.csv"),
    # file.path(DATA_DIR, "ce_em_nn_0.005_1.csv"),
    # file.path(DATA_DIR, "ce_em_nn_0.01_1.csv"),
    # file.path(DATA_DIR, "ce_em_nn_0.02_1.csv"),
    # file.path(DATA_DIR, "ce_em_nn_0.05_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_0.1_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_0.2_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_0.5_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_1_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_2_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_10_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_20_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_50_1.csv"),
    file.path(DATA_DIR, "ce_em_nn_100_1.csv")
)

REGEX_PATTERN <- "^ce_em_nn_([^_]+)_([0-9]+)\\.csv$"

# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

# Plot the last n_turnovers_to_show turnover windows for one file,
# mutation rate vs. update. Returns the plot plus metadata used for naming.
plot_turnover_window <- function(file, n_turnovers_to_show = 100) {
    df <- read_csv(file, show_col_types = FALSE)

    nm <- basename(file)
    m <- regexec(REGEX_PATTERN, nm)
    parts <- regmatches(nm, m)[[1]]
    change_per_update <- as.double(parts[2])

    updates_per_turnover <- 100 / change_per_update
    max_update <- max(df$Update)

    # Go back...
    window_start <- max_update - n_turnovers_to_show * updates_per_turnover

    df_window <- df %>% filter(Update <= window_start)

    p <- ggplot(df_window, aes(x = Update)) +
        geom_line(aes(y = `Average Mutation Rate`, color = "Average")) +
        geom_line(aes(y = `Fittest Organism Mutation Rate`, color = "Fittest")) +
        geom_vline(
            xintercept = seq(max_update, window_start, by = -updates_per_turnover),
            linetype = "dashed", color = "gray50"
        ) +
        labs(
            title = paste("change_per_update =", change_per_update),
            subtitle = paste("dashed lines = turnover boundaries, every", updates_per_turnover, "updates"),
            x = "Update", y = "Mutation rate"
        ) +
        theme_minimal()

        list(plot = p, change_per_update = change_per_update)

    return(p)
}

# ---------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------

# All plots get drawn into one multi-page PDF, one page per file
out_filename <- file.path(OUT_DIR, "turnover_windows.pdf")

pdf(out_filename, width = 8, height = 5)
for (file in FILES_TO_CHECK) {
    p <- plot_turnover_window(file)
    print(p)
}
dev.off()