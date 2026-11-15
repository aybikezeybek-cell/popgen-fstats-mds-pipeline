library(tidyverse)
library(ggrepel)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. OUTGROUP f3 -> (1 - f3) GENETIC DISTANCE MATRIX & CLASSICAL MDS
# ==============================================================================
f3_lines <- grep("^\\s*result:", readLines("outputs/logs/f3_results.log"), value = TRUE)
f3_df <- read.table(text = f3_lines, stringsAsFactors = FALSE)
colnames(f3_df) <- c("tag", "Pop1", "Pop2", "Outgroup", "f3", "SE", "Z", "SNPs")

pops <- sort(unique(c(f3_df$Pop1, f3_df$Pop2)))
dist_mat <- matrix(0, nrow = length(pops), ncol = length(pops),
                   dimnames = list(pops, pops))

for (i in seq_len(nrow(f3_df))) {
  p1 <- f3_df$Pop1[i]
  p2 <- f3_df$Pop2[i]
  d  <- 1 - f3_df$f3[i]
  dist_mat[p1, p2] <- d
  dist_mat[p2, p1] <- d
}

write.csv(round(dist_mat, 6), "outputs/tables/f3_distance_matrix.csv", quote = FALSE)

mds_fit <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE)
var_exp <- round(100 * mds_fit$eig[1:2] / sum(pmax(mds_fit$eig, 0)), 1)

mds_df <- data.frame(
  Population = rownames(mds_fit$points),
  Dim1 = mds_fit$points[, 1],
  Dim2 = mds_fit$points[, 2]
) %>%
  mutate(Region = case_when(
    grepl("^Turkey_", Population) ~ "Anatolia Neolithic",
    grepl("^Jordan_", Population) ~ "Levant Neolithic",
    grepl("^Iran_", Population)   ~ "Zagros / Iran Neolithic",
    TRUE ~ "Other"
  ))

p_mds <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Region, label = Population)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 4.2) +
  geom_text_repel(
    size = 3.5,
    fontface = "bold",
    nudge_x = ifelse(mds_df$Population %in% c("Turkey_Catalhoyuk_MN",
                                              "Turkey_Catalhoyuk_LN"), -0.06, 0.04),
    nudge_y = ifelse(mds_df$Population == "Turkey_N", 0.04, 0),
    box.padding = 0.6,
    point.padding = 0.5,
    force = 5,
    max.overlaps = Inf,
    segment.color = "grey55",
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(-0.2, 0.6, 0.2),
    expand = expansion(mult = c(0.32, 0.20))
  ) +
  scale_y_continuous(
    breaks = seq(-0.2, 0.4, 0.2),
    expand = expansion(mult = c(0.18, 0.15))
  ) +
  scale_color_manual(values = c(
    "Anatolia Neolithic"      = "#1b7837",
    "Levant Neolithic"        = "#2166ac",
    "Zagros / Iran Neolithic" = "#b2182b"
  )) +
  labs(
    title = "Multidimensional Scaling (MDS) of Ancient Near Eastern Populations",
    subtitle = "Pairwise Outgroup f3(Pop1, Pop2; Mbuti) Genetic Distances (1 - f3) | AADR v66",
    x = paste0("Dimension 1 (", var_exp[1], "% variance)"),
    y = paste0("Dimension 2 (", var_exp[2], "% variance)"),
    color = "Geographic Horizon"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("outputs/figures/mds_anatolia_neolithic.png", p_mds,
       width = 9.5, height = 6.0, dpi = 300, bg = "white")

# ==============================================================================
# 2. f4-STATISTICS (qpDstat) -> SUMMARY TABLE & GENE FLOW FOREST PLOT
# ==============================================================================
f4_lines <- grep("^\\s*result:", readLines("outputs/logs/f4_results.log"), value = TRUE)
f4_df <- read.table(text = f4_lines, stringsAsFactors = FALSE)
colnames(f4_df) <- c("tag", "W", "X", "Y", "Z_out", "f4", "Z_score",
                     "BABA", "ABBA", "SNPs")

f4_processed <- f4_df %>%
  mutate(
    SE = abs(f4 / Z_score),
    Test = paste0("f4(", W, ", ", X, "; ", Y, ", Mbuti)"),
    Direction = case_when(
      Z_score > 3  ~ "Excess affinity: (X & Y) [Z > 3]",
      Z_score < -3 ~ "Excess affinity: (W & Y) [Z < -3]",
      TRUE         ~ "Not significant (|Z| <= 3)"
    )
  ) %>%
  arrange(desc(abs(Z_score)))

write.csv(
  f4_processed %>% select(W, X, Y, Outgroup = Z_out, f4, SE, Z_score, BABA, ABBA, SNPs),
  "outputs/tables/f4_summary_results.csv",
  row.names = FALSE,
  quote = FALSE
)

f4_top <- f4_processed %>% slice_head(n = 15)

p_f4 <- ggplot(f4_top, aes(x = f4, y = reorder(Test, f4), color = Direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = f4 - 3 * SE, xmax = f4 + 3 * SE),
                 height = 0.25, linewidth = 0.75) +
  geom_point(size = 3.2) +
  geom_label(
    aes(label = sprintf("Z = %.1f", Z_score)),
    vjust = -0.6,
    size = 3.0,
    fill = "white",
    label.size = 0,
    label.padding = unit(0.12, "lines"),
    show.legend = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.8))) +
  scale_color_manual(values = c(
    "Excess affinity: (X & Y) [Z > 3]"  = "#d95f02",
    "Excess affinity: (W & Y) [Z < -3]" = "#1f78b4",
    "Not significant (|Z| <= 3)"        = "#7f7f7f"
  )) +
  labs(
    title = "Gene Flow Signals via f4-Statistics (qpDstat)",
    subtitle = "Form: f4(W, X; Y, Outgroup) | Point estimates with +/- 3 SE jackknife intervals | AADR v66 1240K",
    x = "f4-statistic",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("outputs/figures/f4_statistics_plot.png", p_f4,
       width = 10.5, height = 6.2, dpi = 300, bg = "white")
