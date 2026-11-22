# ================================
# Final version: Life Cycle Plot
# Export as 600 dpi TIFF / PNG
# ================================

if (!exists("package_root")) stop("Run run_all.R from the package root.")

library(ggplot2)
library(patchwork)

# 1. Data
df <- read.csv(
  file.path(package_root, "Input_Data", "Annual_Publication_Counts_2012_2025.csv")
)
stopifnot(
  identical(names(df), c("year", "publications")),
  nrow(df) == 14L,
  is.numeric(df$year),
  is.numeric(df$publications),
  all(is.finite(df$year)),
  all(is.finite(df$publications)),
  all(df$year == 2012:2025),
  all(df$publications >= 0 & df$publications == floor(df$publications)),
  sum(df$publications) == 770
)

df$cumulative <- cumsum(df$publications)

# 2. Reconstruct the model from recorded Biblioshiny parameters; do not refit.
parameter_data <- read.csv(
  file.path(package_root, "Input_Data", "Lifecycle_Model_Parameters.csv")
)
stopifnot(
  identical(names(parameter_data), c("parameter", "value")),
  nrow(parameter_data) == 3L,
  !anyDuplicated(parameter_data$parameter),
  setequal(parameter_data$parameter, c("K", "peak_year", "peak_publications")),
  is.numeric(parameter_data$value),
  all(is.finite(parameter_data$value)),
  all(parameter_data$value > 0)
)
model_parameters <- setNames(parameter_data$value, parameter_data$parameter)
K <- unname(model_parameters["K"])
peak_year <- unname(model_parameters["peak_year"])
peak_year_display <- round(peak_year)
peak_publications <- unname(model_parameters["peak_publications"])
peak_label <- paste0("Peak year (", peak_year_display, ")")

# 3. Logistic life-cycle model
r <- 4 * peak_publications / K

pred_df <- data.frame(
  year = seq(2012, 2055, length.out = 800)
)

pred_df$cumulative_fit <- K / (1 + exp(-r * (pred_df$year - peak_year)))

z <- exp(-r * (pred_df$year - peak_year))
pred_df$annual_fit <- K * r * z / (1 + z)^2

fit_observed <- subset(pred_df, year <= max(df$year))
fit_forecast <- subset(pred_df, year >= max(df$year))

# 4. Calculate R²
df$annual_fit <- K * r * exp(-r * (df$year - peak_year)) /
  (1 + exp(-r * (df$year - peak_year)))^2

ss_res <- sum((df$publications - df$annual_fit)^2)
ss_tot <- sum((df$publications - mean(df$publications))^2)
r2 <- 1 - ss_res / ss_tot

# ================================
# 5. Left plot
# ================================

p1 <- ggplot() +
  geom_line(
    data = fit_observed,
    aes(x = year, y = annual_fit),
    color = "#1f77b4",
    linewidth = 0.8,
    linetype = "solid"
  ) +
  geom_line(
    data = fit_forecast,
    aes(x = year, y = annual_fit),
    color = "#1f77b4",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(
    data = df,
    aes(x = year, y = publications),
    color = "#1f77b4",
    size = 2
  ) +
  geom_vline(
    xintercept = peak_year,
    color = "red",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  
  # Manual legend
  annotate(
    "rect",
    xmin = 2011.5,
    xmax = 2027.5,
    ymin = 275,
    ymax = 306,
    fill = "white",
    color = "grey85",
    linewidth = 0.4
  ) +
  annotate(
    "segment",
    x = 2012.7,
    xend = 2015.0,
    y = 301,
    yend = 301,
    color = "#1f77b4",
    linewidth = 0.8,
    linetype = "solid"
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 301,
    label = "Logistic fit",
    hjust = 0,
    size = 3
  ) +
  annotate(
    "segment",
    x = 2012.7,
    xend = 2015.0,
    y = 294,
    yend = 294,
    color = "#1f77b4",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 294,
    label = "Forecast",
    hjust = 0,
    size = 3
  ) +
  annotate(
    "point",
    x = 2013.85,
    y = 287,
    color = "#1f77b4",
    size = 2
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 287,
    label = "Observed",
    hjust = 0,
    size = 3
  ) +
  annotate(
    "segment",
    x = 2013.85,
    xend = 2013.85,
    y = 278,
    yend = 283,
    color = "red",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 280.5,
    label = peak_label,
    hjust = 0,
    size = 3
  ) +
  
  labs(
    title = NULL,
    subtitle = paste0(
      "R² = ", sprintf("%.3f", r2),
      " | Peak = ", peak_publications,
      " publications in ", peak_year_display
    ),
    x = "Year",
    y = "Publications (Annual)"
  ) +
  scale_x_continuous(
    breaks = seq(2010, 2050, 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 300, 50)
  ) +
  coord_cartesian(
    xlim = c(2010, 2055),
    ylim = c(0, 310)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ================================
# 6. Right plot
# ================================

p2 <- ggplot() +
  geom_line(
    data = fit_observed,
    aes(x = year, y = cumulative_fit),
    color = "#2ca02c",
    linewidth = 0.8,
    linetype = "solid"
  ) +
  geom_line(
    data = fit_forecast,
    aes(x = year, y = cumulative_fit),
    color = "#2ca02c",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(
    data = df,
    aes(x = year, y = cumulative),
    color = "#2ca02c",
    size = 2
  ) +
  geom_hline(
    yintercept = K * 0.5,
    color = "grey60",
    linetype = "dotted",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = K * 0.9,
    color = "grey60",
    linetype = "dotted",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2013,
    y = K * 0.5 + 80,
    label = "50%",
    color = "grey40",
    size = 3
  ) +
  annotate(
    "text",
    x = 2013,
    y = K * 0.9 + 80,
    label = "90%",
    color = "grey40",
    size = 3
  ) +
  
  # Manual legend
  annotate(
    "rect",
    xmin = 2011.5,
    xmax = 2027.5,
    ymin = 5050,
    ymax = 5500,
    fill = "white",
    color = "grey85",
    linewidth = 0.4
  ) +
  annotate(
    "segment",
    x = 2012.7,
    xend = 2015.0,
    y = 5400,
    yend = 5400,
    color = "#2ca02c",
    linewidth = 0.8,
    linetype = "solid"
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 5400,
    label = "Logistic fit",
    hjust = 0,
    size = 3
  ) +
  annotate(
    "segment",
    x = 2012.7,
    xend = 2015.0,
    y = 5280,
    yend = 5280,
    color = "#2ca02c",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 5280,
    label = "Forecast",
    hjust = 0,
    size = 3
  ) +
  annotate(
    "point",
    x = 2013.85,
    y = 5160,
    color = "#2ca02c",
    size = 2
  ) +
  annotate(
    "text",
    x = 2015.7,
    y = 5160,
    label = "Observed",
    hjust = 0,
    size = 3
  ) +
  
  labs(
    title = NULL,
    subtitle = paste0("Saturation (K) = ", K, " publications"),
    x = "Year",
    y = "Cumulative Publications"
  ) +
  scale_x_continuous(
    breaks = seq(2010, 2050, 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 5000, 1000)
  ) +
  coord_cartesian(
    xlim = c(2010, 2055),
    ylim = c(0, 5600)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ================================
# 7. Combine plots
# ================================

final_plot <- (p1 + p2) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 15,
      face = "bold"
    ),
    plot.tag.position = c(0.1, 0.98)
  )

# 8. Show plot
if (interactive()) print(final_plot)

# 10. Export left plot only
ggsave(
  filename = file.path(package_root, "Figures_TIFF", "Fig3C_Annual_Publication_Lifecycle.tiff"),
  plot = p1,
  width = 8.5,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(package_root, "Figures_PNG", "Fig3C_Annual_Publication_Lifecycle.png"),
  plot = p1,
  width = 8.5,
  height = 7,
  units = "in",
  dpi = 600
)

# 11. Export right plot only
ggsave(
  filename = file.path(package_root, "Figures_TIFF", "Fig3D_Cumulative_Publication_Lifecycle.tiff"),
  plot = p2,
  width = 8.5,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(package_root, "Figures_PNG", "Fig3D_Cumulative_Publication_Lifecycle.png"),
  plot = p2,
  width = 8.5,
  height = 7,
  units = "in",
  dpi = 600
)

# 12. Print output information
cat(
  "R² = ", sprintf("%.3f", r2),
  "\nPeak year = ", peak_year,
  "\nPeak year displayed = ", peak_year_display,
  "\nPeak publications = ", peak_publications,
  "\nSaturation K = ", K,
  "\nSaved in: ", package_root,
  "\n",
  sep = ""
)

# 13. Export observations and the descriptive model reconstruction.
df$cumulative_fit <- K / (1 + exp(-r * (df$year - peak_year)))
write.csv(
  df,
  file.path(package_root, "Source_Data", "Fig3CD_Observed_and_Model_Values.csv"),
  row.names = FALSE
)
write.csv(
  pred_df,
  file.path(package_root, "Source_Data", "Fig3CD_Model_Curves.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    parameter = c("K", "peak_year", "peak_publications", "r", "r2"),
    value = c(K, peak_year, peak_publications, r, r2)
  ),
  file.path(package_root, "Source_Data", "Fig3CD_Model_Statistics.csv"),
  row.names = FALSE
)
