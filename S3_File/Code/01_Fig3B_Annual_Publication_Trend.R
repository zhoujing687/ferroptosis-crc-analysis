# ================================
# Annual Publication Changes
# Final version with one combined legend box
# ================================

if (!exists("package_root")) stop("Run run_all.R from the package root.")

library(ggplot2)

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

# 2. Use x = 0, 1, 2 ... for exponential fitting
df$x <- 0:(nrow(df) - 1)

# 3. Exponential fit: y = a * exp(bx)
fit <- nls(
  publications ~ a * exp(b * x),
  data = df,
  start = list(a = 14, b = 0.17)
)

a <- coef(fit)[["a"]]
b <- coef(fit)[["b"]]

# 4. Generate fitted curve
pred_df <- data.frame(
  x = seq(min(df$x), max(df$x), length.out = 300)
)

pred_df$year <- 2012 + pred_df$x
pred_df$fit <- a * exp(b * pred_df$x)

# 5. Calculate actual R²
df$pred <- predict(fit)

ss_res <- sum((df$publications - df$pred)^2)
ss_tot <- sum((df$publications - mean(df$publications))^2)
r2 <- 1 - ss_res / ss_tot

# 6. Bar colors, close to your original figure
bar_colors <- c(
  "#2c3e63",  # 2012
  "#5f8ca8",  # 2013
  "#a8c9d8",  # 2014
  "#c8ded7",  # 2015
  "#d9dc28",  # 2016
  "#8cc63f",  # 2017
  "#47b5a2",  # 2018
  "#f3a20a",  # 2019
  "#f3701a",  # 2020
  "#ef4c3e",  # 2021
  "#a83244",  # 2022
  "#62006f",  # 2023
  "#7b1fa2",  # 2024
  "#24406e"   # 2025
)

# 7. Legend formula label
fit_label_expr <- paste0(
  "'Exponential Fit: y = ",
  sprintf("%.4f", a),
  "'*italic(e)^{'",
  sprintf("%.4f", b),
  "x'}"
)

# 8. Plot
p <- ggplot(df, aes(x = year, y = publications)) +
  
  # Bars
  geom_col(
    aes(fill = factor(year)),
    width = 0.75,
    show.legend = FALSE
  ) +
  
  # Exponential fit line
  geom_line(
    data = pred_df,
    aes(x = year, y = fit),
    inherit.aes = FALSE,
    color = "orange",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # Value labels
  geom_text(
    aes(label = publications),
    vjust = -0.35,
    size = 3.0,
    color = "black"
  ) +
  
  # R² label, actual calculated value
  annotate(
    "text",
    x = 2022.0,
    y = 115,
    label = paste0("italic(R)^2 == ", sprintf("%.4f", r2)),
    parse = TRUE,
    color = "orange",
    size = 4
  ) +
  
  # One combined legend background box
  annotate(
    "rect",
    xmin = 2011.55,
    xmax = 2016.15,
    ymin = 124.5,
    ymax = 138.5,
    fill = "white",
    color = "grey80",
    linewidth = 0.4
  ) +
  
  # Legend: dashed line sample
  annotate(
    "segment",
    x = 2011.75,
    xend = 2012.18,
    y = 135.2,
    yend = 135.2,
    color = "orange",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # Legend: exponential fit text
  annotate(
    "text",
    x = 2012.35,
    y = 135.2,
    label = fit_label_expr,
    parse = TRUE,
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "black"
  ) +
  
  # Legend: publications square
  annotate(
    "rect",
    xmin = 2011.75,
    xmax = 2012.18,
    ymin = 127.4,
    ymax = 131.0,
    fill = "#2c3e63",
    color = NA
  ) +
  
  # Legend: publications text
  annotate(
    "text",
    x = 2012.35,
    y = 129.2,
    label = "Publications",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "black"
  ) +
  
  # Titles and axis labels
  labs(
    title = NULL,
    x = "Year",
    y = "Number of Publications"
  ) +
  
  # X axis
  scale_x_continuous(
    breaks = df$year,
    labels = df$year,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, 140),
    breaks = seq(0, 140, 20),
    expand = c(0, 0)
  ) +
  
  # Bar colors
  scale_fill_manual(values = bar_colors) +
  
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "plain"
    ),
    
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 12),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 9),
    
    panel.grid.major.y = element_line(
      color = "grey70",
      linetype = "dashed",
      linewidth = 0.35
    ),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.6
    ),
    
    axis.line = element_blank(),
    legend.position = "none"
  )

# 9. Show plot
if (interactive()) print(p)

# 10. Save plot
# Save plot as 600 dpi PNG
ggsave(
  filename = file.path(package_root, "Figures_PNG", "Fig3B_Annual_Publication_Trend.png"),
  plot = p,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600
)

# Save plot as 600 dpi TIFF for journal submission
ggsave(
  filename = file.path(package_root, "Figures_TIFF", "Fig3B_Annual_Publication_Trend.tiff"),
  plot = p,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# 11. Print fitting results
cat(
  "Fitted equation: y = ",
  sprintf("%.4f", a),
  "e^",
  sprintf("%.4f", b),
  "x",
  "\nR-squared = ",
  sprintf("%.4f", r2),
  "\n",
  sep = ""
)

# 12. Export the data and model statistics used in the figure.
write.csv(
  df,
  file.path(package_root, "Source_Data", "Fig3B_Observed_and_Fitted.csv"),
  row.names = FALSE
)
write.csv(
  pred_df,
  file.path(package_root, "Source_Data", "Fig3B_Fitted_Curve.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(parameter = c("a", "b", "r2"), value = c(a, b, r2)),
  file.path(package_root, "Source_Data", "Fig3B_Model_Statistics.csv"),
  row.names = FALSE
)
