# ================================
# Top 10 authors' publication output over time
# Manuscript-style visible-bubble version
# Bubble size = N articles
# Bubble color = TC per Year
# Output = 600 dpi PNG / TIFF
# ================================

# 0. Check and load required packages
if (!exists("package_root")) stop("Run run_all.R from the package root.")

packages <- c("readxl", "ggplot2", "grid", "ragg")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, ". See README.md.")
  }
}

library(readxl)
library(ggplot2)
library(grid)

# ================================
# 1. Read Excel file
# ================================

file_path <- file.path(package_root, "Input_Data", "Author_Production_Over_Time.xlsx")
df_raw <- readxl::read_excel(file_path, sheet = "Data")

cat("Original column names:\n")
print(names(df_raw))

# ================================
# 2. Standardize column names
# ================================

names(df_raw) <- trimws(names(df_raw))
lower_names <- tolower(names(df_raw))

author_col <- names(df_raw)[lower_names %in% c("author", "authors")]
year_col   <- names(df_raw)[lower_names %in% c("year", "years")]
freq_col   <- names(df_raw)[lower_names %in% c("freq", "n_article", "n_articles", "n articles", "articles")]
tcpy_col   <- names(df_raw)[lower_names %in% c("tcpy", "tc_per_year", "tc per year", "tc/year")]

if (length(author_col) == 0) stop("Cannot find Author column.")
if (length(year_col) == 0) stop("Cannot find year column.")
if (length(freq_col) == 0) stop("Cannot find freq / N_articles column.")
if (length(tcpy_col) == 0) stop("Cannot find TCpY / TC_per_year column.")

author_col <- author_col[1]
year_col   <- year_col[1]
freq_col   <- freq_col[1]
tcpy_col   <- tcpy_col[1]

df <- data.frame(
  Author = as.character(df_raw[[author_col]]),
  Year = as.numeric(df_raw[[year_col]]),
  N_articles = as.numeric(df_raw[[freq_col]]),
  TC_per_year = as.numeric(df_raw[[tcpy_col]])
)

stopifnot(
  nrow(df) == 72L,
  !anyNA(df),
  all(nzchar(trimws(df$Author))),
  length(unique(df$Author)) == 10L,
  all(is.finite(df$Year)),
  all(df$Year %in% 2012:2025),
  min(df$Year) == 2012,
  max(df$Year) == 2025,
  all(is.finite(df$N_articles)),
  all(df$N_articles > 0 & df$N_articles == floor(df$N_articles)),
  all(is.finite(df$TC_per_year)),
  all(df$TC_per_year >= 0),
  anyDuplicated(df[c("Author", "Year")]) == 0L
)

cat("\nRows used: ", nrow(df), "\n", sep = "")
cat("Authors: ", length(unique(df$Author)), "\n", sep = "")

# ================================
# 3. Order authors
# ================================

author_total <- aggregate(
  N_articles ~ Author,
  data = df,
  FUN = sum
)

author_total <- author_total[order(author_total$N_articles, decreasing = FALSE), ]

df$Author <- factor(
  df$Author,
  levels = author_total$Author
)

span_df <- aggregate(
  Year ~ Author,
  data = df,
  FUN = function(x) c(min(x), max(x))
)

span_df <- data.frame(
  Author = span_df$Author,
  start_year = span_df$Year[, 1],
  end_year = span_df$Year[, 2]
)

span_df$Author <- factor(
  span_df$Author,
  levels = levels(df$Author)
)

# ================================
# 4. TC per Year groups
# ================================

# Retain the original manuscript's numeric labels for the five color bins.
tc_labels <- c("0", "40", "80", "120", "160")
tc_intervals <- c("[0,20)", "[20,60)", "[60,100)", "[100,140)", "[140,infinity)")

df$TC_group <- cut(
  df$TC_per_year,
  breaks = c(-Inf, 20, 60, 100, 140, Inf),
  labels = tc_labels,
  right = FALSE
)

df$TC_group <- factor(
  df$TC_group,
  levels = tc_labels
)
df$TC_interval <- tc_intervals[as.integer(df$TC_group)]

tc_colors <- c(
  "0"   = "#8fb3d9",
  "40"  = "#6d9dcb",
  "80"  = "#4e86bb",
  "120" = "#2d68a5",
  "160" = "#123f73"
)

# Dummy layer only for complete TC legend
legend_dummy <- data.frame(
  Year = rep(2030, 5),
  Author = factor(rep(levels(df$Author)[1], 5), levels = levels(df$Author)),
  TC_group = factor(
    tc_labels,
    levels = tc_labels
  )
)

# ================================
# 5. Plot
# ================================

p_author <- ggplot() +
  
  geom_segment(
    data = span_df,
    aes(
      x = start_year,
      xend = end_year,
      y = Author,
      yend = Author
    ),
    color = "#d9aaaa",
    linewidth = 0.8,
    alpha = 0.58
  ) +
  
  # Main bubbles: solid circles, safest rendering
  geom_point(
    data = df,
    aes(
      x = Year,
      y = Author,
      size = N_articles,
      colour = TC_group
    ),
    shape = 16,
    alpha = 0.95
  ) +
  
  # Invisible dummy points to force full color legend
  geom_point(
    data = legend_dummy,
    aes(
      x = Year,
      y = Author,
      colour = TC_group
    ),
    shape = 16,
    size = 4.5,
    alpha = 0,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  
  scale_size_area(
    name = "N articles",
    max_size = 9,
    breaks = c(1, 2, 4, 6, 8)
  ) +
  
  scale_colour_manual(
    name = "TC per Year",
    values = tc_colors,
    limits = tc_labels,
    breaks = tc_labels,
    labels = tc_labels,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(2012, 2025, 2),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    xlim = c(2011.3, 2025.7),
    clip = "on"
  ) +
  
  labs(
    title = NULL,
    x = "Year",
    y = "Author"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12),
    
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      size = 10,
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(size = 10),
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.4
    ),
    panel.grid.major.x = element_line(
      color = "grey93",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(
      color = "grey45",
      linewidth = 0.5
    ),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", color = NA),
    legend.spacing.y = unit(0.25, "cm"),
    
    plot.margin = margin(8, 12, 8, 8)
  ) +
  
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        colour = "#2d68a5",
        alpha = 0.95
      )
    ),
    colour = guide_legend(
      order = 2,
      override.aes = list(
        shape = 16,
        size = 4.8,
        alpha = 1
      )
    )
  )

# ================================
# 6. Show plot
# ================================

if (interactive()) print(p_author)

# ================================
# 7. Export files
# ================================

write.csv(
  df,
  file.path(package_root, "Source_Data", "Fig6C_Plotted_Data.csv"),
  row.names = FALSE
)
write.csv(
  author_total,
  file.path(package_root, "Source_Data", "Fig6C_Author_Totals.csv"),
  row.names = FALSE
)

# Keep the supplied manuscript panel and its PNG preview at their original paths.
output_stem <- "Fig6C_Author_Production_Over_Time_Redrawn"
ragg::agg_png(
  filename = file.path(package_root, "Figures_PNG", paste0(output_stem, ".png")),
  width = 10.5,
  height = 5.6,
  units = "in",
  res = 600
)
print(p_author)
invisible(dev.off())

ragg::agg_tiff(
  filename = file.path(package_root, "Figures_TIFF", paste0(output_stem, ".tiff")),
  width = 10.5,
  height = 5.6,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(p_author)
invisible(dev.off())

cat(
  "\nSaved in: ",
  package_root,
  "\n",
  sep = ""
)

