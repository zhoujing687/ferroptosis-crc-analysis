# ================================
# Temporal trends in cumulative publications
# Top 10 journals, 2012–2025
# Legend inside top-left
# ================================

if (!exists("package_root")) stop("Run run_all.R from the package root.")
library(ggplot2)
library(grid)
library(readxl)

# 1. Input data
input_file <- file.path(package_root, "Input_Data", "Journal_Production_Over_Time.xlsx")
source_data <- as.data.frame(read_excel(input_file, sheet = "Data"))
journal_names <- c(
  "CANCERS",
  "INTERNATIONAL JOURNAL OF MOLECULAR SCIENCES",
  "FRONTIERS IN PHARMACOLOGY",
  "FRONTIERS IN ONCOLOGY",
  "PHYTOMEDICINE",
  "SCIENTIFIC REPORTS",
  "PLOS ONE",
  "CANCER LETTERS",
  "LIFE SCIENCES",
  "ONCOTARGET"
)
stopifnot(
  identical(names(source_data), c("Year", journal_names)),
  nrow(source_data) == 14L,
  all(source_data$Year == 2012:2025),
  all(vapply(source_data[journal_names], is.numeric, logical(1)))
)
counts <- as.matrix(source_data[journal_names])
stopifnot(
  all(is.finite(counts)),
  all(counts >= 0 & counts == floor(counts)),
  all(apply(counts, 2, function(x) all(diff(x) >= 0)))
)
# The export already contains cumulative counts; do not accumulate them again.
data_long <- data.frame(
  year = rep(source_data$Year, times = length(journal_names)),
  journal = rep(journal_names, each = nrow(source_data)),
  cumulative_publications = as.vector(counts)
)

# 2. Journal order
data_long$journal <- factor(
  data_long$journal,
  levels = journal_names
)
write.csv(data_long,
          file.path(package_root, "Source_Data", "Fig5D_Plotted_Data.csv"),
          row.names = FALSE)

# 3. Colors
journal_colors <- c(
  "CANCERS" = "#1f77b4",
  "INTERNATIONAL JOURNAL OF MOLECULAR SCIENCES" = "#ff7f0e",
  "FRONTIERS IN PHARMACOLOGY" = "#009E73",
  "FRONTIERS IN ONCOLOGY" = "#CC79A7",
  "PHYTOMEDICINE" = "#E6C300",
  "SCIENTIFIC REPORTS" = "#56B4E9",
  "PLOS ONE" = "#000000",
  "CANCER LETTERS" = "#7F7F7F",
  "LIFE SCIENCES" = "#999999",
  "ONCOTARGET" = "#A6761D"
)

# 4. Plot
p_journals <- ggplot(
  data_long,
  aes(
    x = year,
    y = cumulative_publications,
    color = journal,
    group = journal
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.7) +
  
  # Value labels
  geom_text(
    aes(label = cumulative_publications),
    size = 2.3,
    vjust = -0.6,
    show.legend = FALSE
  ) +
  
  labs(
    title = NULL,
    x = "Year",
    y = "Cumulative number of publications",
    color = NULL
  ) +
  
  scale_x_continuous(
    breaks = 2012:2025,
    labels = 2012:2025
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    limits = c(0, 28),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  
  scale_color_manual(values = journal_colors) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8.5
    ),
    axis.text.y = element_text(size = 9),
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.major.x = element_line(
      color = "grey92",
      linewidth = 0.25
    ),
    
    # Legend inside top-left
    legend.position = c(0.03, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_rect(
      fill = "white",
      color = "grey80",
      linewidth = 0.3
    ),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 7),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    
    plot.margin = margin(8, 10, 8, 10)
  ) +
  
  guides(
    color = guide_legend(
      ncol = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 0.8, size = 2)
    )
  )

# 5. Show plot
if (interactive()) print(p_journals)

# 6. Save as 600 dpi PNG
ggsave(
  filename = file.path(package_root, "Figures_PNG", "Fig5D_Journal_Production_Over_Time.png"),
  plot = p_journals,
  width = 11,
  height = 6,
  units = "in",
  dpi = 600
)

# 7. Save as 600 dpi TIFF
ggsave(
  filename = file.path(package_root, "Figures_TIFF", "Fig5D_Journal_Production_Over_Time.tiff"),
  plot = p_journals,
  width = 11,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# 8. Check file location
cat("Saved Fig 5D PNG and TIFF in the package figure directories.\n")
