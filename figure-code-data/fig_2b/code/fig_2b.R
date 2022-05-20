# load ggplot2
library(ggplot2)
library(PNWColors)

# figure directory
fig_dir <- "/path/to/fig_2b"

# Create Data
data <- data.frame(
  # not covered, PDB, SWISS-MODEL, AF2
  group = as.factor(c("19.7%", "21.7%", "24.5%", "34.1%")),
  value = c(0.197, 0.217, 0.245, 0.341)
)

colors <- c("gray", pnw_palette(name = "Shuksan2", n = 7)[c(2, 4, 6)])

# Basic pie chart
fig_2b <- ggplot(
  data, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", color = "white", size = 2) +
  scale_fill_manual(values = colors) +
  coord_polar("y", start = 0, direction = -1) +
  theme_void() + 
  theme(legend.position = "none")

# save the pie chart to disk
ggsave(
  filename = paste(fig_dir, "fig_2b.svg", sep = "/"),
  plot = fig_2b,
  width = 8,
  height = 8,
  units = "in",
  device = "svg",
)
