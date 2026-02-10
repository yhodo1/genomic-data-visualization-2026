library(ggplot2)

# ---- Load data ----
data <- read.csv("~/Desktop/Visium-IRI-ShamR_matrix.csv.gz")

# ---- Spatial coordinates ----
pos <- data[, c("x", "y")]

# ---- Gene expression: numeric only ----
gexp <- data[, sapply(data, is.numeric)]
gexp <- gexp[, !(colnames(gexp) %in% c("x", "y"))]

# ---- Remove constant/zero-variance genes ----
# (Genes that are all the same across spots, often all zeros)
keep_var <- apply(gexp, 2, sd, na.rm = TRUE) > 0
gexp <- gexp[, keep_var]

# Optional: also remove extremely low-signal genes (helps stability)
# keep_mean <- colMeans(gexp, na.rm = TRUE) > 0
# gexp <- gexp[, keep_mean]

# ---- PCA (let prcomp do scaling safely) ----
pca <- prcomp(gexp, center = TRUE, scale. = TRUE)

# ---- Data for plotting ----
df <- data.frame(
  x = pos$x,
  y = pos$y,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# ---- 2-panel spatial PCA visualization ----
p1 <- ggplot(df, aes(x, y, color = PC1)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial pattern of PC1")

p2 <- ggplot(df, aes(x, y, color = PC2)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial pattern of PC2")

library(patchwork)
p1 + p2


