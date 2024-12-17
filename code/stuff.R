# Define distance bands
distances <- seq(0, 500000, by = 10000)

# Compute Moran's I for each distance band
correlogram <- sp.correlogram(
  neighbours = nb,   # Neighbor list object
  var = data$BDsample_0,  # Variable to analyze
  order = length(distances),  # Number of distance bands
  style = "W",  # Standardized weights
  method = "I",  # Moran's I
  zero.policy = TRUE
)

print(correlogram)

correlogram_df <- data.frame(
  distance = distances,
  morans_I = unlist(correlogram$res)
)

# Plot the correlogram
library(ggplot2)
ggplot(correlogram_df, aes(x = distance, y = morans_I)) +
  geom_line(color = "blue") +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Spatial Correlogram of Bulk Density",
    x = "Distance",
    y = "Moran's I"
  ) +
  theme_minimal()