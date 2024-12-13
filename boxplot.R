# Install required packages
install.packages("ggplot2")
install.packages("ggdist")

# Load the libraries
library(ggplot2)
library(ggdist)

# Example data
set.seed(123)
data <- data.frame(value = rnorm(100, mean = 1.18, sd = 0.5))

# Create the plot
ggplot(data, aes(x = 1, y = value)) +
  ggdist::stat_halfeye( # For KDE and violin-like plot
    adjust = 0.5, width = 0.6, justification = -0.2, .width = 0, point_color = NA
  ) +
  geom_boxplot(width = 0.2, outlier.shape = NA) + # Boxplot
  geom_jitter(width = 0.1, alpha = 0.4, size = 2) + # Points
  geom_text(aes(label = "median: 1.18"), x = 1, y = 1.18, vjust = -1.5) + # Annotation
  theme_minimal() +
  xlab("") +
  ylab("Value") +
  coord_flip()
