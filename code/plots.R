##### PLOTS #####
# plot 
ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = lucas_2018, aes(colour = '2018'), size = 0.5, alpha = 0.7) + # LUCAS 2018 points
  geom_sf(data = lucas_2015, aes(colour = '2015'), size = 0.5, alpha = 0.7) + # LUCAS 2018 points
  #geom_sf(data = lucas_2018_common, aes(colour = '2015&2018'), size = 0.5, alpha = 0.9) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS Soil 2018 Points on Europe Map")

# Plot for 2015 point
plot_2015 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2015, aes(color = "2015"), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS 2015", color = "Year")

# Plot for 2018 point
plot_2018 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2018, aes(color = "2018"), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS 2018", color = "Year")

# plot common points
combined_plot <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2018_common, aes(color = '2015&2018'), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS Combined 2015 & 2018", color = "Year")

# Arrange in a three-column layout
(plot_2015 | plot_2018 | combined_plot)