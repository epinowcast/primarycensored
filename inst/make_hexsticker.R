library(hexSticker)
library(sysfonts)
library(ggplot2)

# font setup
font_add_google("Zilla Slab Highlight", "useme")

# Create data for the plot
set.seed(123)
n_events <- 5
data <- data.frame(
  id = 1:n_events,
  primary = runif(n_events, 0, 10)
)

data$secondary <- data$primary + runif(n_events, 1, 5)
data$observed <- sample(
  c(TRUE, FALSE),
  n_events,
  replace = TRUE,
  prob = c(0.7, 0.3)
)

# Generate random error bar widths for each event
data$primary_error <- runif(n_events, 0.2, 0.8)
data$secondary_error <- runif(n_events, 0.2, 0.8)

# Determine which events cross the observation line
observation_time <- 10
data$crosses_observation <- data$primary <= observation_time & data$secondary > observation_time

# Filter out the 5th and 2nd events
data <- data[data$id != 5, ]

# Create the plot
plot <- ggplot(data, aes(y = id)) +
  # Add linking bar between primary and secondary events (behind other elements)
  geom_segment(data = subset(data, !crosses_observation), aes(x = primary, xend = secondary, yend = id), size = 0.6, color = "#708090") +
  geom_segment(data = subset(data, crosses_observation), aes(x = primary, xend = secondary, yend = id), size = 0.6, color = "#708090", linetype = "dotted") +
  # Add primary events
  geom_point(aes(x = primary), color = "#4682B4", size = 1) +
  # Add secondary events
  geom_point(aes(x = secondary), color = "#20B2AA", size = 1) +
  # Add uncertainty brackets with varying widths
  geom_errorbarh(aes(xmin = primary - primary_error, xmax = primary + primary_error), height = 0.6, color = "#4682B4", size = 0.9) +
  geom_errorbarh(aes(xmin = secondary - secondary_error, xmax = secondary + secondary_error), height = 0.6, color = "#20B2AA", size = 0.9) +
  # Add observation time line
  geom_vline(xintercept = observation_time, linetype = "dashed", color = "#999999", size = 0.3) +
  # Customize the plot
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0.5, n_events - 0.5)) + # Adjusted for removed event
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

# make and save hexsticker
sticker(
  plot,
  package = "primarycensoreddist",
  p_size = 12,
  p_color = "#4682B4",
  s_x = 1,
  s_y = 0.85,
  s_width = 1.3,
  s_height = 0.75,
  h_fill = "#F0F8FF",
  h_color = "#20B2AA",
  filename = file.path("man", "figures", "logo.png"),
  u_color = "#708090",
  u_size = 3.5
)
