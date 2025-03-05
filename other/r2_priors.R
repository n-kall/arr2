library(ggplot2)
library(dplyr)
library(directlabels)

papertheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
  )


# Parameters for Beta distributions with labels
beta_params <- data.frame(
  mean = c(0.5, 0.5, 1/3),
  precision = c(2, 1, 3),
  label = c("beta(0.5, 2)", "beta(0.5, 1)", "beta(0.33, 3)")
)

# Function to calculate Beta density
generate_beta_data <- function(mean, precision, label, x_seq) {
  alpha <- mean * precision
  beta <- (1 - mean) * precision
  data.frame(
    x = x_seq,
    density = dbeta(x_seq, shape1 = alpha, shape2 = beta),
    label = label
  )
}

# Generate data for plotting
x_seq <- seq(0, 1, length.out = 100)
beta_data <- beta_params %>%
  rowwise() %>%
  do(generate_beta_data(.$mean, .$precision, .$label, x_seq)) %>%
  ungroup()

new_colours <- c("beta(0.5, 2)" = "#999999",  # Neutral grey
                 "beta(0.5, 1)" = "#D55E00",  # Orange-red
                                  "beta(0.33, 3)" = "#56B4E9")  # Sky blue


# Plot with ggplot2
p <- ggplot(beta_data, aes(x = x, y = density, color = label)) +
  geom_line(size = 1) +
  labs(
    x = "$R^2$"
  ) +
  guides(colour = "none") +
  papertheme +
  scale_colour_manual(values = new_colours) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())




p <- direct.label(p, method = list(dl.trans(y = y + 0.1, x = x + 0.8), "smart.grid", cex = 0.8))

bayesflow::save_tikz_plot(p, "r2_plot.tex", width = 4, height = 2)
