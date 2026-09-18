library(ggplot2)
library(dplyr)
library(gifski)

set.seed(42)

n_groupes <- 20
n_frames  <- 120

X     <- seq(0, 7, length.out = n_groupes)
Poids <- sample(5:10, n_groupes, replace = TRUE)
X_mid <- mean(X)

beta1_amp <- 1.2   # amplitude du changement de pente

# décalage de phase et amplitude d'oscillation propres à chaque point
phases    <- runif(n_groupes, 0, 2 * pi)
amplitude <- runif(n_groupes, 0.08, 0.15)

# dispersion individuelle fixe (certains points restent plus loin de la courbe)
offset_indiv <- rnorm(n_groupes, mean = 0, sd = 0.10)

dir_out <- "/Users/jmiron/Library/CloudStorage/OneDrive-UniversitéLaval/Quarto/STT-4300/images/frames"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(n_frames)) {

  t <- (i - 1) / n_frames * 2 * pi

  beta1_t <- beta1_amp * cos(t)
  eta_t   <- beta1_t * (X - X_mid)
  p_t     <- plogis(eta_t)

  jitter <- amplitude * sin(t + phases)
  Y_t    <- pmin(pmax(p_t + jitter + offset_indiv, 0.01), 0.99)

  df <- data.frame(X = X, Y = Y_t, Poids = Poids)

  modele <- glm(Y ~ X, data = df, family = binomial(link = "logit"), weights = Poids)

  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille, type = "response")

  p <- ggplot() +
    geom_point(
      data = df,
      aes(x = X, y = Y, size = Poids),
      color = "mediumpurple1",
      alpha = 0.6
    ) +
    geom_line(
      data = grille,
      aes(x = X, y = Y_hat),
      color = "palegreen4",
      linewidth = 1
    ) +
    coord_cartesian(xlim = range(X), ylim = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none")

  ggsave(
    filename = file.path(dir_out, sprintf("frame_%03d.png", i)),
    plot = p, width = 8, height = 5.5, dpi = 300, bg = "transparent"
  )
}

# assemblage en gif
png_files <- sprintf(file.path(dir_out, "frame_%03d.png"), 1:n_frames)

gifski(
  png_files,
  gif_file = file.path(dirname(dir_out), "logit_plot.gif"),
  width = 800, height = 550, delay = 0.08, loop = TRUE
)