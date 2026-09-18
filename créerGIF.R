# ================================================================
# Paramètres communs
# ================================================================
library(ggplot2)
library(gifski)
library(mgcv)   # requis seulement pour le chapitre GAM

base_dir    <- "/Users/jmiron/Library/CloudStorage/OneDrive-UniversitéLaval/GitHub/STT-4300/images"
n_frames    <- 120
largeur_img <- 8
hauteur_img <- 5.5
dpi_img     <- 300
largeur_gif <- 800
hauteur_gif <- 550
delai_gif   <- 0.08

couleur_pts   <- "mediumpurple1"
couleur_ligne <- "palegreen4"

theme_gif <- theme_void() + theme(legend.position = "none")

# ================================================================
# Fonctions utilitaires
# ================================================================

# Crée (si besoin) le dossier de frames d'un chapitre et retourne son chemin
preparer_dossier <- function(nom_dossier) {
  dir_out <- file.path(base_dir, nom_dossier)
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  dir_out
}

# Sauvegarde une frame i dans dir_out avec les dimensions/dpi communs
sauvegarder_frame <- function(plot, dir_out, i) {
  ggsave(
    file.path(dir_out, sprintf("frame_%03d.png", i)),
    plot, width = largeur_img, height = hauteur_img,
    dpi = dpi_img, bg = "transparent"
  )
}

# Assemble les n_frames_gif frames d'un dossier en gif, sauvegardé dans base_dir
creer_gif <- function(dir_out, nom_gif, n_frames_gif = n_frames) {
  png_files <- sprintf(file.path(dir_out, "frame_%03d.png"), seq_len(n_frames_gif))
  gifski(
    png_files, gif_file = file.path(base_dir, nom_gif),
    width = largeur_gif, height = hauteur_gif,
    delay = delai_gif, loop = TRUE
  )
}


# ================================================================
# 02 — Régression logistique
# ================================================================
set.seed(42)
dir_out <- preparer_dossier("frames")

n_groupes <- 20
X     <- seq(0, 7, length.out = n_groupes)
Poids <- sample(5:10, n_groupes, replace = TRUE)
X_mid <- mean(X)
beta1_amp <- 1.2

phases       <- runif(n_groupes, 0, 2 * pi)
amplitude    <- runif(n_groupes, 0.08, 0.15)
offset_indiv <- rnorm(n_groupes, 0, 0.10)

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  beta1_t <- beta1_amp * cos(t)
  eta_t   <- beta1_t * (X - X_mid)
  p_t     <- plogis(eta_t)
  jitter  <- amplitude * sin(t + phases)
  Y_t     <- pmin(pmax(p_t + jitter + offset_indiv, 0.01), 0.99)

  df <- data.frame(X = X, Y = Y_t, Poids = Poids)
  modele <- glm(Y ~ X, data = df, family = binomial(link = "logit"), weights = Poids)
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille, type = "response")

  p <- ggplot() +
    geom_point(data = df, aes(X, Y, size = Poids), color = couleur_pts, alpha = 0.6) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = 1) +
    coord_cartesian(xlim = range(X), ylim = c(0, 1)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "logit_plot.gif")


# ================================================================
# 01 — Régression linéaire
# ================================================================
set.seed(1)
dir_out <- preparer_dossier("frames_lineaire")

n_pts <- 30
X     <- runif(n_pts, 0, 10)
X_mid <- mean(X)
beta1_amp <- 2

phases     <- runif(n_pts, 0, 2 * pi)
amplitude  <- runif(n_pts, 0.3, 0.6)
bruit_fixe <- rnorm(n_pts, 0, 0.5)

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  beta1_t <- beta1_amp * cos(t)

  # bruit perpendiculaire constant, converti en bruit vertical pour compenser la pente
  bruit_perp <- amplitude * sin(t + phases) + bruit_fixe
  facteur    <- sqrt(1 + beta1_t^2)
  Y_t        <- beta1_t * (X - X_mid) + bruit_perp * facteur

  df <- data.frame(X = X, Y = Y_t)
  modele <- lm(Y ~ X, data = df)
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = 0.6, size = 2.5) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = 1) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "lineaire_plot.gif")


# ================================================================
# 03 — Modèles linéaires généralisés (Poisson)
# ================================================================
set.seed(2)
dir_out <- preparer_dossier("frames_poisson")

n_pts <- 25
X     <- seq(0, 6, length.out = n_pts)
X_mid <- mean(X)
beta1_amp <- 0.5

phases_eta <- runif(n_pts, 0, 2 * pi)
amplitude  <- runif(n_pts, 0.05, 0.15)
phases_z1  <- runif(n_pts, 0, 2 * pi)
phases_z2  <- runif(n_pts, 0, 2 * pi)
freq_z2    <- runif(n_pts, 1.5, 2.5)

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  beta1_t  <- beta1_amp * cos(t)
  eta_t    <- 1.5 + beta1_t * (X - X_mid) + amplitude * sin(t + phases_eta)
  lambda_t <- exp(eta_t)

  # signal individuel lisse (variance ~1), converti en bruit d'échelle sqrt(lambda_t)
  z_t <- (sin(t + phases_z1) + 0.5 * sin(freq_z2 * t + phases_z2)) / sqrt(1.25)
  Y_t <- pmax(lambda_t + sqrt(lambda_t) * z_t, 0)

  df <- data.frame(X = X, Y = Y_t)
  modele <- glm(Y ~ X, data = df, family = poisson(link = "log"))
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille, type = "response")

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = 0.6, size = 2.5) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = 1) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "glm_plot.gif")


# ================================================================
# 04 — Prédiction et validation croisée
# ================================================================
set.seed(3)
dir_out <- preparer_dossier("frames_vc")

n_pts       <- 150
k           <- 5
n_frames_vc <- k * 24   # 24 frames de pause par pli

X     <- runif(n_pts, 0, 10)
Y     <- 0.8 * X + rnorm(n_pts, 0, 1.5)
folds <- sample(rep(1:k, length.out = n_pts))

for (i in seq_len(n_frames_vc)) {
  pli_actif <- ((i - 1) %/% 24) %% k + 1

  df <- data.frame(X = X, Y = Y, fold = folds,
                    role = ifelse(folds == pli_actif, "test", "train"))

  modele <- lm(Y ~ X, data = subset(df, role == "train"))
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y, color = role), alpha = 0.7, size = 2.8) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = 1) +
    scale_color_manual(values = c(train = couleur_pts, test = "tomato")) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "cv_plot.gif", n_frames_gif = n_frames_vc)


# ================================================================
# 05 — Sélection et régularisation (LASSO)
# ================================================================
set.seed(4)
dir_out <- preparer_dossier("frames_lasso")

n_coef    <- 5
coef_base <- runif(n_coef, -2, 2)
noms      <- factor(paste0("b", 1:n_coef), levels = paste0("b", 1:n_coef))
labels_coef <- paste0("abs(hat(beta)[", 1:n_coef, "])")   # |beta_j chapeau|, syntaxe plotmath

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  lambda_t <- (cos(t) + 1) / 2   # oscille entre 0 et 1
  coef_t   <- sign(coef_base) * pmax(abs(coef_base) - lambda_t * 1.8, 0)

  df <- data.frame(coef = noms, valeur = abs(coef_t), label = labels_coef)

  p <- ggplot(df, aes(coef, valeur)) +
    geom_col(fill = couleur_pts, alpha = 0.75) +
    geom_text(aes(label = label), parse = TRUE, vjust = -0.4,
              color = couleur_ligne, size = 16, fontface = "bold") +
    geom_hline(yintercept = 0, color = couleur_ligne, linewidth = 1) +
    coord_cartesian(ylim = c(0, 2.4)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "regularisation_plot.gif")


# ================================================================
# 06 — KNN et noyau
# ================================================================
set.seed(5)
dir_out <- preparer_dossier("frames_knn")

n_pts <- 100
X <- runif(n_pts, 0, 10)
Y <- sin(X) + rnorm(n_pts, 0, 0.25)

for (i in seq_len(n_frames)) {
  t  <- (i - 1) / n_frames * 2 * pi
  bw <- 1.1 + 0.95 * cos(t)   # largeur de bande oscillant

  lissage <- ksmooth(X, Y, kernel = "normal", bandwidth = bw,
                      n.points = 200, x.points = seq(0, 10, length.out = 200))

  df_pts <- data.frame(X = X, Y = Y)
  df_l   <- data.frame(X = lissage$x, Y = lissage$y)

  p <- ggplot() +
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts, alpha = 0.6, size = 2.5) +
    geom_line(data = df_l, aes(X, Y), color = couleur_ligne, linewidth = 1) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "knn_noyau_plot.gif")


# ================================================================
# 07 — Splines
# ================================================================
set.seed(6)
dir_out <- preparer_dossier("frames_splines")

n_pts <- 40
X <- sort(runif(n_pts, 0, 10))
Y <- sin(X) + rnorm(n_pts, 0, 0.25)

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  df_t <- max(8.5 + 6.5 * cos(t), 2)   # pas de round() : évite les sauts entre frames

  modele <- smooth.spline(X, Y, df = df_t)
  grille <- seq(min(X), max(X), length.out = 200)
  Y_hat  <- predict(modele, grille)$y

  df_pts <- data.frame(X = X, Y = Y)
  df_l   <- data.frame(X = grille, Y = Y_hat)

  p <- ggplot() +
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts, alpha = 0.6, size = 2.5) +
    geom_line(data = df_l, aes(X, Y), color = couleur_ligne, linewidth = 1) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "splines_plot.gif")


# ================================================================
# 08 — GAM
# ================================================================
set.seed(7)
dir_out <- preparer_dossier("frames_gam")   # corrigé : était "frames_splines" par erreur

n_pts <- 60
X <- sort(runif(n_pts, 0, 10))
bruit_fixe <- rnorm(n_pts, 0, 0.3)   # fixé une fois, hors boucle

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  amp_t <- 1 + 0.8 * cos(t)

  Y <- amp_t * sin(X) + 0.3 * X + bruit_fixe

  df <- data.frame(X = X, Y = Y)
  modele <- gam(Y ~ s(X, k = 10, fx = TRUE), data = df)   # fx = TRUE : degrés de liberté fixes
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = 0.6, size = 2.5) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = 1) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "gam_plot.gif")