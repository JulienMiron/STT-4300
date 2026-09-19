# ================================================================
# Paramètres communs
# ================================================================
library(ggplot2)
library(gifski)
library(mgcv)   # requis seulement pour le chapitre GAM

base_dir    <- "/Users/jmiron/Library/CloudStorage/OneDrive-UniversitéLaval/GitHub/STT-4300/images"
n_frames    <- 150
largeur_img <- 8
hauteur_img <- 5.5
dpi_img     <- 300
largeur_gif <- 800
hauteur_gif <- 550
delai_gif   <- 0.08
n_pts       <- 75 
n_groupes   <- 75 # pour la régression linéaire logistique
n_coef      <- 5 # pour le nombre de paramètres de LASSO



couleur_pts   <- "#9467bd"
couleur_ligne <- "#0E9E6C"
couleur_accent <- "#d2b827"   # 3e couleur, utilisée uniquement pour l'état final de l'index (3 classes)

alpha_pts       <- 1   # transparence des points (hors logit pondéré, voir plus bas)
taille_pts      <- 4   # taille des points (hors logit, dont la taille encode le poids)
epaisseur_ligne <- 1.2     # épaisseur des courbes/lignes ajustées

theme_gif <- theme_void() + theme(legend.position = "none")

# ================================================================
# Fonctions utilitaires
# ================================================================

# Crée (si besoin) le dossier de frames d'un chapitre et retourne son chemin
preparer_dossier <- function(nom_dossier) {
  dir_out <- file.path(base_dir, nom_dossier)
  if (dir.exists(dir_out)) {
    unlink(dir_out, recursive = TRUE)
  }
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
    geom_point(data = df, aes(X, Y, size = Poids), color = couleur_pts, alpha = alpha_pts) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = epaisseur_ligne) +
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
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "lineaire_plot.gif")


# ================================================================
# 03 — Modèles linéaires généralisés (Poisson)
# ================================================================
set.seed(2)
dir_out <- preparer_dossier("frames_poisson")

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
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "glm_plot.gif")


# ================================================================
# 04 — Prédiction et validation croisée
# ================================================================
set.seed(3)
dir_out <- preparer_dossier("frames_vc")

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
    geom_point(data = df, aes(X, Y, color = role), alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = epaisseur_ligne) +
    scale_color_manual(values = c(train = couleur_pts, test = couleur_accent)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "cv_plot.gif", n_frames_gif = n_frames_vc)


# ================================================================
# 05 — Sélection et régularisation (LASSO)
# ================================================================
set.seed(4)
dir_out <- preparer_dossier("frames_lasso")

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
    geom_hline(yintercept = 0, color = couleur_ligne, linewidth = epaisseur_ligne) +
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
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts, alpha = alpha_pts, size = taille_pts) +
    geom_line(data = df_l, aes(X, Y), color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "knn_noyau_plot.gif")


# ================================================================
# 07 — Splines
# ================================================================
set.seed(6)
dir_out <- preparer_dossier("frames_splines")

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
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts, alpha = alpha_pts, size = taille_pts) +
    geom_line(data = df_l, aes(X, Y), color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "splines_plot.gif")


# ================================================================
# 08 — GAM
# ================================================================
set.seed(7)
dir_out <- preparer_dossier("frames_gam")   # corrigé : était "frames_splines" par erreur

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
    geom_point(data = df, aes(X, Y), color = couleur_pts, alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat), color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "gam_plot.gif")


# ================================================================
# Index — nuage de points qui se réorganise (sans ajustement)
# ================================================================
set.seed(9)
dir_out <- preparer_dossier("frames_index")

 
X <- sort(runif(n_pts, 0, 10))
 
# chaque état est un triplet (x, y, couleur) par point, dans le même ordre pour
# tous les états : c'est ce qui permet d'interpoler à la fois la position et la
# couleur de façon continue d'un état à l'autre.
bruit_indiv <- rnorm(n_pts, 0, 1)

x_final <- runif(n_pts, 0, 10)
y_final <- runif(n_pts, -5, 5)
 
etat_aleatoire <- list(x = x_final, y = y_final,
                        col = rep(couleur_pts, n_pts))
etat_lineaire  <- list(x = X, y = 0.9 * (X - mean(X)) + bruit_indiv * 0.4,
                        col = rep(couleur_pts, n_pts))
etat_sinus     <- list(x = X, y = 3 * sin(0.8 * X) + bruit_indiv * 0.35,
                        col = rep(couleur_pts, n_pts))
 
# état final : points aléatoires sur le rectangle, répartis en 3 secteurs angulaires
# (3 vecteurs partant du centre, séparés de 120°, délimitent les classes)

 
centre_x <- 5   # centre du rectangle (domaine x : 0 à 10)
centre_y <- 0   # centre du rectangle (domaine y : -5 à 5)
 
angle          <- atan2(y_final - centre_y, x_final - centre_x)   # entre -pi et pi
angle_positif  <- ifelse(angle < 0, angle + 2 * pi, angle)         # entre 0 et 2*pi
classe_finale  <- floor(angle_positif / (2 * pi / 3)) + 1          # 3 secteurs de 120°
 
couleurs_classes <- c(couleur_pts, couleur_ligne, couleur_accent)
 
etat_grille <- list(x = x_final, y = y_final, col = couleurs_classes[classe_finale])
 
etats   <- list(etat_aleatoire, etat_lineaire, etat_sinus, etat_grille)
n_etats <- length(etats)
 
# angles des 3 vecteurs de séparation (mêmes bornes que celles utilisées pour classe_finale)
angles_separation <- c(0, 2 * pi / 3, 4 * pi / 3)
rayon_separation   <- 6   # longueur des vecteurs une fois complètement affichés
 
ease <- function(frac) (1 - cos(pi * frac)) / 2   # lissage ease-in-out entre deux états
 
# interpolation linéaire de couleurs point par point, dans l'espace RGB
interp_couleur <- function(col_a, col_b, frac) {
  rgb_a <- col2rgb(col_a)
  rgb_b <- col2rgb(col_b)
  rgb_t <- (1 - frac) * rgb_a + frac * rgb_b
  rgb(rgb_t[1, ], rgb_t[2, ], rgb_t[3, ], maxColorValue = 255)
}
 
for (i in seq_len(n_frames_ix)) {
  pos <- (i - 1) / n_frames_ix * n_etats
 
  idx_a <- (floor(pos) %% n_etats) + 1
  idx_b <- (floor(pos) + 1) %% n_etats + 1
  frac  <- ease(pos - floor(pos))
 
  X_t   <- (1 - frac) * etats[[idx_a]]$x   + frac * etats[[idx_b]]$x
  Y_t   <- (1 - frac) * etats[[idx_a]]$y   + frac * etats[[idx_b]]$y
  Col_t <- interp_couleur(etats[[idx_a]]$col, etats[[idx_b]]$col, frac)
 
  # les vecteurs n'apparaissent qu'en entrant dans l'état "grille" (croissance, poids 0→1)
  # ou en le quittant (rétraction, poids 1→0) ; nuls lors des autres transitions
  if (idx_b == n_etats) {
    poids_separation <- frac
  } else if (idx_a == n_etats) {
    poids_separation <- 1 - frac
  } else {
    poids_separation <- 0
  }
 
  df <- data.frame(X = X_t, Y = Y_t, Couleur = Col_t)
 
  df_separation <- data.frame(
    x    = centre_x,
    y    = centre_y,
    xend = centre_x + poids_separation * rayon_separation * cos(angles_separation),
    yend = centre_y + poids_separation * rayon_separation * sin(angles_separation)
  )
 
    # couleur de la flèche : blanche (invisible) quand poids_separation = 0,
  # grise et pleinement visible quand poids_separation = 1
  couleur_separation <- interp_couleur("white", "gray50", poids_separation)
 
  p <- ggplot() +
    geom_segment(data = df_separation, aes(x, y, xend = xend, yend = yend),
                 color = couleur_separation, linewidth = epaisseur_ligne * 0.6,
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    geom_point(data = df, aes(X, Y, color = Couleur), alpha = alpha_pts, size = taille_pts) +
    scale_color_identity() +
    coord_cartesian(xlim = c(0, 10), ylim = c(-5, 5)) +
    theme_gif
 
  sauvegarder_frame(p, dir_out, i)
}
 
creer_gif(dir_out, "index_plot.gif", n_frames_gif = n_frames_ix)