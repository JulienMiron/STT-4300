# ================================================================
# Paramètres communs
# ================================================================
library(ggplot2)
library(gifski)
library(mgcv)   # requis seulement pour le chapitre GAM

base_dir    <- "/Users/jmiron/Documents/GitHub/STT-4300/notes/images"
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

changer <- c(10) 
                #  chapitres à générer 
                #  0 = index, 
                #  1 = régression linéaire, 
                #  2 = Monte Carlo, 
                #  3 = régression logistique, 
                #  4 = GLM Poisson, 
                #  5 = validation croisée, 
                #  6 = LASSO, 
                #  7 = KNN/noyau, 
                #  8 = splines, 
                #  9 = GAM, 
                #  10 = exercices)



couleur_pts    <- "#9467bd"
couleur_ligne  <- "#0E9E6C"
# 3e couleur, utilisée seulement pour l'état final de l'index (3 classes)
couleur_accent <- "tomato"

# transparence des points (hors logit pondéré, voir plus bas)
alpha_pts <- 1
# taille des points (hors logit, dont la taille encode le poids)
taille_pts <- 4
# épaisseur des courbes/lignes ajustées
epaisseur_ligne <- 1.2

theme_gif <- theme_void() + theme(legend.position = "none")

# ================================================================
# Fonctions utilitaires
# ================================================================

# lissage ease-in-out entre deux états (gifs à états : index, exercices)
ease <- function(frac) (1 - cos(pi * frac)) / 2

# interpolation linéaire de couleurs point par point, dans l'espace RGB
interp_couleur <- function(col_a, col_b, frac) {
  rgb_a <- col2rgb(col_a)
  rgb_b <- col2rgb(col_b)
  rgb_t <- (1 - frac) * rgb_a + frac * rgb_b
  rgb(rgb_t[1, ], rgb_t[2, ], rgb_t[3, ], maxColorValue = 255)
}

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

# Assemble les n_frames_gif frames d'un dossier en gif, sauvegardé
# dans base_dir
creer_gif <- function(dir_out, nom_gif, n_frames_gif = n_frames) {
  png_files <- sprintf(
    file.path(dir_out, "frame_%03d.png"), seq_len(n_frames_gif)
  )
  gifski(
    png_files, gif_file = file.path(base_dir, nom_gif),
    width = largeur_gif, height = hauteur_gif,
    delay = delai_gif, loop = TRUE
  )
}





# ================================================================
# 01 — Régression linéaire
# ================================================================
if (1 %in% changer){
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

  # bruit perpendiculaire constant, converti en bruit vertical pour
  # compenser la pente
  bruit_perp <- amplitude * sin(t + phases) + bruit_fixe
  facteur    <- sqrt(1 + beta1_t^2)
  Y_t        <- beta1_t * (X - X_mid) + bruit_perp * facteur

  df <- data.frame(X = X, Y = Y_t)
  modele <- lm(Y ~ X, data = df)
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts,
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "lineaire_plot.gif")
}

# ================================================================
# 02 — Simulations de Monte Carlo (SUGGESTION, à adapter)
# ================================================================
# Idée : illustrer la loi des grands nombres et le rétrécissement de
# l'erreur de Monte Carlo à mesure que B augmente, sur l'exemple de
# couverture du chapitre (theta = p = 0.95). B balaie une seule fois
# B_min → B_max (ease-in-out) sur toute la durée du gif ; en boucle,
# ça donne une succession de balayages B_min → B_max, avec un retour
# net (pas de ping-pong) au redémarrage de chaque boucle.
if (2 %in% changer){
set.seed(8)
dir_out <- preparer_dossier("frames_monte_carlo")

p_vrai <- 0.95
B_min  <- 10
B_max  <- 1000

# 1 = l'intervalle contient la vraie valeur
indicatrices <- rbinom(B_max, 1, p_vrai)
theta_hat    <- cumsum(indicatrices) / seq_len(B_max)
erreur_type  <- sqrt(p_vrai * (1 - p_vrai) / seq_len(B_max))

for (i in seq_len(n_frames)) {
  # t de 0 à pi sur tout le gif : balayage à sens unique B_min → B_max
  t   <- (i - 1) / (n_frames - 1) * pi
  B_t <- round(B_min + (B_max - B_min) * (1 - cos(t)) / 2)

  df     <- data.frame(B = seq_len(B_t), theta = theta_hat[seq_len(B_t)])
  bande  <- data.frame(
    B    = seq_len(B_t),
    bas  = p_vrai - 1.96 * erreur_type[seq_len(B_t)],
    haut = p_vrai + 1.96 * erreur_type[seq_len(B_t)]
  )

  p <- ggplot() +
    geom_line(data = bande, aes(B, bas), color = couleur_pts,
              linewidth = epaisseur_ligne, linetype = "dashed") +
    geom_line(data = bande, aes(B, haut), color = couleur_pts,
              linewidth = epaisseur_ligne, linetype = "dashed") +
    # geom_ribbon(data = bande, aes(B, ymin = bas, ymax = haut),
    #             fill = couleur_pts, alpha = 0.2) +
    geom_hline(yintercept = p_vrai, color = couleur_accent,
               linewidth = epaisseur_ligne) +
    geom_line(data = df, aes(B, theta),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    coord_cartesian(xlim = c(1, B_max),
                    ylim = c(p_vrai - 0.15, p_vrai + 0.15)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "monte_carlo_plot.gif")
}


# ================================================================
# 03 — Régression logistique
# ================================================================
if (3 %in% changer){
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
    modele <- glm(
      Y ~ X, data = df, family = binomial(link = "logit"), weights = Poids
    )
    grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
    grille$Y_hat <- predict(modele, newdata = grille, type = "response")

    p <- ggplot() +
      geom_point(data = df, aes(X, Y, size = Poids),
                 color = couleur_pts, alpha = alpha_pts) +
      geom_line(data = grille, aes(X, Y_hat),
                color = couleur_ligne, linewidth = epaisseur_ligne) +
      coord_cartesian(xlim = range(X), ylim = c(0, 1)) +
      theme_gif

    sauvegarder_frame(p, dir_out, i)
  }

  creer_gif(dir_out, "logit_plot.gif")
}



# ================================================================
# 04 — Modèles linéaires généralisés (Poisson)
# ================================================================
if (4 %in% changer){
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

  # signal individuel lisse (variance ~1), converti en bruit d'échelle
  # sqrt(lambda_t)
  z_t <- (sin(t + phases_z1) + 0.5 * sin(freq_z2 * t + phases_z2)) / sqrt(1.25)
  Y_t <- pmax(lambda_t + sqrt(lambda_t) * z_t, 0)

  df <- data.frame(X = X, Y = Y_t)
  modele <- glm(Y ~ X, data = df, family = poisson(link = "log"))
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille, type = "response")

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts,
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "glm_plot.gif")
}

# ================================================================
# 05 — Prédiction et validation croisée
# ================================================================
if (5 %in% changer){
set.seed(3)
dir_out <- preparer_dossier("frames_vc")
 
k           <- 3
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
    geom_point(data = df, aes(X, Y, color = role),
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    scale_color_manual(values = c(train = couleur_pts, test = couleur_accent)) +
    theme_gif
 
  sauvegarder_frame(p, dir_out, i)
}
 
creer_gif(dir_out, "cv_plot.gif", n_frames_gif = n_frames_vc)
}
# ================================================================
# 06 — Sélection et régularisation (LASSO)
# ================================================================
if (6 %in% changer){
set.seed(4)
dir_out <- preparer_dossier("frames_lasso")

coef_base  <- runif(n_coef, -2, 2)
noms_coef  <- c("lambda", paste0("b", 1:n_coef))
noms       <- factor(noms_coef, levels = noms_coef)
# |beta_j chapeau|, syntaxe plotmath
labels_coef <- c("lambda", paste0("abs(hat(beta)[", 1:n_coef, "])"))

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  lambda_t <- (cos(t) + 1) / 2   # oscille entre 0 et 1
  coef_t   <- sign(coef_base) * pmax(abs(coef_base) - lambda_t * 1.8, 0)
  para <- c(lambda_t, coef_t)

  df <- data.frame(coef = noms, valeur = abs(para), label = labels_coef)

  p <- ggplot(df, aes(coef, valeur)) +
    geom_col(fill = c(couleur_accent, rep(couleur_pts, n_coef)), alpha = 0.75) +
    geom_text(aes(label = label), parse = TRUE, vjust = -0.4,
              color = couleur_ligne, size = 16, fontface = "bold") +
    geom_hline(yintercept = 0, color = couleur_ligne,
               linewidth = epaisseur_ligne) +
    coord_cartesian(ylim = c(0, 2.4)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "regularisation_plot.gif")
}

# ================================================================
# 07 — KNN et noyau
# ================================================================
if (7 %in% changer){
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
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts,
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = df_l, aes(X, Y),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "knn_noyau_plot.gif")
}

# ================================================================
# 08 — Splines
# ================================================================
if (8 %in% changer){
set.seed(6)
dir_out <- preparer_dossier("frames_splines")

X <- sort(runif(n_pts, 0, 10))
Y <- sin(X) + rnorm(n_pts, 0, 0.25)

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  # pas de round() : évite les sauts entre frames
  df_t <- max(8.5 + 6.5 * cos(t), 2)

  modele <- smooth.spline(X, Y, df = df_t)
  grille <- seq(min(X), max(X), length.out = 200)
  Y_hat  <- predict(modele, grille)$y

  df_pts <- data.frame(X = X, Y = Y)
  df_l   <- data.frame(X = grille, Y = Y_hat)

  p <- ggplot() +
    geom_point(data = df_pts, aes(X, Y), color = couleur_pts,
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = df_l, aes(X, Y),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "splines_plot.gif")
}

# ================================================================
# 09 — GAM
# ================================================================
if (9 %in% changer){
set.seed(7)
# corrigé : était "frames_splines" par erreur
dir_out <- preparer_dossier("frames_gam")

X <- sort(runif(n_pts, 0, 10))
bruit_fixe <- rnorm(n_pts, 0, 0.3)   # fixé une fois, hors boucle

for (i in seq_len(n_frames)) {
  t <- (i - 1) / n_frames * 2 * pi
  amp_t <- 1 + 0.8 * cos(t)

  Y <- amp_t * sin(X) + 0.3 * X + bruit_fixe

  df <- data.frame(X = X, Y = Y)
  # fx = TRUE : degrés de liberté fixes
  modele <- gam(Y ~ s(X, k = 10, fx = TRUE), data = df)
  grille <- data.frame(X = seq(min(X), max(X), length.out = 200))
  grille$Y_hat <- predict(modele, newdata = grille)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y), color = couleur_pts,
               alpha = alpha_pts, size = taille_pts) +
    geom_line(data = grille, aes(X, Y_hat),
              color = couleur_ligne, linewidth = epaisseur_ligne) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "gam_plot.gif")
}

# ================================================================
# Exercices — méta-gif, synthèse animée des 8 chapitres (SUGGESTION)
# ================================================================
# Même mécanique que l'index (interpolation d'un nuage de points
# entre états), mais les états rejouent ici la forme caractéristique
# de chaque chapitre (hors Monte Carlo), dans l'ordre du cours : un
# clin d'œil général à "tout ce qu'on pratique dans les exercices".
if (10 %in% changer){
set.seed(11)
dir_out <- preparer_dossier("frames_exercices")

X     <- sort(runif(n_pts, 0, 10))
X_mid <- mean(X)
bruit_indiv <- rnorm(n_pts, 0, 1)

# 01 — régression linéaire : nuage autour d'une droite
etat_lineaire <- list(
  x = X, y = 0.9 * (X - X_mid) + bruit_indiv * 0.4,
  col = rep(couleur_pts, n_pts)
)

# 03 — régression logistique : courbe en S
etat_logistique <- list(
  x = X,
  y = 8 * (plogis(1.3 * (X - X_mid)) - 0.5) + bruit_indiv * 0.25,
  col = rep(couleur_pts, n_pts)
)

# 04 — GLM (Poisson) : croissance exponentielle de comptages
lambda_X <- exp(0.3 * (X - X_mid))
etat_glm <- list(
  x = X,
  y = 4 * (lambda_X - mean(lambda_X)) / diff(range(lambda_X)) +
    bruit_indiv * 0.3,
  col = rep(couleur_pts, n_pts)
)

# 05 — validation croisée : même droite, coloration train/test
role_cv <- sample(c("train", "test"), n_pts, replace = TRUE, prob = c(2, 1))
etat_cv <- list(
  x = X, y = etat_lineaire$y,
  col = ifelse(role_cv == "train", couleur_pts, couleur_accent)
)

# # 06 — sélection/régularisation (LASSO) : coefficients décroissants,
# # plusieurs ramenés à zéro (X sert d'indice de coefficient)
# y_lasso <- rev(sort(abs(rnorm(n_pts, 0, 1.3))))
# y_lasso[round(n_pts * 0.55):n_pts] <- 0
# etat_lasso <- list(
#   x = X, y = y_lasso - 1.5,
#   col = ifelse(y_lasso == 0, couleur_accent, couleur_pts)
# )

# 07 — KNN et noyau : lissage local d'une sinusoïde
etat_knn <- list(
  x = X, y = 3 * sin(0.8 * X) + bruit_indiv * 0.3,
  col = rep(couleur_pts, n_pts)
)

# 08 — splines : sinusoïde déphasée (forme non-linéaire distincte du KNN)
etat_splines <- list(
  x = X, y = 3 * sin(1.3 * X - 1) + bruit_indiv * 0.3,
  col = rep(couleur_pts, n_pts)
)

# 09 — GAM : somme d'une tendance lisse et d'une tendance linéaire
etat_gam <- list(
  x = X,
  y = 2 * sin(0.8 * X) + 0.4 * (X - X_mid) + bruit_indiv * 0.25,
  col = rep(couleur_pts, n_pts)
)

etats <- list(
  etat_lineaire, etat_logistique, etat_glm, etat_cv,
  #etat_lasso, 
  etat_knn, etat_splines, etat_gam
)
n_etats <- length(etats)

for (i in seq_len(n_frames)) {
  pos <- (i - 1) / n_frames * n_etats

  idx_a <- (floor(pos) %% n_etats) + 1
  idx_b <- (floor(pos) + 1) %% n_etats + 1
  frac  <- ease(pos - floor(pos))

  X_t   <- (1 - frac) * etats[[idx_a]]$x   + frac * etats[[idx_b]]$x
  Y_t   <- (1 - frac) * etats[[idx_a]]$y   + frac * etats[[idx_b]]$y
  Col_t <- interp_couleur(etats[[idx_a]]$col, etats[[idx_b]]$col, frac)

  df <- data.frame(X = X_t, Y = Y_t, Couleur = Col_t)

  p <- ggplot() +
    geom_point(data = df, aes(X, Y, color = Couleur),
               alpha = alpha_pts, size = taille_pts) +
    scale_color_identity() +
    coord_cartesian(xlim = c(0, 10), ylim = c(-5, 5)) +
    theme_gif

  sauvegarder_frame(p, dir_out, i)
}

creer_gif(dir_out, "exercices_plot.gif")
}

# ================================================================
# Index — nuage de points qui se réorganise (sans ajustement)
# ================================================================
if (0 %in% changer){
set.seed(9)
dir_out <- preparer_dossier("frames_index")

 
X <- sort(runif(n_pts, 0, 10))
 
# chaque état est un triplet (x, y, couleur) par point, dans le même
# ordre pour tous les états : c'est ce qui permet d'interpoler à la fois
# la position et la couleur de façon continue d'un état à l'autre.
bruit_indiv <- rnorm(n_pts, 0, 1)

x_final <- runif(n_pts, 0, 10)
y_final <- runif(n_pts, -5, 5)
 
etat_aleatoire <- list(x = x_final, y = y_final,
                        col = rep(couleur_pts, n_pts))
etat_lineaire  <- list(x = X, y = 0.9 * (X - mean(X)) + bruit_indiv * 0.4,
                        col = rep(couleur_pts, n_pts))
etat_sinus     <- list(x = X, y = 3 * sin(0.8 * X) + bruit_indiv * 0.35,
                        col = rep(couleur_pts, n_pts))
 
# état final : points aléatoires sur le rectangle, répartis en 3
# secteurs angulaires (3 vecteurs partant du centre, séparés de 120°,
# délimitent les classes)


centre_x <- 5   # centre du rectangle (domaine x : 0 à 10)
centre_y <- 0   # centre du rectangle (domaine y : -5 à 5)

angle         <- atan2(y_final - centre_y, x_final - centre_x)   # -pi à pi
angle_positif <- ifelse(angle < 0, angle + 2 * pi, angle)        # 0 à 2*pi
# 3 secteurs de 120°
classe_finale <- floor(angle_positif / (2 * pi / 3)) + 1

couleurs_classes <- c(couleur_pts, couleur_ligne, couleur_accent)

etat_grille <- list(
  x = x_final, y = y_final, col = couleurs_classes[classe_finale]
)

etats   <- list(etat_aleatoire, etat_lineaire, etat_sinus, etat_grille)
n_etats <- length(etats)

# angles des 3 vecteurs de séparation (mêmes bornes que classe_finale)
angles_separation <- c(0, 2 * pi / 3, 4 * pi / 3)
# longueur des vecteurs une fois complètement affichés
rayon_separation <- 6

for (i in seq_len(n_frames)) {
  pos <- (i - 1) / n_frames * n_etats

  idx_a <- (floor(pos) %% n_etats) + 1
  idx_b <- (floor(pos) + 1) %% n_etats + 1
  frac  <- ease(pos - floor(pos))

  X_t   <- (1 - frac) * etats[[idx_a]]$x   + frac * etats[[idx_b]]$x
  Y_t   <- (1 - frac) * etats[[idx_a]]$y   + frac * etats[[idx_b]]$y
  Col_t <- interp_couleur(etats[[idx_a]]$col, etats[[idx_b]]$col, frac)
 
  # les vecteurs n'apparaissent qu'en entrant dans l'état "grille"
  # (croissance, poids 0→1) ou en le quittant (rétraction, poids 1→0) ;
  # nuls lors des autres transitions
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
    xend = centre_x +
      poids_separation * rayon_separation * cos(angles_separation),
    yend = centre_y +
      poids_separation * rayon_separation * sin(angles_separation)
  )
 
    # couleur de la flèche : blanche (invisible) quand poids_separation = 0,
  # grise et pleinement visible quand poids_separation = 1
  couleur_separation <- interp_couleur("white", "gray50", poids_separation)
 
  p <- ggplot() +
    geom_segment(data = df_separation, aes(x, y, xend = xend, yend = yend),
                 color = couleur_separation, linewidth = epaisseur_ligne * 0.6,
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    geom_point(data = df, aes(X, Y, color = Couleur),
               alpha = alpha_pts, size = taille_pts) +
    scale_color_identity() +
    coord_cartesian(xlim = c(0, 10), ylim = c(-5, 5)) +
    theme_gif
 
  sauvegarder_frame(p, dir_out, i)
}
 
creer_gif(dir_out, "index_plot.gif", n_frames_gif = n_frames)
}