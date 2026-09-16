library(ggplot2)
library(dplyr)

set.seed(42)

# Simulation de données groupées (style turbines)
n_groupes <- 20
df <- data.frame(
  X = seq(0, 7, length.out = n_groupes)
)

# Vrai modèle logistique sous-jacent
beta0_vrai <- -4
beta1_vrai <- 1

df <- df |>
  mutate(
    Poids = sample(5:10, n_groupes, replace = TRUE),      # nombre d'essais par groupe
    eta_vrai = beta0_vrai + beta1_vrai * X,
    p_vrai = plogis(eta_vrai),
    Succes = rbinom(n_groupes, size = Poids, prob = p_vrai),
    Y = Succes / Poids                                       # proportion observée
  )

# Ajustement du modèle logit
modele <- glm(
  cbind(Succes, Poids - Succes) ~ X,
  data = df,
  family = binomial(link = "logit")
)

# Grille de prédiction
grille <- data.frame(
  X = seq(min(df$X), max(df$X), length.out = 200)
)
grille$Y_hat <- predict(modele, newdata = grille, type = "response")

# Graphique
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
  theme_void() +
  theme(legend.position = "none")

p  # affiche dans le Plots pane de Positron (fond blanc là, normal)

setwd("/Users/jmiron/Library/CloudStorage/OneDrive-UniversitéLaval/Quarto/STT-4300/images")

ggsave(
  "logit_plot.png",
  plot = p,
  width = 8, height = 5.5, dpi = 300,
  bg = "transparent"
)