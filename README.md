# Notes de cours — STT-4300

Notes de cours sur la régression linéaire, les modèles linéaires généralisés (régression logistique, régression de Poisson), la validation croisée, la régularisation, les méthodes non-paramétriques (kNN, noyau, splines) et les modèles additifs généralisés (GAM).

Publié comme livre [Quarto](https://quarto.org).

## Développement local

Installer [Quarto](https://quarto.org/docs/get-started/), puis :

```bash
quarto preview        # prévisualisation locale avec rechargement automatique
quarto render         # génère le site final dans _book/
```

## Publication

Le site est publié automatiquement sur GitHub Pages via GitHub Actions à chaque push sur `main` (voir `.github/workflows/publish.yml`). La première fois, active GitHub Pages dans Settings → Pages en sélectionnant la branche `gh-pages` comme source.

## Structure

```
.
├── _quarto.yml               # configuration du livre
├── index.qmd                 # préface
├── 01-regression-lineaire.qmd
├── 02-regression-logistique.qmd
├── 03-modeles-lineaires-generalises.qmd
├── 04-prediction-validation-croisee.qmd
├── 05-selection-regularisation.qmd
├── 06-knn-noyau.qmd
├── 07-splines.qmd
├── 08-gam.qmd
├── styles.css                 # style des blocs Définition/Propriété/Remarque/Exemple
└── .github/workflows/publish.yml
```

## Notes sur la conversion

Ces fichiers ont été convertis automatiquement à partir d'un document LaTeX source. Points à vérifier/possiblement ajuster :

- Les environnements `Définition`, `Propriété`, `Remarque`, `Exemple` sont rendus comme des blocs stylisés (`.env-*` dans `styles.css`), sans numérotation automatique (contrairement à LaTeX/`amsthm`). Si tu veux une numérotation automatique, envisage l'extension Quarto `theorem` ou les [environnements crossref natifs](https://quarto.org/docs/authoring/cross-references.html#theorems-and-proofs) (`#def-`, `#thm-`, etc.).
- Les blocs de code R (` ```r `) affichent le code mais ne l'exécutent pas (les jeux de données comme `turbines`, `quilpie`, `Auto`, `Hitters` ne sont pas fournis). Pour les exécuter réellement, ajoute les packages nécessaires (`GLMsData`, `ISLR`, `MASS`, `mgcv`, `splines`, `KernSmooth`, `boot`, `glmnet`, `FNN`, `clinfun`, `pROC`) et retire les commentaires `#` s'il y a lieu.
- Les renvois internes (`@sec-...`) pointent vers les titres de chapitres/sections ; les renvois vers des remarques spécifiques (ex. la remarque sur l'interprétation du lien logit) sont des liens simples vers l'ancre correspondante.
