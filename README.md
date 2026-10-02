# Site de cours — STT-4300

Site du cours STT-4300 (régression linéaire, modèles linéaires généralisés, validation croisée, régularisation, méthodes non-paramétriques, GAM) : page d'accueil, plan de cours, échéancier, et les notes de cours complètes sous forme de livre [Quarto](https://quarto.org).

Le dépôt contient **deux projets Quarto indépendants** : le site (racine) et le livre (`notes/`), rendus séparément puis fusionnés sous une même arborescence déployée — voir « Structure » ci-dessous.

## Développement local

Installer [Quarto](https://quarto.org/docs/get-started/), puis :

```bash
quarto preview                  # prévisualise le site (racine)
quarto preview notes            # prévisualise le livre
quarto render --to html         # génère le site dans _site/
quarto render notes --to html   # génère le livre dans _site/notes/ (voir notes/_quarto.yml)
```

Pour voir le site complet tel qu'il sera déployé (site + livre fusionnés), il faut lancer les deux rendus dans cet ordre (le rendu du livre doit être la dernière écriture dans `_site/`), puis ouvrir `_site/index.html`.

## Publication

Le site est publié automatiquement sur GitHub Pages via GitHub Actions à chaque push sur `main` (voir `.github/workflows/publish.yml`) : rendu du site, rendu du livre dans `_site/notes/`, puis déploiement de `_site/` sur la branche `gh-pages`. La première fois, active GitHub Pages dans Settings → Pages en sélectionnant la branche `gh-pages` comme source.

## Structure

```
.
├── _quarto.yml                 # configuration du site (racine)
├── index.qmd                   # page d'accueil du site
├── plan-de-cours.qmd
├── echeancier.qmd
├── notebooks-classe/           # notebooks d'exemples en classe (webR, .qmd + .html versionné)
│   └── NN-nom-chapitre-exemples.qmd
├── notes/                      # le livre : projet Quarto indépendant (type: book)
│   ├── _quarto.yml              # configuration du livre, output-dir: ../_site/notes
│   ├── index.qmd                 # préface
│   ├── 01-regression-lineaire.qmd
│   ├── 02-monte-carlo.qmd
│   ├── 03-regression-logistique.qmd
│   ├── 04-modeles-lineaires-generalises.qmd
│   ├── 05-prediction-validation-croisee.qmd
│   ├── 06-selection-regularisation.qmd
│   ├── 07-knn-noyau.qmd
│   ├── 08-splines.qmd
│   ├── 09-gam.qmd
│   ├── exercices-theoriques.qmd   # exercices théoriques (hors ordre des chapitres)
│   ├── exercices-en-classe.qmd    # liste les notebooks-classe/*-exemples.qmd marqués afficher: true
│   └── styles.css                 # style des blocs Définition/Propriété/Remarque/Exemple
└── .github/workflows/publish.yml
```

## Notes sur la conversion

Ces fichiers ont été convertis automatiquement à partir d'un document LaTeX source. Points à vérifier/possiblement ajuster :

- Les environnements `Définition`, `Propriété`, `Remarque`, `Exemple` sont rendus comme des blocs stylisés (`.env-*` dans `notes/styles.css`), sans numérotation automatique (contrairement à LaTeX/`amsthm`). Si tu veux une numérotation automatique, envisage l'extension Quarto `theorem` ou les [environnements crossref natifs](https://quarto.org/docs/authoring/cross-references.html#theorems-and-proofs) (`#def-`, `#thm-`, etc.).
- Les blocs de code R (` ```r `) affichent le code mais ne l'exécutent pas (les jeux de données comme `turbines`, `quilpie`, `Auto`, `Hitters` ne sont pas fournis). Pour les exécuter réellement, ajoute les packages nécessaires (`GLMsData`, `ISLR`, `MASS`, `mgcv`, `splines`, `KernSmooth`, `boot`, `glmnet`, `FNN`, `clinfun`, `pROC`) et retire les commentaires `#` s'il y a lieu.
- Les renvois internes (`@sec-...`) pointent vers les titres de chapitres/sections ; les renvois vers des remarques spécifiques (ex. la remarque sur l'interprétation du lien logit) sont des liens simples vers l'ancre correspondante.
