## Permet à rmarkdown::render() (bouton « Knit ») de trouver pandoc,
## qui n'est pas installé séparément sur ce poste : seul celui fourni
## avec Quarto est disponible (quarto::quarto_render() n'en a pas besoin).
if (Sys.getenv("RSTUDIO_PANDOC") == "") {
  quarto_pandoc_dir <- Sys.which("quarto")
  if (nzchar(quarto_pandoc_dir)) {
    quarto_bin <- dirname(normalizePath(quarto_pandoc_dir))
    candidats <- file.path(quarto_bin, "tools", c("aarch64", "x86_64"))
    candidat <- candidats[file.exists(file.path(candidats, "pandoc"))]
    if (length(candidat) > 0) {
      Sys.setenv(RSTUDIO_PANDOC = candidat[1])
    }
  }
}
