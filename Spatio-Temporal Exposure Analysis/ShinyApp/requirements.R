pkgs <- c(
  "data.table", "tidytable", "dplyr", "purrr", "ggplot2", 
  "gt", "purrr", "shiny", "shinydashboard", "leaflet", "glue",
  "lubridate", "ggthemes", "shinythemes", "Hmisc"
)
dependencies <- c("sp", "checkmate")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  } else {
    require(pkg, character.only = TRUE)
  }
}