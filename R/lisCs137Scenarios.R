#' List available 137Cs fallout scenarios
#'
#' Shows the built-in fallout scenarios included with the package.
#'
#' @return Character vector of scenario names
#' @export
listCs137Scenarios <- function() {
  
  scenario_dir <- system.file(
    "extdata",
    "fallout",
    package = "rCMEM"
  )
  
  files <- list.files(
    scenario_dir,
    pattern = "\\.Rds$",
    full.names = FALSE
  )
  
  sub("\\.rds$", "", files)
}
