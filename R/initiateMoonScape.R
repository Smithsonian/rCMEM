#' Initialize a mineral-only soil profile
#'
#' @param depth Total initial soil depth (cm).
#' @param increment Thickness of each cohort (cm).
#' @param initElv Initial surface elevation (cm).
#' @param species_codes Character vector of species codes.
#' @param meanMineralPackingDensity Mineral packing density (g/cm3).
#' @param meanOmPackingDensity Organic packing density (g/cm3).
#' @param initial_year Calendar year of initial conditions.
#'
#' @return A list containing initial_cohorts,
#'   initial_scenario, and initial_species.
#' @export
initiateMoonScape <- function(
    depth = 100,
    increment = 1,
    initElv = 0,
    species_codes = c("SCAM", "SPPA"),
    meanMineralPackingDensity = 1.5,
    meanOmPackingDensity = 0.1,
    initial_year = 2000
) {
  
  stopifnot(
    length(depth) == 1, is.finite(depth), depth > 0,
    length(increment) == 1, is.finite(increment),
    increment > 0,
    length(initElv) == 1, is.finite(initElv),
    length(initial_year) == 1, is.finite(initial_year),
    meanMineralPackingDensity > 0,
    meanOmPackingDensity > 0,
    length(species_codes) > 0,
    !anyDuplicated(species_codes)
  )
  
  n_cohorts <- round(depth / increment)
  
  if (!isTRUE(all.equal(
    n_cohorts * increment, depth
  ))) {
    stop("depth must be divisible by increment.")
  }
  
  # Cohorts are ordered from surface to bottom.
  initial_cohorts <- data.frame(
    year = rep(initial_year, n_cohorts),
    cohort_index = n_cohorts:1,
    mineral = rep(
      increment * meanMineralPackingDensity,
      n_cohorts
    ),
    fast = rep(0, n_cohorts),
    slow = rep(0, n_cohorts),
    root = rep(0, n_cohorts),
    vol = rep(increment, n_cohorts),
    nonRootVol = seq_len(n_cohorts) * increment,
    cumVol = seq_len(n_cohorts) * increment,
    inputYrs = rep(1, n_cohorts),
    omPackingDensity = rep(
      meanOmPackingDensity, n_cohorts
    ),
    mineralPackingDensity = rep(
      meanMineralPackingDensity, n_cohorts
    )
  )
  
  initial_scenario <- data.frame(
    year = initial_year,
    years_per_iteration = 1,
    surface_elevation = initElv,
    aboveground_biomass = 0,
    belowground_biomass = 0,
    totalRootVolume = 0,
    sediment_delivered = 0,
    rootToShoot = 0,
    rootTurnover = 0,
    abovegroundTurnover = 0,
    rootDepthMax = 0,
    rootPackingDensity = 0.07,
    lambda = 1,
    rootShape = 0
  )
  
  initial_species <- as.data.frame(
    matrix(
      0,
      nrow = 1,
      ncol = length(species_codes),
      dimnames = list(NULL, species_codes)
    )
  )
  
  initial_species <- data.frame(
    year = initial_year,
    initial_species,
    check.names = FALSE
  )
  
  return(list(
    initial_cohorts = initial_cohorts,
    initial_scenario = initial_scenario,
    initial_species = initial_species
  ))
}
