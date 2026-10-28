#' Calculate sediment mass needed to offset a projected accretion deficit
#' @param msl mean sea level vector
#' @param surface_elevation surface elevation vector
#' @param years vector of calendar years
#' @param window window of time that sea level rise and elevation change are classified over
#' @param mineral_density density of added mineral matter in grams per cubic centimetseer
#' 
#' @export 
calculateSedimentAugmentation <- function(
    msl,
    surface_elevation,
    years,
    window = 10,
    mineral_density = 1.5
) {
  
  n <- length(years)
  
  if (n < 2 || years[n] - years[1] < window) {
    return(0)
  }
  
  start <- max(which(years <= years[n] - window))
  
  slr <- (msl[n] - msl[start]) /
    (years[n] - years[start])
  
  sec <- (surface_elevation[n] -
            surface_elevation[start]) /
    (years[n] - years[start])
  
  deficit <- max(0, slr - sec)
  
  # Additional mineral mass per unit area per year
  deficit * mineral_density
}