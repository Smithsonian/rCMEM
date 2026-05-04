#' Calculate Carbon Flux from cohorts and scenario tables
#' 
#' @param cohorts a data frame, output from runCohortMem, tracking mineral and organic mass cohorts over each year of the simulation
#' @param scenario a data frame, output from runCohortMem, tracking annualized variables
#' @param omToOcParams a list, of numerics listed B0, B1 and optionally B3. These define a linear or parabolic relationship between organic matter and organic
#'  
#' @return a list of data frames, including the annualized summaries, and mapped cohorts tracked for every year of the simulation.
#' @export
calculateCarbonFlux <- function(cohorts,
                             scenario,
                             omToOcParams = list(B0=0, B1=0.48)) {
  
  # Calculate C sequestration rate from cohorts table and add it to scenario table

  # To convert OM to OC
  # If parmeter list is 2 long then simple linear correlation
  if (length(omToOcParams) == 2) {
    omToOc <- function(om, B0=omToOcParams$B0, B1=omToOcParams$B1) {return(B0 + om*B1)}
  } else if (length(omToOcParams) == 3) {
    # If parameter list is 3 long, then it's quadratic
    omToOc <- function(om, B0=omToOcParams$B, B1=omToOcParams$B1,
                       B2=omToOcParams$B2) {return(B0 + om*B1 + om^2*B2)}
  } else {
    # If something else then trip an error message
    stop("Invalid number of organic matter to organic carbon conversion parameters,")
  }
  
  # Apparent Carbon Burial Rate
  # Carbon Sequestration Rate
  carbonFluxTab <- cohorts %>%
    # Total organic matter per cohort
    dplyr::mutate(total_om_perCoh = fast + slow + root) %>%
    dplyr::group_by(year) %>%
    # Get the total cumulative observed and sequestered organic matter for the profile
    dplyr::summarise(cumulativeTotalOm = sum(total_om_perCoh),
                     cumulativeSequesteredOm = sum(slow)) %>% 
    # Caluclate fluxes by comparing total at time step i to time step i - 1
    dplyr::mutate(omFlux = cumulativeTotalOm - lag(cumulativeTotalOm),
                  omSequestration = cumulativeSequesteredOm - lag(cumulativeSequesteredOm),
                  # Convert organic matter to organic carbon using function defined in previous step.
                  cFlux = omToOc(om=omFlux),
                  cSequestration = omToOc(om = omSequestration))
  
  # Join flux table to annual time step table
  scenario <- scenario %>% 
    dplyr::left_join(carbonFluxTab, by="year") %>% 
    # Normalize by time
    dplyr::mutate(omFlux = omFlux / years_per_iteration,
                  omSequestration = omSequestration / years_per_iteration,
                  cFlux = cFlux / years_per_iteration,
                  cSequestration = cSequestration / years_per_iteration
                  )

 # Return annual time steps and full set of cohorts
 return(list(scenario = scenario, cohorts = cohorts))
}

