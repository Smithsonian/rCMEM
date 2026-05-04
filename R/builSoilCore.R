#' Simulate a sediment core from a set of cohorts
#' @param cohorts a data frame, output from runCohortMem, tracking mineral and organic mass cohorts over each year of the simulation
#' @param coreYear an integer, year in form YYYY, specify a year to simulate taking a sediment core
#' @param coreDepth an integer, depth, specify a depth to simulate coring to and assume 1 cm sampling intervals
#' @param coreMins a vector of sampling depth minimums to simulate coring subsamples, this is an alternative to depth, and 1cm increments
#' @param coreMaxs a vector of sampling depth maximums to simulate coring subsamples, this is an alternative to depth, and 1cm increments
#' 
#' @return a dataframe with variables simulated from a soil core
#' @export
buildSoilCore <- function(cohorts, coreYear, coreDepth=100, 
                             coreMaxs=1:coreDepth, coreMins=coreMaxs-1
                          ) {
  
  # Filter only to cohorts in the core year
  cohortsInCoreYear <- cohorts[cohorts$year == coreYear,]
  
  if(max(coreMaxs)>max(cohortsInCoreYear$cumVol)) {
    stop("Simulated core is deeper than deposit.")
  }
  
  mineral <- cohortsInCoreYear$mineral
  fast <- cohortsInCoreYear$fast
  slow <- cohortsInCoreYear$slow
  root <- cohortsInCoreYear$root
  vol <- cohortsInCoreYear$vol
  nonRootVol <- cohortsInCoreYear$nonRootVol
  cumVol <- cohortsInCoreYear$cumVol
  inputYrs <- cohortsInCoreYear$inputYrs
  mineralPackingDensity <- cohortsInCoreYear$mineralPackingDensity
  omPackingDensity <- cohortsInCoreYear$omPackingDensity
  
  # Preallocate memory
  
  # Create tracking vectors
  # Include the surface cohort which hasn't been incorporated yet
  core_mineral = c(mineral[1], rep(NA, length(coreMaxs)))
  core_fast = c(fast[1], rep(NA, length(coreMaxs)))
  core_slow = c(slow[1], rep(NA, length(coreMaxs)))
  core_root = c(root[1], rep(NA, length(coreMaxs)))
  core_vol = c(vol[1], rep(NA, length(coreMaxs)))
  core_nonRootVol = c(nonRootVol[1], rep(NA, length(coreMaxs)))
  core_cumVol = c(cumVol[1], rep(NA, length(coreMaxs)))
  core_inputYrs = c(inputYrs[1], rep(NA, length(coreMaxs)))
  core_omPackingDensity = c(omPackingDensity[1], rep(NA, length(coreMaxs)))
  core_mineralPackingDensity = c(mineralPackingDensity[1], rep(NA, length(coreMaxs)))
  core_age_from_surface = c(inputYrs[1], rep(NA, length(coreMaxs)))
  core_calendar_age = c(coreYear-inputYrs[1], rep(NA, length(coreMaxs)))
  
  # Iterate through core increments
  cohort_n_start <- 2
  
  n_cohorts <- nrow(cohortsInCoreYear)
  
  # Express surface cohort as negative volume
  coreMins <- c(-vol[1], coreMins)
  coreMaxs <- c(0, coreMaxs)
  
  for (core_increment in 2:length(coreMaxs)) {
    
    weighted_mineral = 0
    weighted_fast = 0
    weighted_slow = 0
    weighted_root = 0
    weighted_inputYrs = 0
    
    weighted_opd = 0
    weighted_mpd = 0
    
    total_weights = 0

    # Do the rest of the cohorts
    for (cohort in cohort_n_start:n_cohorts) {
      
      cohort_depth_min <- cumVol[cohort-1]
      cohort_depth_max <- cumVol[cohort]
      
      # Is this a zero vol cohort?
      weight_calc_part1 = as.numeric(cohort_depth_max == cohort_depth_min)
      weight_calc_part2 = min(cohort_depth_max, coreMaxs[core_increment], na.rm=T)-max(cohort_depth_min, coreMins[core_increment], na.rm=T)
      
      # IF we pass the bottom of matching depth increments, break the loop so we don't have to do more calculations
      if (weight_calc_part2 < 0) {
        # Update cohort_n_start
        break()
      }  
        
      # Depth weight matrix
      # Is the increment partially in or totallly in?
      depth_weight <- max(weight_calc_part1, max(weight_calc_part2, 0, na.rm=T) / (cohort_depth_max-cohort_depth_min), na.rm=T)
      
      weighted_mineral = weighted_mineral + mineral[cohort] * depth_weight
      weighted_fast = weighted_fast + fast[cohort] * depth_weight
      weighted_slow = weighted_slow + slow[cohort] * depth_weight
      weighted_root = weighted_root + root[cohort] * depth_weight
      weighted_inputYrs = weighted_inputYrs + inputYrs[cohort] * depth_weight
      
      weighted_inputYrs = weighted_inputYrs + inputYrs[cohort] * depth_weight
      
      weighted_opd = weighted_opd + omPackingDensity[cohort] * depth_weight
      weighted_mpd = weighted_mpd + mineralPackingDensity[cohort] * depth_weight
      
      total_weights = total_weights + depth_weight 
      
      cohort_n_start <- cohort
      
    } # End of cohort Iteration loop
    
    core_mineral[core_increment] = weighted_mineral
    core_fast[core_increment] = weighted_fast
    core_slow[core_increment] = weighted_slow
    core_root[core_increment] = weighted_root
    core_vol[core_increment] = coreMaxs[core_increment] - coreMins[core_increment]
    core_cumVol[core_increment] = core_vol[core_increment] + core_cumVol[core_increment-1] 
    core_inputYrs[core_increment] = weighted_inputYrs
    
    core_omPackingDensity[core_increment] = weighted_opd/total_weights
    core_mineralPackingDensity[core_increment] = weighted_mpd/total_weights
    
    core_nonRootVol[core_increment] = core_mineral[core_increment]/core_mineralPackingDensity[core_increment] + (core_fast[core_increment]+core_slow[core_increment])/core_omPackingDensity[core_increment]
    
    core_age_from_surface[core_increment] = core_age_from_surface[core_increment-1] + weighted_inputYrs
    core_calendar_age[core_increment] = core_calendar_age[core_increment-1] - weighted_inputYrs
    
    
  } # end of depht increment iteration loop
  
  # Add to the list of outputs
  return(data.frame(year = rep(coreYear, length(coreMaxs)),
                    cohortIndex = 1:length(coreMaxs),
                    mineral =  core_mineral,
                    fast = core_fast,
                    slow = core_slow,
                    root = core_root,
                    vol = core_vol,
                    cumVol = core_cumVol,
                    inputYrs = core_inputYrs,
                    omPackingDensity = core_omPackingDensity,
                    mineralPackingDensity = core_mineralPackingDensity,
                    nonRootVol = core_nonRootVol,
                    age_from_surface = core_age_from_surface,
                    calendar_age = core_calendar_age
                    )) 
}
