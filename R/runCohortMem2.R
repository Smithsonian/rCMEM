#' Run the cohort MEM model given a sea level rise scenario, initial elevation, a set of plant trait and geopohysical parameters, and optionally, an initial set of soil mass cohorts.
#' 
#' @param run_spinup binary 1 or 0, 1 runs a spinup, 0 skips over a spinup an initiates using a separate initial set of soil cohorts 
#' @param run_scenario bionary 1 or 0, 1 runs a sea-level rise scnario, 0 only runs a spinup
#' @param spinup_iterations integer, the number of spinup iterations to run 
#' @param years_per_spinup_iter integer, the number of years of marsh formation to aggregate into each spinup iteration 
#' @param competition_function binary 1 or 0, 0 indicates a winner take all plant competition scenario, 1 indicates a coexistence with competition scenario 
#' @param scenario_calendar_years a numeric vector of calendar years 
#' @param scenario_iterations integer, the number of iterations in the sea-level change scenario
#' @param years_per_scenario_iter integer, the number of years included in each iteration of the sea-level rise scenario
#' @param initElv numeric, initial elevation of the scenario stage
#' @param initial_cohort_vol for the intial cohort, cohort 1, the initial volume, generically set to 1 m deep
#' @param initial_cohort_years for the intial cohort, cohort 1, the number of years of formation. This is arbitrarily set to a high number
#' @param initial_scenario data frame, a scenario table from a previous rCMEM run
#' @param initial_cohorts data frame, a cohort table from a previous rCMEM run
#' @param initial_species data frame, a species table from a previous rCMEM run
#' @param bMax a numeric vector, for each species a peak biomass in grams per square centimeters
#' @param zVegMin a numeric vector, for each species minimum growing elevation in relative tidal elevation (elevation - msl) / (mhw - msl)
#' @param zVegMax a numeric vector, for each species maximum growing elevation in relative tidal elevation (elevation - msl) / (mhw - msl)
#' @param zVegPeak a numeric vector, for each species peak growing elevation in relative tidal elevation (elevation - msl) / (mhw - msl)
#' @param rootToShoot a numeric vector, for each species, the belowground to aboveground biomass ratio 
#' @param rootTurnover a numeric vector, for each species, the belowground annual turnover rate
#' @param abovegroundTurnover a numeric vector, for each species, the aboveground annual turnover rate
#' @param rootDepthMax a numeric vector, for each species, the maximum rooting depth in centimeters below the surface
#' @param rootPackingDensity a numeric vector, for each species, the root self packing density in grams per square centimeters
#' @param species_codes a character vector, for each species, a unique code or name to track species in output tables
#' @param meanOmPackingDensity numeric, mean soil organic matter self packing density in grams per square centimeters
#' @param meanMineralPackingDensity numeric, mean soil mineral matter self packing density in grams per square centimeters
#' @param omDecayRateFast numeric, annual mass lost of the fast decaying portion of soil organic matter
#' @param omDecayRateSlow numeric, annual mass lost of the slow decaying portion of soil organic matter
#' @param recalcitrantFrac numeric, fraction of soil organic matter in the slow decaying pool
#' @param captureRate numeric, the number of times that a column of sediment will settle out in a full tidal cycle
#' @param suspendedSediment numeric, suspended sediment concentration in grams per cubic centimeter
#' @param rootShape  binary indicating the shape of the rooting zone, 0 for linear, 1 for exponential shape
#' @param relTol numeric, relative tolerance of the non-linear equation solving algorithm
#' @param safe_mode true or false, turns on safe_mode which will stop code and throw errors if some incongrous inputs are made
#' @param flood_frequency numeric vector, indicates the frequence of flood events, must be the same length as the widths of mhwMat and mlwMat
#' @param msl numeric vector, the mean sea levels of each stage of a sea-level change scenario
#' @param mhwMat numeric matrix, the levels of annual average flood events for each stage of the sea level change scenario 
#' @param mlwMat numeric matrix, the levels of annual average ebb events for each stage of the sea level change scenario
#' 
#' @return a list with three data frames, scenario: an annual time series tracking annual inputs and outputs, cohorts: a long form table tracking mass cohorts and soil volumes, and species: species level predictions of 
#' @export
runCohortMem2 <- function(
    run_spinup = 1,
    run_scenario = 1,
    spinup_iterations = 10,
    years_per_spinup_iter = 10,
    competition_function = 0,
    
    # Scenario
    scenario_calendar_years = 1928:2018,
    scenario_iterations = length(scenario_calendar_years),
    
    years_per_scenario_iter = 1,
    
    # Initial conditions
    initElv = -11,
    
    initial_cohort_vol = 100,
    initial_cohort_years = 99999,
    
    initial_scenario = NA,
    initial_cohorts = NA,
    initial_species = NA,
    
    bMax = c(0.0522, 0.0522),
    zVegMin = c(-1.4,-1.4+0.1),
    zVegMax = c(1.8, 1.8+0.1),
    zVegPeak = c(1.5, 1.5+0.1),
    rootToShoot  = c(1.636903, 1.636903),
    rootTurnover = c(0.5, 0.5),
    abovegroundTurnover = c(1, 1),
    rootDepthMax = c(48, 48),
    
    rootPackingDensity = c(0.07, 0.07),
    
    species_codes = c("SCAM", "SPPA"),
    
    meanOmPackingDensity = 0.1,
    meanMineralPackingDensity = 1.5,
    
    omDecayRateFast = 0.99,
    omDecayRateSlow = 0.001,
    recalcitrantFrac = 0.4,
    captureRate = 5,
    suspendedSediment = 6 * 1e-06, # original
    # 0 for linear, 1 for exponential
    rootShape = c(0, 0),
    
    # Other
    relTol = 1e-6,
    
    # Safe mode
    safe_mode = T,
    flood_frequency = c(0.5, 0.46497542, 0.03502458) * 705.79,
    msl,
    mhwMat,
    mlwMat) {
  
  # Counters for later
  nFloods = length(flood_frequency)
  nSpecies = length(bMax)
  
  # If a binary search needs to get used
  maxIter <- ceiling(log(relTol, base = 0.5))
  
  # Run checks
  if (safe_mode == T) {
    
    # Make sure run spinup or run scenario are checked
    if (! (run_spinup == 1 | run_scenario == 1)) {
      stop("Spinup and/or scenario runs must be specified, set one or both to 1.")
    } # end first spinup test
    
    if (run_spinup == 0  & (!is.data.frame(initial_cohorts) | !is.data.frame(initial_scenario) | !is.data.frame(initial_species))) {
      stop("If not running a spinup, then you need to input a cohort table, specify initial aboveground biomass, and initial speices")
    } # end second spinup test
    
    if (run_spinup != 1 &
        !all(names(initial_cohorts) %in% c("year", "cohort_index", "fast", "slow", "mineral", "root", "respired_OM", "years_per_cohort", 
                                           "nonRootVol", "vol", "cumVol", "inputYrs", "omPackingDensity", "mineralPackingDensity"
        ))) {
      stop(
        paste0(
          c("Check to make sure that the inial cohorts file has all of the following:",
            paste(c("year", "cohort_index", "fast", "slow", "mineral", "root", "respired_OM", "years_per_cohort",
                    "nonRootVol", "vol", "cumVol", "inputYrs", "omPackingDensity", "mineralPackingDensity"),
                  collapse = ", ")
          ), collapse = " "
        )
      )
    } # end third spinup test
    
    # Check to see if all the biological parameters have the same number of species
    if (max(length(bMax), 
            length(zVegMin),
            length(zVegMax), 
            length(zVegPeak), 
            length(zVegPeak), 
            length(rootToShoot), 
            length(rootTurnover), 
            length(abovegroundTurnover),
            length(rootPackingDensity),
            length(species_codes)) != nSpecies |
        min(length(bMax), 
            length(zVegMin),
            length(zVegMax), 
            length(zVegPeak), 
            length(zVegPeak), 
            length(rootToShoot), 
            length(rootTurnover), 
            length(abovegroundTurnover),
            length(rootPackingDensity),
            length(species_codes)) != nSpecies
    ) {
      stop("Check plant inputs and make sure that each vector is of equal length.")
    } # end bio-parameters check
    
    if (run_spinup == 1) {
      if (initial_cohort_vol <= max(rootDepthMax)) {
        stop("Iniitial cohort volume isn't enough to handle the specified root depth.") 
      } # end cohort vol max check test
    }  # end spinup check
    
    if (run_spinup == 0) {
      if (max(initial_cohorts$cumVol)<=max(rootDepthMax)) {
        stop("The volume of the cohorts from previous simlation isn't enough to handle the specified root depth.") 
      } # end initial cohort vol test
    } # end spinup check
    
    # Make sure the length of scenario years, msl, MHW and MLW are all the same
    
    # Make sure nrows of MHW and MLW
    
    # And 
    
    # One thing we can't handle right now is a competition function set to coexistance
    # and root shapes that are different
    if (unique(rootShape)>1 & competition_function == 1) {
      stop("The model can't yet handle plant species that coexist, but have different root shapes. That's too hard.")
    } # end root shape check
    
    
  } # end save mode check
  
  # Create some 
  is_spinup <- c(1)
  years_per_iteration <-c(years_per_scenario_iter)
  scenario_i <- c(1)
  
  # Create some indexes
  if (run_spinup == 1) {
    is_spinup <- c(c(1,1), rep(1, spinup_iterations), is_spinup)
    years_per_iteration <-c(c(initial_cohort_years, 1), rep(years_per_spinup_iter, spinup_iterations), years_per_iteration)
    nInitCohorts <- 0
    scenario_i <- c(c(1,1), rep(1, spinup_iterations), scenario_i)
    total_years <- 3 + spinup_iterations
  } else {
    nInitCohorts <- nrow(initial_cohorts)
    total_years <- 1
  }
  
  if (run_scenario == 1) {
    
    is_spinup <- c(is_spinup, rep(0, scenario_iterations-1))
    years_per_iteration <-c(years_per_iteration, rep(years_per_scenario_iter, scenario_iterations-1))
    scenario_i <- c(scenario_i, 2:scenario_iterations)
    
    total_years <- total_years + scenario_iterations-1
    
  } 
  
  # calcualte simulation year 
  simulation_years <- rep(NA, total_years)
  
  simulation_years[total_years] = max(scenario_calendar_years)
  
  for (iter in (length(simulation_years)-1):1) {
    simulation_years[iter] = simulation_years[iter+1] - years_per_iteration[iter+1]  
  }
  
  # Set up vectors and matrices
  {
    cohort_year <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    cohort_index <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    mineral <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    fast <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    slow <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    root <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    nonRootVol <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    # tempVol <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    vol <-  matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    cumVol <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    respired_OM <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    years_per_cohort <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    omPackingDensity <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    mineralPackingDensity <- matrix(nrow=total_years+nInitCohorts, ncol=total_years)
    
    # reference_elevation <- rep(NA, total_years)
    
    # Set up annual time structure
    
    aboveground_biomass_species <- matrix(ncol = nSpecies, nrow = total_years)
    
    surface_elevation <- rep(NA, total_years)
    aboveground_biomass <- rep(NA, total_years)
    belowground_biomass <- rep(NA, total_years)
    totalRootVolume <- rep(NA, total_years)
    sediment_delivered <- rep(NA, total_years)
    rootToShoot_tracked  = rep(NA, total_years)
    rootTurnover_tracked = rep(NA, total_years)
    rootShape_tracked = rep(NA, total_years)
    abovegroundTurnover_tracked = rep(NA, total_years)
    rootDepthMax_tracked = rep(NA, total_years)
    rootPackingDensity_tracked = rep(NA, total_years)
    lambda_tracked <- rep(NA, total_years)
  }
  
  # Pre calculate values
  {
    # Preallocate memory for species
    zVegMin_up = rep(NA, nSpecies)
    zVegMax_low = rep(NA, nSpecies)
    a_up = rep(NA, nSpecies)
    b_up = rep(NA, nSpecies)
    c_up = rep(NA, nSpecies)
    a_low = rep(NA, nSpecies)
    b_low = rep(NA, nSpecies)
    c_low = rep(NA, nSpecies)
    lambda = rep(NA, nSpecies)
    
    # Pre calculate inputs from parameters
    # For the curve applied on the 'upper end' create a new minimum elevation
    # mirroring  upper elevation limit accross the peak biomass elevation
    
    for (species in 1:nSpecies) {
      
      zVegMin_up[species] <- zVegPeak[species]-((zVegMax[species]-zVegPeak[species])) 
      
      # For the curve applied at the 'lower end', same. Create new maximum elevation 
      # tolerance mirroring, lower elevation limit accros the peak biomass elevation.
      zVegMax_low[species] <-zVegPeak[species]+((zVegPeak[species]-zVegMin[species]))
      
      # Solve for the parameters of the upper curve.
      a_up[species] <- -((-zVegMin_up[species] * bMax[species] - zVegMax[species] * bMax[species]) / ((zVegMin_up[species] - zVegPeak[species]) * (-zVegMax[species] + zVegPeak[species])))
      b_up[species] <- -(bMax[species] / ((zVegMin_up[species] - zVegPeak[species]) * (-zVegMax[species] + zVegPeak[species])))
      c_up[species] <- (zVegMin_up[species] * zVegMax[species] * bMax[species]) / ((zVegMin_up[species] - zVegPeak[species]) * (zVegMax[species] - zVegPeak[species]))
      
      # Solve for the parametrs of the lower curve.
      a_low[species] <- -((-zVegMin[species] * bMax[species] - zVegMax_low[species] * bMax[species]) / ((zVegMin[species]- zVegPeak[species]) * (-zVegMax_low[species] + zVegPeak[species])))
      b_low[species] <- -(bMax[species] / ((zVegMin[species] - zVegPeak[species]) * (-zVegMax_low[species] + zVegPeak[species])))
      c_low[species] <- (zVegMin[species] * zVegMax_low[species] * bMax[species]) / ((zVegMin[species] - zVegPeak[species]) * (zVegMax_low[species] - zVegPeak[species]))
      
      # In case  we have exponential root distribution, pre calculate these
      lambda[species] <- -log(0.05) / rootDepthMax[species]
      
    }
    
    # Convert annual decay rate to decay constant
    k_fast<- -log(1-omDecayRateFast)
    k_slow<- -log(1-omDecayRateSlow)
    
  }
  
  # Initialize cohort and year 1 as a 1 m block of mineral sediment
  # Years are rows, depth increments are columns
  if (run_spinup == 1) {
    
    mineral[1,1] <- meanMineralPackingDensity*initial_cohort_vol
    # Set root, OM to zero
    fast[1,1] <- 0
    slow[1,1] <- 0
    root[1,1] <- 0
    nonRootVol[1,1] <- initial_cohort_vol
    # tempVol[1,1] <- 100
    # calculate cumulative volume
    cumVol[1,1] <- initial_cohort_vol
    vol[1,1] <- initial_cohort_vol
    omPackingDensity[1,1] <- meanOmPackingDensity
    mineralPackingDensity[1,1] <- meanMineralPackingDensity
    
    # Set accretion as n years as very high (12,000) IDK between LGM and now?
    years_per_cohort[1,1] <- initial_cohort_years
    # reference_elevation[1] <- 0
    
    cohort_year[1,1] <- simulation_years[1]
    cohort_index[1,1] <- 1
    
    aboveground_biomass_species[1, 1:nSpecies] <- rep(0, nSpecies)
    
    aboveground_biomass[1] <- 0
    belowground_biomass[1] <- 0
    
    surface_elevation[1] <- initElv
    
    totalRootVolume[1] <- 0
    sediment_delivered[1] <- 0
    rootToShoot_tracked[1]  = 0
    rootTurnover_tracked[1] = 0
    abovegroundTurnover_tracked[1] = 0
    rootDepthMax_tracked[1] = 0
    # !!! May have to switch these values, I just added non zeros so they wouldn't break any math
    # Need to double check that that's working
    rootPackingDensity_tracked[1] = 0.07
    lambda_tracked[1] = 1
    rootShape_tracked[1] = 0
    
    
  } else { # If we do not run a spinup, carry over previous cohorts
    
    cohort_year[1:nInitCohorts,1] <- initial_cohorts$year
    cohort_index[1:nInitCohorts,1] <- initial_cohorts$cohort_index
    
    mineral[1:nInitCohorts,1] <- initial_cohorts$mineral
    # Set root, OM to zero
    fast[1:nInitCohorts,1] <- initial_cohorts$fast
    slow[1:nInitCohorts,1] <- initial_cohorts$slow
    root[1:nInitCohorts,1] <- initial_cohorts$root
    nonRootVol[1:nInitCohorts,1] <- initial_cohorts$nonRootVol
    # tempVol[1,1] <- 100
    # calculate cumulative volume
    cumVol[1:nInitCohorts,1] <- initial_cohorts$cumVol
    vol[1:nInitCohorts,1] <- initial_cohorts$vol
    
    omPackingDensity[1:nInitCohorts,1] <- initial_cohorts$omPackingDensity
    mineralPackingDensity[1:nInitCohorts,1] <- initial_cohorts$mineralPackingDensity
    
    # Set accretion as n years as very high (12,000) IDK between LGM and now?
    years_per_cohort[1:nInitCohorts,1] <- initial_cohorts$inputYrs
    
    # All of these inputs to read in a previous sceanario table
    n_previous_iter <- nrow(initial_scenario)
    n_previous_species <- ncol(initial_species)-1
    
    aboveground_biomass_species[1, 1:n_previous_species] <- t(initial_species[n_previous_iter,2:(n_previous_species+1)])
    
    aboveground_biomass[1] <- initial_scenario$aboveground_biomass[n_previous_iter]
    belowground_biomass[1] <- initial_scenario$belowground_biomass[n_previous_iter]
    surface_elevation[1] <- initial_scenario$surface_elevation[n_previous_iter]
    totalRootVolume[1] <- initial_scenario$totalRootVolume[n_previous_iter]
    sediment_delivered[1] <- initial_scenario$sediment_delivered[n_previous_iter]
    rootToShoot_tracked[1]  = initial_scenario$rootToShoot[n_previous_iter]
    rootTurnover_tracked[1] = initial_scenario$rootTurnover[n_previous_iter]
    abovegroundTurnover_tracked[1] = initial_scenario$abovegroundTurnover[n_previous_iter]
    rootDepthMax_tracked[1] = initial_scenario$rootDepthMax[n_previous_iter]
    rootPackingDensity_tracked[1] = initial_scenario$rootPackingDensity[n_previous_iter]
    lambda_tracked[1] = initial_scenario$lambda[n_previous_iter]
    rootShape_tracked[1] = initial_scenario$rootShape[n_previous_iter]
    
  }
  
  year_index <- cohort_index[1,1]:(cohort_index[1,1]+total_years-1)
  
  # Two loops
  for (year in 2:total_years) {
    
    # IF
    # It is either a dynamic step
    # Or it is a spinup and year == 2 
    # We have a list of things to calculate
    cohort_index[1,year] <- year_index[year]
    
    if (is_spinup[year] == 0 | (run_spinup == 1 & year == 2)) {
      
      # Determine reference elevation
      # Reference elevation
      
      # If it is either the switchover from a spinup to dynamic 
      # Or there is no spinup and this is year 2
      if ((run_spinup == 1 & (is_spinup[year] != is_spinup[year-1])) |
          (run_spinup == 0 & year == 2)) { 
        
        # Then calculate reference elevation,
        # Initial elevation = reference elevaation + cumulative volume of the cohort stack
        # So rearranged this is reference elevation = init elevation - cumVol
        reference_elevation <- (surface_elevation[year-1] - cumVol[year_index[year-1], year - 1])
        
      }
      
      # Aboveground biomass
      # Calculate z_star
      weighted_mhw <- 0
      flood_sum <- 0
      for (flood in 1:nFloods) {
        flood_sum = flood_sum + flood_frequency[flood];
      }
      
      for (flood in 1:nFloods) {
        weighted_mhw = weighted_mhw + mhwMat[scenario_i[year],flood] * flood_frequency[flood] / flood_sum;
      }
      
      z_star <- (surface_elevation[year-1] - msl[scenario_i[year]]) / (weighted_mhw - msl[scenario_i[year]])
      
      # Surface biomass addition
      potential_agb <- rep(0, nSpecies)
      max_agb <- 0
      sum_agb <- 0
      
      for (species in 1:nSpecies) {
        
        # If elevation is above the specified peak biomass elevation, apply the upper curve,
        # if it's under apply the lower curve
        agb_above_peak <-  a_up[species]*z_star + b_up[species]*z_star^2 + c_up[species]
        agb_below_peak <-  a_low[species]*z_star + b_low[species]*z_star^2 + c_low[species]
        
        is_above_peak <- as.numeric(z_star > zVegPeak[species])
        
        uncensored_agb <- agb_above_peak*is_above_peak + (1-is_above_peak)*agb_below_peak
        
        # Censor aboveground biomass so it can't be negative
        is_gt_zero <- as.numeric(uncensored_agb > 0)
        
        potential_agb[species] <- is_gt_zero * uncensored_agb
        max_agb <- max(max_agb, potential_agb[species])
        sum_agb <- sum_agb + potential_agb[species] 
        
      }
      
      
      # What are all of the biological parameters we need?
      # total_aboveground_biomass
      
      # If biomass values are
      if (max_agb != 0) {
        
        agb_temp <- 0
        bgb_temp <- 0 
        totalRootVolume_temp <- 0
        rootToShoot_temp <- 0
        bgTurnover_temp <- 0
        agTurnover_temp <- 0
        rootDepth_temp <- 0
        rootPackingDensity_temp <- 0
        lambda_temp <- 0
        rootShape_temp <- 0
        
        for (species in 1:nSpecies) {
          
          # winner take all
          if (competition_function == 0) {
            
            agb_weight <- as.numeric(potential_agb[species] == max_agb)
            
          } else { # coexistence
            
            agb_weight <- potential_agb[species]/sum_agb
            
          } # End if else to calculation weights 
          
          # Calculate weighted averages
          aboveground_biomass_species[year, species] <- max_agb * agb_weight
          
          agb_temp <- agb_temp + max_agb * agb_weight
          bgb_temp <- bgb_temp +  aboveground_biomass_species[year, species] * rootToShoot[species] 
          totalRootVolume_temp <- totalRootVolume_temp +  aboveground_biomass_species[year, species] * rootToShoot[species] / rootPackingDensity[species]
          rootToShoot_temp <- rootToShoot_temp + rootToShoot[species] * agb_weight
          bgTurnover_temp <- bgTurnover_temp + rootTurnover[species] * agb_weight
          agTurnover_temp <- agTurnover_temp + abovegroundTurnover[species] * agb_weight
          rootDepth_temp <- rootDepth_temp + rootDepthMax[species] * agb_weight
          rootPackingDensity_temp <- rootPackingDensity_temp + rootPackingDensity[species] * agb_weight
          lambda_temp = lambda_temp + lambda[species] * agb_weight
          rootShape_temp <- rootShape_temp + rootShape[species] * agb_weight
          
        }
        
        
      } else {
        
        aboveground_biomass_species[year, ] <- rep(0, nSpecies)
        
        agb_temp <- 0
        bgb_temp <- 0 
        totalRootVolume_temp <- 0
        rootToShoot_temp <- rootToShoot[1]
        bgTurnover_temp <- rootTurnover[1]
        agTurnover_temp <- abovegroundTurnover[1]
        rootDepth_temp <- rootDepthMax[1]
        rootPackingDensity_temp <- rootPackingDensity[1]
        lambda_temp <- lambda[1]
        rootShape_temp <- 0
        
      }  # End of check to make sure agb is positive
      
      # Add weighted averages to storage vectors
      aboveground_biomass[year] <- agb_temp
      belowground_biomass[year] <- bgb_temp
      totalRootVolume[year] <- totalRootVolume_temp
      rootToShoot_tracked[year]  = rootToShoot_temp
      rootTurnover_tracked[year] = bgTurnover_temp
      abovegroundTurnover_tracked[year] = agTurnover_temp
      rootDepthMax_tracked[year] = rootDepth_temp
      rootPackingDensity_tracked[year] = rootPackingDensity_temp
      lambda_tracked[year] = lambda_temp
      rootShape_tracked[year] = round(rootShape_temp)
      
      # Things we need for root calculations that only need to be calculated once per profile
      rootWidth <- totalRootVolume[year]*2/(rootDepthMax_tracked[year])
      
      # Calculte mineral type delivered  
      tmp_mineral <- 0 
      
      for (flood in 1:nFloods) {
        
        # Flooding time expressed as fractional time, all, or expressed as a linear function with flood depth,
        # .. we have a trig flood time function we could swap out here in the future
        is_flooded <- as.numeric(surface_elevation[year-1] <= mhwMat[scenario_i[year],flood])
        
        # Flood count is either all or nothing depending on elevation
        flood_count <-flood_frequency[flood] * is_flooded
        
        # Tidal height relative to surface
        flood_depth <- (mhwMat[scenario_i[year],flood]-surface_elevation[year-1])*0.5
        
        flood_time <- min(1, (mhwMat[scenario_i[year],flood]-surface_elevation[year-1])/(mhwMat[scenario_i[year],flood]-mlwMat[scenario_i[year], flood])) * is_flooded
        
        # Calculate fraction of sediment captured
        # if the sediment column IS NOT able to clear
        fraction_captured <- min( 
          # available suspendedSediment is total possible caputre 
          captureRate*flood_time, 
          # if the sediment column IS able to clear
          1)
        
        # Calculate available sediment as a cumulative block of water
        sediment_available<- suspendedSediment * flood_count * flood_depth
        
        # Calculated delivered sediment
        tmp_mineral <- tmp_mineral + sediment_available * fraction_captured
        
      } # End flood loop 
      
      # For each type of flood we're tracking
      sediment_delivered[year] <- tmp_mineral
      
    } else {
      
      # If not first spinup iteration or dynamic step, update with inputs with previous values 
      aboveground_biomass[year] <- aboveground_biomass[year-1]
      belowground_biomass[year] <- belowground_biomass[year-1]
      totalRootVolume[year] <- totalRootVolume[year-1]
      rootToShoot_tracked[year]  = rootToShoot_tracked[year-1]
      rootTurnover_tracked[year] = rootTurnover_tracked[year-1]
      abovegroundTurnover_tracked[year] = abovegroundTurnover_tracked[year-1] 
      rootDepthMax_tracked[year] = rootDepthMax_tracked[year-1]
      rootPackingDensity_tracked[year] = rootPackingDensity_tracked[year-1]
      lambda_tracked[year] =  lambda_tracked[year-1]
      rootShape_tracked[year] = rootShape_tracked[year-1]
      
      sediment_delivered[year] <-  sediment_delivered[year-1] 
      
      aboveground_biomass_species[year, ] <- aboveground_biomass_species[year-1,]
      
      # Species 
      
    } # End of defining surface inputs
    
    # Add surface inputs
    mineral[1, year] <- sediment_delivered[year]  * years_per_iteration[year]
    years_per_cohort[1, year] <-  years_per_iteration[year]
    
    # In the future calc decomp at surface
    total_agb_produced <- 0
    dropped_fast_pool <- 0
    dropped_slow_pool <- 0
    agb_respired <- 0
    
    # We simulate multiple drops but have the option to batch them into a single cohort
    for (drop in (1:years_per_iteration[year])) {
      
      agb_produced <- aboveground_biomass[year-1] * abovegroundTurnover_tracked[year-1]
      
      total_agb_produced <- total_agb_produced + agb_produced
      
      # Decomp
      agb_fast_added <- agb_produced * (1-recalcitrantFrac) * exp(-k_fast * (drop - 1))
      agb_slow_added <- agb_produced * recalcitrantFrac * exp(-k_slow * (drop - 1))
      
      # Update sum total
      dropped_fast_pool <- dropped_fast_pool +  agb_fast_added
      dropped_slow_pool <- dropped_slow_pool + agb_slow_added
      
      # Calculate agb lost
      agb_lost <-  (agb_produced * (1-recalcitrantFrac) * (1-exp(-k_fast * (drop - 1))))  + 
        (agb_produced * recalcitrantFrac * (1-exp(-k_slow * (drop - 1))))
      
      agb_respired <- agb_respired + agb_lost
      
      # Quick check
      # (agb_produced == (agb_fast_added + agb_slow_added + agb_lost))
      
    }
    
    # (total_agb_produced == dropped_fast_pool + dropped_slow_pool + agb_respired)
    
    # Batch any sub years
    fast[1,year] <- dropped_fast_pool
    
    # For slow do the same, alhtough we could make this the same as fast if we add capability for loosing slow pool OM, for example, outside of
    slow[1,year] <- dropped_slow_pool
    
    respired_OM[1,year] <- agb_respired
    
    # Some other stuff
    cumVol[1,year] <- 0
    nonRootVol[1,year] <- 0
    root[1,year] <- 0
    
    omPackingDensity[1,year] <- meanOmPackingDensity
    mineralPackingDensity[1,year] <- meanMineralPackingDensity
    
    vol[1,year] <- (fast[1,year] + slow[1,year]) / omPackingDensity[1,year] + 
      # mineral section of the packing density equation
      mineral[1,year]/ mineralPackingDensity[1,year]
    
    # This is a change, we don't actually need to track this, so I'm only writing to memory when we need it
    tempVol <- rep(NA, year)
    tempVol[1] <- 0
    
    # Pre define this for the case of exponential binary search algorithm
    possibleDepth.arr <- c(0, rep(NA, year-1))
    
    #What is the non rooting volumne down to rooting max
    nonRootVolumeToRootMax <- rootDepthMax_tracked[year] - totalRootVolume[year]
    
    # Check if roots are present
    has_roots <- as.numeric(aboveground_biomass[year] != 0)
    
    cohort_year[1,year] <- simulation_years[year]
    
    cohort_i <- year_index[year]:1
    
    for (cohort in 2:year_index[year]) {
      
      cohort_year[cohort,year] <- simulation_years[year]
      cohort_index[cohort,year] <- cohort_i[cohort]
      
      # Carry over values which don't change
      mineral[cohort, year] <- mineral[cohort-1, year-1]
      years_per_cohort[cohort, year] <- years_per_cohort[cohort-1, year-1]
      
      
      # Age organic cohorts, decompose the previous year's fast pool
      total_bgb_produced <- 0
      turned_over_fast_pool <- fast[cohort-1,year-1] * exp(-k_fast * years_per_iteration[year])
      turned_over_slow_pool <- slow[cohort-1,year-1] * exp(-k_slow * years_per_iteration[year])
      turnover_respiration <- (fast[cohort-1,year-1] * (1-exp(-k_fast * years_per_iteration[year]))) +
        (slow[cohort-1,year-1] * (1-exp(-k_slow * years_per_iteration[year])))
      
      # We simulate single turnover events but have the option to batch them into a single cohort
      # if we are running more than one year per cohort
      # !!! triple check all this. make sure prod = soil om + repiration
      for (turnover in (1:years_per_iteration[year])) {
        
        bgb_produced <- root[cohort-1, year-1] * rootTurnover_tracked[year-1]
        total_bgb_produced <- total_bgb_produced + bgb_produced
        
        # Decomp
        bgb_fast_added <- bgb_produced * (1-recalcitrantFrac) * exp(-k_fast * (turnover - 1))
        bgb_slow_added <- bgb_produced * recalcitrantFrac * exp(-k_slow * (turnover - 1))
        
        # Update sum total
        turned_over_fast_pool <- turned_over_fast_pool +  bgb_fast_added
        turned_over_slow_pool <- turned_over_slow_pool + bgb_slow_added
        
        # Calculate bgb lost
        bgb_lost <-  (bgb_produced * (1-recalcitrantFrac) * (1-exp(-k_fast * (turnover - 1))))  + 
          (bgb_produced * recalcitrantFrac * (1-exp(-k_slow * (turnover - 1))))
        
        turnover_respiration <- turnover_respiration + bgb_lost
        
        # Quick check
        # (bgb_produced == (bgb_fast_added + bgb_slow_added + bgb_lost))
        
      }
      
      # if (total_bgb_produced != turnover_respiration + turned_over_fast_pool + turned_over_slow_pool) {
      #   stop("Total bgb produced is not equal to turnover and carbon loss.")
      # }
      
      # Add fast pool cohort
      fast[cohort,year] <- turned_over_fast_pool
      
      # ... and slow pool.
      slow[cohort,year] <- turned_over_slow_pool
      
      respired_OM[cohort,year] <- turnover_respiration
      
      # Update packing densities
      # ... this may seem extraneous, but we actually need to have separately tracked OM and Min packing densities in case we want to 
      # ... change them manually when simulating real soils and integrating data with the model.
      omPackingDensity[cohort,year] <- meanOmPackingDensity
      mineralPackingDensity[cohort,year] <- meanMineralPackingDensity
      
      # Calculate non root volume, this is cumulative with depth, so we're adding to the next one up
      nonRootVol[cohort,year] <- nonRootVol[cohort-1,year] + 
        # organic section of the packing density equation
        (fast[cohort,year] + slow[cohort,year])/ omPackingDensity[cohort,year] + 
        # mineral section of the packing density equation
        mineral[cohort,year]/ mineralPackingDensity[cohort,year] 
      
      if (has_roots == 1) {
        
        # Updating here
        if (rootShape_tracked[year] == 0) {
          
          # Root width and coefficents simplify equation for calculating new volume with roots
          coef1 <- rootWidth / (2*rootDepthMax_tracked[year])
          coef2 <- 1-rootWidth
          coef3 <- -nonRootVol[cohort,year]
          
          # Calculate new volume with roots; if root volume is present
          if (totalRootVolume[year] != 0) {
            tempVol[cohort] =   (-coef2 + sqrt(coef2^2-4*coef1*coef3))/(2*coef1)
          } else {
            tempVol[cohort] = nonRootVol[cohort, year]
          }
          
          # Add the roots section
          #mass = integral(mass_per_depth, depth)
          # The minimum and maximum layer should not be deeper than root depth max
          # so that the equation returns zero when lower than roots
          layer_bottom <- min(rootDepthMax_tracked[year], tempVol[cohort])
          layer_top <- min(rootDepthMax_tracked[year], tempVol[cohort-1])
          
          # slope and intercept simplify the add roots equation
          slope <- -2 * belowground_biomass[year] / (rootDepthMax_tracked[year]^2)
          intercept <- 2 * belowground_biomass[year] / rootDepthMax_tracked[year]
          
          root_ans <- intercept*(layer_bottom-layer_top) + 
            slope/2*(layer_bottom^2 - layer_top^2)
          
        } else if (rootShape_tracked[year] == 1) {
          
          # The goal of this code chunk is to determine how much of an censored exponential root section
          # ... goes into an unfilled cohort of soil. 
          
          # House keeping first
          # Assign the previous depth as zero
          previousDepth <- possibleDepth.arr[cohort - 1]
          
          # Take care of a couple of special conditions first
          
          # First if the volume of the cohort is very small, lump it in with the previous
          if (nonRootVol[cohort,year] < relTol) {
            
            possibleDepth <- previousDepth
            
            #If we overfill the root zone volume
          } else if (nonRootVolumeToRootMax <= nonRootVol[cohort,year]) {
            
            #Subtract the root zone volume, find the depth beyond and add to the max rooting depth
            possibleDepth <- (nonRootVol[cohort,year]-nonRootVolumeToRootMax) + rootDepthMax_tracked[year]
            
          } else {
            
            # Implement binary search algorithm
            increment <- rootDepthMax_tracked[year] - previousDepth
            possibleDepth <- previousDepth + increment/2
            
            for (ii in 1:maxIter){
              
              effective_depth <- min(possibleDepth, rootDepthMax_tracked[year])
              rootVolume <- totalRootVolume[year] * (1 - exp(-lambda_tracked[year] * effective_depth))/0.95
              
              nonRootVolume_test <- possibleDepth - rootVolume
              
              if ((abs(nonRootVol[cohort,year]-nonRootVolume_test) / nonRootVol[cohort,year]) < relTol) {
                break
              } else {
                increment <- increment / 2
                #Should we go up ...
                if(nonRootVol[cohort,year] > nonRootVolume_test){
                  possibleDepth <- possibleDepth + increment / 2
                }else{ # ... or down?
                  possibleDepth <- possibleDepth - increment / 2
                  
                } #if-else up/down
                
              } #if-else relTol
              
            } #if else statement refining search
            
          } #if else statement for which case to apply
          
          # Calculate root mass  
          layer_top <- min(previousDepth, rootDepthMax_tracked[year])
          layer_bottom <- min(possibleDepth, rootDepthMax_tracked[year])
          
          root_ans <- belowground_biomass[year] * (exp(-lambda_tracked[year] * layer_top) - exp(-lambda_tracked[year] * layer_bottom))/0.95
          
          # Update tracker vector
          possibleDepth.arr[cohort] <- possibleDepth  
          
        } else { # Root shape is invalide
          
          stop("Invalid root shape specified")
          
        }
        
      } else {
        
        # If roots are not present
        root_ans <- 0
        
      } # end of check to see if roots are present
      
      # Calculate the new root mass
      #mass = integral(mass_per_depth, depth)
      root[cohort,year]  <- root_ans
      
      # calculate the volume of the individual cohort
      # organic section of the packing density equation
      
      # Changed this with specific root density
      vol[cohort,year] <- (fast[cohort,year] + slow[cohort,year])/omPackingDensity[cohort,year] + 
        root[cohort,year]/rootPackingDensity_tracked[year] + 
        # mineral section of the packing density equation
        mineral[cohort,year]/mineralPackingDensity[cohort,year]
      
      # Calculate new cumulative volume
      cumVol[cohort,year] <- cumVol[cohort-1,year] + vol[cohort,year]
      
    } 
    
    # Update surface elevation
    if (is_spinup[year] == 1) {
      surface_elevation[year] <- initElv
    } else {
      # If it is not a spinup iteration, Determine surface elevation
      surface_elevation[year] <- (reference_elevation + cumVol[year_index[year], year]) 
    }
    
  }
  
  # nr <- nrow(cumVol)
  # nc <- ncol(cumVol)
  
  # index <- expand.grid(
  #   cohort = seq_len(nr),
  #   time   = seq_len(nc)
  # )
  
  long_df <- data.frame(
    year = as.vector(cohort_year),
    cohort_index = as.vector(cohort_index),
    mineral     = as.vector(mineral),
    fast     = as.vector(fast),
    slow     = as.vector(slow),
    root = as.vector(root),
    vol = as.vector(vol),
    nonRootVol = as.vector(nonRootVol),
    cumVol = as.vector(cumVol),
    inputYrs = as.vector(years_per_cohort),
    omPackingDensity = as.vector(omPackingDensity),
    mineralPackingDensity = as.vector(mineralPackingDensity)
  ) %>% 
    filter(complete.cases(.))
  
  aboveground_biomass_species <- as.data.frame(aboveground_biomass_species)
  names(aboveground_biomass_species) <- species_codes
  aboveground_biomass_species <- bind_cols(data.frame(year = simulation_years),
                                           aboveground_biomass_species)
  
  return(list(cohorts = long_df,
              scenario = data.frame(year = simulation_years,
                                    years_per_iteration = years_per_iteration,
                                    surface_elevation = surface_elevation, 
                                    aboveground_biomass = aboveground_biomass,
                                    belowground_biomass = belowground_biomass,
                                    totalRootVolume = totalRootVolume,
                                    sediment_delivered = sediment_delivered,
                                    rootToShoot = rootToShoot_tracked,
                                    rootTurnover = rootTurnover_tracked,
                                    abovegroundTurnover = abovegroundTurnover_tracked,
                                    rootDepthMax = rootDepthMax_tracked,
                                    rootPackingDensity = rootPackingDensity_tracked,
                                    lambda = lambda_tracked,
                                    rootShape = rootShape_tracked),
              species = aboveground_biomass_species))
  
}
