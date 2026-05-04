#' Simulate radioisotop dates
#' @param core dataframe, output from build soil core
#' @param pb210_influx_Bq_m2_yr numeric, atmospheric flux in Bequerels per kilogram per square meter
#' @param pb210_supported_Bq_kg numeric, Concnetration of supported 210Pb in Bequerels per gram
#' @param cs137_scenario data frame with a year column indicating calendar years, and a second column indicating cs137 influx in Bequerels per square cm 
#' @param cs137_anomaly numeric, fractional anomoly inidicating site offset from regional scenario
#' @param cs137_baseline_Bq_m2_yr numeric, baseline fallout of cs137 without aboveground nuclear wepons testing. Defaults to the minimum value in the 137Cs fallout scenario
#' @param cs137_erosion_loss numeric, bound between 1 and 0, how much of the annual influx of 137Cs is lost to erosion
#' @param pb210_erosion_loss numeric, bound between 1 and 0, how much of the annual influx of 210Pb is lost to erosion
#' @param cs137_mobility_factor integer or zero, indicates how many adjacent depth increments will exchange 137Cs when simulating mobility
#' @param pb210_mobility_factorinteger or zero, indicates how many adjacent depth increments will exchange 210Pb when simulating mobility
#' 
#' @return a data frame with the same structure as core input, but with additional  
#' 
#' @export
#' 
simulateRadioisotopeDates <- function(core, pb210_influx_Bq_m2_yr = 200,
                                      pb210_supported_Bq_kg = 19.8,
                                      cs137_scenario = "cs137_N_non_europe",
                                      cs137_anomaly = 0.6224604,
                                      cs137_baseline_Bq_m2_yr = NA,
                                      cs137_erosion_loss = 0.5,
                                      pb210_erosion_loss = 0.01,
                                      cs137_mobility_factor = 3,
                                      pb210_mobility_factor = 2
                                      ) {
  
  # Constant decay rate
  pb210_decay_rate_yr <- 0.03122285
  cs137_decay_rate_yr <- 0.023
 
  scenario_path <- system.file(
    "extdata", "fallout",
    paste0(cs137_scenario, ".rds"),
    package = "rCMEM"
  )
  
  # Laod up 137Cs curve
  if (file.exists(scenario_path)) {
 
    cs137_scenario_file <- readRDS(scenario_path)
    
    if (! all(c("cs137_influx_Bq_m2_yr", "year") %in% names(cs137_scenario_file))) {
      
      stop("Cs 137 Scenario File Needs to at least two columns labeled cs137_influx_Bq_m2_yr and year")
    }
    
  } else {
    stop("137 Cs scenario file not recognized")
  }
  
  # Conversions and anomalies
  # pb210_influx_Bq_cm2_yr <- pb210_influx_Bq_m2_yr / 10000
  # 
  # pb210_supported_Bq_g <- pb210_supported_Bq_kg / 1000
  # 
  cs137_influx_Bq_m2_yr <- rev(cs137_scenario_file$cs137_influx_Bq_m2_yr)
  

  # cs137_influx_Bq_cm2_yr <- cs137_influx_Bq_m2_yr / 10000
  
  if (is.na(cs137_baseline_Bq_m2_yr)) {
    
    cs137_baseline_Bq_m2_yr <- min(cs137_influx_Bq_m2_yr)
    
  }
  
  cs137_year <- rev(cs137_scenario_file$year)

  # Read in core variables
  # And calculate bulk density
  core_calendar_age <- core$calendar_age
  core_age_from_surface <- core$age_from_surface
  core_fast <- core$fast
  core_slow <- core$slow
  core_mineral <- core$mineral
  core_root <- core$root
  core_vol <- core$vol
  core_cumVol <- core$cumVol
  
  min_core_year <- min(core_calendar_age)
  max_core_year <- core_calendar_age[1]+core_age_from_surface[1]
  
  cs137_year_max <- cs137_year
  cs137_year_min <- cs137_year - 1
  
  min_cs137_year <- min(cs137_year)
  max_cs137_year <- max(cs137_year)
  
  # If core min and max are older or younger than our cs 137 record
  if (max_core_year > max_cs137_year) {
    
    # Then fill in with bacground values
    cs137_influx_Bq_m2_yr <- c(cs137_baseline_Bq_m2_yr, cs137_influx_Bq_m2_yr)
    
    cs137_year_max <- c(max_core_year, cs137_year_max)
    cs137_year_min <- c(max_cs137_year, cs137_year_min)
    
  }
  
  if (min_core_year < min_cs137_year) {
    
    # Then fill in with bacground values
    cs137_influx_Bq_m2_yr <- c(cs137_influx_Bq_m2_yr, cs137_baseline_Bq_m2_yr)
    
    cs137_year_max <- c(cs137_year_max, min_cs137_year)
    cs137_year_min <- c(cs137_year_min, min_core_year)
    
  }
  
  # Take into account erosional losses and anomalies
  pb210_influx_Bq_m2_yr <- pb210_influx_Bq_m2_yr - (pb210_influx_Bq_m2_yr * pb210_erosion_loss)
  
  cs137_influx_Bq_m2_yr <- cs137_influx_Bq_m2_yr + (cs137_influx_Bq_m2_yr * cs137_anomaly)
  
  cs137_influx_Bq_m2_yr <- cs137_influx_Bq_m2_yr - (cs137_influx_Bq_m2_yr * cs137_erosion_loss)
  
  # Iterate through depth increments
  start_year_n <- 1
  
  # Preallocate memory
  core_dbd_g_cc <- rep(NA, nrow(core))

  core_supported_pb210_Bq_kg <- rep(pb210_supported_Bq_kg, nrow(core))
  core_unsupported_ideal_pb210_Bq_kg <- rep(NA, nrow(core))
  core_total_ideal_pb210_Bq_kg <- rep(NA, nrow(core))
  core_cs137_influx_Bq_m2_yr <- rep(NA, nrow(core))
  core_cs137_ideal_Bq_kg <- rep(NA, nrow(core))
  
  for (depth_increment in 1:nrow(core)) {
    
    if (depth_increment == 1) {
      increment_age_min = 0
      increment_calendar_age_max = core_calendar_age[1] + core_age_from_surface[1]
    } else {
      increment_age_min = core_age_from_surface[depth_increment-1]
      increment_calendar_age_max = core_calendar_age[depth_increment-1]
    } # end of if else for determining minimum age
     
    # Calculate increment mass
    core_dbd_g_cc[depth_increment] <- (core_fast[depth_increment] + core_slow[depth_increment] + core_root[depth_increment] + core_mineral[depth_increment]) / core_vol[depth_increment]
    
    # Calculate suppported and unsupported 210Pb
    pb_unsupported_Bq_m2 <- pb210_influx_Bq_m2_yr / pb210_decay_rate_yr * (exp(-pb210_decay_rate_yr*increment_age_min)-exp(-pb210_decay_rate_yr*core_age_from_surface[depth_increment]))
    
    # Iterate through 137Cs calendar years
    # volume is technically cm since we always assume have 1 cm2 width x length
    # g/cm3 * 1,000,000 cm3/1m3 * 1 kg/1,000 g = 1,000 kg/m3 per 1 g/cm3 
    # Divide by 1,000 to go from g/cm3 to 
    
    # Bq/m2 / kg/m2 = Bq/kg
    core_unsupported_ideal_pb210_Bq_kg[depth_increment] <- pb_unsupported_Bq_m2 / (core_dbd_g_cc[depth_increment] * core_vol[depth_increment] * 10)
    
    # Calculate weighted average 137Cs influx for the time span
    core_total_ideal_pb210_Bq_kg[depth_increment] <- core_supported_pb210_Bq_kg[depth_increment] + core_unsupported_ideal_pb210_Bq_kg[depth_increment]
    
    weighted_cs137_influx = 0
    total_weights = 0
    for (calendar_year in start_year_n:length(cs137_year_min)) {
      
      # Is this a zero vol cohort?
      weight_calc_part1 = as.numeric(increment_calendar_age_max == core_calendar_age[calendar_year])
      
      weight_calc_part2 = min(cs137_year_max[calendar_year], increment_calendar_age_max, na.rm=T) -
        max(cs137_year_min[calendar_year], core_calendar_age[depth_increment], na.rm=T)
      
      # IF we pass the bottom of matching depth increments, break the loop so we don't have to do more calculations
      if (weight_calc_part2 < 0) {
        # Update cohort_n_start
        break()
      } # end of break condition
      
      # Depth weight matrix
      # Is the increment partially in or totallly in?
      depth_weight <- max(weight_calc_part1, max(weight_calc_part2, 0, na.rm=T) / (cs137_year_max[calendar_year]-cs137_year_min[calendar_year]), na.rm=T)
      
      weighted_cs137_influx = weighted_cs137_influx + cs137_influx_Bq_m2_yr[calendar_year] * depth_weight
      
      total_weights = total_weights + depth_weight 
      
      start_year_n <- calendar_year
      
      } # end of cs 137 flux iteration
    
    core_cs137_influx_Bq_m2_yr[depth_increment] <- weighted_cs137_influx / total_weights
    
    core_influx_Bq_m2 <- core_cs137_influx_Bq_m2_yr[depth_increment]/ cs137_decay_rate_yr * (exp(-cs137_decay_rate_yr*increment_age_min)-exp(-cs137_decay_rate_yr*core_age_from_surface[depth_increment])) 
    
    core_cs137_ideal_Bq_kg[depth_increment] <- core_influx_Bq_m2 / (core_dbd_g_cc[depth_increment] * core_vol[depth_increment] * 10)
      
    
    } # End of depth increment iteration loop
    
    # Simulate mixing with geometric WA
    # Preallocate memory
    core_unsupported_actual_pb210_Bq_kg <- rep(NA, nrow(core))
    core_total_actual_pb210_Bq_kg <- rep(NA, nrow(core))
    core_cs137_actual_Bq_kg <- rep(NA, nrow(core))
    
    
    if ( any(c(cs137_mobility_factor, pb210_mobility_factor) > 0)) {
      
      n <- nrow(core)
      
      for (depth_increment in 1:n) {
        
        # valid window indices
        dist <- core_cumVol - core_cumVol[depth_increment]
        
        # Gaussian weights
        w_pb210 <- exp(-(dist^2) / (2 * pb210_mobility_factor^2))
        
        # renormalize
        w_pb210 <- w_pb210 / sum(w_pb210)
        
        core_unsupported_actual_pb210_Bq_kg[depth_increment] <- sum(core_unsupported_ideal_pb210_Bq_kg * w_pb210)
        core_total_actual_pb210_Bq_kg[depth_increment] <- core_unsupported_actual_pb210_Bq_kg[depth_increment] + core_supported_pb210_Bq_kg[depth_increment]
        
        w_cs137 <- exp(-(dist^2) / (2 * cs137_mobility_factor^2))
      
        # renormalize
        w_cs137 <- w_cs137 / sum(w_cs137)
        
        core_cs137_actual_Bq_kg[depth_increment] <- sum(core_cs137_ideal_Bq_kg * w_cs137)
        
      } # end of depth iterating wa loop 
      
    } # end of if else mobility factor is greater than 0

  return(
    bind_cols(core,
                   data.frame(unsupported_pb210_ideal_Bq_kg = core_unsupported_ideal_pb210_Bq_kg,
                              supported_pb210_Bq_kg = core_supported_pb210_Bq_kg,
                              total_pb210_ideal_Bq_kg = core_total_ideal_pb210_Bq_kg,
                              cs137_influx_Bq_m2_yr = core_cs137_influx_Bq_m2_yr,
                              cs137_ideal_Bq_kg = core_cs137_ideal_Bq_kg,
                              unsupported_pb210_actual_Bq_kg = core_unsupported_actual_pb210_Bq_kg,
                              total_pb210_actual_Bq_kg = core_total_actual_pb210_Bq_kg,
                              cs137_actual_Bq_kg = core_cs137_actual_Bq_kg
                              )
                   )
         )
  
}
