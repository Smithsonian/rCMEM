# SSC Decline Demo

# Driver
{
  # SLR
  annapolis_gap_filled <- read_csv("inst/extdata/sea_level/annapolis_gap_filled.csv")
  
  loess_model <- loess(msl_cm~year, data = annapolis_gap_filled)
  annapolis_gap_filled$loess_msl <- predict(loess_model, data.frame(year = annapolis_gap_filled$year))
  
  meanSeaLevel = annapolis_gap_filled$loess_msl 
  meanSeaLevelDatum = -1.958671252
  meanHighWaterDatum = 5.521639027 
  meanHighHighWaterDatum = 19.81578852 
  meanHighHighWaterSpringDatum=26.25399172
  lunarNodalAmp = 0
  # suspendedSediment = 22 * 1e-06, # assume 75.8 mg/l in creek and interior is 1/3 to 1/4 creek
  
  meanHighWater <- meanHighWaterDatum-meanSeaLevelDatum
  meanHighHighWater <- meanHighHighWaterDatum-meanSeaLevelDatum
  meanHighHighWaterSpring <- meanHighHighWaterSpringDatum-meanSeaLevelDatum
  
  # suspendedSediment = 10 * 1e-06 # original
  nFloods <- 3
  flood_frequency <- c(0.5, 0.46497542, 0.03502458) * 705.79
  msl <- meanSeaLevel
  mhwMat <- matrix(c(msl+meanHighWater,msl+meanHighHighWater,msl+meanHighHighWaterSpring),
                   ncol = 3)
  mlwMat <- msl - (mhwMat-msl)
}

scenario_calendar_years <- 1928:2018
n_years <- length(scenario_calendar_years)

# Convert mg/L to g/cm3
ssc_scenario <- seq(
  from = 20,
  to = 5,
  length.out = n_years
) * 1e-6

# Run time-varying scenario
run_variable <- runCohortMem2(
  run_spinup = 1,
  run_scenario = 1,
  msl = msl,
  mhwMat = mhwMat,
  mlwMat = mlwMat,
  scenario_calendar_years = scenario_calendar_years,
  suspendedSediment = ssc_scenario
)

# Compare with constant sediment supply
run_constant <- runCohortMem2(
  run_spinup = 1,
  run_scenario = 1,
  msl = msl,
  mhwMat = mhwMat,
  mlwMat = mlwMat,
  scenario_calendar_years = scenario_calendar_years,
  suspendedSediment = 20e-6
)

output_df <- bind_rows(run_variable[[2]] %>% mutate(run_id = "variable SSC"),
                       run_constant[[2]] %>% mutate(run_id = "stable SSC")
                       )

ggplot(output_df, aes(x = year, y = sediment_delivered)) +
  geom_point(aes(color = run_id))

ggplot(output_df, aes(x = year, y = surface_elevation)) +
  geom_point(aes(color = run_id))


