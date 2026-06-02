# Generate some output files

library(tidyverse)
library(rCMEM)

# Inputs
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

a_mem_run <- runCohortMem2(run_spinup = 1, 
                          run_scenario = 1,
                          msl = msl,
                          mhwMat = mhwMat,
                          mlwMat = mlwMat,
                          spinup_iterations = 200, 
                          years_per_spinup_iter = 20,
                          rootShape = c(0),
                          omDecayRateFast = 0.99,
                          omDecayRateSlow = 0.0001,
                          bMax = c(0.0622),
                          rootTurnover = c(0.5),
                          initElv = -11,
                          recalcitrantFrac = 0.5,
                          rootDepthMax = c(40),
                          captureRate = 5,
                          suspendedSediment = 10 * 1e-06,
                          species_codes = "SCAM",
                          zVegMax = 1.8,
                          zVegPeak = 1.5, 
                          zVegMin = - 1.4,
                          abovegroundTurnover = 1.5,
                          rootPackingDensity = 0.07,
                          rootToShoot = 1.7
)

scenario <- a_mem_run[[2]]
cohorts <- a_mem_run[[1]]
species <- a_mem_run[[3]]

a_core <- buildSoilCore(cohorts = a_mem_run[[1]], coreYear = 2018, coreDepth = 200)

a_dates_4 <- simulateRadioisotopeDates(core = a_core, cs137_anomaly = 0.6, 
                                       pb210_mobility_factor = 5, 
                                       cs137_mobility_factor = 5,
                                       cs137_erosion_loss = 0.25, 
                                       pb210_erosion_loss = 0.25)

write_csv(scenario, "vignettes/outputs/scenario_example.csv")
write_csv(cohorts, "vignettes/outputs/cohorts_example.csv")
write_csv(a_core, "vignettes/outputs/core_example.csv")
write_csv(a_dates_4, "vignettes/outputs/dates_example.csv")
write_csv(annapolis_gap_filled, "vignettes/outputs/msl_example.csv")
write_csv(as.data.frame(mhwMat), "vignettes/outputs/mhw_example.csv")
write_csv(as.data.frame(mlwMat), "vignettes/outputs/mlw_example.csv")
write_csv(species, "vignettes/outputs/species_example.csv")



