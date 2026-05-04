# 

library(tidyverse)

# Driver
{
  # SLR
  annapolis_gap_filled <- read_csv("data/sea_level/annapolis_gap_filled.csv")
  
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

scenario_1 <- runCohortMem2(run_spinup = 1, 
                                  run_scenario = 1,
                                  msl = msl,
                                  mhwMat = mhwMat,
                                  mlwMat = mlwMat,
                                  spinup_iterations = 200, 
                                  years_per_spinup_iter = 1,
                                  rootShape = c(1,1))[[2]] %>% 
  mutate(spinup_settings = "1 year iter, 200 iterations")

scenario_5 <- runCohortMem2(run_spinup = 1, 
                                  run_scenario = 1,
                                  msl = msl,
                                  mhwMat = mhwMat,
                                  mlwMat = mlwMat,
                                  spinup_iterations = 40, 
                                  years_per_spinup_iter = 5,
                                  rootShape = c(1,1))[[2]] %>% 
  mutate(spinup_settings = "5 year iter, 40 iterations")

scenario_10 <- runCohortMem2(run_spinup = 1, 
                                   run_scenario = 1,
                                   msl = msl,
                                   mhwMat = mhwMat,
                                   mlwMat = mlwMat,
                                   spinup_iterations = 20, 
                                   years_per_spinup_iter = 10,
                                   rootShape = c(1,1)) [[2]] %>% 
  mutate(spinup_settings = "10 year iter, 20 iterations")


species_10a <- runCohortMem2(run_spinup = 1, 
                                   run_scenario = 1,
                                   msl = msl,
                                   mhwMat = mhwMat,
                                   mlwMat = mlwMat,
                                   spinup_iterations = 20, 
                                   years_per_spinup_iter = 10,
                                   competition_function = 0,
                                   rootShape = c(1,1))[[3]] %>% 
  as.data.frame() %>% 
  mutate(spinup_settings = "Zero-sum")

species_10b <- runCohortMem2(run_spinup = 1, 
                                   run_scenario = 1,
                                   msl = msl,
                                   mhwMat = mhwMat,
                                   mlwMat = mlwMat,
                                   spinup_iterations = 20, 
                                   years_per_spinup_iter = 10,
                                   competition_function = 1,
                                   rootShape = c(1,1)) [[3]] %>% 
  as.data.frame() %>% 
  mutate(spinup_settings = "Co-existance")


scenarios <- scenario_1 %>% bind_rows(scenario_5) %>% bind_rows(scenario_10)

ggplot(scenarios, aes(x = year, y = surface_elevation)) +
  geom_point(aes(color = spinup_settings, shape = spinup_settings)) + 
  scale_shape_manual(values = c(1,3,4))

ggplot(scenarios, aes(x = year, y = aboveground_biomass)) +
  geom_point(aes(color = spinup_settings, shape = spinup_settings)) + 
  scale_shape_manual(values = c(1,3,4))

ggplot(scenarios, aes(x = year, y = belowground_biomass)) +
  geom_point(aes(color = spinup_settings, shape = spinup_settings)) + 
  scale_shape_manual(values = c(1,3,4))

species_codes <- c("SCAM", "SPPA")

species_plot <- species_10b %>% 
  bind_rows(species_10a) %>% 
  gather(key = "species", value = "aboveground_biomass", -c(year, spinup_settings)) %>%
  mutate(species = factor(species, levels = rev(species_codes))) %>% 
  group_by(spinup_settings, year) %>% 
  mutate(cumulative_biomass = cumsum(aboveground_biomass))

ggplot(species_plot, aes(x = year, y = cumulative_biomass)) +
  geom_ribbon(aes(ymax = cumulative_biomass, ymin = 0, fill = species)) +
  facet_wrap(.~spinup_settings)

mem_1 <- runCohortMem2(run_spinup = 1, 
                             run_scenario = 1,
                             msl = msl,
                             mhwMat = mhwMat,
                             mlwMat = mlwMat,
                             scenario_calendar_years = 1928:2018,
                             spinup_iterations = 20, 
                             years_per_spinup_iter = 10) [[2]] %>% 
  mutate(run = "Full run")

mem_1_cohorts <- runCohortMem2(run_spinup = 1, 
                                     run_scenario = 1,
                                     msl = msl,
                                     mhwMat = mhwMat,
                                     mlwMat = mlwMat,
                                     scenario_calendar_years = 1928:2018,
                                     spinup_iterations = 20, 
                                     years_per_spinup_iter = 10) [[1]] %>% 
  mutate(run = "Full run")

mem_2a <- runCohortMem2(run_spinup = 1, 
                              run_scenario = 0,
                              msl = msl[1],
                              mhwMat = t(mhwMat[1,]),
                              mlwMat = t(mlwMat[1,]),
                              scenario_calendar_years = 1928,
                              spinup_iterations = 20, 
                              years_per_spinup_iter = 10)

cohort_2a <- mem_2a[[1]] %>%
  filter(year == 1928)

scenario_2a <- mem_2a[[2]] %>%
  filter(year == 1928)

species_2a <- mem_2a[[3]] %>% 
  filter(year == 1928)

mem_2b <- runCohortMem2(run_spinup = 0, 
                              run_scenario = 1,
                              msl = msl,
                              mhwMat = mhwMat,
                              mlwMat = mlwMat,
                              initial_cohorts = cohort_2a,
                              initial_scenario = scenario_2a,
                              initial_species = species_2a,
)[[2]] %>% 
  mutate(run = "Cohorta and scenario run separately")

mem_2b_cohorts <- runCohortMem2(run_spinup = 0, 
                                      run_scenario = 1,
                                      msl = msl,
                                      mhwMat = mhwMat,
                                      mlwMat = mlwMat,
                                      initial_cohorts = cohort_2a,
                                      initial_scenario = scenario_2a,
                                      initial_species = species_2a,
)[[1]] %>% 
  mutate(run = "Cohorta and scenario run separately")



together <- bind_rows(mem_1,
                      mem_2b)

together_2 <- bind_rows(mem_1_cohorts,
                        mem_2b_cohorts)

ggplot(together, aes(x = year, y = surface_elevation)) +
  geom_point(aes(color = run, shape = run))+
  scale_shape_manual(values = c(1,3,4))

ggplot(together, aes(x = year, y = aboveground_biomass)) +
  geom_point(aes(color = run, shape = run))+
  scale_shape_manual(values = c(1,3,4))

a_mem_run <- runCohortMem2(run_spinup = 1, 
                           run_scenario = 1,
                           msl = msl,
                           mhwMat = mhwMat,
                           mlwMat = mlwMat,
                           spinup_iterations = 200, 
                           years_per_spinup_iter = 20,
                           rootShape = c(0),
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

a_core <- buildSoilCore(cohorts = a_mem_run[[1]], coreYear = 2018, coreDepth = 200)

a_dates <- simulateRadioisotopeDates(core = a_core, cs137_anomaly = 0.6, 
                                     pb210_mobility_factor = 3, 
                                     cs137_mobility_factor = 3,
                                     cs137_erosion_loss = 0.5, 
                                     pb210_erosion_loss = 0.25)
