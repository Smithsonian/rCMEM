
library(tidyverse)
library(rCMEM)


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

# --------------------------------------------------
# 1. Species parameters
# --------------------------------------------------

species_codes <- c("SCAM", "SPPA", "PHAU")

bMax <- c(0.0522, 0.0522, 0.1044)

zVegMin  <- c(-1.4, -1.3, -2.0)
zVegPeak <- c(1.5, 1.8, 1)
zVegMax  <- c(1.8, 1.9, 2.5)

rootDepthMax <- c(40, 40, 60)
rootShape <- c(1, 1, 1)

rootToShoot <- rep(1.636903, 3)
rootTurnover <- rep(0.5, 3)
abovegroundTurnover <- rep(1, 3)
rootPackingDensity <- rep(0.07, 3)

# --------------------------------------------------
# 2. Scenario and species presence
# --------------------------------------------------

scenario_calendar_years <- 1928:2018
n_years <- length(scenario_calendar_years)

# All species permitted throughout
presence_all <- matrix(
  1,
  nrow = n_years,
  ncol = length(species_codes),
  dimnames = list(NULL, species_codes)
)

# PHAU absent for first 10 years
presence_introduction <- presence_all

presence_introduction[
  1:10, "PHAU"
] <- 0

# Control: PHAU never present
presence_control <- presence_all
presence_control[, "PHAU"] <- 0

# Check first years
head(presence_introduction, 12)

# --------------------------------------------------
# 3. Shared model inputs
# --------------------------------------------------

model_args <- list(
  run_spinup = 1,
  run_scenario = 1,
  
  scenario_calendar_years = scenario_calendar_years,
  
  msl = msl,
  mhwMat = mhwMat,
  mlwMat = mlwMat,
  
  spinup_iterations = 20,
  years_per_spinup_iter = 10,
  
  competition_function = 1,
  
  species_codes = species_codes,
  
  bMax = bMax,
  zVegMin = zVegMin,
  zVegPeak = zVegPeak,
  zVegMax = zVegMax,
  
  rootDepthMax = rootDepthMax,
  rootShape = rootShape,
  rootToShoot = rootToShoot,
  rootTurnover = rootTurnover,
  abovegroundTurnover = abovegroundTurnover,
  rootPackingDensity = rootPackingDensity,
  
  initElv = -13,
  captureRate = 5,
  suspendedSediment = 6e-06
)

# --------------------------------------------------
# 4. Run simulations
# --------------------------------------------------

run_control <- do.call(
  runCohortMem2,
  c(model_args, list(
    species_presence = presence_control
  ))
)

run_introduction <- do.call(
  runCohortMem2,
  c(model_args, list(
    species_presence = presence_introduction
  ))
)

# --------------------------------------------------
# 5. Compare species biomass
# --------------------------------------------------

species_results <- bind_rows(
  as_tibble(run_control[[3]]) %>%
    mutate(run = "No PHAU"),
  as_tibble(run_introduction[[3]]) %>%
    mutate(run = "PHAU introduced")
) %>%
  filter(year %in% scenario_calendar_years) %>%
  pivot_longer(
    cols = all_of(species_codes),
    names_to = "species",
    values_to = "biomass"
  )

ggplot(
  species_results,
  aes(x = year, y = biomass, fill = species)
) +
  geom_area() +
  facet_wrap(~run, ncol = 1) +
  geom_vline(
    xintercept = 1938,
    linetype = "dashed"
  ) +
  labs(
    x = "Year",
    y = "Aboveground biomass (g/cm2)",
    fill = "Species"
  ) +
  theme_bw()

# --------------------------------------------------
# 6. Compare elevation and total biomass
# --------------------------------------------------

scenario_results <- bind_rows(
  as_tibble(run_control[[2]]) %>%
    mutate(run = "No PHAU"),
  as_tibble(run_introduction[[2]]) %>%
    mutate(run = "PHAU introduced")
) %>%
  filter(year %in% scenario_calendar_years)

scenario_results %>%
  select(
    year, run,
    surface_elevation,
    aboveground_biomass,
    belowground_biomass
  ) %>%
  pivot_longer(
    cols = c(
      surface_elevation,
      aboveground_biomass,
      belowground_biomass
    )
  ) %>%
  ggplot(aes(year, value, color = run)) +
  geom_line() +
  geom_vline(
    xintercept = 1938,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~name,
    scales = "free_y",
    ncol = 1
  ) +
  labs(x = "Year", y = NULL) +
  theme_bw()

