# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATA ANALYSIS SCRIPT: ALLOMETRIC MODELS, DATA PREPARATION & ENVIRONMENTAL DATA ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Data used           : 1) RangeX_clean_BiomassSeedling_2023_2024.csv, RangeX_clean_BiomassDryweight_2023.csv, RangeX_clean_YearlySize_2021_2023_CHE.csv,
###                       RangeX_clean_InitialSize_2021.csv, Seeds/RangeX_clean_SeedTraits_2022.csv
###                       2) RangeX_clean_EnvHOBO_2021_2023_CHE.csv, RangeX_clean_EnvTMS4_2021_2023_CHE.csv, 2023_CAPHE_CleanData_20240211.csv, terra_final.tif,
###                       DEM.tif, study_region_2024.shp, study_region.shp
### Date last modified  : 27.10.2025
### Purpose             : 1) Create allometric models to calculate biomass for all focal individuals over time. Same for seedlings and initial size. Calculate number of seeds.
###                       Clean up data and save files ready to fit IPMs.
###                       2) Organize environmental data and make some plots for IPM paper. Get species distribution from Caphe data. Make overview plot.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rm(list = ls())

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1) DATA PREPARATION & ALLOMETRIC MODELS ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## PACKAGES etc. ---------------------------------------------------------------

# basic packages
library(dplyr); library(tidylog); library(janitor); library(tidyverse) # data manipulation
library(ggplot2) # test-plotting
library(stringr) # working with regex

# task-specific packages (include short description of what it is used for)
library(car) # get vif to check collinearity
library(lme4) # mixed models
library(ggResidpanel) # check model assumptions


## LOAD & PREPARE DATA SET -----------------------------------------------------

bio_seedling <- read_csv("data/raw/RangeX_clean_SeedlingBiomass_2023_2024_CHE.csv")
bio_adult <- read_csv("data/raw/RangeX_clean_FocalBiomass_2023_CHE.csv")
size_adult <- read_csv("data/raw/RangeX_clean_YearlyDemographics_2021_2023_CHE.csv")
size_initial <- read_csv("data/raw/RangeX_clean_InitialSize_2021_CHE.csv")
seeds <- read_csv("data/raw/RangeX_clean_SeedTraits_2022_CHE.csv")


# define data types
bio_seedling <- bio_seedling %>%
  mutate(across(c(height_vegetative_str, leaf_length1, leaf_length2, leaf_length3, number_leaves, dry_mass), as.numeric),
         across(c(germination_position_ID, germination_seedling_ID, species, functional_group, collector, seed_origin), as.character),
         across(c(date_measurement, date_sowing), ~ as.Date(., format = "%Y-%m-%d")))

bio_adult <- bio_adult %>%
  mutate(across(c(dry_mass), as.numeric),
         across(c(unique_plant_ID, species), as.character),
         across(c(date_collection), ~ as.Date(., format = "%Y-%m-%d")))

size_adult <- size_adult %>%
  mutate(across(c(height_vegetative_str, height_reproductive_str, height_vegetative, leaf_length1, leaf_length2, leaf_length3, number_leaves, number_flowers, mean_inflorescence_size), as.numeric),
         across(c(unique_plant_ID, species, functional_group, collector), as.character),
         across(c(survival, herbivory), as.factor),
         across(c(date_measurement, date_planting), ~ as.Date(., format = "%Y-%m-%d")),
         across(c(height_reproductive, vegetative_width, height_nathan, stem_diameter, leaf_width, petiole_length, number_tillers, number_branches, number_leafclusters), as.character),
         across(c(height_reproductive, vegetative_width, height_nathan, stem_diameter, leaf_width, petiole_length, number_tillers, number_branches, number_leafclusters), as.numeric))

size_initial <- size_initial %>%
  mutate(across(c(height_vegetative_str, leaf_length1, leaf_length2, leaf_length3, number_leaves), as.numeric),
         across(c(unique_plant_ID, species, functional_group, collector), as.character),
         across(c(date_measurement, date_planting), ~ as.Date(., format = "%Y-%m-%d")),
         across(c(vegetative_width), as.character),
         across(c(vegetative_width), as.numeric))

seeds <- seeds %>%
  mutate(across(c(unique_plant_ID, species, counter), as.character),
         across(c(inflorescence_size, no_seeds, seedweight), as.numeric),
         across(c(date_collection), ~ as.Date(., format = "%Y-%m-%d")))

bio_seedling_org <- bio_seedling
bio_adult_org <- bio_adult
size_adult_org <- size_adult
size_initial_org <- size_initial
seeds_org <- seeds

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIX GENERAL PROBLEMS ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# solve hypper problem: in 2023, leaf_length was forgotten to be measured for a subset of bare plots at the high site --> replace with average leaf lengths for other bare plots (both lo and hi) in 2023
hypper_bare_2023 <- size_adult %>%
  filter(year(date_measurement) == 2023 & species == "hypper" & grepl("bare", unique_plant_ID)) %>%
  filter(!if_all(c(leaf_length1, leaf_length2, leaf_length3), is.na) ) %>% # it's either all 3 leaves or none, no individuals with 1 or 2 NAs
  rowwise() %>%
  mutate(leaf_lengths_sorted = list(sort(c(leaf_length1, leaf_length2, leaf_length3), decreasing = TRUE)), # sorted vector
         leaf_max = leaf_lengths_sorted[1],
         leaf_max2 = leaf_lengths_sorted[2],
         leaf_min = leaf_lengths_sorted[3]) %>%
  ungroup() %>%
  dplyr::select(-leaf_lengths_sorted) 

# get average length of longest, second longest and 3rd longest leaf
longest <- mean(hypper_bare_2023$leaf_max)
longest2 <- mean(hypper_bare_2023$leaf_max2)
longest3 <- mean(hypper_bare_2023$leaf_min)

# now replace the NAs
size_adult <- size_adult %>%
  mutate(leaf_length1 = ifelse(year(date_measurement) == 2023 & species == "hypper" & grepl("bare", unique_plant_ID) & is.na(height_vegetative_str) == FALSE & grepl("hi", unique_plant_ID),
                               longest, leaf_length1),
         leaf_length2 = ifelse(year(date_measurement) == 2023 & species == "hypper" & grepl("bare", unique_plant_ID) & is.na(height_vegetative_str) == FALSE & grepl("hi", unique_plant_ID),
                               longest2, leaf_length2),
         leaf_length3 = ifelse(year(date_measurement) == 2023 & species == "hypper" & grepl("bare", unique_plant_ID) & is.na(height_vegetative_str) == FALSE & grepl("hi", unique_plant_ID),
                               longest3, leaf_length3))
  

# check how many of the alive plants have missing values in height_vegetative_str, number_leaves and Leaf_length1
size_adult_na <- size_adult %>%
  filter(survival == 1,
         is.na(height_vegetative_str) | is.na(number_leaves) | is.na(leaf_length1))

# only 30 out of 5400 rows (individuals completely forgotten in one year, single measurements forgotten etc.) --> can be deleted

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMASS CALCULATION ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## PREPARE DATA ----------------------------------------------------------------

### SIZE DATA: PREPARE

# add year column and flower status to all size data frames
size_adult <- size_adult %>%
  mutate(year = year(date_measurement),
         flower_status = case_when(is.na(height_reproductive_str) == FALSE ~ 1, 
                                   is.na(height_reproductive_str) == TRUE & is.na(height_vegetative_str) == TRUE ~ NA,
                                   is.na(height_reproductive_str) == TRUE & is.na(height_vegetative_str) == FALSE ~ 0),
         flower_status = as.factor(flower_status))
size_initial <- size_initial %>%
  mutate(year = year(date_measurement),
         flower_status = as.factor(ifelse(is.na(height_vegetative_str) == FALSE | is.na(number_leaves) == FALSE, "0", NA)))


# calculate mean leaf length, add empty reproductive columns (for initial size), change NA to 0 as number_flower when no reproductive height and no number of flowers
size_adult <- size_adult %>%
  rowwise() %>%
  mutate(mean_leaf_length = ifelse(is.na(leaf_length1) == FALSE | is.na(leaf_length2) == FALSE | is.na(leaf_length3) == FALSE, 
                                   mean(c_across(c(leaf_length1, leaf_length2, leaf_length3)), na.rm = TRUE), NA), # NAs will be deleted before calculation, mean(leaf_length1, leaf_length2, leaf_length3, na.rm = TRUE)
         number_flowers = case_when(survival == 1 & (is.na(height_vegetative_str) == FALSE | is.na(number_leaves) == FALSE) & is.na(number_flowers) == TRUE & is.na(height_reproductive_str) == TRUE ~ 0,
                                    is.na(number_flowers) == FALSE ~ number_flowers,
                                    .default = NA)) %>% 
  ungroup()
size_initial <- size_initial %>%
  rowwise() %>%
  mutate(mean_leaf_length = ifelse(is.na(leaf_length1) == FALSE | is.na(leaf_length2) == FALSE | is.na(leaf_length3) == FALSE, 
                                   mean(c_across(c(leaf_length1, leaf_length2, leaf_length3)), na.rm = TRUE), NA), # NAs will be deleted before calculation, mean(leaf_length1, leaf_length2, leaf_length3, na.rm = TRUE) 
         height_reproductive_str = NA,
         number_flowers = ifelse((is.na(height_vegetative_str) == FALSE | is.na(number_leaves) == FALSE) & is.na(height_reproductive_str) == TRUE, 0, NA)) %>% 
  ungroup()

# check for NAs in metadata
size_adult_na <- size_adult %>%
  filter(if_any(1:7, is.na)) # the 36 individuals which were replaced right before winter in 2021 and therefore don't have size measurements for 2021 --> add year later manually
size_initial_na <- size_initial %>%
  filter(if_any(1:7, is.na)) # none

# check for NAs in relevant columns
size_adult_na <- size_adult %>%
  filter(if_any(c(number_leaves, mean_leaf_length, height_vegetative_str, number_flowers), ~ is.na(.)) & survival == 1) # 41 (different reasons)
size_initial_na <- size_initial %>%
  filter(if_any(c(number_leaves, mean_leaf_length, height_vegetative_str), ~ is.na(.))) # 13, some completely empty rows, some forgotten leaf counts
  
  
### BIOMASS DATA: PREPARE

# add year column to bio data frames
bio_adult <- bio_adult %>%
  mutate(year = year(date_collection))
bio_seedling <- bio_seedling %>%
  mutate(year = year(date_measurement))


# for adult: filter only the size measurements from 2023
size_adult_2023 <- size_adult %>%
  filter(year(date_measurement) == 2023)

# add 2023 size measurements on biomass data frame
bio_adult_2023 <- full_join(bio_adult, size_adult, by = c("unique_plant_ID", "species", "year")) # 3 rows only in bio
anti_join(bio_adult, size_adult, by = c("unique_plant_ID", "species", "year"))

# it's the three rows with no collection date and therefore also no dry mass --> delete

bio_adult_2023 <- left_join(bio_adult, size_adult, by = c("unique_plant_ID", "species", "year")) %>% # the 412 rows only in size_adult_2023 were dead in 2023 and not harvested
  filter(!is.na(date_collection))

# are there any NAs in significant columns? 
bio_adult_2023_na <- bio_adult_2023 %>%
  filter(if_any(c(number_leaves, mean_leaf_length, height_vegetative_str, number_flowers), ~ is.na(.))) 

# CHE.lo.ambi.bare.wf.02.09.2 was forgotten to measure
# CHE.lo.ambi.bare.wf.10.10.2 dito
# CHE.lo.ambi.vege.wf.08.15.1 leaf measurements were forgotten

# rest has a reproductive height, but flowers were not counted 

# delete rows with NA for model fitting
bio_adult_2023 <- bio_adult_2023 %>%
  filter(if_all(c(number_leaves, mean_leaf_length, height_vegetative_str, number_flowers), ~ !is.na(.))) 


# for seedlings: create variable mean_leaf_length, add reproductive height, flower status, number of flowers
bio_seedling <- bio_seedling %>%
  rowwise() %>%
  mutate(height_reproductive_str = NA,
         number_flowers = ifelse(is.na(height_vegetative_str) == FALSE & is.na(height_reproductive_str) == TRUE, 0, NA),
         flower_status = as.factor(ifelse(number_flowers == 0 & is.na(height_reproductive_str) == TRUE & is.na(height_vegetative_str) == FALSE, "0", "1")),
         mean_leaf_length = mean(c_across(c(leaf_length1, leaf_length2, leaf_length3)), na.rm = TRUE)) %>% # NAs will be deleted before calculation, mean(leaf_length1, leaf_length2, leaf_length3, na.rm = TRUE) 
  ungroup()

# are there any NAs in significant columns? 
bio_seedling_na <- bio_seedling %>%
  filter(if_any(c(number_leaves, mean_leaf_length, height_vegetative_str), ~ is.na(.))) # no

# check how many entries
table(bio_seedling$species) # between 19 (hypper, scacol) and 46 (cenjac)

### BIOMASS DATA: COMBINE

# delete unnecessary column for later adding up seedling and adult bio
bio_adult_2023_ready <- bio_adult_2023 %>%
  dplyr::select(unique_plant_ID, species, dry_mass, date_measurement, height_vegetative_str, height_reproductive_str, number_leaves, mean_leaf_length, number_flowers, flower_status)

bio_seedling_ready <- bio_seedling %>%
  dplyr::select(germination_seedling_ID, species, dry_mass, date_measurement, height_vegetative_str, height_reproductive_str, number_leaves, mean_leaf_length, number_flowers, flower_status) %>%
  rename("unique_plant_ID" = "germination_seedling_ID")

# add together
bio_ready <- bind_rows(bio_adult_2023_ready, bio_seedling_ready)

# model flowering and non-flowering species separately: check how many individuals per level
table(bio_ready$species, bio_ready$flower_status) # looks alright, min. of 16 (daucar flowering)

### LOGGING ALL

# exclude rows with any NA in any of the variables, plus log all the variables
# log transform all model variables, exclude rows with missing relevant variables
bio_ready_log <- bio_ready %>%
  mutate(log_dry_mass = log(dry_mass),
         log_number_leaves = log(number_leaves),
         log_mean_leaf_length = log(mean_leaf_length),
         log_height_vegetative_str = log(height_vegetative_str),
         log_height_reproductive_str = log(height_reproductive_str),
         log_number_flowers = ifelse(number_flowers == 0, 0, log(number_flowers))) %>%
  filter(if_else(flower_status == 0,
                 !is.na(height_vegetative_str) & !is.na(number_leaves) & !is.na(mean_leaf_length) & !is.na(dry_mass), # rows removed: individuals from seed experiment 2 which were to small to weight
                 !is.na(height_vegetative_str) & !is.na(number_leaves) & !is.na(mean_leaf_length) & !is.na(number_flowers) & !is.na(height_reproductive_str) & !is.na(dry_mass)))

# log transform all model variables (no exclusions if possible)
size_adult_log <- size_adult %>%
  mutate(log_number_leaves = log(number_leaves),
         log_mean_leaf_length = log(mean_leaf_length),
         log_height_vegetative_str = log(height_vegetative_str),
         log_height_reproductive_str = log(height_reproductive_str),
         log_number_flowers = ifelse(number_flowers == 0, 0, log(number_flowers)))


# log transform all model variables (no exclusions if possible)
size_initial_log <- size_initial %>%
  mutate(log_number_leaves = log(number_leaves),
         log_mean_leaf_length = log(mean_leaf_length),
         log_height_vegetative_str = log(height_vegetative_str),
         log_height_reproductive_str = log(height_reproductive_str),
         log_number_flowers = ifelse(number_flowers == 0, 0, log(number_flowers)))


## FIT MODELS & PREDICT --------------------------------------------------------


### DEFINE FUNCTION

# make a function to fit models and predict

predict_biomass <- function(dat1, dat2) { # where dat1 = data to fit model on and dat2 = new data to do predictions on
  # prepare empty tibbles
  model_allspecies <- tibble("species" = character(), "flower_status" = numeric(), "model" = list(), "R2" = numeric()) # store model summaries
  variables_allspecies <- tibble("species" = character(), "flower_status" = numeric(), "intercept" = numeric(), "log_number_leaves" = numeric(), 
                                 "log_mean_leaf_length" = numeric(), "log_height_vegetative_str" = numeric(), "log_height_reproductive_str" = numeric(), "log_number_flowers" = numeric()) # store model variables
  all_predictions <- tibble() # store predictions
  
  # define unique species for looping
  species <- unique(dat1$species)
  
  for (spec in species) {
    
    dat_spec <- dat1[dat1$species == spec, ]
    
    
    # fit different models for flowering and non-flowering species
    for (j in 1:2) {
      
      flow <- j-1
      
      dat_spec_flowstat <- dat_spec[dat_spec$flower_status == flow, ]
      
      if (flow == 0) {
        # leaf variables separate
        mod_spec_flow <- lm(
          log_dry_mass ~ log_number_leaves + log_mean_leaf_length + log_height_vegetative_str,
          data = dat_spec_flowstat)
        
        ## leaf variables combined
        #mod_spec_flow <- lm(
        #  log_dry_mass ~ log_leaf_var + log_height_vegetative_str,
        #  data = dat_spec_flowstat)
        
        coef_spec_flow <- coef(mod_spec_flow)
        
        # save results
        model_allspecies <- model_allspecies %>%
          add_row(species = spec, 
                  flower_status = flow,
                  model = list(mod_spec_flow),
                  R2 = summary(mod_spec_flow)$r.squared)
        
        variables_allspecies <- variables_allspecies %>%
          add_row(species = spec,
                  flower_status = flow,
                  intercept = coef_spec_flow["(Intercept)"],
                  log_number_leaves = coef_spec_flow["log_number_leaves"],
                  log_mean_leaf_length = coef_spec_flow["log_mean_leaf_length"],
                  log_height_vegetative_str = coef_spec_flow["log_height_vegetative_str"])
        
      } else {
        # leaf variables separate
        mod_spec_flow <- lm(
          log_dry_mass ~ log_number_leaves + log_mean_leaf_length + log_height_vegetative_str + log_height_reproductive_str + log_number_flowers,
          data = dat_spec_flowstat)
        
        ## leaf variables combined
        #mod_spec_flow <- lm(
        #  log_dry_mass ~ log_leaf_var + log_height_vegetative_str + log_height_reproductive_str + log_number_flowers,
        #  data = dat_spec_flowstat)
        
        coef_spec_flow <- coef(mod_spec_flow)
        
        # save results
        model_allspecies <- model_allspecies %>%
          add_row(species = spec, 
                  flower_status = flow,
                  model = list(mod_spec_flow),
                  R2 = summary(mod_spec_flow)$r.squared)
        
        variables_allspecies <- variables_allspecies %>%
          add_row(species = spec,
                  flower_status = flow,
                  intercept = coef_spec_flow["(Intercept)"],
                  log_number_leaves = coef_spec_flow["log_number_leaves"],
                  log_mean_leaf_length = coef_spec_flow["log_mean_leaf_length"],
                  log_height_vegetative_str = coef_spec_flow["log_height_vegetative_str"],
                  log_height_reproductive_str = coef_spec_flow["log_height_reproductive_str"],
                  log_number_flowers = coef_spec_flow["log_number_flowers"])
      }
      
      if (flow %in% unique(dat2$flower_status[dat2$species == spec])) { # only predict if the flower status is present (i.e. not for flower_status == 1 in intial size data)
        
        # predict
        pred_data <- dat2[dat2$species == spec & dat2$flower_status == flow,] %>%
          janitor::remove_empty(which = "rows")
        
        pred_data <- pred_data %>%
          mutate(fit = exp(predict(mod_spec_flow, pred_data))) %>%
          dplyr::select(unique_plant_ID, date_measurement, species, flower_status, fit)
        
        # add predictions
        all_predictions <- bind_rows(all_predictions, pred_data)
        
      }
    }
  }
  
  # return the results as a list
  list(models = model_allspecies, variables = variables_allspecies, predictions = all_predictions)
}

# for the model fitting: there's an outlier in daucar and plamed --> those are weird measurements (CHE.lo.ambi.vege.wf.02.12.1, CHE.lo.ambi.bare.wf.03.29.1), filter out and refit
bio_ready_log_noout <- bio_ready_log %>%
  filter(!unique_plant_ID == "CHE.lo.ambi.vege.wf.02.12.1", ! unique_plant_ID == "CHE.lo.ambi.bare.wf.03.29.1")

### RUN FUNCTION: ADULT SIZE

# conditions: all variables logged (explanatory and response variables), leaf variables separately, plamed and daucar outliers removed for model fitting (but not for predictions)

# run function for adult size data
results_YS_noout <- predict_biomass(dat1 = bio_ready_log_noout, dat2 = size_adult_log)

# extract results (noout = no outlier)
model_summaries_YS_noout <- results_YS_noout$models # TABLE S3
model_variables_YS_noout <- results_YS_noout$variables
predictions_YS_noout <- results_YS_noout$predictions

# get average R2 for focal individual models
mean(model_summaries_YS_noout$R2) # 0.8787244
sd(model_summaries_YS_noout$R2) # 0.08625693

# get sample size of the groups
summary <- bio_ready_log_noout %>%
  group_by(species, flower_status) %>%
  summarize(count = n())


# add predictions back on
size_adult_log_pred <- left_join(size_adult_log, predictions_YS_noout, by = c("unique_plant_ID", "date_measurement", "species", "flower_status"))

# which rows are only in size_adult_ect. and have no predictions?
b <- size_adult_log[size_adult_log$survival == 0, ]
length(b$unique_plant_ID) # 631 individuals are dead...
a <- size_adult_log[is.na(size_adult_log$height_vegetative_str) == TRUE, ]
length(a$unique_plant_ID) 


c <- anti_join(a, b, by = "unique_plant_ID") 

# 19 were not found in one year or forgotten - they are also all NA, but survival = 1 (added later due to illogical time series)
# 4 are plants from 2023 with only reproductive height (rest already dry/ dead)


# check how many rows have no fit, but are alive
size_adult_log_pred_na <- size_adult_log_pred %>%
  filter(if_any(c(fit), ~ is.na(.)) & survival == 1) # 41 rows: 23 as seen above (forgotten/ only rep. height), rest missing no. leaves, leaf_length or similar)


### RUN FUNCTION: INITIAL SIZE

# conditions: all variables logged (explanatory and response variables), leaf variables separately, plamed and daucar outliers removed for model fitting (but not for predictions)

# run function for initial size data
results_IS_noout <- predict_biomass(dat1 = bio_ready_log_noout, dat2 = size_initial_log)

# extarct results
model_summaries_IS_noout <- results_IS_noout$models # they are identical to model_summaries_YS as the same data is used to fit the models
model_variables_IS_noout <- results_IS_noout$variables
predictions_IS_noout <- results_IS_noout$predictions

# get average R2 for focal individual models (also identical)
mean(model_summaries_IS_noout$R2) # 0.8787244
sd(model_summaries_IS_noout$R2)

# add predictions back on
size_initial_log_pred <- left_join(size_initial_log, predictions_IS_noout, by = c("unique_plant_ID", "date_measurement", "species", "flower_status"))

# which rows are only in size_initial_ect. and have no predictions?
a <- size_initial_log[is.na(size_initial_log$height_vegetative_str) == TRUE, ]
length(a$unique_plant_ID) # 3 rows completely NA (forgotten?)



## CHECK: MEASURED VS. FITTED --------------------------------------------------

# run modelling and prediction for complete bio data set (not just add predictions from size_adult on, then the ones for the seedlings are missing)

results_bio <- predict_biomass(dat1 = bio_ready_log_noout, dat2 = bio_ready_log_noout)

# extract results
model_summaries_bio <- results_bio$models
model_variables_bio <- results_bio$variables
predictions_bio <- results_bio$predictions

# rename fit column
predictions_bio <- predictions_bio %>%
  rename("fit_bio" = "fit")


# add predictions back on
bio_ready_log_pred <- left_join(bio_ready_log, predictions_bio, by = c("unique_plant_ID", "date_measurement", "species", "flower_status"))


# plot
ggplot(data = bio_ready_log_pred[is.na(bio_ready_log_pred$fit_bio) == FALSE,], aes(x = dry_mass, y = fit_bio, col = flower_status)) +
  geom_point() +
  facet_wrap(~species, scales = "free") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_smooth(method = "lm") +
  ggtitle("original (all log, but backtransformed)") + # column "fit" is backtransformed in predict_biomass function
  theme_bw() 
  


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SEED CALCULATION ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# add treatment columns to seed data
seeds <- seeds %>%
  mutate(site = str_sub(unique_plant_ID, 5, 6),
         treat_warm = str_sub(unique_plant_ID, 8, 11),
         treat_comp = str_sub(unique_plant_ID, 13, 16),
         plot_ID = paste(site, treat_warm, treat_comp, sep = "."))



# calculate number of seeds produced of each individual every year
# for species with no inflorescence measurements: no seeds = flowers * average number of seeds per flower in treatment --> or ignore
# for species with inflorescence measurements: calculate no seeds for this inflorescence length, then * number of flowers

# check out which species have consistent inflorescence measurements
size_adult %>%
  group_by(species, year) %>%
  filter(!is.na(mean_inflorescence_size)) %>%
  dplyr::summarize(mean_inflorescence_size = n())

# values more or less complete for plamed, salpra, daucar and brapin
# count individuals with flowers and inflorescence length vs. those missing the inflorescence size
size_adult %>%
  group_by(species, year) %>%
  filter(species %in% c("daucar", "brapin", "salpra", "plamed"), 
         is.na(number_flowers) == FALSE & is.na(mean_inflorescence_size) == TRUE) %>%
  dplyr::summarize(mean_inflorescence_size = n())

# only a few cases where the individual was flowering, but inflroescence length not measured: 
# 2 daucar in 2022, 3 in 2023
# 2 plamed in 2021, 2 in 2023, 2 in 2023
# 2 salpra in 2022, 3 in 2023

# calculate mean no seeds per species and treatment: see how many individuals counted per species & treatment
seed_counts <- seeds %>%
  group_by(species, site, treat_warm, treat_comp) %>%
  dplyr::summarize(seed_counts = n())

# 20 combinations with less than 10 samples --> leave out warming treatment
seed_counts <- seeds %>%
  group_by(species, site, treat_comp) %>%
  dplyr::summarize(seed_counts = n())

# still 9 combinations with less than 10 samples (all hi vege) --> check out mean to see whether it's ok to leave out site as well
seed_counts <- seeds %>%
  group_by(species, site, treat_comp) %>%
  dplyr::summarize(seed_counts = n(),
            mean_seeds = mean(no_seeds))


# for now, ignore inflorescence length
# for silvul hi, silvul lo, brapin lo, broere hi, broere lo, hypper hi and plamed hi, combine on site level (<5 samples)
seed_counts_sitetreat <- seeds %>%
  group_by(species, site, treat_comp) %>%
  dplyr::summarize(seed_counts = n(),
            mean_seeds = mean(no_seeds))
seed_counts_site <- seeds %>%
  group_by(species, site) %>%
  dplyr::summarize(seed_counts = n(),
            mean_seeds = mean(no_seeds))

# add the seed counts on
size_adult_log_pred <- size_adult_log_pred %>%
  mutate(site = str_sub(unique_plant_ID, 5, 6),
         treat_warm = str_sub(unique_plant_ID, 8, 11),
         treat_comp = str_sub(unique_plant_ID, 13, 16)) %>%
  left_join(seed_counts_sitetreat, by = c("site", "treat_comp", "species")) %>% # rows only in x: species which didn't flower in certain treatments at all (e.g. salpra in vege etc.)
  dplyr::select(-seed_counts) %>%
  rename("mean_seeds_sitetreat" = "mean_seeds") %>%
  left_join(seed_counts_site, by = c("site", "species")) %>%
  dplyr::select(-seed_counts) %>%
  rename("mean_seeds_site" = "mean_seeds")

# now define what to use for each species: generally use the site - treat_comp calculation, but for the combinations above, use site only, then calculate number of seeds
size_adult_log_pred <- size_adult_log_pred %>%
  mutate(mean_seeds = case_when(species %in% c("silvul", "broere", "brapin", "plamed", "salpra") ~ mean_seeds_site,
                                species == "hypper" & site == "hi" ~ mean_seeds_site,
                                species == "scacol" & site == "lo" ~ mean_seeds_site,
                                .default = mean_seeds_sitetreat),
         number_seeds = mean_seeds * number_flowers)


# check whether there are any flowering species (flower_status == 1) with no number_seeds
check <- size_adult_log_pred %>%
  filter(flower_status == 1 & is.na(number_seeds))

# 22 individuals have flower_status == 1, but no number_seeds --> reasons:
# reproductive height was measured, flower count forgotten
# hi ambi vege brapin/ lo ambi vege salpra/ lo ambi vege scacol & lo ambi vege plamed: wasn't collected --> use mean_seeds_site 

# 12 rows still NA for number_seeds after corrections: no number_flowers, so ignore


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# COMBINE ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### SIZE ADULT

# now combine calculated biomass predictions, 2023 biomass measurements and calculated seed numbers to final data frame for IPMs (delete columns not needed)
size_adult_log_pred <- size_adult_log_pred  %>%
  left_join(bio_adult, by = c("unique_plant_ID", "species", "year")) 

# why are there 3 rows only in bio_adult but not size_adult?
anti_join(bio_adult, size_adult_log_pred[size_adult_log_pred$year == 2023,], by = c("unique_plant_ID", "species", "year"))

# they have no collection date nor a dry mass (use predictions later on)

# check whether there are any individuals with no dry_mass but measured variables
size_adult_log_pred_na <- size_adult_log_pred %>%
  filter(year == 2023 & is.na(dry_mass) == TRUE & is.na(height_vegetative_str) == FALSE) # yes, four --> use fit later on

# check out 2023 (where there's both a dry mass as well as predictions)
size_adult_log_pred_2023 <- dplyr::filter(size_adult_log_pred, lubridate::year(date_collection) == 2023)# looks good


# now select wanted columns and combine real dry mass and predicted dry mass
size_adult_final <- size_adult_log_pred %>%
  rowwise() %>%
  mutate(biomass_comb = ifelse(is.na(dry_mass) == FALSE, dry_mass, fit)) %>%
  ungroup() %>%
  dplyr::select(c(1:7, flower_status, biomass_comb, number_seeds, height_vegetative_str, number_leaves)) # or without height_vegetative_str, number_leaves


### SIZE INITIAL

# select wanted columns
size_initial_final <- size_initial_log_pred %>%
  rename("biomass" = "fit") %>%
  dplyr::select(c(1:7, flower_status, biomass, height_vegetative_str, number_leaves), -ind_number) # or without height_vegetative_str, number_leaves


### FINAL NA CHECK

# NA in metadata
size_adult_final_na <- size_adult_final %>%
  filter(if_any(1:7, is.na)) # 36 without any measurement date & collector
size_initial_final_na <- size_initial_final %>%
  filter(if_any(1:6, is.na)) # none

# check for NAs in relevant columns
size_adult_final_na <- size_adult_final %>%
  filter(if_any(c(biomass_comb), ~ is.na(.)) & survival == 1) # 33 individuals which are alive but no dry mass (prob. because some value was forgotten etc.)
size_adult_final_na <- size_adult_final %>%
  filter(if_any(c(number_seeds), ~ is.na(.)) & flower_status == 1) # 12 individuals which have flower status 1 but no number of seeds (prob. because no flower count)


size_initial_final_na <- size_initial_final %>%
  filter(if_any(c(biomass), ~ is.na(.))) # 13, some completely empty rows, some forgotten leaf counts and therefore no biomass predictions

### make wide
size_adult_final <- size_adult_final %>%
  mutate(year = year(date_measurement))

size_adult_final[is.na(size_adult_final$year) == TRUE,] # the 36 individuals with no measurements in 2021 because late planting --> add 2021 as year manually

size_adult_final <- size_adult_final %>%
  mutate(year = ifelse(is.na(year) == TRUE, 2021, year))

size_adult_final[is.na(size_adult_final$year) == TRUE,] # fixed


size_adult_final_wide <- size_adult_final %>%
  dplyr::select(-date_measurement, -collector) %>%  # ignore these columns
  pivot_wider(id_cols = c(unique_plant_ID, species, functional_group, date_planting),  # columns to keep the same
              names_from = year,  # use year column to create the new column names
              values_from = c(survival, flower_status, biomass_comb, number_seeds, height_vegetative_str, number_leaves))  # columns to pivot (or without height_vegetative_str, number_leaves)



# explanation of NAs
# survival = 1, no biomass NA in corresponding year: all good
# survival = 1, biomass NA in corresponding year: measurement forgotten or plant died
# survival = 0: everything should be NA in corresponding years

# flower status = 1, >0 seeds in corresponding year: all good
# flower status = 1, NA seeds in corresponding year: reproductive height was measured, but flower count forgotten
# flower status = 0, 0 seeds in corresponding year: no flowering

# also why are there 1836 rows??
size_adult_final_wide %>%
  dplyr::select(unique_plant_ID) %>%
  group_by(unique_plant_ID) %>% 
  filter(n()>1)

# 0 duplicates - so it's the 36 individuals replaced after the first yearly size measurements (they will have NA survival in 2021)


# count how many full rows there are (i.e. species living until the end, never forgotten any measurements)
complete <- size_adult_final_wide[complete.cases(size_adult_final_wide) == TRUE,] # 1334 individuals

complete_died <- size_adult_final_wide %>%
  filter(
    case_when( # conditions for 2021
      survival_2021 == 1 ~ !is.na(flower_status_2021) & !is.na(biomass_comb_2021) & !is.na(number_seeds_2021),
      survival_2021 == 0 ~ TRUE, # allow NAs if survival is 0
      TRUE ~ TRUE), # if survival is NA 
    case_when( # conditions for 2022
      survival_2022 == 1 ~ !is.na(flower_status_2022) & !is.na(biomass_comb_2022) & !is.na(number_seeds_2022),
      survival_2022 == 0 ~ TRUE,
      TRUE ~ TRUE),
    case_when(# conditions for 2023
      survival_2023 == 1 ~ !is.na(flower_status_2023) & !is.na(biomass_comb_2023) & !is.na(number_seeds_2023),
      survival_2023 == 0 ~ TRUE,
      TRUE ~ TRUE))

# 1761 rows are alright (i.e. they're either complete (surviving until the end) or they died along the way or they were planted after the 2021 yearly measurements) --> no variables missing

# get the rows with missing variables
incomplete <- anti_join(size_adult_final_wide, complete_died, by = colnames(size_adult_final_wide)) %>%
  mutate(site = str_sub(unique_plant_ID, 5, 6),
         treat_warm = str_sub(unique_plant_ID, 8, 11),
         treat_comp = str_sub(unique_plant_ID, 13, 16))
  

# 39 rows: get the different species and treatments
summary_table <- incomplete %>%
  group_by(species, site, treat_warm, treat_comp) %>%
  summarise(count = n_distinct(unique_plant_ID), .groups = "drop")

# fairly evenly distributed between 1 and 3 --> check the ones with 5 or more
# daucar lo ambi bare
# plamed lo ambi bare

a <- size_adult_final %>%
  filter(grepl("lo.ambi.bare", unique_plant_ID) == TRUE & species == "daucar") # 7 missing biomass calcualtions in 2022
a_a <- size_adult_log_pred %>%
  filter(grepl("lo.ambi.bare", unique_plant_ID) == TRUE & species == "daucar" & year == 2022)  
# it's because the number of leaves is missing (they were already dry as data collection as too late in 2022)

b <- size_adult_final %>%
  filter(grepl("lo.ambi.bare", unique_plant_ID) == TRUE & species == "plamed") # 4 missing biomass calcualtions in 2023, 1 missing seed number in 2023
b_b <- size_adult_log_pred %>%
  filter(grepl("lo.ambi.bare", unique_plant_ID) == TRUE & species == "plamed" & year == 2023)  
# number of leaves, leaf length and veg. height are missing because already dry for the 4 missing biomasses
# the missing seeds have a weighted biomass (CHE.lo.ambi.bare.wf.10.10.2), but was never recorded in the yearly size (no number of leaves etc. - nothing)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SAVE DATA ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

write_csv(size_adult_final_wide, "data/derived/RangeX_ready_YSwide_IPM.csv")
write_csv(size_adult_final, "data/derived/RangeX_ready_YS_IPM.csv")
write_csv(size_initial_final, "data/derived/RangeX_ready_IS_IPM.csv")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2) ENVIRONMENTAL DATA ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## PACKAGES & FUNCTIONS -------------------------------------------------------------------

# task-specific packages
library(suncalc)
library(terra)
library(rnaturalearth)
library(sf)
library(lwgeom) # smooth europe outline
library(rnaturalearth)
library(smoothr)

# small functions

# function to get the opposite of %in%
`%nin%` <- Negate(`%in%`)

# function to extract legend from plot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# naming etc.
species_names <- c("brapin" = "Brachypodium\npinnatum", "broere" = "Bromus\nerectus", "daucar" = "Daucus\ncarota", "hypper" = "Hypericum\nperforatum", "medlup" = "Medicago\nlupulina", 
                   "plamed" = "Plantago\nmedia", "silvul" = "Silene\nvulgaris", "scacol" = "Scabiosa\ncolumbaria", "cenjac" = "Centaurea\njacea", "salpra" = "Salvia\npratensis")



## LOAD DATA -------------------------------------------------------------------

# load environmental data
dat_hobo <- read_csv("/Users/mac/Desktop/ETH_Phd+/Projects/RangeX/RangeX_Data/6_DataClean/RangeX_clean_EnvHOBO_2021_2023_CHE.csv") # cleaned data paper hobo data from RangeX experiment 2021 - 2023
dat_tms4 <- read_csv("/Users/mac/Desktop/ETH_Phd+/Projects/RangeX/RangeX_Data/6_DataClean/RangeX_clean_EnvTMS4_2021_2023_CHE.csv")  # cleaned data paper TMS4 data from RangeX experiment 2021 - 2023

# read vegetation height data
dat_height <- read_csv("data/raw/RangeX_clean_VegHeight_2023_CHE.csv")

# load Caphe data
dat_caphe <- read_csv("/Users/mac/Desktop/ETH_Phd+/Projects/Others/Caphe/2024_CAPHE_CleanData_20250306.csv") # this data is not freely available
dat_alti <- read_csv("data/raw/SpeciesAltitudes.csv") # Caphe altitudes provided by Mikko Tiusanen on 21.03.2024


# load macroclimate data (not provided --> )
dat_clim <- terra::rast("data/raw/terra_final.tif")
dat_clim_int <- terra::rast("data/raw/terra_final_integrated.tif")
calanda_1 <- st_read("data/raw/study_region_2024.shp") # shapefile of study region 2024 (received from Billur Jan 2025)
calanda_2 <- st_read("data/raw/study_region.shp") # shapefile of study region 2023 (received from Billur Jan 2025)

dem <- rast("/Volumes/ExtDesktop/Downloaded_Data/for_Evelin/data/topography/DEM.tif") # altitudes compiled by Billur from https://earthexplorer.usgs.gov/ (received Jan 2025)

# NOTE: The climate data is 1x1 km resolution of mean annual temperature. It is based on downscaled TerraClimate data, with the downscaling based on CHELSAcruts min and max
# temperatures between 2000 - 2016. Downscaling was done using the factor change methodology and using the metods described in Iseli et al. (2025),  https://doi.org/10.1111/1365-2745.70114
# The integrated data is the running mean over 5 years to smooth out yearly variation.

# Topography and Caphe data are not included in the githib repository. See "Data availability" for more information.


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MICROCLIMATE  ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### PREPARE DATA ---------------------------------------------------------------

### HOBO

# separate site, treatments etc., set time zone to Europe/Zurich
dat_hobo <- dat_hobo %>%
  mutate(treat_warm = str_sub(unique_plot_ID, 8, 11),
         treat_comp = str_sub(unique_plot_ID, 13, 16),
         site = str_sub(unique_plot_ID, 5, 6),
         unique_block_id = str_sub(unique_plot_ID, -2),
         date_time = with_tz(date_time, tzone = "Europe/Zurich"),
         date = date(date_time),
         lat = if_else(site == "lo", 46.869000, 46.888000),
         lon = if_else(site == "lo", 9.490000, 9.489000)) 

# get sunrise and sunset times for Calanda
get_sunrise <- dat_hobo %>%
  dplyr::select(date, lat, lon) %>%
  group_by(date, lat, lon) %>%
  distinct()


sunrise_sunset <- get_sunrise %>%
  mutate(sunrise = getSunlightTimes(date, lat, lon, tz = "Europe/Zurich")$sunrise,
         sunset = getSunlightTimes(date, lat, lon, tz = "Europe/Zurich")$sunset)

# add on to hobo data
dat_hobo2 <- dat_hobo %>%
  left_join(sunrise_sunset, by = c("date", "lat", "lon")) %>%
  filter(date_time >= sunrise & date_time <= sunset)

# calculate mean per day
daily_sun_mean <- dat_hobo2 %>%
  group_by(site, treat_comp, treat_warm, date) %>%
  summarize(daily_sun_mean = mean(temperature)) %>%
  ungroup() %>%
  mutate(site_warm = paste(site, treat_warm, sep = "."))

# calculate mean per day, but combined over years
daily_sun_mean_yearly <- dat_hobo2 %>%
  mutate(date_day = format(as.Date(date), "%m-%d")) %>%
  group_by(site, treat_comp, treat_warm, date_day) %>%
  summarize(daily_sun_mean_yearly = mean(temperature)) %>%
  ungroup() %>%
  mutate(site_warm = paste(site, treat_warm, sep = "."))

# calculate difference between warmed and not warmed at high site
daily_sun_mean_diff <- daily_sun_mean %>%
  filter(site == "hi") %>%
  dplyr::select(-site_warm) %>%
  group_by(treat_comp) %>%
  pivot_wider(names_from = treat_warm, 
              values_from = daily_sun_mean,
              names_prefix = "treat_warm_") %>%
  mutate(temp_diff = treat_warm_warm - treat_warm_ambi,
         year = year(date)) %>%
  ungroup() 

# when were the OTCs up?
dat_otc <- data.frame("year" = year(c("2021-06-29", "2021-04-26", "2023-05-11")), "OTC_up" = ymd(c("2021-06-29", "2021-04-26", "2023-05-11")), "OTC_down" = ymd(c("2021-10-16", "2021-10-10", "2023-10-10")))

# filter data for only when OTCs were up for average warming due to OTCs
daily_sun_mean_diff_otc <- daily_sun_mean_diff %>%
  left_join(dat_otc, by = "year") %>%
  filter(date >= OTC_up & date <= OTC_down) %>%
  dplyr::select(-OTC_up, -OTC_down)

# calculate mean and sd for each year
mean_diff_year <- daily_sun_mean_diff_otc %>%
  group_by(treat_comp, year) %>%
  summarize(mean = mean(temp_diff),
            sd = sd(temp_diff))

mean_diff_total <- daily_sun_mean_diff_otc %>%
  group_by(treat_comp) %>%
  summarize(mean = mean(temp_diff),
            sd = sd(temp_diff)) 

daily_sun_mean_diff_otc %>%
  summarize(mean = mean(temp_diff),
            sd = sd(temp_diff))

# 1.03 +- 0.703

# how many loggers where in the end?
hobo_summary <- dat_hobo2 %>%
  dplyr::select(site, unique_block_id, treat_comp, treat_warm) %>%
  distinct() # 28


### TMS4

# separate site, treatments etc., set time zone to Europe/Zurich
dat_tms4 <- dat_tms4 %>%
  mutate(treat_warm = str_sub(unique_plot_ID, 8, 11),
         treat_comp = str_sub(unique_plot_ID, 13, 16),
         site = str_sub(unique_plot_ID, 5, 6),
         unique_block_id = str_sub(unique_plot_ID, -2),
         date_time = with_tz(date_time, tzone = "Europe/Zurich"),
         date = date(date_time),
         lat = if_else(site == "lo", 46.869000, 46.888000),
         lon = if_else(site == "lo", 9.490000, 9.489000)) 

# sunrise data already saved in get_sunrise and sunrise_sunset --> but get again as some additional dates for tms4
get_sunrise_tms4 <- dat_tms4 %>%
  dplyr::select(date, lat, lon) %>%
  group_by(date, lat, lon) %>%
  distinct()

sunrise_sunset_tms4 <- get_sunrise_tms4 %>%
  mutate(sunrise = getSunlightTimes(date, lat, lon, tz = "Europe/Zurich")$sunrise,
         sunset = getSunlightTimes(date, lat, lon, tz = "Europe/Zurich")$sunset)

# add on to tms4 data
dat_tms4.2 <- dat_tms4 %>%
  left_join(sunrise_sunset, by = c("date", "lat", "lon")) %>%
  filter(date_time >= sunrise & date_time <= sunset)

# calculate mean per day
daily_sun_mean_tms4 <- dat_tms4.2 %>%
  group_by(site, treat_comp, treat_warm, date) %>%
  summarize(daily_sun_mean_T1 = mean(TMS_T1),
            daily_sun_mean_T2 = mean(TMS_T2),
            daily_sun_mean_T3 = mean(TMS_T3),
            daily_sun_mean_moist = mean(TMS_moist)) %>%
  ungroup() %>%
  mutate(site_warm = paste(site, treat_warm, sep = "."))

# calculate difference between warmed and not warmed at high site
daily_sun_mean_diff_tms4 <- daily_sun_mean_tms4 %>%
  filter(site == "hi") %>%
  dplyr::select(-site_warm) %>%
  group_by(treat_comp) %>%
  pivot_wider(names_from = treat_warm, 
              values_from = c(daily_sun_mean_T1, daily_sun_mean_T2, daily_sun_mean_T3, daily_sun_mean_moist),
              names_prefix = "treat_warm_") %>%
  mutate(temp_diff_T1 = daily_sun_mean_T1_treat_warm_warm - daily_sun_mean_T1_treat_warm_ambi,
         temp_diff_T2 = daily_sun_mean_T2_treat_warm_warm - daily_sun_mean_T2_treat_warm_ambi,
         temp_diff_T3 = daily_sun_mean_T3_treat_warm_warm - daily_sun_mean_T3_treat_warm_ambi,
         temp_diff_moist = daily_sun_mean_moist_treat_warm_warm - daily_sun_mean_moist_treat_warm_ambi,
         year = year(date)) %>%
  ungroup() 

# when were the OTCs up? --> see at hobo code

# filter data for only when OTCs were up for average warming due to OTCs
daily_sun_mean_diff_otc_tms4 <- daily_sun_mean_diff_tms4 %>%
  left_join(dat_otc, by = "year") %>%
  filter(date >= OTC_up & date <= OTC_down) %>%
  dplyr::select(-OTC_up, -OTC_down)

# calculate mean and sd for each year
mean_diff_year_tms4 <- daily_sun_mean_diff_otc_tms4 %>%
  group_by(treat_comp, year) %>%
  summarize(mean_T1 = mean(temp_diff_T1),
            sd_T1 = sd(temp_diff_T1),
            mean_T2 = mean(temp_diff_T2),
            sd_T2 = sd(temp_diff_T2),
            mean_T3 = mean(temp_diff_T3),
            sd_T3 = sd(temp_diff_T3),
            mean_moist = mean(temp_diff_moist),
            sd_moist = sd(temp_diff_moist))

mean_diff_total_tms4 <- daily_sun_mean_diff_otc_tms4 %>%
  group_by(treat_comp) %>%
  summarize(mean_T1 = mean(temp_diff_T1),
            sd_T1 = sd(temp_diff_T1),
            mean_T2 = mean(temp_diff_T2),
            sd_T2 = sd(temp_diff_T2),
            mean_T3 = mean(temp_diff_T3),
            sd_T3 = sd(temp_diff_T3),
            mean_moist = mean(temp_diff_moist),
            sd_moist = sd(temp_diff_moist)) 

daily_sun_mean_diff_otc_tms4 %>%
  summarize(mean_T1 = mean(temp_diff_T1),
            sd_T1 = sd(temp_diff_T1),
            mean_T2 = mean(temp_diff_T2),
            sd_T2 = sd(temp_diff_T2),
            mean_T3 = mean(temp_diff_T3),
            sd_T3 = sd(temp_diff_T3),
            mean_moist = mean(temp_diff_moist),
            sd_moist = sd(temp_diff_moist))

# calculate mean per day but for complete day, not just sunlight hours
dat_tms4_daily_mean <- dat_tms4 %>%
  group_by(site, treat_comp, treat_warm, date) %>%
  summarize(daily_mean_T1 = mean(TMS_T1),
            daily_mean_T2 = mean(TMS_T2),
            daily_mean_T3 = mean(TMS_T3),
            daily_mean_moist = mean(TMS_moist)) %>%
  ungroup() %>%
  mutate(site_warm = paste(site, treat_warm, sep = "."))


# test whether in general, soil moisture between OTCs up and down is reduced (only sunlight hours)
daily_sun_mean_tms4_hi <- daily_sun_mean_tms4 %>%
  filter(site == "hi") %>%
  mutate(year = year(date)) %>%
  left_join(dat_otc, by = "year") %>%
  filter(date > OTC_up & date < OTC_down)

t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi",]$daily_sun_mean_moist)

mean(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm",]$daily_sun_mean_moist)
mean(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi",]$daily_sun_mean_moist)

t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm" & daily_sun_mean_tms4_hi$treat_comp == "bare",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi" & daily_sun_mean_tms4_hi$treat_comp == "bare",]$daily_sun_mean_moist)
t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm" & daily_sun_mean_tms4_hi$treat_comp == "vege",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi" & daily_sun_mean_tms4_hi$treat_comp == "vege",]$daily_sun_mean_moist)

# test whether in general, soil moisture between OTCs up and down is reduced (whole day)
dat_tms4_daily_mean_hi <- dat_tms4_daily_mean %>%
  filter(site == "hi") %>%
  mutate(year = year(date)) %>%
  left_join(dat_otc, by = "year") %>%
  filter(date > OTC_up & date < OTC_down)

t.test(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm",]$daily_mean_moist,
       dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi",]$daily_mean_moist)

mean(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm",]$daily_mean_moist)
mean(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi",]$daily_mean_moist)

t.test(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm" & dat_tms4_daily_mean_hi$treat_comp == "bare",]$daily_mean_moist,
       dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi" & dat_tms4_daily_mean_hi$treat_comp == "bare",]$daily_mean_moist)
t.test(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist,
       dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist)

mean(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist)
mean(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist)

sd(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "warm" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist)
sd(dat_tms4_daily_mean_hi[dat_tms4_daily_mean_hi$treat_warm == "ambi" & dat_tms4_daily_mean_hi$treat_comp == "vege",]$daily_mean_moist)

# how many loggers per treatment are included in cleaned data? (TABLE S6)
tms_summary <- dat_tms4 %>%
  mutate(year = year(date)) %>%
  left_join(dat_otc, by = "year") %>%
  filter(date > OTC_up & date < OTC_down) %>% # only use data between OTCs up and down
  group_by(site, treat_warm, treat_comp) %>%
  summarize(no_datapoints = n(),
            mean_moist = mean(VWC),
            sd_moist = sd(VWC))

no_loggers <- dat_tms4 %>%
  mutate(year = year(date)) %>%
  left_join(dat_otc, by = "year") %>%
  filter(date > OTC_up & date < OTC_down) %>% # only use data between OTCs up and down
  dplyr::select(site, treat_warm, treat_comp, unique_block_id) %>%
  group_by(site, treat_warm, treat_comp, unique_block_id) %>%
  distinct() %>%
  ungroup() %>%
  group_by(site, treat_warm, treat_comp) %>%
  summarize(no_loggers = n())

## PLOT ------------------------------------------------------------------------

### HOBO

# daily mean for different treatments (only 1 year)
daily_sun_mean_yearly$date_day <- as.Date(daily_sun_mean_yearly$date_day, format = "%m-%d")
daily_sun_mean_yearly_plot <- daily_sun_mean_yearly %>%
  filter(format(date_day, "%m-%d") > "05-01" & format(date_day, "%m-%d") < "10-01")

# FIGURE S2
png("plots/FigS_EnvDat_20251027.png", width = 17, height = 17, units="cm", res=800)

ggplot(dat = daily_sun_mean_yearly_plot, aes(x = date_day, y = daily_sun_mean_yearly, col = site_warm, linetype = treat_comp, group = treat_comp)) +
  #geom_point()
  geom_line() +
  labs(x = "Month", y = "Mean Daytime Temperature") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  #facet_wrap(~site, scales = "free") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face="italic", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  #stat_summary(aes(group=focalorigin), fun=mean, geom="line", size = 3)
  scale_colour_manual(values=c("#4D531B", "#4D531B", "#244B71"),
                      #breaks = c("#4D531B", "#244B71"),
                      labels = c("beyond-range site", "beyond-range site", "within-range site")) +
  scale_linetype_manual(labels = c("bare soil",
                                   "vegetation"),
                        values = c("solid", "longdash")) +
  scale_x_date(breaks = seq(from = min(daily_sun_mean_yearly$date_day), to = max(daily_sun_mean_yearly$date_day), by = "1 month"),  # Set breaks for every month
               labels = scales::date_format("%b")) +
  guides(linetype = "none", colour = "none") +
  facet_wrap(~ site, labeller = labeller(site = c("hi" = "beyond-range site", "lo" = "within-range site")), nrow = 2, strip.position = "right")

dev.off()


# legend
hobo_with_legend <- ggplot(dat = daily_sun_mean_yearly_plot, aes(x = date_day, y = daily_sun_mean_yearly, col = site_warm, linetype = treat_comp, group = treat_comp)) +
  #geom_point()
  geom_line(linewidth = 1.3) +
  labs(x = "Month", y = "Mean Daily Temperature") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  #facet_wrap(~site, scales = "free") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face="italic", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14)) +
  #stat_summary(aes(group=focalorigin), fun=mean, geom="line", size = 3)
  scale_colour_manual(values=c("#4D531B", "#4D531B", "#244B71"),
                      #breaks = c("#4D531B", "#244B71"),
                      labels = c("beyond-range site", "beyond-range site", "within-range site")) +
  scale_linetype_manual(labels = c("bare soil",
                                   "vegetation"),
                        values = c("solid", "dotted")) +
  scale_x_date(breaks = seq(from = min(daily_sun_mean_yearly$date_day), to = max(daily_sun_mean_yearly$date_day), by = "1 month"),  # Set breaks for every month
               labels = scales::date_format("%b")) +
  #guides(linetype = "none", colour = "none") +
  facet_wrap(~ site, labeller = labeller(site = c("hi" = "beyond-range site", "lo" = "within-range site")), nrow = 2)


legend <- g_legend(hobo_with_legend)

ggsave("plots/FigS_EnvDat_legend_20251027.png", legend, width = 9, height = 1, dpi = 300)


### TMS4

# test whether in general, soil moisture between OTCs up and down is reduced
daily_sun_mean_tms4_hi <- daily_sun_mean_tms4 %>%
  filter(site == "hi") %>%
  mutate(year = year(date)) %>%
  left_join(dat_otc, by = "year") %>%
  filter(date > OTC_up & date < OTC_down)

t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi",]$daily_sun_mean_moist)

mean(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm",]$daily_sun_mean_moist)
mean(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi",]$daily_sun_mean_moist)

t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm" & daily_sun_mean_tms4_hi$treat_comp == "bare",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi" & daily_sun_mean_tms4_hi$treat_comp == "bare",]$daily_sun_mean_moist)
t.test(daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "warm" & daily_sun_mean_tms4_hi$treat_comp == "vege",]$daily_sun_mean_moist,
       daily_sun_mean_tms4_hi[daily_sun_mean_tms4_hi$treat_warm == "ambi" & daily_sun_mean_tms4_hi$treat_comp == "vege",]$daily_sun_mean_moist)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SPECIES DISTRIBUTION ON CALANDA  ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# how many plots in Caphe?
length(unique(dat_caphe$plot_id))

max(dat_caphe$altitude_max, na.rm = TRUE)
min(dat_caphe$altitude_min, na.rm = TRUE)


# extract focal species
dat_focal <- dat_caphe %>%
  filter(grepl(c("Salvia p|Medicago lu|Brachypodium pi|Bromus er|Daucus car|Silene vul|Scabiosa col|Centaurea jac|Plantago med|Hypericum per"), taxon_CH),
         !(grepl(c("Silene vulgaris subsp.|Centaurea jacea subsp. a"), taxon_CH)))

unique(dat_focal$taxon_CH)

dat_focal <- dat_focal %>%
  rename("taxon_name" = "taxon_CH")

# check out Sielene
unique(dat_caphe[grepl("Sil", dat_caphe$taxon_CH),]$taxon_CH) # two Silene vulgaris varieties

# change names (combine Brachypodium pinnatum (L.) P. Beauv. and Brachypodium pinnatum aggr. as well as Centaurea jacea L. and Centaurea jacea aggr.)
dat_focal <- dat_focal %>%
  mutate(species = case_when(grepl("Brachypodium pi", taxon_name) ~ "brapin",
                             grepl("Plantago me", taxon_name) ~ "plamed",
                             grepl("Medicago lup", taxon_name) ~ "medlup",
                             grepl("Bromus erec", taxon_name) ~ "broere",
                             grepl("Centaurea jac", taxon_name) ~ "cenjac",
                             grepl("Daucus caro", taxon_name) ~ "daucar",
                             grepl("Silene vulga", taxon_name) ~ "silvul",
                             grepl("Salvia prate", taxon_name) ~ "salpra",
                             grepl("Scabiosa columb", taxon_name) ~ "scacol",
                             grepl("Hypericum perf", taxon_name) ~ "hypper"))



# calculate 0.9 quantile for upper range edge (TABLE S1)
dat_quant <- dat_focal %>%
  group_by(species) %>%
  summarize(quant0.9_elev = round(quantile(altitude_max, probs = 0.9, na.rm = TRUE), 0),
            max_elev = round(max(altitude_max, na.rm = TRUE), 0),
            quant0.1_elev = round(quantile(altitude_max, probs = 0.1, na.rm = TRUE), 0),
            median_elev = median(altitude_max, na.rm = TRUE),
            min_elev = round(min(altitude_max, na.rm = TRUE), 0))

# look at silvul to see distribution
dat_focal %>%
  filter(species == "silvul" & altitude_max > 2000) %>%
  summarize(count = n()) # only 6 occurrences above 2000 masl (out of the 52 total occurrences)

# how many occurrences in total for all species?
no_occ <- dat_focal %>%
  group_by(species) %>%
  summarize(count = n()) # between 25 (scacol) and 168 (brapin)

# what's the overall highest and lowest caphe plot?

min(dat_caphe$altitude_max, na.rm = TRUE)
max(dat_caphe$altitude_max, na.rm = TRUE)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MACROCLIMATE  ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### prepare the data -----------------------------------------------------------

# some infos
names(dat_clim)

names(dat_clim_int)

# get Calanda limits
calanda_1 <- st_transform(calanda_1, crs = "EPSG:4326")
calanda_2 <- st_transform(calanda_2, crs = "EPSG:4326")

common_cols <- intersect(names(calanda_1), names(calanda_2))
calanda_1 <- calanda_1[, common_cols]
calanda_2 <- calanda_2[, common_cols]
calanda_combined <- rbind(calanda_1, calanda_2)
calanda_dissolved <- st_union(calanda_combined)
calanda_mask <- st_cast(st_make_valid(calanda_dissolved), "POLYGON")
plot(st_geometry(calanda_mask), border = "black", lwd = 2, main = "Merged Polygon without Internal Lines")

# get climate data and reproject
cropped_clim <- crop(dat_clim, ext(vect(calanda_mask))) # yearly values
cropped_clim_df <- 
  as.data.frame(cropped_clim, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = mean_1995:mean_2021) %>%
  mutate(year = as.integer(str_extract(name, "\\d{4}")))
cropped_clim_int <- crop(dat_clim_int, ext(vect(calanda_mask))) # integrated values
cropped_clim_int_df <- 
  as.data.frame(cropped_clim_int, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = Y_2000_5_mean:Y_2021_5_mean) %>%
  mutate(year = str_extract(name, "_\\d{4}_"),
         year = as.integer(str_extract(year, "\\d{4}")))


# prepare altitude file (make it coarser to match the climate data)
cropped_dem <- crop(dem, ext(vect(calanda_mask)))
cropped_dem_df <- 
  as.data.frame(cropped_dem, xy = TRUE, na.rm = TRUE)
# apply a focal mean filter with a window size of 3x3
smoothed_dem <- terra::focal(cropped_dem, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
agg_dem <- aggregate(cropped_dem, fact = 100, fun = mean)

# get altitudes
# convert cropped_clim_df to a SpatVector for extraction
clim_coords <- vect(cropped_clim_df[, c("x", "y")], geom = c("x", "y"))
clim_coords_int <- vect(cropped_clim_int_df[, c("x", "y")], geom = c("x", "y"))

# extract the elevation values for each climate pixel
elevation_values <- extract(cropped_dem, clim_coords)
elevation_values_int <- extract(cropped_dem, clim_coords_int)

# add the elevations back on
#cropped_clim_df$elevation <- elevation_values$focal_mean
cropped_clim_df$elevation <- elevation_values$n43_e004_1arc_v3
cropped_clim_int_df$elevation <- elevation_values_int$n43_e004_1arc_v3

unique(cropped_clim_df$elevation)

### model ----------------------------------------------------------------------

# 1) model 2021 separately
mat_2021_df <- cropped_clim_int_df %>%
  filter(year == 2021)

# fit linear model
model_2021 <- lm(value ~ elevation, data = mat_2021_df)
summary(model_2021)

# predict the MAT (temperature) at 2000 masl
predicted_temp_2021 <- predict(model_2021, newdata = data.frame(elevation = 2000))
predicted_temp_2021 # temp. at 2000 masl in 2021

# predict the MAT (temperature) at 1391 masl (highest occurrence of broere)
predicted_temp_2021_broere <- predict(model_2021, newdata = data.frame(elevation = 1391))
predicted_temp_2021_broere

# predict the MAT (temperature) at all highest occurrences of all species
max_elev <- dat_quant$max_elev
predicted_temp_2021_maxelev <- predict(model_2021, newdata = data.frame(elevation = max_elev))
mean(predicted_temp_2021_maxelev)
sd(predicted_temp_2021_maxelev)

# calculate difference between temp. at 2000 masl and temp at other elevations, re-arrange species
dat_quant <- dat_quant %>%
  mutate(predicted_temp = predict(model_2021, newdata = data.frame(elevation = dat_quant$max_elev)),
         temp_diff = predicted_temp - predicted_temp_2021,
         species_ordered = factor(species, levels = c("silvul", "medlup", "plamed", "brapin", "hypper", "scacol",
                                                      "cenjac", "salpra", "daucar", "broere")))

mean(dat_quant$temp_diff)
sd(dat_quant$temp_diff)
max(dat_quant$temp_diff)
min(dat_quant$temp_diff)

# 2) model with interactions for all years

# center elevation at 2000 m so intercepts is directly the MAT at 2000 masl
cropped_clim_int_df_centered <- cropped_clim_int_df %>%
  mutate(year = factor(year),
         elev_c = elevation - 2000)

# pooled model: per-year intercepts and per-year slopes
mod_all <- lm(value ~ elev_c * year, data = cropped_clim_int_df_centered)

# predict MATs at 2000 masl at in different years
pred_values <- tibble(elev_c = 0, # 2000 m
                      year   = factor(c(2000:2021), levels = levels(cropped_clim_int_df_centered$year))) 

pred_values <- bind_cols(pred_values, as.data.frame(predict(mod_all, newdata = pred_values, interval = "confidence"))) %>%
  transmute(year = as.integer(as.character(year)), mat_hat = fit, ci_low = lwr, ci_high = upr)

ggplot(data = pred_values, aes(x = year, y = mat_hat)) +
  geom_point()

# post-2010 linear trend (°C/decade)
trend_2010p <- lm(mat_hat ~ year, data = filter(pred_values, year >= 2010))
slope_decade <- coef(trend_2010p)["year"] * 10
ci_decade    <- confint(trend_2010p)["year", ] * 10
slope_decade; ci_decade  # headline number + 95% CI


### plot -----------------------------------------------------------------------

### PLOT RANGE FOCALS (FIGURE 1)

# manually define axis breaks
elevation_breaks <- seq(0, max(dat_quant$quant0.9_elev), by = 500)
# add lower limits
label_data <- no_occ %>%
  mutate(y_position = 420, 
         label = paste(count))

sp_order <- c("daucar", "medlup", "cenjac", "hypper", "plamed", "broere", "salpra", "scacol", "silvul", "brapin")
dat_quant$species <- factor(dat_quant$species, levels = sp_order)

png("plots/Fig1_focals_20251210.png", width = 15, height = 10, units="cm", res=800) # wide: 20, not wide: 13

ggplot(data = dat_quant) +
  geom_linerange(aes(x = species, y = quant0.9_elev, ymin = quant0.1_elev, ymax = quant0.9_elev), 
                 linewidth = 15, col = "#DEDFE4") +
  coord_cartesian(ylim = c(550, 2600)) +
  geom_abline(intercept = 2000, slope = 0, col = "#5B7D6D", size = 2, linetype = "dashed") +
  geom_abline(intercept = 1400, slope = 0, col = "#C23A2C", size = 2, linetype = "dashed") +
  #geom_abline(intercept = elevation_at_predicted_temp_2011, slope = 0, col = "darkred", size = 1.3) +
  #geom_abline(intercept = elevation_at_predicted_temp_2001, slope = 0, col = "#C1C2C7", size = 1) +
  geom_point(aes(x = species_ordered, y = median_elev), size = 4) +
  geom_point(aes(x = species_ordered, y = max_elev), col = "white", size = 3) +
  geom_point(aes(x = species_ordered, y = max_elev)) +
  #geom_point(aes(x = species, y = min_elev)) +
  #geom_errorbar(aes(x = species, ymin = quant0.1_elev, ymax = quant0.9_elev), width = 0.2) +
  theme_bw() +
  labs(y = "Elevation (m a.s.l.)") +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_y_continuous("Elevation (m a.s.l.)", breaks = elevation_breaks#, 
                     #sec.axis = sec_axis(~ (. - a)/b, name = "Temperature difference to\nMAT at 2000 m a.s.l. (°C)", breaks = temp_diff_breaks)
  ) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, face = "italic"),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank()) +
  geom_text(data = label_data, 
            aes(x = species, y = y_position, label = label), 
            size = 4.5, hjust = 0.5, vjust = -1, fontface = "italic") +
  scale_y_continuous("Elevation (m a.s.l.)", 
                     breaks = seq(500, 2500, by = 500)#,   # Specify breaks at 500, 1000, 1500, 2000, 2500
                     #sec.axis = sec_axis(~ (. - a)/b, name = "Temperature difference to\nMAT at 2000 m a.s.l. (°C)", breaks = temp_diff_breaks)
  )

dev.off()

### PLOT SWITZERLAND

# get shapefile
switzerland <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name == "Switzerland")

png("plots/Fig1_Switzerland_20250107.png", width = 8, height = 5, units="cm", res=800)

ggplot() +
  geom_sf(data = switzerland, fill = "white", color = "black", size = 10, linewidth = 1) + 
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())

dev.off()

### PLOT EUROPE

# get shapefile
europe <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

sf_use_s2(FALSE)

# land outline (low detail)
land <- ne_download(scale = "small", type = "land",
                    category = "physical", returnclass = "sf") %>%
  st_zm(drop = TRUE) %>%
  st_make_valid()

# crop to Europe
bbox <- st_as_sfc(st_bbox(c(xmin = -15, ymin = 34, xmax = 35, ymax = 72),
                          crs = st_crs(land)))
europe <- st_intersection(land, bbox) %>%
  st_union() # dissolve borders for one silhouette

# add vertices (densify) then smooth (round)
europe_dense  <- smoothr::densify(europe, max_distance = 0.5)   # adds points ~0.5° apart
europe_smooth <- smoothr::smooth(europe_dense, method = "chaikin")

# simplfied Switzerland
switzerland_small <- ne_countries(scale = "small", returnclass = "sf") %>% 
  filter(name == "Switzerland")

# plot

png("plots/Fig1_Europe_20251015.png", width = 5, height = 8, units="cm", res=800)

ggplot() +
  geom_sf(data = europe_smooth, fill = "#DEDFE4", color = "#DEDFE4", linewidth = 0.1) +
  geom_sf(data = switzerland,   fill = "black", color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(-15, 35), ylim = c(34, 72), expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())

dev.off()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VEGETATION HEIGHT  ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dat_height <- dat_height %>%
  mutate(treat_warm = str_sub(unique_plot_ID, 8, 11),
         treat_comp = str_sub(unique_plot_ID, 13, 16),
         site = str_sub(unique_plot_ID, 5, 6)) %>%
  filter(grepl("wf", unique_plot_ID) == TRUE)


dat_height_summary <- dat_height %>%
  group_by(site) %>%
  summarize(mean_height = mean(veg_height),
            sd_height = sd(veg_height))

t.test(dat_height[dat_height$site == "lo",]$veg_height, 
       dat_height[dat_height$site == "hi",]$veg_height)














         

