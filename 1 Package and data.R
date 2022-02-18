# Packages ----
library(Rmisc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(PerformanceAnalytics)
library(dunn.test)
library(BiodiversityR)
# if problem occurs with BiodiversityR package, try to install the tool XQuartz
library(codyn)
library(ggpubr)
library(iNEXT)
library(openxlsx)
library(leaps)
library(geosphere)

# Settings ----
lvl_landuse <- 
  c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow")

# Get data ----
## Quadrat info ----
# socio-economic data
socio_economic <- read.xlsx("RawData/GIS Socio_economic.xlsx") %>%
  select(Qua_id, land_price, PopTTotal) %>% 
  rename_with(tolower)

# land cover proportion data
land_cover <- read.xlsx("RawData/GIS Quadrat_land_cover.xlsx") %>% 
  select(Quadrat_ID, Detailed_land_cover, Shape_Area) %>% 
  rename_with(tolower) %>% 
  rename(qua_id = quadrat_id) %>% 
  group_by(qua_id, detailed_land_cover) %>% 
  summarise(area = sum(shape_area)/400) %>% 
  pivot_wider(
    id_cols = qua_id, names_from = detailed_land_cover, values_from = area)

# other quadrat info
qua_data <- read.xlsx("RawData/KUP Quadrat_info.xlsx", "Qua_info") %>% 
  select(Qua_id, Landuse_class, Lat, Long) %>% 
  rename_with(tolower) %>% 
  rename(landuse = landuse_class) %>% 
  tibble() %>% 
  mutate(dist_ctr = distm(
    x = matrix(c(long, lat), ncol = 2), y = c(135.76928, 35.00213))[, 1]) %>% 
  full_join(socio_economic, by = "qua_id") %>% 
  full_join(land_cover, by = "qua_id")

## Plant info ----
# information of all the plants species: provenance and taxonomy 
all_plant_info <- 
  read.csv("RawData/KUP Plant_info.csv", stringsAsFactors = FALSE) %>% 
  rename_with(tolower) %>% 
  tibble()

## Plant data ----
# data of all the plants: taxonomy, attributes, abundance, etc.
all_plant_data <- 
  read.csv("RawData/KUP Plant_data.csv", stringsAsFactors = FALSE) %>%
  rename_with(tolower) %>% 
  rename(qua_id = plot_id) %>% 
  left_join(all_plant_info, by = "species_lt") %>%
  left_join(qua_data,by = "qua_id") %>% 
  tibble()
qua_plant_div <- all_plant_data %>% 
  mutate(presence = 1) %>% 
  select(qua_id, species_lt, presence) %>% 
  rename_with(tolower) %>% 
  pivot_wider(names_from = species_lt, values_from = presence, 
              values_fn = sum, values_fill = 0) %>% 
  mutate(richness = apply(.[2:ncol(.)]>0, 1, sum)) %>% 
  left_join(qua_data, by = "qua_id")

# data of trees and shrubs and diversity table of trees and shrubs 
tree_data <- subset(all_plant_data, tree_shrub == "tree") %>% 
  mutate(Area = NULL)
shrub_data <- subset(all_plant_data, tree_shrub == "shrub") %>% 
  mutate(stem = NULL) %>% 
  mutate(area = area/10000)

# delete unused objects 
rm(socio_economic, land_cover)
