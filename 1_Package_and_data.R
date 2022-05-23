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
  c("ResLow", "ResHigh", "ResOther", "Ind", "ComNbr", "Com")

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
  select(Qua_id, KES_plot_id, Landuse_class, Lat, Long) %>% 
  rename_with(tolower) %>% 
  rename(landuse = landuse_class, 
         kes_quaid = kes_plot_id) %>% 
  tibble() %>% 
  mutate(dist_ctr = distm(
    x = matrix(c(long, lat), ncol = 2), y = c(135.76928, 35.00213))[, 1]) %>% 
  full_join(socio_economic, by = "qua_id") %>% 
  full_join(land_cover, by = "qua_id") %>% 
  rename_with(~ gsub("/", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("-", "_", .x, fixed = TRUE))
qua_data[is.na(qua_data)] <- 0

## Species info ----
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
all_qua_div <- all_plant_data %>% 
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

# Data proc ----
## KUP data ----
# 函数：构建群落宽数据集
fun_comm <- function(x, nq_colgroup, nq_colabund) {
  output <- x %>%
    select({{nq_colgroup}}, species_lt, {{nq_colabund}}) %>%
    pivot_wider(names_from = species_lt, values_from = {{nq_colabund}},
                values_fn = sum, values_fill = 0)
  return(output)
}

# 函数：分组计算属性汇总统计值和多样性指数
fun_div <- function(x, x_comm, nq_colabund, nq_colgroup) {
  # 内置函数：分组汇总统计各项属性
  funin_attrcalc <- function(coltar, tarvalue) {
    x_sub <- x
    x_sub["tarornot"] <- x_sub[coltar] == tarvalue
    x_sub <- x_sub %>% group_by({{nq_colgroup}}) %>% 
      summarise(
        perc = sum({{nq_colabund}} * tarornot) / sum({{nq_colabund}})) %>%
      ungroup() %>% 
      select({{nq_colgroup}}, perc)
    names(x_sub)[2] <- paste0("perc_", tarvalue)
    return(x_sub)
  }
  
  perc_planted <- funin_attrcalc("pla_spo", "planted")
  perc_nonpot <- funin_attrcalc("pla_spo", "non_pot")
  perc_private <- funin_attrcalc("pla_spo", "private")
  perc_nonstreet <- funin_attrcalc("pla_spo", "non_street")
  perc_native <- funin_attrcalc("pla_spo", "native")
  
  output <- x_comm %>%
    mutate(abundance = rowSums(.[3:ncol(.)]),
           richness = apply(.[2:ncol(.)]>0, 1, sum),
           shannon = diversity(.[2:ncol(.)], index = "shannon"),
           simpson = diversity(.[2:ncol(.)], index = "simpson"),
           evenness = shannon / log(richness)) %>%
    select({{nq_colgroup}}, 
           abundance, richness, shannon, simpson, evenness) %>%
    left_join(perc_planted) %>%
    left_join(perc_nonpot) %>%
    left_join(perc_private) %>%
    left_join(perc_nonstreet) %>%
    left_join(perc_native)
  
  return(output)
}

tree_lu_comm <- fun_comm(tree_data, landuse, stem)
tree_lu_div <- fun_div(tree_data, tree_lu_comm, stem, landuse)

tree_qua_comm <- fun_comm(tree_data, qua_id, stem)
tree_qua_div <- fun_div(tree_data, tree_qua_comm, stem, qua_id)
tree_qua_all <- tree_qua_comm %>% 
  left_join(tree_qua_div) %>% 
  left_join(qua_data)

shrub_lu_comm <- fun_comm(shrub_data, landuse, area)
shrub_lu_div <- fun_div(shrub_data, shrub_lu_comm, area, landuse)

shrub_qua_comm <- fun_comm(shrub_data, qua_id, area)
shrub_qua_div <- fun_div(shrub_data, shrub_qua_comm, area, qua_id)
shrub_qua_all <- shrub_qua_comm %>% 
  left_join(shrub_qua_div) %>% 
  left_join(qua_data)

# some other variables 
number_plant_species <- length(unique(all_plant_data$species_lt))
number_tree_species <- length(unique(tree_data$species_lt))
number_shrub_species <- length(unique(shrub_data$species_lt))

## KES data ----
# factor level of annual ES
lvl_es_annual <- c("carbon_seq", 
               "no2_removal", "o3_removal", "pm25_removal", "so2_removal",
               "avo_runoff")

# info species code 
info_species_code <- read.csv("RawData/iTree_species_list.csv") %>% 
  mutate(species = paste(Genus, Species.Name)) %>% 
  rename(species_code = SppCode, 
         genus = Genus, 
         species_name = Species.Name, 
         common_name = Common.Name) %>%
  select(species_code, species, common_name) %>% 
  subset(!duplicated(species_code)) %>% 
  subset(!duplicated(species)) %>% 
  subset(!duplicated(common_name))
info_species_code$species[info_species_code$species == " "] <- NA

# i-Tree input data: from Dr. Hirabayashi
itree_input <- read.csv("RawData/KES i-Tree_input.csv") %>% 
  select(ID, PlotId, TreeId, FieldLandUse, TreeStatus, 
         Species, TreeHeightTotal, CrownWidth1, CrownWidth2, CrownLightExposure,
         PercentCrownMissing, PercentImperviousBelow, PercentShrubBelow) %>% 
  rename(res_tree_id = ID, 
         kes_quaid = PlotId, 
         in_tree_id = TreeId, 
         spo_pla = TreeStatus, 
         species_code = Species, 
         height = TreeHeightTotal, 
         crown_width_ew = CrownWidth1, 
         crown_width_ns = CrownWidth2, 
         light_expo = CrownLightExposure, 
         per_crow_mis = PercentCrownMissing, 
         per_impervious_below = PercentImperviousBelow, 
         per_shrub_below = PercentShrubBelow
  ) %>% 
  left_join(qua_data, by = "kes_quaid") %>% 
  left_join(info_species_code, by = "species_code")

# individual ES data from Excel file
# The Excel file is from Access data
tree_ind_es <- read.xlsx("RawData/KES Tree_es_output.xlsx", sheet = "Trees") %>% 
  rename(res_tree_id = "TreeID", 
         dbh = "DBH.(CM)", 
         lai = "LEAF.AREA.INDEX", 
         carbon_storage = "CARBON.STORAGE.(KG)", 
         carbon_seq = `GROSS.CARBON.SEQ.(KG/YR)`, 
         biomass = `FUTURE.BIOMASS.(KG)`,
         co_removal = "CO.Removal.(g)",
         no2_removal = "NO2.Removal.(g)", 
         o3_removal = "O3.Removal.(g)", 
         pm25_removal = "PM25.Removal.(g)", 
         so2_removal = "SO2.Removal.(g)", 
         compensatory_value = "TREE.VALUE.(Yen)",          
         no2_value = "NO2.Value.($)", 
         o3_value = "O3.Value.($)", 
         pm25_value = "PM25.Value.($)",          
         so2_value = "SO2.Value.($)", 
         avo_runoff = "Avoided.Runoff.(m3)") %>% 
  left_join(itree_input, by = "res_tree_id") %>% 
  mutate(compensatory_value = compensatory_value/100, 
         carbon_storage_value = 188/1000*carbon_storage, 
         carbon_seq_value = 188/1000*carbon_seq, 
         avo_runoff_value = 2.36*avo_runoff, 
         landuse = factor(
           landuse, 
           levels = c("ResLow", "ResHigh", "ResOther", "Ind", "ComNbr", "Com"))
  ) %>% 
  group_by(res_tree_id) %>% 
  mutate(lvl_es_annual_value = sum(carbon_seq_value, 
                               no2_value, o3_value, pm25_value, so2_value,  
                               avo_runoff_value)) %>% 
  ungroup() %>% 
  mutate(
    dbh_class = case_when(
      dbh <= 15 ~ "0 < DBH ≤ 15", 
      dbh <= 30 ~ "15 < DBH ≤ 30", 
      dbh > 30 ~ "DBH > 30"
    ), 
    lai_class = case_when(
      lai <= 3 ~ "0 < LAI ≤ 3", 
      lai <= 6 ~ "3 < LAI ≤ 6", 
      lai <= 9 ~ "6 < LAI ≤ 9", 
      lai > 9 ~ "LAI > 9" 
    )
  ) %>% 
  select(qua_id, res_tree_id, kes_quaid, in_tree_id, species_code, species, 
         common_name, 
         spo_pla, dbh, dbh_class, height, crown_width_ew, crown_width_ns,
         per_crow_mis, light_expo, per_shrub_below, per_impervious_below, 
         lai, lai_class, biomass, 
         landuse, 
         carbon_storage, carbon_seq, 
         no2_removal, o3_removal, pm25_removal, so2_removal, co_removal, 
         avo_runoff, 
         compensatory_value, 
         carbon_storage_value, carbon_seq_value, 
         no2_value, o3_value, pm25_value, so2_value,  
         avo_runoff_value, 
         lvl_es_annual_value)

# Quadrat data
tree_qua_es <- tree_ind_es %>% 
  select(qua_id, kes_quaid, dbh, lai, biomass, 
         carbon_storage, carbon_seq, 
         no2_removal, o3_removal, pm25_removal, so2_removal, co_removal, 
         avo_runoff, 
         compensatory_value, 
         carbon_storage_value, carbon_seq_value, 
         no2_value, o3_value, pm25_value, so2_value, 
         avo_runoff_value, 
         lvl_es_annual_value) %>% 
  group_by(qua_id) %>% 
  summarise(across(!starts_with("qua_id"), sum), 
            treenum = n()) %>% 
  ungroup() %>%  
  left_join(qua_data[c("qua_id", "landuse")], by = "qua_id")

## BEF data ----
# 将生物量数据加入乔木的全数据中
tree_qua_all <- tree_qua_all %>% 
  left_join(
    select(tree_qua_es, -kes_quaid, -landuse), by = "qua_id")

# levels of richness and abundance 
tree_qua_all$richness_level <- "High.rich"
tree_qua_all$richness_level[which(
  tree_qua_all$richness < quantile(tree_qua_all$richness, 0.66))] <- "Mid.rich"
tree_qua_all$richness_level[which(
  tree_qua_all$richness < quantile(tree_qua_all$richness, 0.33))] <- "Low.rich"
tree_qua_all$richness_level <- 
  factor(tree_qua_all$richness_level, 
         levels = c("Low.rich", "Mid.rich", "High.rich"))

tree_qua_all$abundance_level <- "High.abund"
tree_qua_all$abundance_level[which(
  tree_qua_all$abundance < quantile(tree_qua_all$abundance, 0.66))] <- "Mid.abund"
tree_qua_all$abundance_level[which(
  tree_qua_all$abundance < quantile(tree_qua_all$abundance, 0.33))] <- "Low.abund"
tree_qua_all$abundance_level <- 
  factor(tree_qua_all$abundance_level, 
         levels = c("Low.abund", "Mid.abund", "High.abund"))
