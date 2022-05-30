# Package ----
source("1_package_and_data.R")

library(car)
library(gplots)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RODBC)
library(tidyr)
library(openxlsx)
library(dunn.test)

# Get data ----
# factor level of annual ES
kLvlEsAnnual <- 
  c("carbon_seq", 
    "no2_removal", "o3_removal", "pm25_removal", "so2_removal",
    "avo_runoff")

# info species code 
info.spe.code <- read.csv("RawData/iTree_species_list.csv") %>% 
  as_tibble() %>% 
  mutate(species = paste(Genus, Species.Name)) %>% 
  rename(species_code = SppCode, 
         genus = Genus, 
         species_name = Species.Name, 
         common_name = Common.Name) %>%
  select(species_code, species, common_name) %>% 
  subset(!duplicated(species_code)) %>% 
  subset(!duplicated(species)) %>% 
  subset(!duplicated(common_name)) %>% 
  # replace the null value of species with NA
  mutate(species = ifelse(species == " ", NA, species))

# i-Tree input data: from Dr. Hirabayashi
itree.input <- read.csv("RawData/KES i-Tree_input.csv") %>% 
  as_tibble() %>% 
  rename_with(tolower) %>% 
  select(id, plotid, treeid, species, 
         treeheighttotal, crownwidth1, crownwidth2, 
         crownlightexposure, percentcrownmissing) %>% 
  rename(res_tree_id = id, 
         kes_qua_id = plotid, 
         in_tree_id = treeid, 
         species_code = species, 
         height = treeheighttotal, 
         light_expo = crownlightexposure, 
         pct_crow_mis = percentcrownmissing
  ) %>% 
  # calculate average crown width of two directions 
  mutate(crown_width = (crownwidth1 + crownwidth2) / 2) %>% 
  select(-crownwidth1, -crownwidth2) %>% 
  # add scientific name of species
  left_join(select(info.spe.code, species_code, species), 
            by = "species_code") %>%
  select(-species_code) %>% 
  # add information of land use
  left_join(select(qua.data, kes_qua_id, land_use), by = "kes_qua_id")

# individual ES data from Excel file, which is from Access data
ind.es <- read.xlsx("RawData/KES Tree_es_output.xlsx", sheet = "Trees") %>% 
  as_tibble() %>% 
  rename_with(~ tolower(gsub(".", "_", .x, fixed = TRUE))) %>% 
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("/", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("__", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("_$", "", .x, fixed = TRUE)) %>% 
  select(-spjapanese, -crownwidthew, -crownwidthns, -baseht, -spcode, 
         -height_m, -ground_area_m2, -tree_condition, -leaf_area_m2, 
         -leaf_biomass_kg, -future_dbh_cm, -future_height_m, 
         -biomass_adjustment, -leaf_type) %>% 
  rename(res_tree_id = treeid, 
         dbh = dbh_cm, 
         lai = leaf_area_index, 
         carbon_storage = carbon_storage_kg, 
         carbon_seq = gross_carbon_seq_kg_yr, 
         biomass = future_biomass_kg,
         co_removal = co_removal_g,
         no2_removal = no2_removal_g, 
         o3_removal = o3_removal_g, 
         pm25_removal = pm25_removal_g, 
         so2_removal = so2_removal_g, 
         compensatory_value = tree_value_yen,          
         no2_value = no2_value, 
         o3_value = o3_value, 
         pm25_value = pm25_value,          
         so2_value = so2_value, 
         avo_runoff = avoided_runoff_m3) %>% 
  left_join(itree.input, by = "res_tree_id") %>% 
  mutate(
    compensatory_value = compensatory_value/100, 
    carbon_storage_value = 188/1000*carbon_storage, 
    carbon_seq_value = 188/1000*carbon_seq, 
    avo_runoff_value = 2.36*avo_runoff, 
    land_use = factor(land_use, levels = kLvlLanduse)
  ) %>% 
  group_by(res_tree_id) %>% 
  mutate(annual_es_value = sum(across(kLvlEsAnnual))) %>% 
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
      lai > 6 ~ "LAI > 9"
    )
  )

# Quadrat data
qua.es <- ind.es %>% 
  select(kes_qua_id, dbh, lai, biomass, 
         carbon_storage, carbon_seq, 
         no2_removal, o3_removal, pm25_removal, so2_removal, co_removal, 
         avo_runoff, 
         compensatory_value, 
         carbon_storage_value, carbon_seq_value, 
         no2_value, o3_value, pm25_value, so2_value, 
         avo_runoff_value, 
         annual_es_value) %>% 
  group_by(kes_qua_id) %>% 
  summarise(across(!starts_with("kes_qua_id"), sum), 
            treenum = n()) %>% 
  ungroup() %>%  
  left_join(select(qua.data, kes_qua_id, land_use), by = "kes_qua_id")

# Analysis ----
## DBH structure ----
# structure across land use types: calculatd based on land use scale
ggplot(ind.es) + 
  geom_bar(aes(land_use, fill = dbh_class), 
           color = "black", position = "fill") + 
  scale_fill_manual(
    values = c("#ffffff", "#b8b8b8", "#242424")) + 
  labs(x = "Land use class", y = "Proportion", fill = "DBH class") + 
  theme_bw()

## LAI structure ----
# structure across land use types: calculated based on land use scale
ggplot(ind.es) + 
  geom_bar(aes(land_use, fill = lai_class), 
           color = "black", position = "fill") + 
  scale_fill_manual(
    values = c("#ffffff", "#b8b8b8", "#242424")) + 
  labs(x = "Land use class", y = "Proportion", fill = "LAI class") + 
  theme_bw()

## Qua avg ESs value ~ land use ----
qua.val.sum <- qua.es %>% 
  select(land_use, carbon_seq_value, 
         no2_value, o3_value, pm25_value, so2_value, 
         avo_runoff_value) %>% 
  pivot_longer(cols = c(carbon_seq_value, 
                        no2_value, o3_value, pm25_value, so2_value, 
                        avo_runoff_value), 
               names_to = "es", values_to = "es_value") %>%
  group_by(land_use, es) %>% 
  summarise(es_value = mean(es_value)) %>% 
  ungroup() %>%
  mutate(es = factor(
    es, levels = c("carbon_seq_value", 
                   "no2_value", "o3_value", "pm25_value", "so2_value", 
                   "avo_runoff_value")))
ggplot(qua.val.sum, aes(land_use, es_value)) + 
  geom_bar(aes(fill = es), stat = "identity", position = "stack") + 
  scale_fill_discrete(
    limits = c("carbon_seq_value", "no2_value", 
               "o3_value", "pm25_value", 
               "so2_value", "avo_runoff_value"), 
    labels = c("Carbon sequestration", "NO2 removal", 
               "O3 removal", "PM2.5 removal", 
               "SO2 removal", "Avoided runoff")) + 
  labs(x = "Land use", y = "Quadrat ecosystem service values (dollars / year)", 
       fill = "Ecosystem service") + 
  theme_bw()
# export Excel file for plot in article 
qua.val.sum %>% 
  pivot_wider(id_cols = land_use, names_from = es, values_from = es_value) %>% 
  mutate(land_use = factor(land_use, levels = kLvlLanduse)) %>% 
  arrange(land_use) %>% 
  write.xlsx("ProcData/2EcosystemService/Average_value_es.xlsx")
# proportion of the es values 
ggplot(qua.val.sum, aes(land_use, es_value)) + 
  geom_bar(aes(fill = es), stat = "identity", position = "fill")

## Quadrat and individual ESs ~ land use ---- 
# test the assumptions for statistical analysis
for (i in kLvlLanduse) {
  cat("\n", i)
  apply(as.data.frame(qua.es[which(qua.es$land_use == i), kLvlEsAnnual]), 2, 
        function(x) {shapiro.test(x)$p.value > 0.05}) %>% 
    sum() %>% 
    print()
}
for (i in kLvlLanduse) {
  cat("\n", i)
  apply(as.data.frame(ind.es[which(ind.es$land_use == i), kLvlEsAnnual]), 2, 
        function(x) {shapiro.test(x)$p.value > 0.05}) %>% 
    sum() %>% 
    print()
}

# function for non-parametric statistical analysis 
Comparison <- function(var_es, name_depend_var, name_independ_var) {
  var_es <- as.data.frame(var_es)
  if (length(unique(var_es$species)) == 1) {
    var_species_name <- var_es$species[1]
    print(var_species_name)
  } else {
    var_species_name <- ""}
  
  # store the statistics: chi-square and p value
  chi <- vector("numeric")
  pvalue <- vector("numeric")
  
  # prepare the canvas for boxplots 
  par(mfrow = c(floor(sqrt(length(name_depend_var))),
                ceiling(sqrt(length(name_depend_var)))))
  
  # get statistics, post hoc comparison, and plots 
  for (var_loop_colname in name_depend_var) {
    print(var_loop_colname)
    
    var_loop_fit <- 
      kruskal.test(var_es[, var_loop_colname] ~ var_es[, name_independ_var])
    chi <- c(chi, var_loop_fit$statistic)
    var_loop_aov_pvalue <- var_loop_fit$p.value
    pvalue <- c(pvalue, var_loop_aov_pvalue)
    
    if (var_loop_aov_pvalue < 0.05) {
      var_loop_dunn <- 
        dunn.test(var_es[, var_loop_colname], var_es[, name_independ_var], 
                  method = "Bonferroni")
      cat("\n")
      boxplot(var_es[, var_loop_colname] ~ var_es[, name_independ_var], 
                ylab = var_loop_colname, xlab = "", las = 2, 
                main = paste0(var_species_name, "\n", 
                              "anova p-value = ", round(var_loop_aov_pvalue,2)))
    } else {
      boxplot(var_es[, var_loop_colname] ~ var_es[, name_independ_var], 
                ylab = var_loop_colname, xlab = "", las = 2, 
                main = paste0(var_species_name, "\n", 
                              "anova p-value = ", round(var_loop_aov_pvalue,2)),
                border = "grey")
    }
  }
  cat("\n\n")
  par(mfrow = c(1,1))
  
  # data frame of the statistics 
  output <- data.frame(depend_var = name_depend_var, chi = chi, pvalue = pvalue)
  output$mark <- ""
  output$mark[which(output$pvalue < 0.05)] <- "*"
  output$mark[which(output$pvalue < 0.01)] <- "**"
  output$mark[which(output$pvalue < 0.001)] <- "***"
  
  return(output)
}

# quadrat ES ~ land use
Comparison(qua.es, kLvlEsAnnual, "land_use")
# individual ES ~ land use
Comparison(ind.es, kLvlEsAnnual, "land_use")

## Qua structure ~ land use ----
# test assumption of normality 
for (i in kLvlLanduse) {
  cat(
    i, 
    subset(qua.es, land_use == i) %>%  
      as.data.frame() %>% 
      select(dbh, lai, treenum) %>% 
      apply(., 2, function(x) {shapiro.test(x)$p.value > 0.05}) %>% 
      sum(), 
    "\n"
  )
}
for (i in kLvlLanduse) {
  cat(
    i, 
    subset(ind.es, land_use == i) %>%  
      as.data.frame() %>% 
      select(dbh, lai) %>% 
      apply(., 2, function(x) {shapiro.test(x)$p.value > 0.05}) %>% 
      sum(), 
    "\n"
  )
}

# Kruskal-Wallis rank sum test
Comparison(qua.es, c("dbh", "lai", "treenum"), "land_use")
Comparison(ind.es, c("dbh", "lai"), "land_use")

## Carbon seq ~ carbon storage ----
# at quadrat scale 
ggplot(qua.es) + 
  geom_point(aes(carbon_storage, carbon_seq, color = land_use), alpha = 0.5) 
summary(lm(carbon_seq ~ carbon_storage, data = qua.es))

# at individual scale
ggplot(ind.es) + 
  geom_point(aes(carbon_storage, carbon_seq, color = species), alpha = 0.5) + 
  guides(color = "none")
summary(lm(carbon_seq ~ carbon_storage, data = ind.es))

# Species-specific analysis ----
# function to select target species and land use 
SubBySpecies <- function(var_es, name_gp, name_subgp, num_sample, num_subgp) {
  var_es <- as.data.frame(var_es)
  # each pair with sample size larger than or equal to 2
  var_gp_subgp_ct <- table(var_es[,name_gp], var_es[,name_subgp]) %>% 
    as.data.frame() %>% 
    subset(Freq >= num_sample)
  # each group with larger than or equal to 3 subgroups
  var_gp_ct <- var_gp_subgp_ct %>% 
    group_by(Var1) %>% 
    summarise(n = n()) %>% 
    subset(n >= num_subgp)
  # knock out pairs of group-subgroup with not enough sample size  
  var_gp_subgp_tar <- subset(var_gp_subgp_ct, Var1 %in% var_gp_ct$Var1)
  # get required subset from original var_es data frame
  subset(var_es, paste0(var_es[, name_gp], var_es[, name_subgp]) %in% 
           paste0(var_gp_subgp_tar$Var1, var_gp_subgp_tar$Var2))
}
# individual ES ~ land use 
SubBySpecies(ind.es, "species", "land_use", 3, 4) %>% 
  split(.$species) %>% 
  lapply(Comparison, kLvlEsAnnual, "land_use") %>% 
  write.xlsx("ProcData/2EcoService/Species_specific_comparison.xlsx")
