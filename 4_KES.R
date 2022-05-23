library(car)
library(gplots)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RODBC)
library(tidyr)
library(openxlsx)
library(dunn.test)

# Analysis ----
## DBH structure ----
# structure across land use types: calculatd based on land use scale
ggplot(tree_ind_es) + 
  geom_bar(aes(landuse, fill = dbh_class), color = "grey", position = "fill") + 
  labs(x = "Land use class", y = "Proportion", fill = "DBH class") + 
  theme_bw()y

## LAI structure ----
# structure across land use types: calculated based on land use scale
ggplot(tree_ind_es) + 
  geom_bar(aes(landuse, fill = lai_class), color = "grey", position = "fill") + 
  labs(x = "Land use class", y = "Proportion", fill = "LAI class") + 
  theme_bw()

## Qua avg ESs value ~ land use ----
quavalue_summary <- tree_qua_es %>% 
  select(landuse, carbon_seq_value, 
         no2_value, o3_value, pm25_value, so2_value, 
         avo_runoff_value) %>% 
  pivot_longer(cols = c(carbon_seq_value, 
                        no2_value, o3_value, pm25_value, so2_value, 
                        avo_runoff_value), 
               names_to = "es", values_to = "es_value") %>%
  group_by(landuse, es) %>% 
  summarise(es_value = mean(es_value)) %>% 
  ungroup() %>%
  mutate(es = factor(
    es, levels = c("carbon_seq_value", 
                   "no2_value", "o3_value", "pm25_value", "so2_value", 
                   "avo_runoff_value")))
ggplot(quavalue_summary, aes(landuse, es_value)) + 
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
ggplot(quavalue_summary, aes(landuse, es_value)) + 
  geom_bar(aes(fill = es), stat = "identity", position = "fill")

## Quadrat and individual ESs ~ land use ---- 
# test the assumptions for statistical analysis
apply(as.data.frame(tree_qua_es[es_annual]), 2, 
      function(x) {shapiro.test(x)$p.value > 0.05})
apply(as.data.frame(tree_ind_es[es_annual]), 2, 
      function(x) {shapiro.test(x)$p.value > 0.05})

# function for non-parametric statistical analysis 
fun_comparison <- function(var_es, name_depend_var, name_independ_var) {
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
      var_loop_tukey <- 
        dunn.test(var_es[, var_loop_colname], var_es[, name_independ_var], 
                  method = "Bonferroni")
      cat("Tukey result: \n")
      # print(subset(as.data.frame(var_loop_tukey[[1]]), `p adj` < 0.05))
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
  output <- data.frame(chi = chi, pvalue = pvalue)
  output$mark <- ""
  output$mark[which(output$pvalue < 0.05)] <- "*"
  output$mark[which(output$pvalue < 0.01)] <- "**"
  output$mark[which(output$pvalue < 0.001)] <- "***"
  
  return(output)
}

# quadrat ES ~ land use
fun_comparison(tree_qua_es, es_annual, "landuse")
# individual ES ~ land use
fun_comparison(tree_ind_es, es_annual, "landuse")

## Qua structure ~ land use ----
# test assumption of normality 
apply(as.data.frame(tree_qua_es[c("dbh", "lai", "treenum")]), 2, 
      function(x) {shapiro.test(x)$p.value > 0.05})
apply(as.data.frame(tree_ind_es[c("dbh", "lai")]), 2, 
      function(x) {shapiro.test(x)$p.value > 0.05})

# Kruskal-Wallis rank sum test
fun_comparison(tree_qua_es, c("dbh", "lai", "treenum"), "landuse")
fun_comparison(tree_ind_es, c("dbh", "lai"), "landuse")


## Carbon seq ~ carbon storage ----
# at quadrat scale 
ggplot(tree_qua_es) + 
  geom_point(aes(carbon_storage, carbon_seq, color = landuse), alpha = 0.5) 
summary(lm(carbon_seq ~ carbon_storage, data = tree_qua_es))

# at individual scale
ggplot(tree_ind_es) + 
  geom_point(aes(carbon_storage, carbon_seq, color = species), alpha = 0.5) + 
  guides(color = "none")
summary(lm(carbon_seq ~ carbon_storage, data = tree_ind_es))

# Species-specific analysis ----
# function to select target species and land use 
func_var_sub <- function(var_es, name_gp, name_subgp, num_sample, num_subgp) {
  var_es <- as.data.frame(var_es)
  # each pair with sample size larger than or equal to 2
  var_gp_subgp_ct <- table(var_es[,name_gp], var_es[,name_subgp]) %>% 
    as.data.frame() %>% 
    subset(Freq >= num_sample)
  # each group with larger than or equal to 3 subgroups
  var_gp_ct <- var_gp_subgp_ct %>% group_by(Var1) %>% summarise(n = n()) %>% 
    subset(n >= num_subgp)
  # knock out pairs of group-subgroup with not enough sample size  
  var_gp_subgp_tar <- subset(var_gp_subgp_ct, Var1 %in% var_gp_ct$Var1)
  # get required subset from original var_es data frame
  subset(var_es, paste0(var_es[, name_gp], var_es[, name_subgp]) %in% 
           paste0(var_gp_subgp_tar$Var1, var_gp_subgp_tar$Var2))
}
# individual ES ~ land use 
func_var_sub(tree_ind_es, "species", "landuse", 3, 4) %>% 
  split(.$species) %>% 
  lapply(fun_comparison, es_annual, "landuse") %>% 
  write.xlsx("Output_Species-specific_comparison.xlsx")

