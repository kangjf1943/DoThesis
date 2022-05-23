library(priceTools)
# if the package priceTools doesn't exist, install it with: 
# devtools::install_github("ctkremer/priceTools")
library(dplyr)
library(ggplot2)
library(openxlsx)
library(patchwork)

# Linear model ----
# 建立空列表用于存储线性模型结果
lm_sfac_bef <- vector("list", 2)
names(lm_sfac_bef) <- c("biomass_richness", "biomass_abundance")

## Single factor ----
# biomass~richness
lm_sfac_bef$biomass_richness <- lm(biomass ~ richness, data = tree_qua_all)
summary(lm_sfac_bef$biomass_richness)
# biomass~abundance
lm_sfac_bef$biomass_abundance <- lm(biomass ~ abundance, data = tree_qua_all)
summary(lm_sfac_bef$biomass_abundance)
# 作图
plt_lm_sfac_bef <- vector("list", 2)
names(plt_lm_sfac_bef) <- c("biomass_richness", "biomass_abundance")

plt_lm_sfac_bef$biomass_richness <- 
  ggplot(tree_qua_all, aes(richness, biomass/1000)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  labs(x = "Richness", y = "Biomass (t)") + 
  theme_bw()
plt_lm_sfac_bef$biomass_abundance <- 
  ggplot(tree_qua_all, aes(abundance, biomass/1000)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  labs(x = "Abundance", y = "Biomass (t)") + 
  theme_bw()
plt_lm_sfac_bef$biomass_richness + 
  plt_lm_sfac_bef$biomass_abundance

## Double factors ----
fun_2flm <- function(x, name_group_fac, name_var) {
  x_ls <- split(x, f = x[, name_group_fac])
  
  # 用于储存等级和p值结果
  lvl <- vector()
  p <- vector()
  for (i in x_ls) {
    fit <- lm(i$biomass ~ i[[name_var]])
    res <- summary(fit)
    lvl <- c(lvl, as.character(unique(i[, name_group_fac])))
    p <- c(p, res$coefficients[, 4][2])
  }
  
  # 整理输出结果数据框
  output <- data.frame(
    lvl = lvl, 
    p = p
  )
  return(output)
}

res2flm_richness <- fun_2flm(tree_qua_all, "abundance_level", "richness")
res2flm_abundance <- fun_2flm(tree_qua_all, "richness_level", "abundance")
res2flm_richness
res2flm_abundance

# plots
plt_lm_mfac_bef <- vector("list", 2)
names(plt_lm_mfac_bef) <- c("biomass_richness", "biomass_abundance")

plt_lm_mfac_bef$biomass_richness <- 
  ggplot(tree_qua_all, aes(richness, biomass/1000)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  facet_grid(.~ as.factor(abundance_level)) + 
  labs(title = "(a)", x = "Richness", y = "Biomass (t)") + 
  theme_bw()
plt_lm_mfac_bef$biomass_abundance <- 
  ggplot(tree_qua_all, aes(abundance, biomass/100)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  facet_grid(.~ as.factor(richness_level)) + 
  labs(title = "(b)", x = "Abundance", y = "Biomass (t)") + 
  theme_bw()
Reduce("/", plt_lm_mfac_bef)

## Interact with land use ----
res2flm_richness_landuse <- fun_2flm(tree_qua_all, "land_use", "richness")
res2flm_abundance_landuse <- fun_2flm(tree_qua_all, "land_use", "abundance")
write.csv(res2flm_richness_landuse, "AnaData/res2flm_richness_landuse.csv")
write.csv(res2flm_abundance_landuse, "AnaData/res2flm_abundance_landuse.csv")
res2flm_richness_landuse
res2flm_abundance_landuse

# plots
plt_lm_landuse_bef <- vector("list", 2)
names(plt_lm_landuse_bef) <- c("biomass_richness", "biomass_abundance")

plt_lm_landuse_bef$biomass_richness <- 
  ggplot(tree_qua_all, aes(richness, biomass/1000)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  facet_grid(.~ as.factor(landuse)) + 
  labs(title = "(a)", x = "Richness", y = "Biomass (t)") + 
  theme_bw()
plt_lm_landuse_bef$biomass_abundance <- 
  ggplot(tree_qua_all, aes(abundance, biomass/1000)) + 
  geom_point(color = "black", alpha = 0.3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  facet_grid(.~ as.factor(landuse)) + 
  labs(title = "(b)", x = "Abundance", y = "Biomass (t)") + 
  theme_bw()
Reduce("/", plt_lm_landuse_bef)

# Price equation ----
# 待办：有一组功能相同的配对
fun_pairprice <- function(x) {
  x <- subset(x, select = c("qua_id", "species_code", "biomass"))
  x <- x %>% group_by(.dots = c("qua_id"))
  pairprice <- pairwise.price(x, species="species_code", func="biomass")
  
  pairprice$SI <- pairprice$SIE.G + pairprice$SIE.L
  pairprice <- 
    subset(pairprice, x.func > y.func, 
           select = c("SR", "SI", "CDE", "x.func", "y.func"))
  pairprice$EF <- pairprice$SR + pairprice$SI + pairprice$CDE
  pairprice$SR_std <- pairprice$SR / pairprice$x.func
  pairprice$SI_std <- pairprice$SI / pairprice$x.func
  pairprice$CDE_std <- pairprice$CDE / pairprice$x.func
  pairprice$EF_std <- pairprice$EF / pairprice$x.func
  pairprice_lng <- rbind(
    tibble(effect = rep('SR', nrow(pairprice)), 
           value = pairprice$SR_std), 
    tibble(effect = rep('SI', nrow(pairprice)), 
           value = pairprice$SI_std), 
    tibble(effect = rep('CDE', nrow(pairprice)), 
           value = pairprice$CDE_std), 
    tibble(effect = rep('EF', nrow(pairprice)), 
           value = pairprice$EF_std)
  )
  return(pairprice_lng)
}

# 生成空列表以存储图片
plt_price <- vector("list", 2)
names(plt_price) <- c("price_all", "price_landuse")

# for all the trees
res_pp_all <- fun_pairprice(tree_ind_es)
res_pp_all$effect <- 
  factor(res_pp_all$effect, levels = c("EF", "SR", "SI", "CDE"))
plt_price$price_all <- 
  ggplot(res_pp_all) + geom_violin(aes(effect, value)) + 
  labs(x = "Effect", y = "Effect value") + 
  theme_bw()

# interact with land use
res_pp_landuse <- lapply(
  split(tree_ind_es, tree_ind_es$landuse), 
  fun_pairprice
)
res_pp_landuse <- lapply(
  res_pp_landuse, 
  function(x) {
    x$effect <- factor(x$effect, levels = c("EF", "SR", "SI", "CDE"))
    return(x)
  }
)

for (i in names(res_pp_landuse)) {
  res_pp_landuse[[i]]$landuse <- i
}
res_pp_landuse_df <- Reduce(rbind, res_pp_landuse) %>% 
  mutate(effect = factor(effect, levels = c("EF", "SR", "SI", "CDE")))

plt_price$price_landuse <- 
  ggplot(res_pp_landuse_df) + 
  geom_violin(aes(effect, value)) + 
  facet_wrap(~landuse) + 
  labs(x = "Effect", y = "Effect value") + 
  theme_bw()

plot_ls <- lapply(res_pp_landuse, function(x) {
  ggplot(x) + geom_violin(aes(effect, value)) + 
    labs(x = "Effect", y = "Effect value") + 
    theme_bw()
})
for (i in names(res_pp_landuse)) {
  plot_ls[[i]] <- plot_ls[[i]] + labs(title = i)
}
Reduce("+", plot_ls[1:3]) /
  Reduce("+", plot_ls[4:6])

# Biomass across species ----
# whether rare species generally have smaller biomass?
tree_ind_es$species_code <- 
  factor(tree_ind_es$species_code, 
         levels = names(table(tree_ind_es$species_code))[
           order(table(tree_ind_es$species_code))])
ggplot(tree_ind_es) + 
  geom_bar(aes(species_code))
ggplot(tree_ind_es) + 
  geom_boxplot(aes(species_code, biomass/1000)) + 
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(x = "Species rank from low to high by abundance", y = "Biomass (t)")
