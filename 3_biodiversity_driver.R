# Statement ----
# pairwise Pearson correlations between biodiversity indexes and influence factors for plants, trees, and shrubs respectively 
# get "best linear regression models" between biodiversity indexes and influence factors with all possible regression method 

# Package ----
source("1_package_and_data.R")
source("2_biodiversity.R")

library(leaps)

# Data ----
kLvlDrvFac <- 
  c("dist_ctr", "land_price", "pop_tot", 
    "agriculture", "water_wetland", "park", "temple_shrine", 
    "cemetery", "vacant", "golf_course", 
    "residential", "multi_family residential", "institutional", 
    "transportation", "commercial_neighbor", "commercial_industrial")

# Analysis ----
## Pearson correlation between biodiversity and indexes ----
Cortest <- function(x, col_dep) {
  # 建立列表以存储相关分析结果
  cor_res <- vector("list", length = 3)
  
  # 各自变量和因变量的相关性作图和计算
  for (i in kLvlDrvFac) {
    # 如果样本量大于10则计算，否则直接将结果定为NA
    if(length(na.omit(x[[i]])) > 10) {
      # plot(x[[i]], x[[col_dep]], xlab = i, ylab = col_dep)
      
      res <- cor.test(x[[i]], x[[col_dep]])
      cor_res[[1]] <- c(cor_res[[1]], i)
      cor_res[[2]] <- c(cor_res[[2]], res$estimate)
      cor_res[[3]] <- c(cor_res[[3]], res$p.value)
    } else {
      cor_res[[1]] <- c(cor_res[[1]], i)
      cor_res[[2]] <- c(cor_res[[2]], NA)
      cor_res[[3]] <- c(cor_res[[3]], NA)
    }
  }
  
  # 将结果转化为数据框
  cor_res <- tibble(
    variable = cor_res[[1]], 
    correlation = cor_res[[2]], 
    p = cor_res[[3]]
  )
  # 生成显著性标记数据列
  cor_res$p_mark <- ""
  cor_res$p_mark[which(cor_res$p < 0.05)] <- "*"
  cor_res$p_mark[which(cor_res$p < 0.01)] <- "**"
  cor_res$p_mark[which(cor_res$p < 0.001)] <- "***"
  
  cor_res$result <- 
    paste0(round(cor_res$correlation, digits = 2), cor_res$p_mark)
  return(cor_res)
}

# 建立列表以储存皮尔森相关性分析结果
prsn.res <- vector("list", 9)
names(prsn.res) <- paste(
  rep(c("plant", "tree", "shrub"), each = 3), 
  rep(c("richness", "abundance", "evenness"), 3), sep = "_")

# 全体植物皮尔森分析结果
prsn.res$plant_richness <- 
  Cortest(qua.div, "richness")[c("variable", "result")] %>% 
  rename(plant_richness = result)

# 乔木皮尔森分析结果
for (i in c("richness", "abundance", "evenness")) {
  prsn.res[[paste0("tree_", i)]] <- 
    Cortest(qua.tree.comm.div, i)[c("variable", "result")]
  names(prsn.res[[paste0("tree_", i)]])[2] <- paste0("tree_", i)
}

# 灌木皮尔森分析结果
for (i in c("richness", "abundance", "evenness")) {
  prsn.res[[paste0("shrub_", i)]] <- 
    Cortest(qua.shrub.comm.div, i)[c("variable", "result")]
  names(prsn.res[[paste0("shrub_", i)]])[2] <- paste0("shrub_", i)
}

# 合并各项结果
Merge <- function(x, y) {
  output <- merge(x, y, by = "variable")
  return(output)
}
prsn.res <- Reduce(Merge, prsn.res[c(1, 4:9)])

# check and output the results 
prsn.res
write.xlsx(prsn.res, "ProcData/1Biodiversity/Prsn_plant_indexes.xlsx")

## All possible regression for best models ----
# 备选变量：仅考虑上一步皮尔森检测中数据量大于10个的变量
GetRegSubset <- function(x, var_response) {
  leaps <- regsubsets(
    richness ~ dist_ctr + land_price + pop_tot + agriculture + residential +
      multi_family_residential + commercial_industrial + institutional + park + 
      transportation, data = x)
  plot(leaps, scale = "adjr2", main = var_response)
}

# 函数：生成结果数据框
LmRes2Df <- function(x) {
  ressum <- summary(x)
  ressum <- ressum$coefficients
  # 提取变量名、估计值和p值
  output <- data.frame(
    variable = names(ressum[-1, 1]), 
    estimate = ressum[-1, 1], 
    p = ressum[-1, 4]
  )
  
  rownames(output) <- NULL
  output$p_mark <- ""
  output$p_mark[which(output$p < 0.05)] <- "*"
  output$p_mark[which(output$p < 0.01)] <- "**"
  output$p_mark[which(output$p < 0.001)] <- "***"
  
  output$result <- 
    paste0(round(output$estimate, digits = 2), output$p_mark)
  output <- output[c("variable", "estimate", "result")]
  return(output)
}

# 基于选入变量建模
# 生成空列表以存储各配对最优线性模型
fitbest <- vector("list", 9)
names(fitbest) <- names(prsn.res)

# 生成各数据集各指标最佳模型所选入的变量
# 基于图像选择最佳模型自变量
GetRegSubset(qua.div, "richness")
fitbest$plant_richness <- 
  lm(richness ~ dist_ctr + pop_tot + 
       agriculture + residential + multi_family_residential + 
       commercial_industrial, data = qua.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(plant_richness = result)

GetRegSubset(qua.tree.comm.div, "richness")
fitbest$tree_richness <- 
  lm(richness ~ land_price + pop_tot + 
       agriculture + multi_family_residential + 
       commercial_industrial + institutional + transportation, 
     data = qua.tree.comm.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>% 
  rename(tree_richness = result)

GetRegSubset(qua.shrub.comm.div, "richness")
fitbest$shrub_richness <- 
  lm(richness ~ residential + multi_family_residential + commercial_industrial, 
     data = qua.shrub.comm.div) %>%
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(shrub_richness = result)

GetRegSubset(qua.tree.comm.div, "abundance")
fitbest$tree_abundance <- 
  lm(richness ~ land_price + pop_tot + 
       agriculture + multi_family_residential + 
       commercial_industrial + institutional + transportation, 
     data = qua.tree.comm.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(tree_abundance  = result)

GetRegSubset(qua.shrub.comm.div, "abundance")
fitbest$shrub_abundance <- 
  lm(abundance ~ dist_ctr + land_price + pop_tot +agriculture + residential + 
       multi_family_residential + 
       commercial_industrial + park, 
     data = qua.shrub.comm.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(shrub_abund = result)

GetRegSubset(qua.tree.comm.div, "evenness")
fitbest$tree_evenness <- 
  lm(evenness ~ dist_ctr + pop_tot +agriculture + residential + 
       multi_family_residential + institutional + 
       park +  transportation, 
     data = qua.tree.comm.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(tree_even = result)

GetRegSubset(qua.shrub.comm.div, "evenness")
fitbest$shrub_evenness <- 
  lm(evenness ~ land_price + pop_tot +agriculture + residential + 
       multi_family_residential + commercial_industrial + institutional + 
       institutional +  transportation, 
     data = qua.shrub.comm.div) %>% 
  LmRes2Df() %>% 
  select(variable, result) %>%
  rename(shrub_even = result)

Merge <- function(x, y) {
  output <- merge(x, y, by = "variable", all = TRUE)
  return(output)
}
lmdf <- fitbest[!sapply(fitbest, is.null)] %>% 
  Reduce(Merge, x = .) %>% 
  tibble()

# check and output the results
lmdf
write.xlsx(lmdf, "ProcData/1Biodiversity/Best_lm_models.xlsx")
