# Package ----
library(leaps)

# Data ----
lvl_drv_fac <- c("dist_ctr", "land_price", "popttotal", 
             "residential", "transportation", "temple/shrine", 
             "multi-family residential", "agriculture", 
             "commercial-neighbor", "water/wetland", "park", 
             "cemetery", "vacant", "commercial/industrial", "institutional", 
             "golf course")

# Pearson correlation between variables and quadrat diversity index ----
fun_cortest <- function(x, col_dep) {
  # 建立列表以存储相关分析结果
  cor_res <- vector("list", length = 3)
  
  # 各自变量和因变量的相关性作图和计算
  for (i in lvl_drv_fac) {
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
prsn_res <- vector("list", 9)
names(prsn_res) <- paste(
  rep(c("plant", "tree", "shrub"), each = 3), 
  rep(c("richness", "abundance", "evenness"), 3), sep = "_")

# 全体植物皮尔森分析结果
prsn_res$plant_richness <- 
  fun_cortest(all_qua_div, "richness")[c("variable", "result")] %>% 
  rename(plant_richness = result)

# 乔木皮尔森分析结果
for (i in c("richness", "abundance", "evenness")) {
  prsn_res[[paste0("tree_", i)]] <- 
    fun_cortest(tree_qua_all, i)[c("variable", "result")]
  names(prsn_res[[paste0("tree_", i)]])[2] <- paste0("tree_", i)
}

# 灌木皮尔森分析结果
for (i in c("richness", "abundance", "evenness")) {
  prsn_res[[paste0("shrub_", i)]] <- 
    fun_cortest(shrub_qua_all, i)[c("variable", "result")]
  names(prsn_res[[paste0("shrub_", i)]])[2] <- paste0("shrub_", i)
}

# 合并各项结果
fun_merge <- function(x, y) {
  output <- merge(x, y, by = "variable")
  return(output)
}
prsn_res_df <- Reduce(fun_merge, prsn_res[c(1, 4:9)])

# 输出结果
write.xlsx(prsn_res_df, "Out prsn_allplant_richness.xlsx")

# 全子集回归分析 ----
# 备选变量：仅考虑上一步皮尔森检测中数据量大于10个的变量
fun_regsubset <- function(x, var_response) {
  leaps <- regsubsets(
    richness ~ dist_ctr + land_price + popttotal + agriculture + residential +
      multi_family_residential + commercial_industrial + institutional + park + 
      transportation, data = x)
  plot(leaps, scale = "adjr2", main = var_response)
}

# 函数：生成结果数据框
fun_lmdf <- function(x) {
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
names(fitbest) <- names(prsn_res)

# 生成各数据集各指标最佳模型所选入的变量
# 基于图像选择最佳模型自变量
fun_regsubset(all_qua_div, "richness")
fitbest$plant_richness <- 
  lm(richness ~ dist_ctr + popttotal + 
       agriculture + residential + multi_family_residential + 
       commercial_industrial, data = all_qua_div) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(plant_richness = result)

fun_regsubset(tree_qua_all, "richness")
fitbest$tree_richness <- 
  lm(richness ~ land_price + popttotal + 
       agriculture + multi_family_residential + 
       commercial_industrial + institutional + transportation, 
     data = tree_qua_all) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>% 
  rename(tree_richness = result)

fun_regsubset(shrub_qua_all, "richness")
fitbest$shrub_richness <- 
  lm(richness ~ residential + multi_family_residential + commercial_industrial, 
     data = shrub_qua_all) %>%
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(shrub_richness = result)

fun_regsubset(tree_qua_all, "abundance")
fitbest$tree_abundance <- 
  lm(richness ~ land_price + popttotal + 
       agriculture + multi_family_residential + 
       commercial_industrial + institutional + transportation, 
     data = tree_qua_all) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(tree_abundance  = result)

fun_regsubset(shrub_qua_all, "abundance")
fitbest$shrub_abundance <- 
  lm(abundance ~ dist_ctr + land_price + popttotal +agriculture + residential + 
       multi_family_residential + 
       commercial_industrial + park, 
     data = shrub_qua_all) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(shrub_abund = result)

fun_regsubset(tree_qua_all, "evenness")
fitbest$tree_evenness <- 
  lm(evenness ~ dist_ctr + popttotal +agriculture + residential + 
       multi_family_residential + institutional + 
       park +  transportation, 
     data = tree_qua_all) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(tree_even = result)

fun_regsubset(shrub_qua_all, "evenness")
fitbest$shrub_evenness <- 
  lm(evenness ~ land_price + popttotal +agriculture + residential + 
       multi_family_residential + commercial_industrial + institutional + 
       institutional +  transportation, 
     data = shrub_qua_all) %>% 
  fun_lmdf() %>% 
  select(variable, result) %>%
  rename(shrub_even = result)

fun_merge <- function(x, y) {
  output <- merge(x, y, by = "variable", all = TRUE)
  return(output)
}
lmdf <- Reduce(fun_merge, fitbest[!sapply(fitbest, is.null)]) %>% 
  tibble()

write.xlsx(lmdf, "Out Best lm model.xlsx")

