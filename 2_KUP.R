# Functions ----
# function: get "wide data" for community
FunComm <- function(x, nq_colgroup, nq_colabund) {
  x %>%
    select({{nq_colgroup}}, species, {{nq_colabund}}) %>%
    pivot_wider(names_from = species, values_from = {{nq_colabund}},
                values_fn = sum, values_fill = 0) %>% 
    return()
}

# function: summarise attribute and calculate diversity indexes
FunDiv <- function(x, x_comm, nq_colabund, nq_colgroup) {
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

# func to plot species accumulation curve and extrapolation 
# func parameters: x, raw data; y, size of extrapolation; z, extent of x axis; method, plot the curve separately by land use (method = "land_use") or at city level (method = "city")
FunAccum <- function(x, y, z, method) {
  # func to derive list of incidence data of certain land use
  # func para: xtemp, raw data; ytemp, name of land use
  fun_temp <- function(xtemp, ytemp) {
    xtemp %>% filter(land_use %in% ytemp) %>% 
      group_by(species) %>% 
      summarise(n = length(unique(qua_id))) %>% 
      # arrange(desc(n)) %>% 
      .$n %>% 
      c(length(unique(
        xtemp[which(xtemp[["land_use"]] %in% ytemp), ][["qua_id"]])), .)
  }
  if (method == "city") {
    incidence <- list(fun_temp(x, kLvlLanduse))
    names(incidence) <- "city"
  } else {
    incidence <- list(
      fun_temp(x, "Com"), 
      fun_temp(x, "ComNbr"),
      fun_temp(x, "Ind"), 
      fun_temp(x, "ResOther"), 
      fun_temp(x, "ResHigh"), 
      fun_temp(x, "ResLow")
    )
    names(incidence) <- 
      c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow")
  }
  
  accum <- iNEXT(incidence, q = 0, datatype = "incidence_freq", 
                 size = seq(1:y), se = FALSE)
  if (method == "city") {
    accum$iNextEst$city$land_use <- "city"
    accum <- accum$iNextEst$city
  } else {
    accum$iNextEst$Com$land_use <- "Com"
    accum$iNextEst$"ComNbr"$land_use <- "ComNbr"
    accum$iNextEst$Ind$land_use <- "Ind"
    accum$iNextEst$"ResOther"$land_use <- "ResOther"
    accum$iNextEst$"ResHigh"$land_use <- "ResHigh"
    accum$iNextEst$"ResLow"$land_use <- "ResLow"
    accum <- Reduce(rbind, accum$iNextEst)
    accum$land_use <- factor(
      accum$land_use, 
      levels = c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow"))
  }
  accum$method[accum$method == "interpolated"] <- "observed"
  
  hline <- accum %>% group_by(land_use) %>% 
    summarise(asymtote = max(qD))
  print(arrange(hline, desc(asymtote)))
  
  if (method == "city") {
    plotdata <- ggplot(accum) + 
      geom_line(aes(t, qD, color = land_use, linetype = method), size = 1) +
      geom_hline(data = hline, aes(yintercept = asymtote, color = land_use), 
                 linetype = 2, size = 1) + 
      scale_linetype_discrete(limits = c("observed", "extrapolated"))
  } else {
    plotdata <- ggplot(accum[which(accum$method == "observed"), ]) + 
      geom_line(aes(t, qD, color = land_use), size = 1) +
      geom_hline(data = hline, aes(yintercept = asymtote, color = land_use), 
                 linetype = 2, size = 1)
  }
  plotdata + 
    labs(x = "Number of quadrats", y = "Number of species") + 
    scale_x_continuous(limits = c(0, z)) + 
    theme_bw()
}

# func to generate top species with abundance data
# func para: x, raw data; y, level of classification; z, column of abundance data; k, calculation of total abundance; n, number of top "n"
FunTop <- function(x, y, z, k, n) {
  top_plant <- x %>% group_by(get(y)) %>% 
    dplyr::summarise(number = sum(get(z)), prop = number/k) %>% 
    arrange(desc(number)) %>% 
    head(n)
  names(top_plant)[1] = c(y)
  print(top_plant)
}

# func to test if top species are contained in top families
# func para: x, raw data of species; y, raw data of families 
FunContain <- function(x, y) {
  x %>% left_join(plant.info, by = "species") %>% 
    select(species, family) %>% 
    left_join(y, by = "family") %>% 
    select(species, family)
}

FunRankData <- function(x, name_spe) {
  x %>% select({{name_spe}}) %>%
    subset(select = colSums(.) != 0) %>% 
    as.data.frame() %>% 
    rankabundance() %>% 
    as.data.frame() %>%
    mutate(species = rownames(.)) %>% 
    as_tibble()
}

# func to plot rank-abundance curve
# func para: x, raw data; title, plot title; method, plot at city level (method = "city") or at land use level (method = "land_use)
FunRankPlot <- function(x, title, method) {
  if (method == "city") {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + labs(title = title)
  } else {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + 
      facet_wrap(~land_use, nrow = 1) + labs(title = title)
  }
  plotdata + 
    geom_point(aes(color = nt_ex), alpha = 0.3, size = 2) + 
    labs(x = "Scaled rank of species", y = "Relative abundance (%)") +
    scale_color_discrete("Provenance") +
    theme(strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
}

# Get data ----
## Tree data ----
lu.tree.comm <- FunComm(tree.data, land_use, stem)
lu.tree.div <- FunDiv(tree.data, lu.tree.comm, stem, land_use)

qua.tree.comm <- FunComm(tree.data, qua_id, stem)
qua.tree.div <- FunDiv(tree.data, qua.tree.comm, stem, qua_id)
qua.tree.comm.div <- qua.tree.comm %>% 
  left_join(qua.tree.div) %>% 
  left_join(qua.data)

## Shrub data ----
lu.shrub.comm <- FunComm(shrub_data, land_use, area)
lu.shrub.div <- FunDiv(shrub_data, lu.shrub.comm, area, land_use)

qua.shrub.comm <- FunComm(shrub_data, qua_id, area)
qua.shrub.div <- FunDiv(shrub_data, qua.shrub.comm, area, qua_id)
qua.shrub.comm.div <- qua.shrub.comm %>% 
  left_join(qua.shrub.div) %>% 
  left_join(qua.data)

# Constant ----
KNumPlantSpecies <- length(unique(plant.data$species))
KNumTreeSpecies <- length(unique(tree.data$species))
KNumShrubSpecies <- length(unique(shrub_data$species))

KLvlIndex <- c("abundance", "richness", "shannon", "simpson", "evenness")
KTreeSpecies <- unique(tree.data$species)
KShrubSpecies <- unique(shrub_data$species)

# Analysis ----
## City level ----
# Number of species
cat("\n", "total species:", length(unique(plant.data$species)), "\n", 
    "total genera:", length(unique(plant.data$genus)), "\n", 
    "total families:", length(unique(plant.data$family)), "\n")

# species accumulation curve 
# extrapolation up to double the reference sample size
FunAccum(plant.data, 348, 350, method = "city") + 
  theme(legend.position = "none")

## Top taxa ----
# top species families of all plants by species number
plant.info %>% group_by(family) %>% 
  dplyr::summarise(Num_spe = n(), Prop = n()/nrow(plant.info)) %>% 
  arrange(desc(Prop))

# abundance of trees and shrubs
cat("\n", "number of trees:", nrow(tree.data), "in", 
    nrow(qua.tree.div), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$area), "m2 in", 
    nrow(qua.shrub.div), "plot")

# top species families of trees and shrubs by abundance
tree.top.species <- FunTop(
  tree.data, "species", "stem", sum(tree.data$stem), 10)
tree.top.family <- FunTop(
  tree.data, "family", "stem", sum(tree.data$stem), 10)
FunContain(tree.top.species, tree.top.family)

shrub.top.species <- FunTop(
  shrub_data, "species", "area", sum(shrub_data$area), 10)
shrub.top.family <- FunTop(
  shrub_data, "family", "area", sum(shrub_data$area), 10)
intersect(tree.top.family$family, shrub.top.family$family)
FunContain(shrub.top.species, shrub.top.family)

# attributes of trees and shrubs
# the number of exotic vs. native by number of species
table(plant.info$nt_ex)/nrow(plant.info)

# the attributes of trees and shrubs
for (i in c("pla_spo", "pub_pri", "nt_ex")) {
  print(tapply(tree.data$stem, tree.data[,i], sum)/sum(tree.data$stem), 
        digits = 2)
  cat("\n")
}
for (i in c("pla_spo", "pub_pri", "nt_ex")) {
  print(tapply(shrub_data$area, shrub_data[,i], sum)/sum(shrub_data$area), 
        digits = 2)
  cat("\n")
}

### Distribution of species abundance ----
# rank abundance curves for trees and shrubs
city.tree.rank <- FunRankData(qua.tree.comm, KTreeSpecies) %>% 
  left_join(select(plant.info, c("species", "nt_ex")), by = "species")
city.shrub.rank <- FunRankData(qua.shrub.comm, KShrubSpecies) %>% 
  left_join(select(plant.info, c("species", "nt_ex")), by = "species")

ggarrange(plotlist = list(
  FunRankPlot(city.tree.rank, title = "(a)", method = "city"), 
  FunRankPlot(city.shrub.rank, title = "(b)", method = "city")
), nrow = 1, common.legend = TRUE)

# the top 3 species regarding abundance
subset(city.tree.rank[,c("rank", "species", "nt_ex")], rank <= 3)
subset(city.shrub.rank[,c("rank", "species", "nt_ex")], rank <= 3)

# EQ evenness index and plot
city.tree.rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>% 
  arrange(EQ)
city.shrub.rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>%
  arrange(EQ)

# Land use level ----
## Species accumulation curve ----
FunAccum(plant.data, 348, 50, method = "land_use") + 
  labs(color = "Land use") + 
  scale_color_manual(values = c("#FF0000", "#FF7800", "#DF73FF", 
                                "#BFBF30", "#6BE400", "#00733E"))

## Distribution of species abundance ----
lu.tree.rank <- 
  left_join(qua.tree.comm, qua.data[c("qua_id", "land_use")], by = "qua_id") %>% 
  ddply("land_use", name_spe = KTreeSpecies, FunRankData) %>% 
  left_join(select(plant.info, c("species", "nt_ex")), by = "species")
lu.shrub.rank <- 
  left_join(qua.shrub.comm, qua.data[c("qua_id", "land_use")], by = "qua_id") %>%
  ddply("land_use", name_spe = KShrubSpecies, FunRankData) %>% 
  left_join(select(plant.info, c("species", "nt_ex")), by = "species")
ggarrange(FunRankPlot(lu.tree.rank, "(a)", method = "land_use"),
          FunRankPlot(lu.shrub.rank, "(b)", method = "land_use"), 
          nrow = 2, common.legend = TRUE, legend = "bottom")

# the top 3 species regarding abundance
subset(lu.tree.rank[,c("land_use", "rank", "species", "nt_ex")], rank <= 3)
subset(lu.shrub.rank[,c("land_use", "rank", "species", "nt_ex")], rank <= 3)

# calculate the EQ evenness index and plot
community_structure(
  lu.tree.rank, time.var = "land_use", 
  abundance.var = "abundance", metric = "EQ") %>% 
  arrange(desc(EQ))
community_structure(
  lu.shrub.rank, time.var = "land_use", 
  abundance.var = "abundance", metric = "EQ") %>%
  arrange(desc(EQ))

## Species composition ----
# Bray-Curtis dissimilarity of pairs of land use for trees and shrubs
lu.tree.comm %>% 
  select(-land_use) %>% 
  vegdist() %>% 
  as.matrix() %>% 
  round(digits = 2) %>% 
  as_tibble() %>% 
  rename_all(~kLvlLanduse) %>% 
  mutate(pairs = kLvlLanduse) %>% 
  relocate(pairs)

lu.shrub.comm %>% 
  select(-land_use) %>% 
  vegdist() %>% 
  as.matrix() %>% 
  round(digits = 2) %>% 
  as_tibble() %>% 
  rename_all(~kLvlLanduse) %>% 
  mutate(pairs = kLvlLanduse) %>% 
  relocate(pairs)

# plant occupancy of species for different land use types
FunOccupRate <- function(x,y) {
  apply(x[,2:(y+1)], 2, function(k)sum(ifelse(k>0,1,0))/length(k))
}
FunOccupDf <- function(x){
  data.frame("Com" = names(head(sort(x[,"Com"],decreasing = TRUE),10)), 
             "ComNbr" = names(head(sort(x[,"ComNbr"],decreasing = TRUE),10)),
             "ResLow" = names(head(sort(x[,"ResLow"],decreasing = TRUE),10)),
             "ResHigh" = names(head(sort(x[,"ResHigh"],decreasing = TRUE),10)), 
             "ResOther" = names(head(sort(x[,"ResOther"],decreasing = TRUE),10)), 
             "Ind" = names(head(sort(x[,"Ind"],decreasing = TRUE),10)))
}

plant.occup <- ddply(
  qua.div, .(land_use), y = KNumPlantSpecies, FunOccupRate) %>% 
  .[,-1] %>% t() 
colnames(plant.occup) = kLvlLanduse
(plant.occup.top <- FunOccupDf(plant.occup))

# shared species over land use types
Reduce(intersect, list(plant.occup.top[,1], plant.occup.top[,2], 
                       plant.occup.top[,3], plant.occup.top[,4], 
                       plant.occup.top[,5], plant.occup.top[,6]))

# shared top species between land use types
ShareProp <- function(occup_top_data) {
  share_prop <- as.data.frame(matrix(numeric(0),ncol=3, nrow = 36))
  colnames(share_prop) <- c("land_use_1", "land_use_2", "prop")
  k <- 0
  for (i in c(1:6)) {
    for (j in c(1:6)) {
      k <- k+1
      share_prop$land_use_1[k] <- kLvlLanduse[i]
      share_prop$land_use_2[k] <- kLvlLanduse[j]
      share_prop$prop[k] <- 
        (length(intersect(occup_top_data[,i],occup_top_data[,j])))/
        (length(union(occup_top_data[,i],occup_top_data[,j])))
    }
  }
  ggplot(share_prop, aes(land_use_1, land_use_2, fill = prop)) + 
    geom_tile() + geom_text(aes(label = round(prop*100))) +
    scale_fill_gradient2(high = "red", low = "blue", 
                         midpoint = 0.4, limits = c(0,0.8)) +
    theme(axis.text.x = element_text(angle = 90))
}
ShareProp(plant.occup.top)

# unique ubiquitous species in certain land use types
plant.occup.top %>% pivot_longer(
  cols = all_of(kLvlLanduse),  
  names_to = "land_use", values_to = "Species") %>%
  .[which(
    .$Species %in% 
      # unique ubiquitous species - present once only
      names(table(as.character(as.matrix(plant.occup.top)))[
        table(as.character(as.matrix(plant.occup.top))) == 1])),] %>% 
  arrange(land_use)

# Quadrat level ----
## richness ~ land use for all plants ----
ggplot(qua.div) + 
  geom_boxplot(aes(land_use, richness)) + 
  labs(x = "Land use type", y = "Quadrat richness") + 
  geom_text(data = data.frame(
    land_use = kLvlLanduse, 
    Label = c(rep("a", 4), "ab", "b")
  ), aes(x = land_use, y = Inf, label = Label), vjust = 2) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_bw()
kruskal.test(qua.div$richness ~ qua.div$land_use)
dunn.test(x = qua.div$richness, g = qua.div$land_use)

## Indexes ~ land use for trees and shrubs ----
# Kruskal-Wallis test & box plot for trees 
# tree diversity longer and shrub diversity longer data set
ConsLong <- function(x) {
  subset(x, select = c("abundance", "evenness", "land_use")) %>% 
    pivot_longer(cols = c("abundance", "evenness"), 
                 names_to = "Index", values_to = "Index_value") %>% 
    mutate(Index = factor(Index, levels = c("abundance", "evenness")), 
           land_use = factor(land_use, levels = kLvlLanduse), 
           Attr = c("Land use type")) %>% 
    na.omit()
}
qua.tree.div_long <- qua.tree.div %>% 
  left_join(qua.data, by = "qua_id") %>% 
  ConsLong()
qua.shrub.div_long <- qua.shrub.div %>% 
  left_join(qua.data, by = "qua_id") %>% 
  ConsLong()

# get p-values for box plots
GetBoxPvalue <- function(x) {
  y <- data.frame(Index = c("abundance", "evenness"), 
                  Pvalue = NA, Label = NA) %>% 
    mutate(Index = factor(Index, levels = c("abundance", "evenness")))
  j <- 0
  for (i in c("abundance", "evenness")) {
    j <- j+1
    y$Pvalue[j] <- round(kruskal.test(
      x[[i]] ~ x$land_use)$p.value,digits = 3)
  }
  y$Label <- case_when(
    y$Pvalue >= 0.05 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), sep = ""), 
    y$Pvalue < 0.05 & y$Pvalue >= 0.01 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "*", sep = ""), 
    y$Pvalue < 0.01 & y$Pvalue >= 0.001 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "**", sep = ""),
    y$Pvalue < 0.001 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "***", sep = "")
  )
  data.frame(y)
}
qua.tree.box.pvalue <- qua.tree.div %>% 
  left_join(qua.data, by = "qua_id") %>% 
  GetBoxPvalue()
qua.shrub.box.pvalue <- qua.shrub.div %>% 
  left_join(qua.data, by = "qua_id") %>% 
  GetBoxPvalue()

# get box plots
GetBoxPlot <- function(x, y, z) {
  ggplot(x, aes(land_use, Index_value)) + 
    geom_boxplot() + 
    facet_grid(Index ~ Attr, scales = "free", 
               space = "free_x", switch = "both") + 
    scale_y_continuous(expand = expansion(mult = c(0.05,0.3))) +
    geom_text(data = y, aes(x =Inf, y = Inf, label = Label), 
              size=3.5, hjust = 1.05, vjust = 1.5) +
    theme(axis.text.x = element_text(angle = 90, size = 12)) + 
    labs(title = z, x = NULL, y = NULL) + 
    theme_bw()
}
ggarrange(GetBoxPlot(qua.tree.div_long, qua.tree.box.pvalue, "(a)"), 
          GetBoxPlot(qua.shrub.div_long, qua.shrub.box.pvalue, "(b)"))

# pairwise dunn test of indexes ~ land use type
GetDunn <- function(x, taxa, index) {
  dunn_result <- dunn.test(x[[index]], x$land_use, 
                           table = FALSE, kw = FALSE)
  x <- data.frame(
    "taxa" = taxa, 
    "index" = index, 
    "comparison" = dunn_result$comparisons, 
    "p" = dunn_result$P.adjusted
  ) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
}
dunn.df.1 <- rbind(GetDunn(qua.tree.comm.div, "tree", "abundance"), 
                   GetDunn(qua.tree.comm.div, "tree", "evenness"),
                   GetDunn(qua.shrub.comm.div, "shrub", "abundance"), 
                   GetDunn(qua.shrub.comm.div, "shrub", "evenness")
) 
dunn.df.2 <- dunn.df.1[, c("taxa", "index", "comparison_2", "comparison_1", "p")]
names(dunn.df.2) <- c("taxa", "index", "comparison_1", "comparison_2", "p")
dunn.df <- rbind(dunn.df.1, dunn.df.2)%>% 
  mutate(index = factor(index, levels = KLvlIndex), 
         comparison_1 = factor(comparison_1, kLvlLanduse), 
         comparison_2 = factor(comparison_2, kLvlLanduse))

# plot the pairwise test results
PlotDunn <- function(x, title) {
  ggplot(x, aes(comparison_1, comparison_2, fill = p)) +
    geom_tile() + 
    geom_text(aes(label = round(p*100)), size = 2.5) +
    scale_fill_gradient2(high = "blue", low = "red", 
                         midpoint = 0.05, limits = c(0, 0.05)) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    xlab(NULL) + ylab(NULL) + guides(fill = "none") + 
    facet_wrap(~ index, scales = "free", nrow = 1) +
    labs(title = title)
}
ggarrange(PlotDunn(subset(dunn.df, taxa == "tree"), "Tree"), 
          PlotDunn(subset(dunn.df, taxa == "shrub"), "Shrub"), 
          nrow = 2)

## Species composition ----
# use non-metric multidimensional scaling
set.seed(1234)
# nMDS calculation for tree
tree.mds.selected <- subset(qua.tree.comm.div, abundance > 1)
tree.mds.meta <- tree.mds.selected %>% 
  select(2:(KNumTreeSpecies+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree.mds.meta$stress
stressplot(tree.mds.meta)
tree.mds.selected <- cbind(tree.mds.selected, tree.mds.meta$points)

# nMDS calculation for shrub
shrub.mds.selected <- qua.shrub.comm.div %>% filter(abundance > 5)
shrub.mds.meta <- shrub.mds.selected %>% 
  select(2:(KNumShrubSpecies+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub.mds.meta$stress
stressplot(shrub.mds.meta)
shrub.mds.selected <- cbind(shrub.mds.selected, shrub.mds.meta$points)

# ANOSIM of trees and shrubs as labels for the nMDS plots
qua.tree.anosim <- 
  anosim(tree.mds.selected[KTreeSpecies], 
         tree.mds.selected$land_use)
qua.shrub.anosim <- 
  anosim(shrub.mds.selected[KShrubSpecies], 
         shrub.mds.selected$land_use)

# get hull for nMDS plots of trees and shrubs
FindHull <- function(x) {x[chull(x$MDS1, x$MDS2), ]}
qua.tree.hull <- ddply(tree.mds.selected, "land_use", FindHull)
qua.shrub.hull <- ddply(shrub.mds.selected, "land_use", FindHull)

# nMDS plots for trees and shrubs by land use types
PlotNmds <- function(mds_selected, hull, plot_title, mds_meta, anosim) {
  ggplot(mds_selected, aes(MDS1, MDS2, color = land_use)) + 
    geom_point(size=3) +
    geom_polygon(data = hull, alpha = 0, aes(fill=land_use), size=1) +
    labs(title = plot_title, color = "Land use type", fill = "Land use type", 
         subtitle = paste("stress=", sprintf("%.3f",mds_meta$stress),
                          ", R=", sprintf("%.3f",anosim$statistic),
                          ", p=", sprintf("%.3f",anosim$signif),sep = "")) +
    theme(axis.text = element_text(size = 12), 
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 12)) +
    theme_bw()
}
ggarrange(
  PlotNmds(tree.mds.selected, qua.tree.hull, "(a)", 
                tree.mds.meta, qua.tree.anosim), 
  PlotNmds(shrub.mds.selected, qua.shrub.hull, "(b)", 
                shrub.mds.meta, qua.shrub.anosim),
  common.legend = T, legend = "right"
)

# pairwise result of ANOSIM of trees and shrubs by land use
anosim.pair <- combn(kLvlLanduse, 2)
PairAnosim <- function(x, y) {
  result <- NULL
  for (i in 1:ncol(anosim.pair)) {
    set.seed(1234)
    mds_selected_sub <- 
      subset(x, land_use %in% c(anosim.pair[1,i], anosim.pair[2,i]))
    result <- c(result, anosim(mds_selected_sub[2:(y+1)], 
                               mds_selected_sub$land_use)$signif)
  }
  result
}
data.frame("Comp_1" = anosim.pair[1,], "Comp_2" = anosim.pair[2,], 
           "p" = PairAnosim(tree.mds.selected, KNumTreeSpecies)
) %>% subset(p < 0.05)
data.frame("Comp_1" = anosim.pair[1,], "Comp_2" = anosim.pair[2,], 
           "p" = PairAnosim(shrub.mds.selected, KNumShrubSpecies)
) %>% subset(p < 0.05)
rm(qua.tree.anosim, qua.tree.hull, tree.mds.meta, tree.mds.selected, 
   qua.shrub.anosim, qua.shrub.hull, shrub.mds.meta, shrub.mds.selected, 
   anosim.pair, PlotNmds, FindHull, PairAnosim)

# Cor among the indexes ----
chart.Correlation(subset(qua.tree.div, select = KLvlIndex))
chart.Correlation(subset(qua.shrub.div, select = KLvlIndex))

# Data for discussion ----
# means of quadrat abundance and richness for trees
qua.tree.comm.div %>% group_by(land_use) %>% 
  dplyr::summarise(abundance = mean(abundance), richness = mean(richness))

# means of quadrat abundance and richness for trees
qua.div %>% group_by(land_use) %>% 
  dplyr::summarise(richness = mean(richness))
