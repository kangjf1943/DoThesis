# Settings ----
lvl_bd_index <- c("abundance", "richness", "shannon", "simpson", "evenness")
tree_spe <- unique(tree_data$species_lt)
shrub_spe <- unique(shrub_data$species_lt)

# Functions ----
# func to plot species accumulation curve and extrapolation 
# func parameters: x, raw data; y, size of extrapolation; z, extent of x axis; method, plot the curve separately by land use (method = "land_use") or at city level (method = "city")
fun_accum <- function(x, y, z, method) {
  # func to derive list of incidence data of certain land use
  # func para: xtemp, raw data; ytemp, name of land use
  fun_temp <- function(xtemp, ytemp) {
    xtemp %>% filter(landuse %in% ytemp) %>% 
      group_by(species_lt) %>% 
      summarise(n = length(unique(qua_id))) %>% 
      # arrange(desc(n)) %>% 
      .$n %>% 
      c(length(unique(
        xtemp[which(xtemp[["landuse"]] %in% ytemp), ][["qua_id"]])), .)
  }
  if (method == "city") {
    incidence <- list(fun_temp(x, lvl_landuse))
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
fun_top <- function(x, y, z, k, n) {
  top_plant <- x %>% group_by(get(y)) %>% 
    dplyr::summarise(number = sum(get(z)), prop = number/k) %>% 
    arrange(desc(number)) %>% 
    head(n)
  names(top_plant)[1] = c(y)
  print(top_plant)
}

# func to test if top species are contained in top families
# func para: x, raw data of species; y, raw data of families 
fun_contain <- function(x, y) {
  x %>% left_join(all_plant_info, by = "species_lt") %>% 
    select(species_lt, family) %>% 
    left_join(y, by = "family") %>% 
    select(species_lt, family)
}

fun_rank_data <- function(x, name_spe) {
  x %>% select({{name_spe}}) %>%
    subset(select = colSums(.) != 0) %>% 
    as.data.frame() %>% 
    rankabundance() %>% 
    as.data.frame() %>%
    mutate(species_lt = rownames(.))
}

# func to plot rank-abundance curve
# func para: x, raw data; title, plot title; method, plot at city level (method = "city") or at land use level (method = "land_use)
fun_rank_plot <- function(x, title, method) {
  if (method == "city") {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + labs(title = title)
  } else {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + 
      facet_wrap(~landuse, nrow = 1) + labs(title = title)
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

# Analysis ----
## City level ----
# Number of species
cat("total species:", length(unique(all_plant_data$species_lt)), "\n", 
    "total genera:", length(unique(all_plant_data$genus)), "\n", 
    "total families:", length(unique(all_plant_data$family)), "\n")

# species accumulation curve 
# extrapolation up to double the reference sample size
fun_accum(all_plant_data, 348, 350, method = "city") + 
  theme(legend.position = "none")

## Top taxa ----
# top species families of all plants by species number
all_plant_info %>% group_by(family) %>% 
  dplyr::summarise(Num_spe = n(), Prop = n()/nrow(all_plant_info)) %>% 
  arrange(desc(Prop))

# abundance of trees and shrubs
cat("number of trees:", nrow(tree_data), "in", 
    nrow(tree_qua_div), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$area), "m2 in", 
    nrow(shrub_qua_div), "plot")

# top species families of trees and shrubs by abundance
tree_top_species <- fun_top(
  tree_data, "species_lt", "stem", sum(tree_data$stem), 10)
tree_top_family <- fun_top(
  tree_data, "family", "stem", sum(tree_data$stem), 10)
fun_contain(tree_top_species, tree_top_family)

shrub_top_species <- fun_top(
  shrub_data, "species_lt", "area", sum(shrub_data$area), 10)
shrub_top_family <- fun_top(
  shrub_data, "family", "area", sum(shrub_data$area), 10)
intersect(tree_top_family$family, shrub_top_family$family)
fun_contain(shrub_top_species, shrub_top_family)

rm(tree_top_species, tree_top_family, 
   shrub_top_species, shrub_top_family, 
   fun_top, fun_contain)

# attributes of trees and shrubs
# the number of exotic vs. native by number of species
table(all_plant_info$nt_ex)/nrow(all_plant_info)

# the attributes of trees and shrubs
for (i in c("pla_spo", "pub_pri", "nt_ex")) {
  print(tapply(tree_data$stem, tree_data[,i], sum)/sum(tree_data$stem), 
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
city_tree_rank <- fun_rank_data(tree_qua_comm, tree_spe) %>% 
  left_join(select(all_plant_info, c("species_lt", "nt_ex")), by = "species_lt")
city_shrub_rank <- fun_rank_data(shrub_qua_comm, shrub_spe) %>% 
  left_join(select(all_plant_info, c("species_lt", "nt_ex")), by = "species_lt")

ggarrange(plotlist = list(
  fun_rank_plot(city_tree_rank, title = "(a)", method = "city"), 
  fun_rank_plot(city_shrub_rank, title = "(b)", method = "city")
), nrow = 1, common.legend = TRUE)

# the top 3 species regarding abundance
subset(city_tree_rank[,c("rank", "species_lt", "nt_ex")], rank <= 3)
subset(city_shrub_rank[,c("rank", "species_lt", "nt_ex")], rank <= 3)

# EQ evenness index and plot
city_tree_rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>% 
  arrange(EQ)
city_shrub_rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>%
  arrange(EQ)

# Land use level ----
## Species accumulation curve ----
fun_accum(all_plant_data, 348, 50, method = "landuse") + 
  labs(color = "Land use") + 
  scale_color_manual(values = c("#FF0000", "#FF7800", "#DF73FF", 
                                "#BFBF30", "#6BE400", "#00733E"))

## Distribution of species abundance ----
tree_lu_rank <- 
  left_join(tree_qua_comm, qua_data[c("qua_id", "landuse")], by = "qua_id") %>% 
  ddply("landuse", name_spe = tree_spe, fun_rank_data) %>% 
  left_join(select(all_plant_info, c("species_lt", "nt_ex")), by = "species_lt")
shrub_lu_rank <- 
  left_join(shrub_qua_comm, qua_data[c("qua_id", "landuse")], by = "qua_id") %>%
  ddply("landuse", name_spe = shrub_spe, fun_rank_data) %>% 
  left_join(select(all_plant_info, c("species_lt", "nt_ex")), by = "species_lt")
ggarrange(fun_rank_plot(tree_lu_rank, "(a)", method = "land_use"),
          fun_rank_plot(shrub_lu_rank, "(b)", method = "land_use"), 
          nrow = 2, common.legend = TRUE, legend = "bottom")
# the top 3 species regarding abundance
subset(tree_lu_rank[,c("landuse", "rank", "species_lt", "nt_ex")], rank <= 3)
subset(shrub_lu_rank[,c("landuse", "rank", "species_lt", "nt_ex")], rank <= 3)
# calculate the EQ evenness index and plot
community_structure(
  tree_lu_rank, time.var = "landuse", 
  abundance.var = "abundance", metric = "EQ") %>% 
  arrange(desc(EQ))
community_structure(
  shrub_lu_rank, time.var = "landuse", 
  abundance.var = "abundance", metric = "EQ") %>%
  arrange(desc(EQ))

## Species composition ----
# Bray-Curtis dissimilarity of pairs of land use for trees and shrubs
tree_lu_comm %>% 
  select(-landuse) %>% 
  vegdist() %>% 
  as.matrix() %>% 
  round(digits = 2) %>% 
  as_tibble() %>% 
  rename_all(~lvl_landuse) %>% 
  mutate(pairs = lvl_landuse) %>% 
  relocate(pairs)

shrub_lu_comm %>% 
  select(-landuse) %>% 
  vegdist() %>% 
  as.matrix() %>% 
  round(digits = 2) %>% 
  as_tibble() %>% 
  rename_all(~lvl_landuse) %>% 
  mutate(pairs = lvl_landuse) %>% 
  relocate(pairs)

# plant occupancy of species for different land use types
fun_occup_rate <- function(x,y) {
  apply(x[,2:(y+1)], 2, function(k)sum(ifelse(k>0,1,0))/length(k))
}
fun_occup_df <- function(x){
  data.frame("Com" = names(head(sort(x[,"Com"],decreasing = TRUE),10)), 
             "ComNbr" = names(head(sort(x[,"ComNbr"],decreasing = TRUE),10)),
             "ResLow" = names(head(sort(x[,"ResLow"],decreasing = TRUE),10)),
             "ResHigh" = names(head(sort(x[,"ResHigh"],decreasing = TRUE),10)), 
             "ResOther" = names(head(sort(x[,"ResOther"],decreasing = TRUE),10)), 
             "Ind" = names(head(sort(x[,"Ind"],decreasing = TRUE),10)))
}

plant_occup <- ddply(
  all_qua_div, .(landuse), y = number_plant_species, fun_occup_rate) %>% 
  .[,-1] %>% t() 
colnames(plant_occup) = lvl_landuse
(plant_occup_top <- fun_occup_df(plant_occup))

# shared species over land use types
Reduce(intersect, list(plant_occup_top[,1], plant_occup_top[,2], 
                       plant_occup_top[,3], plant_occup_top[,4], 
                       plant_occup_top[,5], plant_occup_top[,6]))

# shared top species between land use types
fun_share_prop <- function(occup_top_data) {
  share_prop <- as.data.frame(matrix(numeric(0),ncol=3, nrow = 36))
  colnames(share_prop) <- c("land_use_1", "land_use_2", "prop")
  k <- 0
  for (i in c(1:6)) {
    for (j in c(1:6)) {
      k <- k+1
      share_prop$land_use_1[k] <- lvl_landuse[i]
      share_prop$land_use_2[k] <- lvl_landuse[j]
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
fun_share_prop(plant_occup_top)

# unique ubiquitous species in certain land use types
plant_occup_top %>% pivot_longer(
  cols = all_of(lvl_landuse),  
  names_to = "landuse", values_to = "Species") %>%
  .[which(
    .$Species %in% 
      # unique ubiquitous species - present once only
      names(table(as.character(as.matrix(plant_occup_top)))[
        table(as.character(as.matrix(plant_occup_top))) == 1])),] %>% 
  arrange(landuse)

# Quadrat level ----
## richness ~ land use for all plants ----
ggplot(all_qua_div) + 
  geom_boxplot(aes(landuse, richness)) + 
  labs(x = "Land use type", y = "Quadrat richness") + 
  geom_text(data = data.frame(
    landuse = lvl_landuse, 
    Label = c(rep("a", 4), "ab", "b")
  ), aes(x = landuse, y = Inf, label = Label), vjust = 2) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_bw()
kruskal.test(all_qua_div$richness ~ all_qua_div$landuse)
dunn.test(x = all_qua_div$richness, g = all_qua_div$landuse)

## Indexes ~ land use for trees and shrubs ----
# Kruskal-Wallis test & box plot for trees 
# tree diversity longer and shrub diversity longer data set
fun_cons_long <- function(x) {
  subset(x, select = c("abundance", "evenness", "landuse")) %>% 
    pivot_longer(cols = c("abundance", "evenness"), 
                 names_to = "Index", values_to = "Index_value") %>% 
    mutate(Index = factor(Index, levels = c("abundance", "evenness")), 
           landuse = factor(landuse, levels = lvl_landuse), 
           Attr = c("Land use type")) %>% 
    na.omit()
}
tree_qua_div_long <- tree_qua_div %>% 
  left_join(qua_data, by = "qua_id") %>% 
  fun_cons_long()
shrub_qua_div_long <- shrub_qua_div %>% 
  left_join(qua_data, by = "qua_id") %>% 
  fun_cons_long()

# get p-values for box plots
fun_get_pvalue <- function(x) {
  y <- data.frame(Index = c("abundance", "evenness"), 
                  Pvalue = NA, Label = NA) %>% 
    mutate(Index = factor(Index, levels = c("abundance", "evenness")))
  j <- 0
  for (i in c("abundance", "evenness")) {
    j <- j+1
    y$Pvalue[j] <- round(kruskal.test(
      x[[i]] ~ x$landuse)$p.value,digits = 3)
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
tree_box_pvalue <- tree_qua_div %>% 
  left_join(qua_data, by = "qua_id") %>% 
  fun_get_pvalue()
shrub_box_pvalue <- shrub_qua_div %>% 
  left_join(qua_data, by = "qua_id") %>% 
  fun_get_pvalue()

# get box plots
fun_box_plot <- function(x, y, z) {
  ggplot(x, aes(landuse, Index_value)) + 
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
ggarrange(fun_box_plot(tree_qua_div_long, tree_box_pvalue, "(a)"), 
          fun_box_plot(shrub_qua_div_long, shrub_box_pvalue, "(b)"))

# pairwise dunn test of indexes ~ land use type
fun_dunn <- function(x, taxa, index) {
  dunn_result <- dunn.test(x[[index]], x$landuse, 
                           table = FALSE, kw = FALSE)
  x <- data.frame(
    "taxa" = taxa, 
    "index" = index, 
    "comparison" = dunn_result$comparisons, 
    "p" = dunn_result$P.adjusted
  ) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
}
dunn_df_1 <- rbind(fun_dunn(tree_qua_all, "tree", "abundance"), 
                   fun_dunn(tree_qua_all, "tree", "evenness"),
                   fun_dunn(shrub_qua_all, "shrub", "abundance"), 
                   fun_dunn(shrub_qua_all, "shrub", "evenness")
) 
dunn_df_2 <- dunn_df_1[, c("taxa", "index", "comparison_2", "comparison_1", "p")]
names(dunn_df_2) <- c("taxa", "index", "comparison_1", "comparison_2", "p")
dunn_df <- rbind(dunn_df_1, dunn_df_2)%>% 
  mutate(index = factor(index, levels = lvl_bd_index), 
         comparison_1 = factor(comparison_1, lvl_landuse), 
         comparison_2 = factor(comparison_2, lvl_landuse))
# plot the pairwise test results
fun_dunn_plot <- function(x, title) {
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
ggarrange(fun_dunn_plot(subset(dunn_df, taxa == "tree"), "Tree"), 
          fun_dunn_plot(subset(dunn_df, taxa == "shrub"), "Shrub"), 
          nrow = 2)

## Species composition ----
# use non-metric multidimensional scaling
set.seed(1234)
# nMDS calculation for tree
tree_mds_selected <- subset(tree_qua_all, abundance > 1)
tree_mds_meta <- tree_mds_selected %>% 
  select(2:(number_tree_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_meta$stress
stressplot(tree_mds_meta)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_meta$points)

# nMDS calculation for shrub
shrub_mds_selected <- shrub_qua_all %>% filter(abundance > 5)
shrub_mds_meta <- shrub_mds_selected %>% 
  select(2:(number_shrub_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_meta$stress
stressplot(shrub_mds_meta)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_meta$points)

# ANOSIM of trees and shrubs as labels for the nMDS plots
tree_anosim <- 
  anosim(tree_mds_selected[tree_spe], 
         tree_mds_selected$landuse)
shrub_anosim <- 
  anosim(shrub_mds_selected[shrub_spe], 
         shrub_mds_selected$landuse)

# get hull for nMDS plots of trees and shrubs
fun_find_hull <- function(x) {x[chull(x$MDS1, x$MDS2), ]}
tree_hulls <- ddply(tree_mds_selected, "landuse", fun_find_hull)
shrub_hulls <- ddply(shrub_mds_selected, "landuse", fun_find_hull)

# nMDS plots for trees and shrubs by land use types
fun_nmds_plot <- function(mds_selected, hull, plot_title, mds_meta, anosim) {
  ggplot(mds_selected, aes(MDS1, MDS2, color = landuse)) + 
    geom_point(size=3) +
    geom_polygon(data = hull, alpha = 0, aes(fill=landuse), size=1) +
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
  fun_nmds_plot(tree_mds_selected, tree_hulls, "(a)", 
                tree_mds_meta, tree_anosim), 
  fun_nmds_plot(shrub_mds_selected, shrub_hulls, "(b)", 
                shrub_mds_meta, shrub_anosim),
  common.legend = T, legend = "right"
)

# pairwise result of ANOSIM of trees and shrubs by land use
anosim_pairs <- combn(lvl_landuse, 2)
fun_anosim_pairs <- function(x, y) {
  result <- NULL
  for (i in 1:ncol(anosim_pairs)) {
    set.seed(1234)
    mds_selected_sub <- 
      subset(x, landuse %in% c(anosim_pairs[1,i], anosim_pairs[2,i]))
    result <- c(result, anosim(mds_selected_sub[2:(y+1)], 
                               mds_selected_sub$landuse)$signif)
  }
  result
}
data.frame("Comp_1" = anosim_pairs[1,], "Comp_2" = anosim_pairs[2,], 
           "p" = fun_anosim_pairs(tree_mds_selected, number_tree_species)
) %>% subset(p < 0.05)
data.frame("Comp_1" = anosim_pairs[1,], "Comp_2" = anosim_pairs[2,], 
           "p" = fun_anosim_pairs(shrub_mds_selected, number_shrub_species)
) %>% subset(p < 0.05)
rm(tree_anosim, tree_hulls, tree_mds_meta, tree_mds_selected, 
   shrub_anosim, shrub_hulls, shrub_mds_meta, shrub_mds_selected, 
   anosim_pairs, fun_nmds_plot, fun_find_hull, fun_anosim_pairs)

# Cor among the indexes ----
chart.Correlation(subset(tree_qua_div, select = lvl_bd_index))
chart.Correlation(subset(shrub_qua_div, select = lvl_bd_index))


# Data for discussion ----
# means of quadrat abundance and richness for trees
tree_qua_all %>% group_by(landuse) %>% 
  dplyr::summarise(abundance = mean(abundance), richness = mean(richness))

# means of quadrat abundance and richness for trees
all_qua_div %>% group_by(landuse) %>% 
  dplyr::summarise(richness = mean(richness))

