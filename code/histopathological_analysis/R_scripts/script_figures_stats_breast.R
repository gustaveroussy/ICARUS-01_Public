library(DirichletReg)
library(ggpubr)
library(lemon)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(lme4)
library(lmerTest)
library(PairedData)
library(rstatix)

# Read data, change path accordingly
df_cluster_perc <- read.csv("./BREAST/clustering/New_cluster_analysis/cluster_percentages_cell_feats_only.csv")
n_clusters <- ncol(df_cluster_perc) - 1

# Prepare compositional data for Dirichlet regression
comp_data <- DR_data(df_cluster_perc[,1:n_clusters] / 100)

# Convert Response.to.treatment to a factor
df_cluster_perc$Response.to.treatment <- as.factor(df_cluster_perc$Response.to.treatment)

# Perform Dirichlet regression
dir_model <- DirichReg(comp_data ~ Response.to.treatment , df_cluster_perc, model="common")
summary(dir_model)

# Extract coefficients and confidence intervals
coef_matrix <- coef(dir_model)
conf_int <- confint(dir_model)$ci
se_matrix <- confint(dir_model)$se
cluster_names <- names(coef_matrix)
cluster_indices <- seq_along(cluster_names)

# Initialize data frame for plotting
plot_data <- data.frame(
  Component = character(),
  Predictor = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric()
)

# Calculate odds ratios and confidence intervals
index <- 2
for (idx in cluster_indices) {
  coef_values <- coef_matrix[[cluster_names[idx]]][2]
  odds_ratios <- exp(coef_values)
  
  conf_values <- conf_int[[1]][[idx]]
  lower_ci <- exp(conf_values[, 1])
  upper_ci <- exp(conf_values[, 2])
  
  for (predictor in names(coef_values)) {
    z_value <- coef_values[predictor] / se_matrix[index]
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    index <- index + 2
    plot_data <- rbind(
      plot_data,
      data.frame(
        Component = paste0("Cluster ", as.character(idx-1)),
        Predictor = predictor,
        OR = odds_ratios[predictor],
        CI_lower = lower_ci[predictor],
        CI_upper = upper_ci[predictor],
        p_value = p_value
      )
    )
  }
}

# Adjust component factor levels for plotting
plot_data$Component <- factor(plot_data$Component, levels = rev(unique(plot_data$Component)))

# Plot odds ratios with confidence intervals
ggplot(plot_data, aes(y = Component, x = OR, color = Predictor)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4, color = "blue", size=1) +
  geom_point(size = 3, shape = 21, fill = "black", colour="black") +
  geom_text(aes(label = sprintf("%.3f", OR)), vjust = -0.5, hjust = -0.1, size = 5, color = "black") +
  geom_text(aes(label = sprintf("p = %.3f", p_value)), vjust = 1.5, hjust = -0.1, size = 5, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  #labs(x = "Odds Ratio", y = "Component", title = "Odds Ratios with 95% Confidence Intervals") +
  theme_minimal(base_size = 8) +
  xlab("Odds ratios") +
  theme(legend.position = "none", axis.text.y = element_text(size = 12), axis.title.y = element_blank(),  axis.title.x = element_text(size = 12)) 

# Helper function to extract cluster ID labels
extract_id_label <- function(x) {
  id <- sub(".*?(\\d+).*", "\\1", x)
  paste("Cluster", id)
}

# Read and plot cluster percentages
df_cluster_perc_melted <- read.csv("./BREAST/clustering/New_cluster_analysis/cluster_percentages_cell_feats_only_melted.csv")
box_plot <-ggboxplot(df_cluster_perc_melted, x = "Cluster", y = "Percentage",
                     color = "Response.to.treatment", palette = "npg",
                     add = "point", xlab="", ylab="Cluster percentage") + labs(color='Objective response to treatment :') +
  scale_x_discrete(labels = extract_id_label) 
box_plot

# Perform Wilcoxon test and adjust p-values
stat.test <- df_cluster_perc_melted %>%
  group_by(Cluster) %>%
  wilcox_test(Percentage ~ Response.to.treatment)
stat.test$p.adj = format(round(p.adjust(stat.test$p, method = "BH", n = length(stat.test$p)), 4), nsmall=4)
stat.test <- stat.test %>% 
  add_xy_position(x = "Cluster", dodge = 0.8)

# Read and plot DAB variance
df_intracell_distr <- read.csv("./BREAST/nuclei_IHC/Cell_staining_distribution/cytoplasm_dab_avg_var.csv")
box_plot <-ggboxplot(df_intracell_distr, x="Response.to.treatment", y = "cytoplasm_dab_avg_var",
                     color = "Response.to.treatment", palette = "npg",
                     add = "point", xlab="", ylab="Average DAB variance OD") + labs(color='Objective response to treatment :') + stat_compare_means(label = "p.signif", label.x = 1.5, label.y = 0.075)
box_plot

# Read and reshape data for feature visualization
df_feats <- read.csv("./BREAST/clustering/New_cluster_analysis/features_cells_feats_only.csv")
names(df_feats) <- gsub("_N_cells.*", "", names(df_feats))

# Uncomment to add star 
#df_feats <- df_feats %>%
#  mutate(cluster = paste0(cluster, "*"))

# Long format data for feature visualization
df_cell_feats <- df_feats %>%
  pivot_longer(cols = -cluster, names_to = "feature", values_to = "value")

# Uncomment when using cell features with DAB
#df_cell_feats <- df_feats %>% select(-DAB_avg) %>%
#  pivot_longer(cols = -cluster, names_to = "feature", values_to = "value")

# Plot boxplots for features across clusters
ggplot(df_cell_feats, aes(x = factor(cluster), y = value, fill = factor(cluster))) +
  geom_boxplot() +
  labs(x = "Cluster id", y = "Number of cells") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~ feature, scales = "free") +  # scales="free_y" allows each feature to have its own y-scale
  theme_minimal() +
  theme(legend.position = "none")

# Plot DAB average across clusters
ggplot(df_feats, aes(x = factor(cluster), y = DAB_avg, fill = factor(cluster))) +
  geom_boxplot() +
  labs(x = "Cluster id", y = "Average DAB optical density") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "none")

# Comparison between number of patches belonging to each cluster with the number of ERG+ cells

df_ERG = read.csv("./BREAST_ERG/dataset_ERG_cluster_occ_with_dab.csv")
df_ERG <- na.omit(df_ERG)

# Plot the correlation between the number of patches belonging to cluster 0 and the number of ERG+ cells using Spearman correlation

ggscatter(df_ERG, x = "cluster_0", y = "Nb.ERG..cells", 
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          xlab = "Number of patches belonging to cluster 0", 
          ylab = "Number of ERG cells"
) +
  stat_cor(method = "spearman",                     # Use Spearman correlation
           label.x = 3,                             # Position on x-axis
           label.y = max(df_ERG$`Nb.ERG..cells`, na.rm = TRUE)-200, # Position on y-axis
           p.accuracy = 0.001,                      # P-value accuracy
           r.accuracy = 0.01,                       # Correlation accuracy
           label.sep = ", ",                        # Separator between rho and p-value
           aes(label = paste("italic(rho)~`=`~", ..r.., 
                             "*`,`~italic(p)~`=`~", 
                             ..p.., 
                             "*", 
                             ifelse(..p.. < 0.0001, "' ****'", 
                                    ifelse(..p.. < 0.001, "' ***'", 
                                           ifelse(..p.. < 0.01, "' **'", 
                                                  ifelse(..p.. < 0.05, "' *'", "' ns'"))))))) + theme(legend.position = "none",
                                                                                                      axis.text.x = element_text(size = 12), 
                                                                                                      axis.text.y = element_text(size = 12), 
                                                                                                      axis.title.y = element_text(size = 12),  
                                                                                                      axis.title.x = element_text(size = 12))

# Compute Spearman correlation between the number of patches belonging to each cluster and the number of ERG+ cells

corr_cluster_ERG <- data.frame(
  Variable = character(),
  Correlation = numeric(),
  P.Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each variable
for (var in names(df_ERG)[1:8]) {  # Exclude the target variable
  test <- cor.test(df_ERG$`Nb.ERG..cells`, df_ERG[[var]], use = "complete.obs", method="spearman", exact = TRUE)
  corr_cluster_ERG <- rbind(corr_cluster_ERG, data.frame(
    Variable = var,
    Correlation = test$estimate,
    P.Value = test$p.value
  ))
}

corr_cluster_ERG$P.Adj = format(round(p.adjust(corr_cluster_ERG$P.Value, method = "BH", n = length(corr_cluster_ERG$P.Value)), 4), nsmall=4)

# Display results
corr_cluster_ERG                                                                                         

# Comparison between the percentage of cluster 0 in the HER3-Dxd positive and negative groups at C1D3

df_HER3_Dxd = read.csv("./BREAST_on_treatment/dataset_HER3_positivity_cluster_perc_with_dab.csv")
df_HER3_Dxd = na.omit(df_HER3_Dxd)
df_HER3_Dxd = df_HER3_Dxd %>% filter(Timepoint %in% c("C1D3"))

# Identify two groups based on the percentage of intact-Dxd in PANK+ cells
df_HER3_Dxd$group <- ifelse(df_HER3_Dxd$X.intact.Dxd.in.PANK. >= 5, "Intact-Dxd >= 5%",  "Intact-Dxd < 5%")
df_HER3_Dxd$group <- factor(df_HER3_Dxd$group, levels = c("Intact-Dxd < 5%", "Intact-Dxd >= 5%"))

# Boxplot comparing the percentage of cluster 0 in the two groups

ggboxplot(df_HER3_Dxd, x = "group", y = "cluster_0",
          color = "group", palette = "npg",
          add = "point", xlab="", ylab="Percentage of cluster 0") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Intact-Dxd >= 5%", "Intact-Dxd < 5%")),
                     label = "p.", # Displays significance stars (e.g., **)
                     label.y = max(df_HER3_Dxd$cluster_0, na.rm = TRUE) + 0.2) + # Adjust label position
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12),  
        axis.title.x = element_text(size = 12))

# Comparison between the percentage of cluster 0 at baseline and EOT

df_baseline = read.csv("./BREAST/clustering/New_cluster_analysis/cluster_percentages_with_DAB_with_info_patients.csv")
df_eot = read.csv("./BREAST_EOT/cluster_percentages_with_dab_with_info_patients_EOT.csv")

# Filter data to keep only common IDs
common_ids <- intersect(df_baseline$Additional.Specimen.ID, df_eot$Additional.Specimen.ID)

df_baseline_filtered <- df_baseline %>%
  filter(Additional.Specimen.ID %in% common_ids) %>%
  arrange(match(Additional.Specimen.ID, common_ids)) %>%
  mutate(timepoint = "Baseline")

df_eot_filtered <- df_eot %>%
  filter(Additional.Specimen.ID %in% common_ids) %>%
  arrange(match(Additional.Specimen.ID, common_ids)) %>%
  mutate(timepoint = "EOT")

df_comparison_timepoints <- bind_rows(df_baseline_filtered[, c("cluster_0", "timepoint")], df_eot_filtered[, c("cluster_0", "timepoint")])

# Boxplot comparing the percentage of cluster 0 at baseline and EOT

ggboxplot(df_comparison_timepoints, x = "timepoint", y = "cluster_0", 
          color = "timepoint", palette = c("#059748", "#faa43a"),
          order = c("Baseline", "EOT"),
          ylab = "Percentage of cluster 0", xlab = "Timepoint") +
  stat_compare_means(method = "wilcox.test",
                     paired=TRUE,
                     comparisons = list(c("Baseline", "EOT")),
                     label = "p.", # Displays significance stars (e.g., **)
                     label.y = max(df_comparison_timepoints$cluster_0, na.rm = TRUE) + 0.2)+
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12),  
        axis.title.x = element_text(size = 12))