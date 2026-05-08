### Libraries
library(readxl)
library(writexl)
library(devtools)
library(leaflet)
library(scatterpie)
library(webr)
library(tibble)
library(vegan)
library(plyr)
library(pairwiseAdonis)
library(ape)  # For phylogenetic tree



# STATISTICAL ANALYSES

## ------------------------ Alpha diversity ------------------------ ##
# ---- Richness 
df <- read_excel("/path/to/folder/file", sheet = "richness")
df_mat <- as.data.frame(df) # Make the first column the row names
rownames(df_mat) <- df_mat$`Etiquetas de fila` 
df_mat <- df_mat[, -1] # Remove that non-numeric column
df_mat[] <- lapply(df_mat, as.numeric) # Ensure all remaining columns are numeric


# --- Observed richness (number of phyla per sample)
richness <- specnumber(df_mat)  # counts non-zero phyla per row (sample)

# --- Shannon diversity index
shannon <- diversity(df_mat, index = "shannon")

# --- Chao1 index (only can be calculated from counts of species, not from relative abundances)
df_chao <- read_excel("/path/to/folder/file", sheet = "chao")
df_mat <- as.data.frame(df_chao) # Make the first column the row names
rownames(df_mat) <- df_mat$`Etiquetas de fila` 
df_mat <- df_mat[, -1] # Remove that non-numeric column
df_mat[] <- lapply(df_mat, as.numeric) # Ensure all remaining columns are numeric

chao1 <- estimateR(df_mat)["S.chao1", ] 

# --- Hill number q=1 (effective number of taxa)
hill_q1 <- exp(shannon)

# --- Simpson and inverse Simpson diversity
simpson <- diversity(df_mat, index = "simpson")
inv_simpson <- 1 / simpson  # Hill q=2
Pielou <- shannon / log(richness)


# Combine into one table
alpha_div <- data.frame(
  Sample = rownames(df_mat),
  Observed_Richness = richness,
  Shannon = shannon,
  Chao1 = chao1,
  Hill_q1 = hill_q1,
  Inverse_Simpson = inv_simpson,
  Pielou = Pielou
)

# head(alpha_div)

merged_df <- merge(alpha_div, metadata[, c("Sample", "Habitat", "Vent_Classification", "Site_28")],
                   by = "Sample", all.x = TRUE)

# Function to summarize all alpha metrics
summarize_alpha <- function(df, group_var) {
  df %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      Mean_ObsRichness = mean(Observed_Richness, na.rm = TRUE),
      Min_ObsRichness  = min(Observed_Richness, na.rm = TRUE),
      Max_ObsRichness  = max(Observed_Richness, na.rm = TRUE),
      Range_ObsRichness = max(Observed_Richness, na.rm = TRUE) - min(Observed_Richness, na.rm = TRUE),
      Mean_Shannon   = mean(Shannon, na.rm = TRUE),
      Min_Shannon    = min(Shannon, na.rm = TRUE),
      Max_Shannon    = max(Shannon, na.rm = TRUE),
      Range_Shannon  = max(Shannon, na.rm = TRUE) - min(Shannon, na.rm = TRUE),
      Mean_Chao1   = mean(Chao1, na.rm = TRUE),
      Min_Chao1    = min(Chao1, na.rm = TRUE),
      Max_Chao1    = max(Chao1, na.rm = TRUE),
      Range_Chao1  = max(Chao1, na.rm = TRUE) - min(Chao1, na.rm = TRUE),
      Mean_Hill_q1   = mean(Hill_q1, na.rm = TRUE),
      Min_Hill_q1    = min(Hill_q1, na.rm = TRUE),
      Max_Hill_q1    = max(Hill_q1, na.rm = TRUE),
      Range_Hill_q1  = max(Hill_q1, na.rm = TRUE) - min(Hill_q1, na.rm = TRUE),
      Mean_Inverse_Simpson   = mean(Inverse_Simpson, na.rm = TRUE),
      Min_Inverse_Simpson    = min(Inverse_Simpson, na.rm = TRUE),
      Max_Inverse_Simpson    = max(Inverse_Simpson, na.rm = TRUE),
      Range_Inverse_Simpson  = max(Inverse_Simpson, na.rm = TRUE) - min(Inverse_Simpson, na.rm = TRUE),
      Mean_Pielou   = mean(Pielou, na.rm = TRUE),
      Min_Pielou    = min(Pielou, na.rm = TRUE),
      Max_Pielou    = max(Pielou, na.rm = TRUE),
      Range_Pielou  = max(Pielou, na.rm = TRUE) - min(Pielou, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

By_Habitat <- summarize_alpha(merged_df, "Habitat")
By_Vent_Classification <- summarize_alpha(merged_df, "Vent_Classification")
By_Site_28 <- summarize_alpha(merged_df, "Site_28")

# Create a list of the summary tables
summary_list <- list(
  By_Habitat = By_Habitat,
  By_Vent_Classification = By_Vent_Classification,
  By_Site_28 = By_Site_28
)

# Write to Excel
# write_xlsx(summary_list, path = "Alpha_Richness_Summary.xlsx")

## Kruskal walis
df_clean <- df_long %>%
  filter(is.finite(Shannon),
    is.finite(Observed_Richness),
    is.finite(Chao1),
    is.finite(Hill_q1),
    is.finite(Inverse_Simpson),
    is.finite(Pielou))

test_results <- df_clean %>%
  pivot_longer(
    cols = c(Observed_Richness, Shannon, Chao1, Hill_q1, Inverse_Simpson, Pielou),
    names_to = "Metric",
    values_to = "Value") %>%
  group_by(Grouping, Metric) %>%
  group_modify(~ {n_groups <- n_distinct(.x$Category)
    if (n_groups == 2) {wilcox_test(.x, Value ~ Category)} else {kruskal_test(.x, Value ~ Category)}}) %>% ungroup()

# Adjust p-values (multiple testing correction)
test_results <- test_results %>%
  mutate(p_adj = p.adjust(p, method = "BH"))




## ------------------------ Beta diversity ------------------------ ##
otu_mat <- read_excel("/path/to/folder/file", sheet = "otu_table")
tax_mat <- read_excel("/path/to/folder/file", sheet = "tax_table")
sample_metadata <- read_excel("/path/to/folder/file", sheet = "sample_data")
tree_file <- "/path/to/folder/file"
imported_tree <- read.tree(tree_file)

# Prepare otu_table
otu_mat <- as.data.frame(otu_mat)
rownames(otu_mat) <- otu_mat[,1]   # Use first column as row names # Assign row names using the first column
otu_mat <- otu_mat[,-1]            # Remove the first column after setting row names
dim(otu_mat)
otu_mat <- na.omit(otu_mat)       # Remove rows with NA values
dim(otu_mat)

# Prepare tax_table
tax_mat <- as.data.frame(tax_mat)
rownames(tax_mat) <- tax_mat[,1] 
tax_mat <- tax_mat[,-1]
dim(tax_mat)
tax_mat <- as.matrix(tax_mat)

# Prepare sample_data
sample_metadata <- as.data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata[,1]   # Use first column as row names # Assign row names using the first column
sample_metadata <- sample_metadata[,-1]            # Remove the first column after setting row names
dim(sample_metadata)

# Transform into a Phyloseq object
otu_table <- otu_table(otu_mat, taxa_are_rows = TRUE)
tax_table <- tax_table(tax_mat)
sample_data <- sample_data(sample_metadata)
phy_tree <- phy_tree(imported_tree) 
View(otu_table)
View(tax_table)
View(sample_data)
View(phy_tree)
length(taxa_names(otu_table))
length(taxa_names(tax_table))
length(taxa_names(phy_tree))

set.seed(4400)

# Create Phyloseq object
ps <- phyloseq(otu_table, tax_table, sample_data, phy_tree) 

# Weighted UniFrac distance matrix
wunifrac_dist <- distance(ps, method = "wunifrac")
# View distance matrix
wunifrac_dist
View(as.matrix(wunifrac_dist))


# Convert distance matrix to a usable format
sample_metadata <- data.frame(sample_data(ps))

### Calculate PERMANOVA
adonis_result <- adonis2(wunifrac_dist ~ Habitat, data = sample_metadata)
adonis_result
# PERMANOVA interpretation
# A low p-value (p < 0.05) suggests significant differences between groups.
# The R² value shows how much variance is explained by the factor.

### Calculate PAIRWISE PERMANOVA
pairwise_result <- pairwise.adonis2(wunifrac_dist ~ Habitat, data = sample_metadata)
pairwise_result

### Calculate BETADISPERS - checks for homogeneity of group variances (an assumption of PERMANOVA)
beta_disp <- betadisper(wunifrac_dist, sample_metadata$Habitat)
beta_disp_test <- anova(beta_disp) # Test for significance
beta_disp_test
plot(beta_disp, main = "Beta Dispersion by Sample Type") 
# Add group ellipses for visual clarity
ordispider(beta_disp, sample_metadata$Habitat, col = as.factor(sample_metadata$Habitat))
ordihull(beta_disp, sample_metadata$Habitat, col = as.factor(sample_metadata$Habitat), draw = "polygon", label = TRUE)
boxplot(beta_disp, xlab = "Sample Type", ylab = "Distance to Centroid", main = "Beta Dispersion")
# Extract scores (PCoA coordinates)
scores_df <- as.data.frame(scores(beta_disp$vectors))
scores_df$Habitat <- sample_metadata$Habitat
eig_vals <- beta_disp$eig # Extract eigenvalues from the betadisper object
var_exp <- eig_vals / sum(eig_vals[eig_vals > 0]) * 100  # Calculate percentage of variation explained # Use only positive eigenvalues
# Round for labeling
pc1_label <- paste0("PCoA 1 (", round(var_exp[1], 1), "%)")
pc2_label <- paste0("PCoA 2 (", round(var_exp[2], 1), "%)")
# Plot 
betadisp<-ggplot(scores_df, aes(x = PCoA1, y = PCoA2, color = Habitat)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(name="Sample Type", values = Habitat_coloring) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +
  labs(title = "Beta Dispersion (PCoA)",
       x = pc1_label,
       y = pc2_label) +
  theme_bw() +
  theme(legend.position = "right");betadisp
# Betadisper tests whether the dispersion (variance) of community composition 
# within groups is significantly different across your Sample_Type groups. 
# This is important because PERMANOVA assumes equal dispersion — if that's 
# violated, it might falsely detect group differences.
# It doesn't invalidate PERMANOVA entirely, but it raises the possibility 
# that differences in group spread (dispersion) are driving the PERMANOVA result, 
# not actual differences in centroid location (i.e., community structure).

### Calculate TURKEY post-hoc pairwise - explore which groups differ in their dispersion
TukeyHSD(beta_disp)




## ------------------------ Correlations ------------------------ ##
df <- read_excel("/path/to/folder/file", sheet = "Sequences_hgcA_abundance")
df <- column_to_rownames(df, var = colnames(df)[1])
df_long <- df %>%
  rownames_to_column(var = "OTU") %>% 
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance")
metadata <- read_excel("/path/to/folder/file", sheet = "Metadata")
df_long <- df_long %>%
  left_join(metadata, by = "Sample")
head(df_long)

# Shapiro test
df_long <- df_long %>%
  filter(!is.na(pH_Numbers), !is.na(Abundance))

shapiro.test(as.numeric(df_long$pH_Numbers))
shapiro.test(df_long$Abundance)
df_long$Abundance_log <- log(df_long$Abundance)
shapiro.test(df_long$Abundance_log)

# Pearson
cor_results_p <- df_long %>%
  summarise(cor_test = list(cor.test(as.numeric(pH_Numbers), Abundance, method = "pearson"))) %>%
  mutate(correlation = sapply(cor_test, function(x) x$estimate),
         p_value = sapply(cor_test, function(x) x$p.value)) %>%
  select(correlation, p_value); print(cor_results_p)
print(as.data.frame(cor_results_p))

# Spearman
cor_results_s <- df_long %>%
  summarise(cor_test = list(cor.test(as.numeric(pH_Numbers), Abundance, method = "spearman")),
            .groups = "drop") %>%
  mutate(correlation = sapply(cor_test, function(x) x$estimate),
         p_value = sapply(cor_test, function(x) x$p.value)) %>%
  select(correlation, p_value); print(cor_results_s)
print(as.data.frame(cor_results_s))
