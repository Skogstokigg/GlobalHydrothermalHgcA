### Libraries
library(reshape2)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggbreak)
library(writexl)
library(colorspace)
library(RColorBrewer)
library(devtools)
library(leaflet)
library(scatterpie)
library(webr)
library(tibble)


# PRINCIPAL FIGURES
## ------------------------ Fig. 1 - Pie charts ------------------------ ##
pie_all <- read_excel("/path/to/folder/file", sheet = "pie_chart")

# Compute percentages
pie_all$fraction <- pie_all$N_Samples_Strict / sum(pie_all$N_Samples_Strict)*100
# Compute the cumulative percentages (top of each rectangle)
pie_all$ymax <- cumsum(pie_all$fraction)
# Compute the bottom of each rectangle
pie_all$ymin <- c(0, head(pie_all$ymax, n=-1))
# Compute label position
pie_all$labelPosition <- (pie_all$ymax + pie_all$ymin)/2
# Compute a good label
pie_all$label <- paste0(pie_all$Habitat) 

habitat_colour <- c("Seawater" = "#56B4E9", # #56B4E9,
                    "Plume_fluid" = "#E41A1C",
                    "Vent_fluid" = "darkgrey",
                    "Chimney" = "black",
                    "Biofilm" = "#77DD77",
                    "Sediment" = "#FD8504")  #D55E00
names(habitat_colour) <- c("Seawater","Plume_fluid","Vent_fluid",
                           "Chimney","Biofilm","Sediment")
pie_all$label <- factor(pie_all$Habitat, levels = c("Seawater","Plume_fluid", "Vent_fluid", "Chimney", "Biofilm", "Sediment"))
ggplot(pie_all, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  scale_fill_manual(name="Vent Type", values=habitat_colour) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "left") +
  guides(fill=guide_legend(ncol = 1))


## ------------------------ Fig. 2A - Map ------------------------ ##
abundances_strict <- read_excel("/path/to/folder/file", sheet = "S2_HgcA+_Sites_Abundance")
#Load world map data
world_map <- map_data("world")
# Create a transparent theme object
transparent_theme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))
# Map code
mybreaks_normalised_rpob_representation <- c(0.01, 0.1, 1, 1.2) 
map_final <- ggplot(world_map, aes(x = long, y = lat, group = group, colour = habitat_colour))+
  geom_polygon(fill="grey93", colour = "grey93")+
  coord_cartesian(xlim=c(-180,180),ylim=c(-90,90)) +
  theme_light()+xlab("Longitude")+ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(), panel.grid.major = element_blank())+ 
  geom_point(data=abundances_strict[1:28,], 
             aes(y=as.numeric(Latitude),x=as.numeric(Longitude),group=Sites_Code,
                 size=rpoB_HgcA_Strict,
                 alpha=0.5,colour=habitat_colour), 
             stroke=TRUE)+ 
  scale_size_continuous(range = c(4, 25),
                        breaks = mybreaks_normalised_rpob_representation) +
  # scale_size_area(max_size = 20)+
  scale_color_manual(values = vent_classif_coloring)+
  guides(fill=guide_legend(title=c("hgcA Norm. Abundance"), override.aes = list(size = 6))); map_final 

#14.13x6.8

## ------------------------ Fig. 2B - Boxplot by site ------------------------ ##
vent_classif_coloring <- c(
  "Type 1A" = "#000075",  
  "Type 1B" = "#008080",  
  "Type 2A" = "#FFDB6D",  
  "Type 2B" = "#F58231", 
  "Type 3" = "#E6194B",  
  "Type 4" = "#77DD77") 

df <- read_excel("path/to/folder/file", sheet = "Sequences_hgcA_abundance")
metadata <- read_excel("path/to/folder/file", sheet = "Metadata")

# Set row names (assuming the first column contains taxa names)
df <- column_to_rownames(df, var = colnames(df)[1])
# Don't change row names, just keep taxa as a column
df_long <- df %>%
  rownames_to_column(var = "OTU") %>% 
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance")
# Check the structure of the transformed data
head(df_long)

df_long <- all_Sequences %>%
  left_join(metadata, by = "Sample")
head(df_long)
df_long$Site_28 <- as.factor(df_long$Site_28)

# Order from geographical longitud
df_long$Site_28 <- factor(df_long$Site_28, levels = c("VFR_Lau Basin", "Hawaii","JdF", "Guaymas", "Pescadero", "9º50'N", 
                                                      "Cayman", "Snake Pit", "Rainbow", "SMAR_2", "GAKR_2", "SWIR_3", 
                                                      "SOKI", "OKI", "NWPO", "Manus", "Prony", "Brother volcano"))

total_abundance <- df_long %>%
  group_by(Site_28) %>%
  summarise(total_log2_abundance = log2(sum(Abundance) + 1))

Sequences_boxplot_dispersion <- ggplot(df_long, 
                                  aes(x = Site_28, 
                                      y = log2(Abundance+1),
                                      stat = "identity", fill = Vent_Classification)) + 
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.3) +
  scale_fill_manual(values = c(  
    "Type 1A" = "#46F0F0",  
    "Type 1B" = "#008080",  
    "Type 2A" = "#FFDB6D",  
    "Type 2B" = "#F58231",    
    "Type 3" = "#E6194B",    
    "Type 4" = "#77DD77")) +
  labs(x = "Site", y = "Log2(hgcA Abundance + 1)") +
  guides(fill=guide_legend(title="Taxa")) +  
  theme_classic(base_size = 20) + 
  theme(legend.position = "none", #"bottom", 
        strip.text = element_text(size = 30),  # Customize strip text size if needed
        axis.text.x = element_text(angle = 45,hjust = 1,size=20),
        axis.text.y = element_text(angle = 0,hjust = 1,size=20),
        strip.text.y = element_text(size = 30, color = "black")) +
  scale_y_break(c(0.5, 0.7)); Sequences_boxplot_dispersion 

# Papier 10x20


## ------------------------ Fig. 3A - Boxplot by habitat ------------------------ ##
df <- read_excel("/path/to/folder/data", sheet = "Sequences_hgcA_abundance")
metadata <- read_excel("/path/to/folder/file", sheet = "Metadata")
# Set row names (assuming the first column contains taxa names)
df <- column_to_rownames(df, var = colnames(df)[1])
# Don't change row names, just keep taxa as a column
df_long <- df %>%
  rownames_to_column(var = "OTU") %>% 
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance")
# Check the structure of the transformed data
head(df_long)
df_long <- all_Sequences %>%
  left_join(metadata, by = "Sample")
head(df_long)

Sequences_boxplot_dispersion <- ggplot(df_long, 
                                  aes(x = Sample_Type, 
                                      y = log2(Abundance+1),
                                      stat = "identity", fill = Sample_Type)) + 
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.3) +
  scale_fill_manual(values = c("Vent_fluid" = "darkgrey", "Biofilm" = "#77DD77", "Sediment" = "#D55E00", "Chimney" = "black")) +
  labs(x = "Hydrothermal Habitat", y = "Log2(hgcA Abundance + 1)") +
  guides(fill=guide_legend(title="Taxa")) +  
  theme_classic(base_size = 20) + 
  theme(legend.position = "none", #"bottom", 
        strip.text = element_text(size = 30),  
        axis.text.x = element_text(angle = 0,hjust = 0.5,size=30),
        axis.text.y = element_text(angle = 0,hjust = 1,size=30),
        strip.text.y = element_text(size = 30, color = "black")); Sequences_boxplot_dispersion



## ------------------------ Fig. 3B - Boxplot by vent type ----------------------- ##
df <- read_excel("/path/to/folder/file", sheet = "Sequences_hgcA_abundance")
metadata <- read_excel("/path/to/folder/file", sheet = "Metadata")

df <- column_to_rownames(df, var = colnames(df)[1])
# Don't change row names, just keep taxa as a column
df_long <- df %>%
  rownames_to_column(var = "OTU") %>% 
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance")
# Check the structure of the transformed data
head(df_long)

df_long <- all_Sequences %>%
  left_join(metadata, by = "Sample")
head(df_long)


Sequences_boxplot_dispersion <- ggplot(df_long, 
                                  aes(x = Vent_Classification, 
                                      y = log2(Abundance+1),
                                      stat = "identity", fill = Vent_Classification)) + 
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.3) +
  scale_fill_manual(values = c(  
    "Type 1A" = "#46F0F0", #"#000075" # Blues & Cyans 
    "Type 1B" = "#008080",  
    "Type 2" = "#FFDB6D",  # Yellows  
    "Type 2B" = "#F58231",  # Oranges  
    "Type 3" = "#E6194B",  # Reds  
    "Type 4" = "#77DD77")) +
  labs(x = "Geological Vent Type", y = "Log2(hgcA Abundance + 1)") +
  guides(fill=guide_legend(title="Taxa")) +  #"Phylum_Class_[Family/Genus]")) +
  theme_classic(base_size = 20) + #coord_flip() +
  theme(legend.position = "none", #"bottom", 
        strip.text = element_text(size = 30),  # Customize strip text size if needed
        axis.text.x = element_text(angle = 0,hjust = 0.5,size=30),
        axis.text.y = element_text(angle = 0,hjust = 1,size=30),
        strip.text.y = element_text(size = 30, color = "black")) +
  scale_y_break(c(0.45, 0.7)); Sequences_boxplot_dispersion

# 15 x 10 (publish)




## ------------------------ Fig. 4 - Weighted UniFrac (example with habitat (A)) ------------------------ ##
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

# Observe Phyloseq object
print(ps)             # Overview of the phyloseq object
otu_table(ps)         # View OTU table
tax_table(ps)         # View taxonomy table
sample_data(ps)       # View sample metadata

# Weighted UniFrac distance matrix
wunifrac_dist <- distance(ps, method = "wunifrac")
# View distance matrix
wunifrac_dist
View(as.matrix(wunifrac_dist))

# Compute PCoA
ordination <- ordinate(ps, method = "PCoA", distance = "wunifrac")

Habitat_coloring<-c("Vent_fluid" = "blue",  
                        "Biofilm" = "#77DD77", 
                        "Sediment" = "#D55E00", 
                        "Chimney" = "black")
names(Habitat_coloring)<-c("Vent_fluid","Biofilm","Sediment","Chimney")
Habitat_shapes<-c("Vent_fluid" = 16, # circle
                      "Biofilm" = 18, # rombo
                      "Sediment" = 15, # square
                      "Chimney" = 17) # triangle
names(Habitat_shapes)<-c("Vent_fluid","Biofilm","Sediment","Chimney")
# Extract scores (PCoA coordinates)
scores_df <- as.data.frame(scores(ordination$vectors))
scores_df$Habitat <- sample_data(ps)$Habitat
var_exp <- ordination$values$Relative_eig[1:2] * 100
# Round for labeling
pc1_label <- paste0("PCoA 1 (", round(var_exp[1], 1), "%)")
pc2_label <- paste0("PCoA 2 (", round(var_exp[2], 1), "%)")


weighted_pcoa<-ggplot(scores_df, aes(x = Axis.1, y = Axis.2, color = Habitat, shape = Habitat)) +
  geom_point(size = 7, alpha = 0.7) + 
  scale_color_manual(name="Sample Type", values = Habitat_coloring) +
  scale_shape_manual(name="Sample Type shape", values = Habitat_shapes) +
  labs(title = "Weighted UniFrac (PCoA)",
       x = pc1_label,
       y = pc2_label) +
  theme_bw() +
  theme(legend.position = "right");weighted_pcoa




## ------------------------ Fig. 5 - Correlation ------------------------ ##
df <- read_excel("/path/to/folder/file", sheet = "Sequences_hgcA_abundance")
df <- column_to_rownames(df, var = colnames(df)[1])
df_long <- df %>%
  rownames_to_column(var = "OTU") %>% 
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance")
metadata <- read_excel("/path/to/folder/file", sheet = "Metadata")
df_long <- df_long %>%
  left_join(metadata, by = "Sample")
head(df_long)

# PCoA1 vs. pH
ggplot(df_long, aes(x = as.numeric(pH_Numbers), y = Axis.1, colour = Vent_Classification)) +
  geom_smooth(method = "lm", color = "blue") +
  geom_point(size = 5, alpha = 0.7) + 
  scale_colour_manual(values = vent_classif_coloring) +
  theme_bw() +
  labs(y = "PCoA1", x = "pH") +
  scale_y_log10(oob = scales::squish_infinite) +
  coord_cartesian(xlim = range(as.numeric(df_long$pH_Numbers)))

# Abundance vs. pH
ggplot(df_long, aes(x = as.numeric(pH_Numbers), y = Abundance, colour = Vent_Classification)) +
  geom_smooth(method = "lm", color = "blue") +
  geom_point(size = 5, alpha = 0.8) + 
  scale_colour_manual(values = vent_classif_coloring) +
  theme_bw() +
  labs(y = "Normalized hgcA abundance", x = "pH", color = "Vent Type") +
  scale_y_log10(oob = scales::squish_infinite)




## ------------------------ Fig. 6 - Barplots HgcA abundance ------------------------ ##

sorted_color_palette <- c(
  "#E31A1C", "#D6604D", "#B2182B", "#FB9A99", "#FDBF6F", "#FF7F00", "#F0E442", "#FFD92F", "#FFFF99","#F2E802",
  "#E69F00", "#6A3D9A", "#B2DF8A", "#33A02C", "#1B7837", "#009E73", "#56B4E9",
  "#A6CEE3", "#4393C3", "#0072B2", "#999999", "grey")
extended_palette_dpc <- colorRampPalette(sorted_color_palette)(39) 
Sequences_752_strict <- read_excel("/path/to/folder/file", sheet = "Sequences_hgcA")

## By Habitat (Fig. 6A)
# Calculate abundances in %
all_Sequences <- Sequences_752_strict[1:752,]
all_Sequences <- all_Sequences %>%
  group_by(Sample_Type) %>%
  mutate(Normalised_Abundance_percent = (Normalised_Abundance_hgcA / sum(Normalised_Abundance_hgcA)) * 100) %>%
  ungroup()

# Prepare data
stack_data <- all_Sequences %>%
  group_by(Sample_Type, dp, dpc) %>%
  summarise(value = sum(Normalised_Abundance_percent), .groups = "drop") %>%
  group_by(Sample_Type) %>%
  mutate(row_order = rev(row_number())) %>%
  arrange(Sample_Type, row_order) %>%
  mutate(
    ymin = cumsum(lag(value, default = 0)),
    ymax = ymin + value
  ) %>%
  select(-row_order) %>%
  ungroup()
stack_data <- as.data.frame(stack_data)

# Plot
Sequences_bar_plot <- ggplot(all_Sequences) + 
  geom_bar(aes(x = Sample_Type, y = Normalised_Abundance_percent, fill = dpc),
           stat = "identity", position = "stack", width = 0.95) +
  scale_fill_manual(values = extended_palette_dpc) + 
  labs(title = "Abundance of Taxa (as % within each Habitat)", 
       x = "Habitat", y = "Relative Abundance (%)") +
  guides(fill = guide_legend(title = "Taxa")) +  
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        strip.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 25),
        axis.text.y = element_text(size = 30),
        strip.text.y = element_text(size = 30, color = "black")) +
  geom_segment(data = stack_data,
               aes(x = as.numeric(factor(Sample_Type)) - 0.475,
                   xend = as.numeric(factor(Sample_Type)) + 0.475,
                   y = ymax, yend = ymax),
               color = "white", linewidth = 0.7, alpha = 0.3); Sequences_bar_plot



## By Vent Type (Fig. 6B)
# Calculate abundances in % 
all_Sequences <- Sequences_752_strict[1:752,]
all_Sequences <- all_Sequences %>%
  group_by(Vent_Type) %>%
  mutate(Normalised_Abundance_percent = (Normalised_Abundance_hgcA / sum(Normalised_Abundance_hgcA)) * 100) %>%
  ungroup()

# Prepare data
stack_data <- all_Sequences %>%
  group_by(Vent_Type, dp, dpc) %>%
  summarise(value = sum(Normalised_Abundance_percent), .groups = "drop") %>%
  group_by(Vent_Type) %>%
  mutate(row_order = rev(row_number())) %>%
  arrange(Vent_Type, row_order) %>%
  mutate(
    ymin = cumsum(lag(value, default = 0)),
    ymax = ymin + value
  ) %>%
  select(-row_order) %>%
  ungroup()
stack_data <- as.data.frame(stack_data)

# Plot
Sequences_bar_plot <- ggplot(all_Sequences) + 
  geom_bar(aes(x = Vent_Type, y = Normalised_Abundance_percent, fill = dpc),
           stat = "identity", position = "stack", width = 0.95) +
  scale_fill_manual(values = extended_palette_dpc) +  
  labs(title = "Abundance of Taxa (as % within each Vent Type)", 
       x = "Vent Type", y = "Relative Abundance (%)") +
  guides(fill = guide_legend(title = "Taxa")) +  
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        strip.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 25),
        axis.text.y = element_text(size = 30),
        strip.text.y = element_text(size = 30, color = "black")) +
  geom_segment(data = stack_data,
               aes(x = as.numeric(factor(Vent_Type)) - 0.475,
                   xend = as.numeric(factor(Vent_Type)) + 0.475,
                   y = ymax, yend = ymax),
               color = "white", linewidth = 0.7, alpha = 0.3); Sequences_bar_plot



## By site (Fig. 6C)
# Calculate abundances in % for each site
all_Sequences <- Sequences_752_strict[1:752,]
all_Sequences <- all_Sequences %>%
  group_by(Site_28) %>%
  mutate(Normalised_Abundance_percent = (Normalised_Abundance_hgcA / sum(Normalised_Abundance_hgcA)) * 100) %>%
  ungroup()
all_Sequences$Site_28 <- factor(all_Sequences$Site_28, levels = c("VFR_Lau Basin", "Hawaii",  "JdF", "Guaymas", "Pescadero", "9º50'N", 
                                                                  "Cayman", "Snake Pit", "Rainbow", "SMAR_2", "GAKR_2", "SWIR_3", 
                                                                  "SOKI", "OKI", "NWPO", "Manus", "Prony", "Brother Volcano"))

# Prepare data
stack_data <- all_Sequences %>%
  group_by(Site_28, dp, dpc) %>%
  summarise(value = sum(Normalised_Abundance_percent), .groups = "drop") %>%
  group_by(Site_28) %>%
  mutate(row_order = rev(row_number())) %>%
  arrange(Site_28, row_order) %>%
  mutate(
    ymin = cumsum(lag(value, default = 0)),
    ymax = ymin + value
  ) %>%
  select(-row_order) %>%
  ungroup()
stack_data <- as.data.frame(stack_data)

# Plot
Sequences_bar_plot <- ggplot(all_Sequences) + 
  geom_bar(aes(x = Site_28, y = Normalised_Abundance_percent, fill = dpc),
           stat = "identity", position = "stack", width = 0.95) +
  scale_fill_manual(values = extended_palette_dpc) + # extended_palette_dp_strict
  labs(title = "Abundance of Taxa (as % within each Site)", 
       x = "Site", y = "Relative Abundance (%)") +
  guides(fill = guide_legend(title = "Taxa")) +  
  theme_classic(base_size = 20) + 
  theme(legend.position = "right",
        strip.text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
        axis.text.y = element_text(size = 30),
        strip.text.y = element_text(size = 30, color = "black")) +
  geom_segment(data = stack_data,
               aes(x = as.numeric(factor(Site_28)) - 0.475,
                   xend = as.numeric(factor(Site_28)) + 0.475,
                   y = ymax, yend = ymax),
               color = "white", linewidth = 0.7, alpha = 0.6); Sequences_bar_plot






