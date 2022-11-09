## Initialization
library(STvEA)
library(AdjacencyScore)
library(RANN)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(patchwork)
library(monocle3)
library(ggplot2)
library(gridExtra)
library(stringr)
library(flowCore)
library(batchelor)
library(igraph)
library(ggraph)
library(mixOmics)
library(Cairo)
source("Additional CODEX Functions.R")


# Prior to this point, export the .fcs files for each imaged CODEX region and concatenate into a single file (all_regions.fcs) using flowJo, Excel, or equivalent software.
# Additionally, define stain gating in the native CODEX software (multiplex analysis viewer; MAV) and export as a single file (gate.csv)


## Pre-Processing Raw Data in STvEA (Govek et al. Science Advances 2019.)
# Read in channel names from .fcs file
codex_names <- ReadNamesFCS("all_regions.fcs")[["channels"]]
# Store names of blank, non-protein, and DAPI "channels"
is_blank <- grepl("Blank",codex_names)
is_protein <- !(codex_names %in% c("cell_id","region","tile_num","x","y","z","x_tile","y_tile","size","tsne_x","tsne_y","homogeneity"))
is_DAPI <- grepl("DAPI",codex_names)
# Read in dataset with information defined above
stvea_object <- ReadDataFCS("all_regions.fcs", is_protein, is_blank, protein_names = codex_names[is_protein])
# Read in MAV-defined gating
gate <- read.delim(file = "gate.csv", header = TRUE, sep = ",")
gate <- as.logical(unlist(gate[,2]))
# Filter out stained debris using method from STvEA paper and MAV-defined gating
stvea_object <- FilterCODEX(stvea_object, inclusion = gate)
# Normalize staining intensities to Gaussian distribution
stvea_object <- CleanCODEX(stvea_object, model="gaussian", num_cores=4)


## Batch Correction in mnnCorrect (Haghverdi et al. Nature Biotechnology 2018.)
# Set batch correction parameters
b <- 25
# Perform batch correction between image regions using mutual nearest neighbors
temp_object <- mnnCorrect(t(stvea_object@codex_clean), batch = stvea_object@codex_region, k = b, cos.norm.in = FALSE, cos.norm.out = FALSE, subset.row = !grepl("DAPI",codex_names[is_protein]) & !grepl("Blank",codex_names[is_protein]), correct.all = TRUE)
stvea_object@codex_clean <- t(temp_object@assays@data@listData$corrected)
rm(temp_object)
gc()


## Generating UMAP Manifold
# Generate counts matrix for testing dataset
seurat_matrix <- stvea_object@codex_clean[,!grepl("Blank",colnames(stvea_object@codex_clean))]
rownames(seurat_matrix) <- paste(stvea_object@codex_region, stvea_object@codex_cell_ID, sep = "_")
seurat_matrix <- t(seurat_matrix)
# Create Seurat object for testing dataset
# Set UMAP parameters of interest
m <- "cosine"; i <- 500; j <- 0.1; k <- 1; l <- 100
# Create Seurat object and run minimum required preprocessing
seurat_object <- CreateSeuratObject(counts = seurat_matrix)
rm(seurat_matrix)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = rownames(seurat_object), npcs = 5)
# UMAP processing (timed)
umap_start_time <- Sys.time()
seurat_object <- RunUMAP(seurat_object, umap.method = "umap-learn", return.model = TRUE, dims = 1:5, metric = m, n.neighbors = i, min.dist = j, repulsion.strength = k, negative.sample.rate = l)
umap_end_time <- Sys.time()
pushover(paste0("PDAC UMAP Complete (",m,", ",i,", ",j,", ",k,", ",l,")"))
# Generate UMAP clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)
seurat_object <- FindClusters(seurat_object, resolution = 0.5) # Use 0.1 if FindNeighbors used UMAP, 0.5 if FindNeighbors used PCA (default)
# Plot UMAP manifold
DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + NoLegend()
gc()


## Transferring Information Back to STvEA
# Populate STvEA object with Seurat UMAP, clusters, etc.
stvea_object@codex_emb <- data.frame(row.names = rownames(stvea_object@codex_protein))
stvea_object@codex_emb$V1 <- seurat_object[["umap"]]@cell.embeddings[,1]
stvea_object@codex_emb$V2 <- seurat_object[["umap"]]@cell.embeddings[,2]
stvea_object@codex_clusters <- seurat_object$seurat_clusters
# Use concatenated cell ID as in Seurat
stvea_object@codex_cell_ID <- paste(stvea_object@codex_region,stvea_object@codex_cell_ID, sep = "_")
# Post-manifold threshold: remove region-specific artifacts caused by heterogeneity in imaging process
stvea_object <- RemoveArtifacts(stvea_object, size_threshold = 0.1, outlier_threshold = 0.4)
# Apply filter above to Seurat object
seurat_object <- seurat_object[,stvea_object@codex_cell_ID]
# Plot manifold in STvEA
PlotClusterCODEXemb(stvea_object)


## Cell Type Annotations
# Import cell type annotations into STvEA
loaded_annotations <- read.delim(file = "cluster_annotations.csv", header = TRUE, row.names=1, sep = ",")
loaded_annotations <- t(loaded_annotations)
stvea_object@codex_clusters <- as.factor(loaded_annotations[,as.character(stvea_object@codex_clusters)])
# Add cell type annotations to Seurat
seurat_object$annotated_cell_types <- as.character(stvea_object@codex_clusters)


## Additional Sub-Clustering
# Get current cluster information
codex_clusters <- as.character(stvea_object@codex_clusters)
names(codex_clusters) <- stvea_object@codex_cell_ID
# Separate Helper and Cytotoxic T Lymphocytes by CD8 expression
codex_clusters[stvea_object@codex_clusters == "T Lymphocytes" & stvea_object@codex_clean[,"CD8"] < 0.62] <- "Helper T Lymphocytes"
codex_clusters[stvea_object@codex_clusters == "T Lymphocytes" & stvea_object@codex_clean[,"CD8"] >= 0.62] <- "Cytotoxic T Lymphocytes" 
# Separate Macrophages subtypes by HLA-DR expression
codex_clusters[stvea_object@codex_clusters == "Macrophages" & stvea_object@codex_clean[,"HLA-DR"] < 0.62] <- "Macrophages 1"
codex_clusters[stvea_object@codex_clusters == "Macrophages" & stvea_object@codex_clean[,"HLA-DR"] >= 0.62] <- "Macrophages 2"
# Transfer new clusters back to STvEA
stvea_object@codex_clusters <- as.factor(codex_clusters)
seurat_object$annotated_cell_types <- codex_clusters
rm(codex_clusters,int_score_summary,fibroblast_names)


## Set Color Scheme for Plots
# Extract color scheme from Seurat UMAP plot
cluster_colors <- DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = "annotated_cell_types")
cluster_colors <- ggplot_build(cluster_colors)
cluster_colors <- data.frame(colours = unique(cluster_colors$data[[1]]["colour"]), label = cluster_colors$plot$scales$scales[[1]]$get_labels()[unlist(unique(cluster_colors$data[[1]]["group"]))]) 
cluster_colors <- cluster_colors[order(cluster_colors[,"label"]),]
cluster_colors <- setNames(cluster_colors$colour, cluster_colors$label)


## Cell Spatial Plots
# Plot cells on original section (spatial plot)
PlotExprCODEXspatial(stvea_object, type = "cluster", region = "1")


## Processing of Patient Metadata
# Load in patient metadata
loaded_metadata <- read.delim(file = "patient_metadata.csv", header = TRUE, row.names=1, sep = ",")
loaded_metadata <- t(loaded_metadata)
# Set blank values to "NA" for compatibility
loaded_metadata[loaded_metadata == ""] <- NA
# Add new columns for processed metadata categories
temp_coords <- matrix(ncol=length(loaded_metadata["Patient.ID",]), nrow = 2)
colnames(temp_coords) <- loaded_metadata["Patient.ID",]
rownames(temp_coords) <- c("Processed.DFS","Survivorship.Category")
colnames(loaded_metadata) <- loaded_metadata["Patient.ID",]
loaded_metadata <- rbind(loaded_metadata,temp_coords)
rm(temp_coords)
# Populate processed metadata categories
for (patient in 1:length(loaded_metadata['DFS',])) # If disease free survival is not listed and patient is deceased, disease free survival is by definition equivalent to overall survival
{
  if (is.na(loaded_metadata['DFS',patient]) & !(is.na(loaded_metadata['Alive.or.Dead',patient])) & loaded_metadata['Alive.or.Dead',patient] == 'Dead')
  {
    loaded_metadata['Processed.DFS',patient] <- loaded_metadata['OS',patient]
  }
  else
  {
    loaded_metadata['Processed.DFS',patient] <- loaded_metadata['DFS',patient]
  }
}
# Create categories for initial ML classification (downstream of CODEX analysis)
for (patient in 1:length(loaded_metadata['Alive.or.Dead',])) # Set survivorship category based on alive or dead status and overall survival
{
  if (!(is.na(loaded_metadata['Alive.or.Dead',patient])) & loaded_metadata['Alive.or.Dead',patient] == "Alive" & !(is.na(loaded_metadata['OS',patient])) & as.numeric(loaded_metadata['OS',patient]) < 730)
  {
    loaded_metadata['Survivorship.Category',patient] <- NA
  } else if (!(is.na(loaded_metadata['Alive.or.Dead',patient])) & loaded_metadata['Alive.or.Dead',patient] == "Alive" & !(is.na(loaded_metadata['OS',patient])) & as.numeric(loaded_metadata['OS',patient]) >= 730)
  {
    loaded_metadata['Survivorship.Category',patient] <- "Good.Survivorship"
  } else if (!(is.na(loaded_metadata['Alive.or.Dead',patient])) & loaded_metadata['Alive.or.Dead',patient] == "Dead" & !(is.na(loaded_metadata['OS',patient])) & as.numeric(loaded_metadata['OS',patient]) < 730)
  {
    loaded_metadata['Survivorship.Category',patient] <- "Poor.Survivorship"
  } else if (!(is.na(loaded_metadata['Alive.or.Dead',patient])) & loaded_metadata['Alive.or.Dead',patient] == "Dead" & !(is.na(loaded_metadata['OS',patient])) & as.numeric(loaded_metadata['OS',patient]) >= 730)
  {
    loaded_metadata['Survivorship.Category',patient] <- "Good.Survivorship"
  } else
  {
    loaded_metadata['Survivorship.Category',patient] <- NA
  }
}
# Define discrete patient groups for later testing
patient_groups <- c(unique(loaded_metadata["ECM.Category",]),unique(loaded_metadata["Recurrence.Category",]),unique(loaded_metadata["Survivorship.Category",]),unique(loaded_metadata["Neoadjuvant.Category",]),unique(loaded_metadata["Adjuvant.Category",]),unique(loaded_metadata["Lymphovascular.Category",]),unique(loaded_metadata["Perineural.Category",]))
patient_groups <- patient_groups[!is.na(patient_groups)]


## Pre-Analysis
# Store region and cluster filters for usage in later calculation
region_filters <- list()
cluster_filters <- list()
for (region in as.character(unique(stvea_object@codex_region)))
{
  region_filters[[region]] <- stvea_object@codex_region == region
}
for (cluster in as.character(unique(stvea_object@codex_clusters)))
{
  cluster_filters[[cluster]] <- stvea_object@codex_clusters == cluster
}
# Specify discrete group comparisons
group_comparisons <- c("Y.Neoadjuvant vs. N.Neoadjuvant")
# Specify continuous categories of interest
continuous_comparisons <- c("OS","Processed.DFS")
# Specify ordinal categories of interest
ordinal_comparisons <- NULL
# Specify survival-positive and survival-negative categories
survival_positive_categories <- c("OS","Processed.DFS")
survival_negative_categories <- NULL


## Cell-Level Analysis
# Quantify cell types in each specimen
cluster_proportions <- matrix(nrow = length(unique(stvea_object@codex_clusters)), ncol = length(unique(stvea_object@codex_region)))
rownames(cluster_proportions) <- as.character(unique(stvea_object@codex_clusters))
colnames(cluster_proportions) <- as.character(unique(stvea_object@codex_region))
for (region_to_count in colnames(cluster_proportions))
{
  for (cluster_to_count in rownames(cluster_proportions))
  {
    cluster_proportions[cluster_to_count,region_to_count] <- sum(region_filters[[region_to_count]] & cluster_filters[[cluster_to_count]])/sum(region_filters[[region_to_count]])
  }
}


## Spatial Analysis - Region Level
# Generate all possible cell pairs
cell_types <- sort(as.character(unique(stvea_object@codex_clusters)))
cell_pairs <- t(combn(cell_types,2))
for (cell in cell_types) 
{
  cell_pairs <- rbind(cell_pairs, c(cell,cell))
}
# Create empty matrix to store for excessively small regions 
empty_matrix <- matrix(nrow = nrow(cell_pairs), ncol = 3)
rownames(empty_matrix) <- paste(cell_pairs[,1], cell_pairs[,2], sep = " - ")
colnames(empty_matrix) <- c("f1","f2","score") 
empty_matrix[,"f1"] <- cell_pairs[,1]
empty_matrix[,"f2"] <- cell_pairs[,2]
empty_matrix[rownames(empty_matrix),"score"] <- 0 
# AUTOMATED: Generate interaction scores 
k_interaction_scores <- list()
k_interaction_network <- list()
for (interaction_k in c(20))
{
  print(paste0("BEGINNING CALCULATIONS FOR k=",interaction_k,"."))
  # Generate interaction scores and network for each region 
  interaction_scores <- list()
  interaction_network <- list()
  for (region in as.character(unique(stvea_object@codex_region)))
  {
    print(paste0("Calculating interaction scores for region: ",region,"."))
    n_cells <- sum(region_filters[[region]])
    # If too few cells, populate with empty matrix
    if (n_cells < interaction_k)
    {
      interaction_scores[[region]] <- empty_matrix
    }
    else
    {
      # Generate adjacency matrix based on cell locations
      adjacency_matrix <- knn_graph(stvea_object@codex_spatial[region_filters[[region]],], k=interaction_k)
      # Generate feature matrix of assigned cell types
      cell_matrix <- t(sapply(cell_types, function(x) (stvea_object@codex_clusters[region_filters[[region]]]==x)*1))
      row.names(cell_matrix) <- cell_types
      # Calculate interaction scores as dot-product of paired cell feature vectors with adjacency matrix (f1*j*f2)
      interaction_scores[[region]] <- CalculateInteractionScores(adjacency_matrix = adjacency_matrix, cell_matrix = cell_matrix, cell_pairs = cell_pairs, calculate_individual_scores = FALSE, normalize = TRUE)
      # Generate interaction matrix and resulting network for plotting
      interaction_score_matrix <- AdjScoreHeatmapMatrix(interaction_scores[[region]])
      class(interaction_score_matrix) <- "numeric"
      interaction_network[[region]] <- graph_from_adjacency_matrix(interaction_score_matrix, weighted = TRUE)
    }
  }
  # Store interaction scores and networks
  k_interaction_scores[[as.character(interaction_k)]] <- interaction_scores
  k_interaction_network[[as.character(interaction_k)]] <- interaction_network
  rm(adjacency_matrix,cell_matrix,interaction_score_matrix)
}
rm(interaction_scores,interaction_network)


## Spatial Analysis - Group Level
# AUTOMATED: calculate average interaction scores for each group and produce resulting interaction network
k_group_interaction_scores <- list()
k_group_interaction_scores_with_zeroes <- list()
k_group_interaction_network <- list()
for (interaction_k in c(20))
{
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  # Calculate average interaction scores for each group and produce resulting interaction network
  group_interaction_scores <- list()
  group_interaction_scores_with_zeroes <- list()
  group_interaction_network <- list()
  for (group in patient_groups)
  {
    # Define group-specific regions and variable for storing group-specific interaction scores
    metadata_category <- paste0(strsplit(group, split = "\\.")[[1]][2],".Category")
    group_regions <- names(which(loaded_metadata[metadata_category,] == group))
    group_regions <- group_regions[group_regions %in% as.character(unique(stvea_object@codex_region))]
    group_interaction_scores[[group]] <- matrix(nrow = nrow(cell_pairs), ncol=3)
    colnames(group_interaction_scores[[group]]) <- c("f1","f2","score") 
    group_interaction_scores[[group]][,"f1"] <- cell_pairs[,1]
    group_interaction_scores[[group]][,"f2"] <- cell_pairs[,2]
    rownames(group_interaction_scores[[group]]) <- paste(group_interaction_scores[[group]][,"f1"], group_interaction_scores[[group]][,"f2"], sep = " - ")
    # For each possible cell pair, calculate the average interaction score across all group-specific regions
    for (pair in 1:nrow(cell_pairs))
    {
      score_sum <- 0
      for (region in group_regions)
      {
        matching_row_1 <- (interaction_scores[[region]][,"f1"] == cell_pairs[pair,1]) & (interaction_scores[[region]][,"f2"] == cell_pairs[pair,2])
        matching_row_2 <- (interaction_scores[[region]][,"f2"] == cell_pairs[pair,1]) & (interaction_scores[[region]][,"f1"] == cell_pairs[pair,2])
        if (sum(matching_row_1) == 1)
        {
          if (!is.na(interaction_scores[[region]][matching_row_1,"score"]))
          {
            score_sum <- score_sum + as.numeric(interaction_scores[[region]][matching_row_1,"score"])
          }
        }
        else if (sum(matching_row_2) == 1)
        {
          if (!is.na(interaction_scores[[region]][matching_row_2,"score"]))
          {
            score_sum <- score_sum + as.numeric(interaction_scores[[region]][matching_row_2,"score"])
          }
        }
      }
      if (score_sum != 0)
      {
        group_interaction_scores[[group]][pair,"score"] <- score_sum/length(group_regions)
      }
    }
    # Save full interaction score variable (including zero/NA values) for later usage in differential analysis
    group_interaction_scores_with_zeroes[[group]] <- group_interaction_scores[[group]]
    group_interaction_scores_with_zeroes[[group]][is.na(group_interaction_scores_with_zeroes[[group]][,"score"]),"score"] <- 0
    # Remove non-existing interactions from group-averaged interaction score variable
    group_interaction_scores[[group]] <- group_interaction_scores[[group]][!is.na(group_interaction_scores[[group]][,"score"]),]
    # Generate interaction matrix and resulting network for plotting
    group_interaction_score_matrix <- AdjScoreHeatmapMatrix(group_interaction_scores[[group]])
    class(group_interaction_score_matrix) <- "numeric"
    group_interaction_network[[group]] <- graph_from_adjacency_matrix(group_interaction_score_matrix, weighted = TRUE)
  }
  # Store group interaction scores and networks
  k_group_interaction_scores[[as.character(interaction_k)]] <- group_interaction_scores
  k_group_interaction_scores_with_zeroes[[as.character(interaction_k)]] <- group_interaction_scores_with_zeroes
  k_group_interaction_network[[as.character(interaction_k)]] <- group_interaction_network
  rm(group_interaction_score_matrix)
}
rm(group_interaction_scores,group_interaction_scores_with_zeroes,group_interaction_network)


## Spatial Analysis - Differential Level
# AUTOMATED: calculate differential interaction scores and produce resulting differential interaction network
k_differential_interaction_scores <- list()
k_differential_interaction_network <- list()
k_differential_interaction_color <- list()
for (interaction_k in c(20))
{
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  group_interaction_scores_with_zeroes <- k_group_interaction_scores_with_zeroes[[as.character(interaction_k)]]
  # For discrete groups, calculate differential interaction scores and generate resulting differential interaction networks
  differential_interaction_scores <- list()
  differential_interaction_network <- list()
  differential_interaction_color <- list()
  for (comparison in group_comparisons)
  {
    # Calculate differences in interaction scores between groups
    group_1 <- str_split(comparison," vs. ")[[1]][1]
    group_2 <- str_split(comparison," vs. ")[[1]][2]
    differential_interaction_scores[[comparison]] <- group_interaction_scores_with_zeroes[[group_1]]
    differential_interaction_scores[[comparison]][,"score"] <- as.numeric(group_interaction_scores_with_zeroes[[group_2]][,"score"]) - as.numeric(group_interaction_scores_with_zeroes[[group_1]][,"score"])
    # Generate differential interaction score matrix
    differential_interaction_score_matrix <- AdjScoreHeatmapMatrix(differential_interaction_scores[[comparison]])
    class(differential_interaction_score_matrix) <- "numeric"
    # Convert differential interaction score matrix to absolute values and separate color variables for plotting
    differential_interaction_color[[comparison]] <- differential_interaction_score_matrix
    differential_interaction_color[[comparison]][differential_interaction_score_matrix < 0] <- "Negative"
    differential_interaction_color[[comparison]][differential_interaction_score_matrix > 0] <- "Positive"
    differential_interaction_score_matrix <- abs(differential_interaction_score_matrix)
    # Generate differential interaction network
    differential_interaction_network[[comparison]] <- graph_from_adjacency_matrix(differential_interaction_score_matrix, weighted = TRUE)
  }
  rm(differential_interaction_score_matrix)
  # For continuous/ordinal categories, calculate differential interaction scores (Pearson/Spearman correlations) and generate resulting differential interaction networks
  for (comparison in c(continuous_comparisons,ordinal_comparisons))
  {
    # Define variable for storing Pearson correlations 
    differential_interaction_scores[[comparison]] <- matrix(nrow = nrow(cell_pairs), ncol=3)
    colnames(differential_interaction_scores[[comparison]]) <- c("f1","f2","score") 
    differential_interaction_scores[[comparison]][,"f1"] <- cell_pairs[,1]
    differential_interaction_scores[[comparison]][,"f2"] <- cell_pairs[,2]
    rownames(differential_interaction_scores[[comparison]]) <- paste(differential_interaction_scores[[comparison]][,"f1"], differential_interaction_scores[[comparison]][,"f2"], sep = " - ")
    # For each possible cell pair, calculate the Pearson correlation with the metadata category of interest
    for (pair in 1:nrow(cell_pairs))
    {
      score_values <- vector(mode = "numeric", length = length(unique(stvea_object@codex_region)))
      names(score_values) <- as.character(unique(stvea_object@codex_region))
      # For particular cell pair, extract each region's score value
      for (region in names(score_values))
      {
        matching_row_1 <- (interaction_scores[[region]][,"f1"] == cell_pairs[pair,1]) & (interaction_scores[[region]][,"f2"] == cell_pairs[pair,2])
        matching_row_2 <- (interaction_scores[[region]][,"f2"] == cell_pairs[pair,1]) & (interaction_scores[[region]][,"f1"] == cell_pairs[pair,2])
        if (sum(matching_row_1) == 1)
        {
          if (!is.na(interaction_scores[[region]][matching_row_1,"score"]))
          {
            score_values[region] <- interaction_scores[[region]][matching_row_1,"score"]
          }
        }
        else if (sum(matching_row_2) == 1)
        {
          if (!is.na(interaction_scores[[region]][matching_row_2,"score"]))
          {
            score_values[region] <- interaction_scores[[region]][matching_row_2,"score"]
          }
        }
      }
      # Calculate and store Pearson/Spearman correlation 
      if (comparison %in% continuous_comparisons)
      {
        correlation <- cor(x = as.numeric(loaded_metadata[comparison,names(score_values)]), y = as.numeric(score_values), use = "na.or.complete")
      }
      if (comparison %in% ordinal_comparisons)
      {
        correlation <- cor(x = as.numeric(loaded_metadata[comparison,names(score_values)]), y = as.numeric(score_values), use = "na.or.complete", method = "spearman")
      }
      if (!is.null(correlation))
      {
        differential_interaction_scores[[comparison]][pair,"score"] <- correlation
      }
      if(is.na(correlation))
      {
        differential_interaction_scores[[comparison]][pair,"score"] <- 0
      }
    }
    # Generate differential interaction score matrix
    differential_interaction_score_matrix <- AdjScoreHeatmapMatrix(differential_interaction_scores[[comparison]])
    class(differential_interaction_score_matrix) <- "numeric"
    # Convert differential interaction score matrix to absolute values and separate color variables for plotting
    # Note colors are set so that positive correlations produce blue lines and negative correlations produce red lines
    # This is done by manipulating the alphabetical order of the positive/negative designations
    differential_interaction_color[[comparison]] <- differential_interaction_score_matrix
    if (comparison %in% survival_positive_categories)
    {
      differential_interaction_color[[comparison]][differential_interaction_score_matrix < 0] <- "Survival-Negative"
      differential_interaction_color[[comparison]][differential_interaction_score_matrix > 0] <- "Positive"
    }
    if (comparison %in% survival_negative_categories)
    {
      differential_interaction_color[[comparison]][differential_interaction_score_matrix < 0] <- "Negative"
      differential_interaction_color[[comparison]][differential_interaction_score_matrix > 0] <- "Positive"
    }
    differential_interaction_score_matrix <- abs(differential_interaction_score_matrix)
    # Generate differential interaction network
    differential_interaction_network[[comparison]] <- graph_from_adjacency_matrix(differential_interaction_score_matrix, weighted = TRUE)
  }
  # Store differential interaction scores, networks, and colors
  k_differential_interaction_scores[[as.character(interaction_k)]] <- differential_interaction_scores
  k_differential_interaction_network[[as.character(interaction_k)]] <- differential_interaction_network
  k_differential_interaction_color[[as.character(interaction_k)]] <- differential_interaction_color
  rm(differential_interaction_score_matrix)
}
rm(differential_interaction_scores,differential_interaction_network,differential_interaction_color)
# Plot differential interaction network
interaction_k <- 20 # Specify k of interest
comparison_of_interest <- "OS" # Specify comparison of interest (can be double-discrete or single-continuous)
# Populate variables for plotting
differential_interaction_network <- k_differential_interaction_network[[as.character(interaction_k)]]
differential_interaction_color <- k_differential_interaction_color[[as.character(interaction_k)]]
# Plot interaction network for specific comparison
ggraph(differential_interaction_network[[comparison_of_interest]], layout = "linear", circular = TRUE) + geom_edge_link(aes(edge_width=E(differential_interaction_network[[comparison_of_interest]])$weight, edge_colour = as.factor(as.vector(differential_interaction_color[[comparison_of_interest]][differential_interaction_color[[comparison_of_interest]] != 0])))) + geom_node_point( color="black", size=8) + geom_node_text( aes(label=name), repel = TRUE, size=8, color="#681984") + theme_void() + theme(legend.position="none", plot.margin=unit(rep(1,4), "cm")) + scale_edge_colour_manual(values = c("#317cca48","#C43d3d48"))
# Plot interaction network with self-loops for specific comparison
ggraph(differential_interaction_network[[comparison_of_interest]], layout = "linear", circular = TRUE) + geom_edge_link(aes(edge_width=E(differential_interaction_network[[comparison_of_interest]])$weight, edge_colour = as.factor(as.vector(differential_interaction_color[[comparison_of_interest]][differential_interaction_color[[comparison_of_interest]] != 0])))) + geom_node_point( color="black", size=8) + geom_node_text( aes(label=name), repel = TRUE, size=8, color="#681984") + geom_edge_loop(aes(strength = 0.25, direction = 90, edge_width=E(differential_interaction_network[[comparison_of_interest]])$weight, edge_colour = as.factor(as.vector(differential_interaction_color[[comparison_of_interest]][differential_interaction_color[[comparison_of_interest]] != 0])))) + theme_void() + theme(legend.position="none", plot.margin=unit(rep(1,4), "cm")) + scale_edge_colour_manual(values = c("#317cca48","#C43d3d48"))


## Spatial Analysis - Individual Interactions
k_interaction_score_overview <- list()
for (interaction_k in c(20))
{  
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  # Extract individual interaction scores per region
  interaction_score_overview <- matrix(nrow = nrow(cell_pairs), ncol = length(unique(stvea_object@codex_region)))
  rownames(interaction_score_overview) <- paste(cell_pairs[,1],cell_pairs[,2], sep = " - ")
  colnames(interaction_score_overview) <- as.character(unique(stvea_object@codex_region))
  for (interaction in rownames(interaction_score_overview))
  {
    for (region in colnames(interaction_score_overview))
    {
      interaction_score_overview[interaction,region] <- as.numeric(interaction_scores[[region]][interaction,"score"])
    }
  }
  k_interaction_score_overview[[as.character(interaction_k)]] <- interaction_score_overview
}


## Initialize Variables for Subsample-Level Analysis
# Set number of desired subsamples per region and calculate intervals by which to split x and y coordinates
x_subsample <- 10
y_subsample <- 10
x_interval <- max(stvea_object@codex_spatial[,"x"])/x_subsample
y_interval <- max(stvea_object@codex_spatial[,"y"])/y_subsample
# Generate spatial filters for later usage
x_spatial_filters <- list()
y_spatial_filters <- list()
for (x in 1:x_subsample)
{
  x_spatial_filters[[as.character(x)]] <- (stvea_object@codex_spatial[,"x"] > (x-1)*x_interval) & (stvea_object@codex_spatial[,"x"] <= x*x_interval)
}
for (y in 1:y_subsample)
{
  y_spatial_filters[[as.character(y)]] <- (stvea_object@codex_spatial[,"y"] > (y-1)*y_interval) & (stvea_object@codex_spatial[,"y"] <= y*y_interval)
}
# Generate subsampled region names
subsample_names <- vector(mode = "character" , length = 0)
for (region in unique(stvea_object@codex_region))
{
  for (x in 1:x_subsample)
  {
    temp_names <- paste(region, x, 1:y_subsample, sep = "_")
    subsample_names <- c(subsample_names,temp_names)
  }
}
# Create vector to store size (number of cells) and number of unique cell types in each subsample
subsample_sizes <- vector(mode = "numeric", length = length(subsample_names))
names(subsample_sizes) <- subsample_names
subsample_cell_types <- vector(mode = "numeric", length = length(subsample_names))
names(subsample_cell_types) <- subsample_names


## Subsample Assignment
# Create subsample assignment variable
assigned_subsamples <- vector(mode = "character", length = length(stvea_object@codex_region))
names(assigned_subsamples) <- paste(stvea_object@codex_region,stvea_object@codex_cell_ID, sep = "_")
# Assign subsamples
for (subsample in subsample_names)
{
  # Extract region name, x position, and y position from subsample name
  region <- strsplit(subsample,"_")[[1]][1]
  x_position <- strsplit(subsample,"_")[[1]][2]
  y_position <- strsplit(subsample,"_")[[1]][3]
  # Identify cells in region and within dimensions of subsample 
  subsample_filter <- x_spatial_filters[[x_position]] & y_spatial_filters[[y_position]] & region_filters[[region]]
  # Assign subsample to appropriate cells  
  assigned_subsamples[subsample_filter] <- subsample
  # Store size (number of cells) and number of unique cell types in the subsample
  subsample_sizes[subsample] <- sum(subsample_filter)
  subsample_cell_types[subsample] <- length(unique(stvea_object@codex_clusters[subsample_filter]))
}
rm(region,x_position,y_position,subsample_filter)
# Cull subsamples by cell number and cell type requirements (usually, set cell number requirement to desired k for interaction analysis)
subsample_names <- subsample_names[subsample_sizes >= 20 & subsample_cell_types >= 2]


## Subsampled Matrix Plots
pt_size <- 5
# Generate subsampled matrix plots
for (subsample in subsample_names)
{
  # Create temporary STvEA object and cull to current subsample
  subsample_stvea_object <- stvea_object
  subsample_stvea_object@codex_spatial <- subsample_stvea_object@codex_spatial[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_clean <- subsample_stvea_object@codex_clean[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_protein <- subsample_stvea_object@codex_protein[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_clusters <- subsample_stvea_object@codex_clusters[assigned_subsamples == subsample]
  # Export spatial plot of subsample with stains of interest
  PlotExprCODEXspatial(subsample_stvea_object, type = "protein", name = c("COLI","COLIV"), high_color = "red", high_color2 = "green") + geom_point(size = pt_size) + labs(title = "", subtitle = "") + theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"))
  if (x_subsample == 1 & y_subsample == 1) # If not subsampling, export region name rather than subsample name
  {
    ggsave(paste0(strsplit(subsample, split = "_1_1")[[1]][1],".png"), width = 5, height = 5)
  } else
  {
    #ggsave(paste0(subsample,".png"), width = 5, height = 5)
    ggsave(paste0(subsample,".png"), type = "cairo-png", width = 5, height = 5, bg = "black")
  }
}
rm(subsample_stvea_object)


## Subsampled Matrix Plots (COLI Only)
pt_size <- 5
# Generate subsampled matrix plots
for (subsample in subsample_names)
{
  # Create temporary STvEA object and cull to current subsample
  subsample_stvea_object <- stvea_object
  subsample_stvea_object@codex_spatial <- subsample_stvea_object@codex_spatial[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_clean <- subsample_stvea_object@codex_clean[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_protein <- subsample_stvea_object@codex_protein[assigned_subsamples == subsample,]
  subsample_stvea_object@codex_clusters <- subsample_stvea_object@codex_clusters[assigned_subsamples == subsample]
  # Export spatial plot of subsample with stains of interest
  PlotExprCODEXspatial(subsample_stvea_object, type = "protein", name = "COLI", high_color = "#0a2363", low_color = "#DCDEEA") + geom_point(size = pt_size) + labs(title = "", subtitle = "") 
  if (x_subsample == 1 & y_subsample == 1) # If not subsampling, export region name rather than subsample name
  {
    ggsave(paste0(strsplit(subsample, split = "_1_1")[[1]][1],".png"), width = 5, height = 5, bg = "white")
  } else
  {
    ggsave(paste0(subsample,".png"), type = "cairo-png", width = 5, height = 5, bg = "white")
  }
}
rm(subsample_stvea_object)


## Spatial Analysis - Subsample Level
# Set k parameter for interaction calculations
interaction_k <- 20
# Generate interaction scores and network for each region 
subsample_interaction_scores <- list()
for (subsample in subsample_names)
{
  # Limit spatial and cluster information to current subsample
  spatial_info <- stvea_object@codex_spatial[assigned_subsamples == subsample,]
  cluster_info <- stvea_object@codex_clusters[assigned_subsamples == subsample]
  # Generate adjacency matrix based on cell locations
  adjacency_matrix <- knn_graph(spatial_info, k=interaction_k)
  # Generate "feature" matrix of assigned cell types
  cell_matrix <- t(sapply(cell_types, function(x) (cluster_info==x)*1))
  row.names(cell_matrix) <- cell_types
  # Calculate interaction scores with hypergeometric null distribution for binary data 
  subsample_interaction_scores[[subsample]] <- CalculateInteractionScores(adjacency_matrix = adjacency_matrix, cell_matrix = cell_matrix, cell_pairs = cell_pairs, calculate_individual_scores = FALSE, normalize = TRUE)
  # Add interactions as row names
  rownames(subsample_interaction_scores[[subsample]]) <- paste(subsample_interaction_scores[[subsample]][,"f1"],subsample_interaction_scores[[subsample]][,"f2"], sep = " - ")
}
rm(adjacency_matrix,cell_matrix)


## Mapping Col I Ultrastructure to Trichrome Ultrastructure
# Prior to this point, run the full ultrastructural analysis pipeline on the
# exported collagen I matrix plots and save the cell_data_set object as 
# test_us_object for mapping. Also import us_object from trichrome analysis.
# Generate nn index in original reference dataset and save transform model
nn_object <- us_object
nn_object <- reduce_dimension(cds = nn_object, reduction_method = "UMAP", preprocess_method="PCA", build_nn_index = TRUE) 
us_object@reduce_dim_aux$UMAP$nn_index <- nn_object@reduce_dim_aux$UMAP$nn_index
us_object@reduce_dim_aux$UMAP$model <- nn_object@reduce_dim_aux$UMAP$model
rm(nn_object)
save_transform_models(us_object, 'reference_model')
# Load reference transform model into query dataset
test_us_object <- load_transform_models(test_us_object, "reference_model")
# Run PCA on query dataset
test_us_object <- preprocess_cds(cds = test_us_object, method = "PCA", num_dim=5)
# Apply original transform model to query dataset
test_us_object <- preprocess_transform(test_us_object)
test_us_object <- reduce_dimension_transform(test_us_object)
# Map pseudotime values to new cells using nearest neighbors
nn_indices <- as.vector(nn2(data = reducedDim(us_object, "UMAP"), query = reducedDim(test_us_object, "UMAP"), k = 1)$nn.idx)
test_us_object$mapped_pseudotime <- us_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][as.numeric(nn_indices)]
names(colData(test_us_object)$mapped_pseudotime) <- colnames(loaded_counts)
# Plot query cells on manifold
plot_cells(cds = test_us_object, color_cells_by = "mapped_pseudotime", cell_size = 1.5, show_trajectory_graph = TRUE, trajectory_graph_segment_size = 1.5, label_branch_points = FALSE, graph_label_size = 4, group_label_size = 0)


## ECM - Cell Interaction Correlations
# Extract cell interaction scores for subsamples
subsample_interaction_overview <- matrix(nrow = length(subsample_names), ncol = nrow(subsample_interaction_scores[[1]]))
rownames(subsample_interaction_overview) <- subsample_names
colnames(subsample_interaction_overview) <- rownames(subsample_interaction_scores[[1]])
for (subsample in rownames(subsample_interaction_overview))
{
  for (interaction in colnames(subsample_interaction_overview))
  {
    subsample_interaction_overview[subsample,interaction] <- subsample_interaction_scores[[subsample]][interaction,"score"]
  }
}
# Calculate correlations between ECM features (incl. mapped pseudotime) and individual cell interactions
us_int_correlations <- matrix(nrow = nrow(test_us_object)+1, ncol = ncol(subsample_interaction_overview))
rownames(us_int_correlations) <- c("Mapped.Pseudotime",rownames(test_us_object))
colnames(us_int_correlations) <- colnames(subsample_interaction_overview)
for (interaction in colnames(us_int_correlations))
{
  pearson <- cor(x = as.numeric(subsample_interaction_overview[rownames(subsample_interaction_overview),interaction]), y = as.numeric(colData(test_us_object)$mapped_pseudotime[rownames(subsample_interaction_overview)]), use = "na.or.complete")
  if (!is.null(pearson))
  {
    us_int_correlations["Mapped.Pseudotime",interaction] <- pearson
  }
  for (ecm_feature in rownames(us_int_correlations)[rownames(us_int_correlations) != "Mapped.Pseudotime"])
  {
    pearson <- cor(x = as.numeric(subsample_interaction_overview[rownames(subsample_interaction_overview),interaction]), y = as.numeric(test_us_object@assays@data$counts[ecm_feature,rownames(subsample_interaction_overview)]), use = "na.or.complete")
    if (!is.null(pearson))
    {
      us_int_correlations[ecm_feature,interaction] <- pearson
    }
  }
}
# Calculate p-values for correlations
us_int_pvalues <- matrix(nrow = nrow(test_us_object)+1, ncol = ncol(subsample_interaction_overview))
rownames(us_int_pvalues) <- c("Mapped.Pseudotime",rownames(test_us_object))
colnames(us_int_pvalues) <- colnames(subsample_interaction_overview)
for (interaction in colnames(us_int_pvalues))
{
  us_int_pvalues["Mapped.Pseudotime",interaction] <- cor.test(x = as.numeric(subsample_interaction_overview[rownames(subsample_interaction_overview),interaction]), y = as.numeric(colData(test_us_object)$mapped_pseudotime[rownames(subsample_interaction_overview)]), alternative = "two.sided", method = "pearson")$p.value
  for (ecm_feature in rownames(us_int_pvalues)[rownames(us_int_pvalues) != "Mapped.Pseudotime"])
  {
    us_int_pvalues[ecm_feature,interaction] <- cor.test(x = as.numeric(subsample_interaction_overview[rownames(subsample_interaction_overview),interaction]), y = as.numeric(test_us_object@assays@data$counts[ecm_feature,rownames(subsample_interaction_overview)]), alternative = "two.sided", method = "pearson")$p.value
  }
}
for (ecm_feature in rownames(us_int_pvalues))
{
  us_int_pvalues[ecm_feature,] <- p.adjust(us_int_pvalues[ecm_feature,], method = "BH")
}
 
 
## ECM Pseudotime-Individual US Parameter Correlations
# Calculate correlations between mapped pseudotime and individual ultrastructure parameters
us_pseudotime_correlations <- matrix(nrow = nrow(test_us_object), ncol = 1)
rownames(us_pseudotime_correlations) <- rownames(test_us_object)
colnames(us_pseudotime_correlations) <- c("Mapped.Pseudotime")
for (ecm_feature in rownames(us_pseudotime_correlations))
{
  pearson <- cor(x = as.numeric(test_us_object@assays@data$counts[ecm_feature,colnames(test_us_object@assays@data$counts)]), y = as.numeric(colData(test_us_object)$mapped_pseudotime[colnames(test_us_object@assays@data$counts)]), use = "na.or.complete")
  if (!is.null(pearson))
  {
    us_pseudotime_correlations[ecm_feature,"Mapped.Pseudotime"] <- pearson
  }
}
us_pseudotime_correlations <- as.matrix(sort(us_pseudotime_correlations[,"Mapped.Pseudotime"]))
colnames(us_pseudotime_correlations) <- "Mapped.Pseudotime"
# Calculate p-values for correlations
us_pseudotime_pvalues <- matrix(nrow = nrow(test_us_object), ncol = 1)
rownames(us_pseudotime_pvalues) <- rownames(test_us_object)
colnames(us_pseudotime_pvalues) <- c("Mapped.Pseudotime")
for (ecm_feature in rownames(us_pseudotime_pvalues))
{
  us_pseudotime_pvalues[ecm_feature,"Mapped.Pseudotime"] <- cor.test(x = as.numeric(test_us_object@assays@data$counts[ecm_feature,colnames(test_us_object@assays@data$counts)]), y = as.numeric(colData(test_us_object)$mapped_pseudotime[colnames(test_us_object@assays@data$counts)]), alternative = "two.sided", method = "pearson")$p.value
}
us_pseudotime_pvalues[,"Mapped.Pseudotime"] <- p.adjust(us_pseudotime_pvalues[,"Mapped.Pseudotime"], method = "BH")


## Loading Plots of Interactions and US Parameters Correlated w/ Pseudotime
# Generate dataframes for plotting
us_int_df <- data.frame(row.names = names(sort(us_int_correlations["Mapped.Pseudotime",])))
us_int_df$Mapped.Pseudotime <- sort(us_int_correlations["Mapped.Pseudotime",]) 
us_int_df$interactions <- rownames(us_int_df)
us_pseudotime_df <- data.frame(us_pseudotime_correlations)
us_pseudotime_df$parameters <- rownames(us_pseudotime_df)
# Plot loadings of top n positively and negatively correlated interactions w/ pseudotime (log10 scale)
n <- 5
ggplot(data = us_int_df[c(1:n,(nrow(us_int_df)-n+1):nrow(us_int_df)),]) + geom_segment(aes(x=0, y=0, xend=sign(.data[["Mapped.Pseudotime"]]), yend=abs(.data[["Mapped.Pseudotime"]]), colour = abs(.data[["Mapped.Pseudotime"]])), arrow = arrow(length = unit(1, "picas"))) + annotate("text", x = sign(us_int_df[c(1:n,(nrow(us_int_df)-n+1):nrow(us_int_df)),]$Mapped.Pseudotime), y = (abs(us_int_df[c(1:n,(nrow(us_int_df)-n+1):nrow(us_int_df)),]$Mapped.Pseudotime)), label = us_int_df[c(1:n,(nrow(us_int_df)-n+1):nrow(us_int_df)),]$interactions) + theme_bw() + scale_y_continuous(trans = "log10") + scale_color_viridis(option = "plasma") + guides(colour=guide_colorbar(ticks.colour = NA, title = "Strength (Abs(Pearson))")) + scale_x_continuous(limits = c(-2,2))
ggsave(paste0("Top ",n," Interactions Pseudotime Loadings Plot.png"),width=7,height=5)
# Plot top n ultrastructure parameters positively and negatively correlated with mapped pseudotime (log10 scale)
n <- 5
ggplot(data = us_pseudotime_df[c(1:n,(nrow(us_pseudotime_df)-n+1):nrow(us_pseudotime_df)),]) + geom_segment(aes(x=0, y=0, xend=sign(.data[["Mapped.Pseudotime"]]), yend=abs(.data[["Mapped.Pseudotime"]]), colour = abs(.data[["Mapped.Pseudotime"]])), arrow = arrow(length = unit(1, "picas"))) + annotate("text", x = sign(us_pseudotime_df[c(1:n,(nrow(us_pseudotime_df)-n+1):nrow(us_pseudotime_df)),]$Mapped.Pseudotime), y = (abs(us_pseudotime_df[c(1:n,(nrow(us_pseudotime_df)-n+1):nrow(us_pseudotime_df)),]$Mapped.Pseudotime)), label = us_pseudotime_df[c(1:n,(nrow(us_pseudotime_df)-n+1):nrow(us_pseudotime_df)),]$parameters) + theme_bw() + scale_y_continuous(trans = "log10") + scale_color_viridis(option = "plasma") + guides(colour=guide_colorbar(ticks.colour = NA, title = "Strength (Abs(Pearson))")) + scale_x_continuous(limits = c(-2,2))
ggsave(paste0("Top ",n," US Parameters Pseudotime Loadings Plot.png"),width=7,height=5)


## Anchor Transfer for Blinded Testing Dataset
# Prior to this point, initialize and batch correct the blinded testing 
# dataset (stvea_test) as performed previously for stvea_object.
# Generate counts matrix for testing dataset
test_matrix <- stvea_test@codex_clean[,!grepl("Blank",colnames(stvea_test@codex_clean))]
rownames(test_matrix) <- paste(stvea_test@codex_region, stvea_test@codex_cell_ID, sep = "_")
test_matrix <- t(test_matrix)
# Create Seurat object and run minimum required preprocessing
seurat_test <- CreateSeuratObject(counts = test_matrix)
rm(test_matrix)
seurat_test <- ScaleData(seurat_test, features = rownames(seurat_test))
seurat_test <- RunPCA(seurat_test, features = rownames(seurat_test), npcs = 5)
# Find transfer anchors by projecting testing dataset on PCA space of training dataset
anchors <- FindTransferAnchors(reference = seurat_object, query = seurat_test, features = rownames(seurat_test), reference.reduction = "pca", dims = 1:5)
# Classify testing cells by PCA-based transfer 
pca_classifications <- TransferData(anchorset = anchors, refdata = seurat_object$annotated_cell_types, dims = 1:5)
# Add cell type classifications to testing data
seurat_test$annotated_cell_types <- pca_classifications$predicted.id
stvea_test@codex_clusters <- as.factor(seurat_test$annotated_cell_types)

