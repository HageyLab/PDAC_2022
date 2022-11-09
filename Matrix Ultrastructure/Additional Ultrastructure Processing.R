## Initialization
library(dplyr)
library(patchwork)
library(DDRTree)
library(ggplot2)
library(monocle3)
library(gridExtra)


## Loading in MATLAB-Processed Ultrastructure Data
# Load data (.csv should have one row per image and one column per matrix parameter, plus headers)
loaded_counts <- read.delim(file = "pdac.csv", header = TRUE, row.names=1, sep = ",")
loaded_counts <- t(loaded_counts)
us_object <- new_cell_data_set(loaded_counts)


## Processing of Patient Metadata
# Load in patient metadata
loaded_metadata <- read.delim(file = "patient_metadata.csv", header = TRUE, row.names=1, sep = ",")
loaded_metadata <- t(loaded_metadata)
# Load image to patient key (.csv should have one row per image and a column with patient ID, plus headers)
loadedkey <- read.delim(file = "image_to_patient.csv", header = TRUE, row.names=1, sep = ",")
loadedkey <- t(loadedkey)
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
# Create defined survivorship category for classification
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
# Initialize patient metadata categories
colData(us_object)$patient_ID <- NA
colData(us_object)$cohort <- NA
colData(us_object)$age <- NA
colData(us_object)$sex <- NA
colData(us_object)$alive_or_dead <- NA
colData(us_object)$overall_survival <- NA
colData(us_object)$disease_free_survival <- NA
colData(us_object)$diagnosis <- NA
colData(us_object)$t <- NA
colData(us_object)$n <- NA
colData(us_object)$m <- NA
colData(us_object)$ajcc_stage <- NA
colData(us_object)$grade <- NA
colData(us_object)$tumor_size <- NA
colData(us_object)$lymph_nodes <- NA
colData(us_object)$preop_ca199 <- NA
colData(us_object)$neoadjuvant_chemo <- NA
colData(us_object)$include <- NA
colData(us_object)$survivorship_category <- NA
names(colData(us_object)$patient_ID) <- colnames(loaded_counts)
names(colData(us_object)$cohort) <- colnames(loaded_counts)
names(colData(us_object)$age) <- colnames(loaded_counts)
names(colData(us_object)$sex) <- colnames(loaded_counts)
names(colData(us_object)$alive_or_dead) <- colnames(loaded_counts)
names(colData(us_object)$overall_survival) <- colnames(loaded_counts)
names(colData(us_object)$disease_free_survival) <- colnames(loaded_counts)
names(colData(us_object)$diagnosis) <- colnames(loaded_counts)
names(colData(us_object)$t) <- colnames(loaded_counts)
names(colData(us_object)$n) <- colnames(loaded_counts)
names(colData(us_object)$m) <- colnames(loaded_counts)
names(colData(us_object)$ajcc_stage) <- colnames(loaded_counts)
names(colData(us_object)$grade) <- colnames(loaded_counts)
names(colData(us_object)$tumor_size) <- colnames(loaded_counts)
names(colData(us_object)$lymph_nodes) <- colnames(loaded_counts)
names(colData(us_object)$preop_ca199) <- colnames(loaded_counts)
names(colData(us_object)$neoadjuvant_chemo) <- colnames(loaded_counts)
names(colData(us_object)$include) <- colnames(loaded_counts)
names(colData(us_object)$survivorship_category) <- colnames(loaded_counts)
# Populate patient metadata based on patient ID lookup
colData(us_object)$patient_ID <- as.character(loaded_metadata['Patient.ID',as.character(loadedkey['Patient.ID',])])
colData(us_object)$cohort <- as.character(loaded_metadata['Cohort',as.character(loadedkey['Patient.ID',])])
colData(us_object)$age <- as.numeric(loaded_metadata['Age',as.character(loadedkey['Patient.ID',])])
colData(us_object)$sex <- as.character(loaded_metadata['Sex',as.character(loadedkey['Patient.ID',])])
colData(us_object)$alive_or_dead <- as.character(loaded_metadata['Alive.or.Dead',as.character(loadedkey['Patient.ID',])])
colData(us_object)$overall_survival <- as.numeric(loaded_metadata['OS',as.character(loadedkey['Patient.ID',])])
colData(us_object)$disease_free_survival <- as.numeric(loaded_metadata['Processed.DFS',as.character(loadedkey['Patient.ID',])])
colData(us_object)$diagnosis <- as.character(loaded_metadata['Diagnosis',as.character(loadedkey['Patient.ID',])])
colData(us_object)$t <- as.character(loaded_metadata['T',as.character(loadedkey['Patient.ID',])])
colData(us_object)$n <- as.character(loaded_metadata['N',as.character(loadedkey['Patient.ID',])])
colData(us_object)$m <- as.character(loaded_metadata['M',as.character(loadedkey['Patient.ID',])])
colData(us_object)$ajcc_stage <- as.character(loaded_metadata['AJCC.Stage',as.character(loadedkey['Patient.ID',])])
colData(us_object)$grade <- as.character(loaded_metadata['Grade',as.character(loadedkey['Patient.ID',])])
colData(us_object)$tumor_size <- as.numeric(loaded_metadata['Tumor.size',as.character(loadedkey['Patient.ID',])])
colData(us_object)$lymph_nodes <- as.numeric(loaded_metadata['Lymph.nodes',as.character(loadedkey['Patient.ID',])])
colData(us_object)$preop_ca199 <- as.numeric(loaded_metadata['Preop.CA199',as.character(loadedkey['Patient.ID',])])
colData(us_object)$neoadjuvant_chemo <- as.character(loaded_metadata['Neoadjuvant.Chemo',as.character(loadedkey['Patient.ID',])])
colData(us_object)$include <- as.logical(as.numeric(loaded_metadata['Include',as.character(loadedkey['Patient.ID',])]))
colData(us_object)$survivorship_category <- as.character(loaded_metadata['Survivorship.Category',as.character(loadedkey['Patient.ID',])])
# Store ultrastructural parameter names 
rowData(us_object)$gene_name <- rownames(us_object)
rowData(us_object)$gene_short_name <- rowData(us_object)$gene_name


## UMAP and Trajectory Analysis (Uwot, DDRTree via Monocle Wrappers)
# Pre-UMAP PCA
us_object <- preprocess_cds(cds = us_object, method = "PCA", num_dim=5)
# Run/plot UMAP (Uwot via Monocle)
us_object <- reduce_dimension(cds = us_object, reduction_method = "UMAP", preprocess_method="PCA") 
us_object <- cluster_cells(cds = us_object, reduction_method = "UMAP", resolution = 5e-5) 
plot_cells(us_object, cell_size = 1)
# Build trajectories (DDRTree via Monocle)
us_object <- learn_graph(us_object, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 30))
# Order cells in pseudotime 
us_object <- order_cells(us_object, reduction_method = "UMAP")
# Plot pseudotime trajectories
plot_cells(cds = us_object, color_cells_by = "pseudotime", cell_size = 1.5, show_trajectory_graph = TRUE, trajectory_graph_segment_size = 1.5, label_branch_points = FALSE, graph_label_size = 4, group_label_size = 0)


## Storing Processed Metrics
# Store pseudotime for each datapoint
pseudotime <- us_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
# Calculate summary metrics
temp_coords <- matrix(ncol=length(loaded_metadata["Patient.ID",]), nrow = 3)
colnames(temp_coords) <- loaded_metadata["Patient.ID",]
rownames(temp_coords) <- c("Median.Pseudotime","Dispersion.Index","Bin")
colnames(loaded_metadata) <- loaded_metadata["Patient.ID",]
loaded_metadata <- rbind(loaded_metadata,temp_coords)
for (namedpatient in loaded_metadata["Patient.ID",])
{
  loaded_metadata["Median.Pseudotime",namedpatient] <- median(pseudotime[colData(us_object)$patient_ID == namedpatient])
  if (is.na(loaded_metadata["Median.Pseudotime",namedpatient]))
  {
    loaded_metadata["Bin",namedpatient] <- NA
  } else if (as.numeric(loaded_metadata["Median.Pseudotime",namedpatient]) < 24)
  {
    loaded_metadata["Bin",namedpatient] <- "L"
  } else if (as.numeric(loaded_metadata["Median.Pseudotime",namedpatient]) >= 24)
  {
    loaded_metadata["Bin",namedpatient] <- "R"
  }
  ecm_variance <- vector(mode = "numeric", length = 147)
  names(ecm_variance) <- rownames(loaded_counts)
  for (parameter in names(ecm_variance))
  {
    ecm_variance[parameter] <- var(rescale(loaded_counts[parameter,colData(us_object)$patient_ID == namedpatient]), na.rm = TRUE)
  }
  loaded_metadata["Dispersion.Index",namedpatient] <- sum(ecm_variance, na.rm = TRUE)
}

