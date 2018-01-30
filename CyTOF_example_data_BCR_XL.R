##########################################################################################
# Script to load CyTOF example data set: BCR-XL
# 
# The original 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and has
# previously been used for benchmark evaluations by Bruggner et al. (2014) and Nowicka et
# al. (2017).
# 
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section 'PBMC Example 1'):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
# 
# Cell population labels are reproduced from Nowicka et al. (2017), where they were
# generated using a strategy of expert-guided manual merging of automatically generated
# clusters from the FlowSOM algorithm. Code to reproduce the cell population labels is
# available in the script 'cell_population_labels_BCR_XL.R'.
# 
# Lukas Weber, January 2018
##########################################################################################


library(flowCore)
library(SummarizedExperiment)
library(Rtsne)
library(ggplot2)


DIR_RAW_DATA <- "../raw_data"
DIR_LABELS <- "../population_IDs"
DIR_DATA_OUT <- "../data"
DIR_PLOTS <- "../plots"



# -------------
# Download data
# -------------

url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
fcs_filename <- "PBMC8_fcs_files.zip"

download.file(file.path(url, fcs_filename), destfile = file.path(DIR_RAW_DATA, fcs_filename))
unzip(file.path(DIR_RAW_DATA, fcs_filename), exdir = DIR_RAW_DATA)



# ---------
# Filenames
# ---------

# .fcs files
files <- list.files(DIR_RAW_DATA, pattern = "\\.fcs$", full.names = TRUE)
files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

files_all <- c(files_BCRXL, files_ref)

# cell population labels
files_labels <- list.files(DIR_LABELS, pattern = "\\.csv$", full.names = TRUE)
files_labels_BCRXL <- files_labels[grep("patient[1-8]_BCR-XL\\.csv$", files_labels)]
files_labels_ref <- files_labels[grep("patient[1-8]_Reference\\.csv$", files_labels)]

files_labels_all <- c(files_labels_BCRXL, files_labels_ref)



# ---------
# Load data
# ---------

data <- lapply(files_all, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# sample IDs
sample_IDs <- gsub("\\.fcs$", "", gsub("^PBMC8_30min_", "", basename(files_all)))
sample_IDs

names(data) <- sample_IDs

# group IDs
group_IDs <- factor(gsub("^patient[0-9+]_", "", sample_IDs), levels = c("Reference", "BCR-XL"))
group_IDs

# patient IDs
patient_IDs <- factor(gsub("_.*$", "", sample_IDs))
patient_IDs

# cell population labels
labels <- lapply(files_labels_all, read.csv)



# --------------------------
# Protein marker information
# --------------------------

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

all(sapply(seq_along(data), function(i) all(colnames(data[[i]]) == colnames(data[[1]]))))

marker_names <- colnames(data[[1]])
marker_names <- gsub("\\(.*$", "", marker_names)

is_marker <- is_celltype_marker <- is_state_marker <- rep(FALSE, length(marker_names))

is_marker[cols_markers] <- TRUE
is_celltype_marker[cols_lineage] <- TRUE
is_state_marker[cols_func] <- TRUE

marker_info <- data.frame(marker_names, is_marker, is_celltype_marker, is_state_marker)
marker_info



# --------------
# Transform data
# --------------

# apply 'arcsinh' transform with 'cofactor' = 5; this is a standard transformation
# commonly applied to mass cytometry data (see Bendall et al., 2011, Supp. Fig. S2)

data <- lapply(data, function(d) {
  asinh(d / 5)
})



# ---------------------------
# Create SummarizedExperiment
# ---------------------------

# number of cells
n_cells <- sapply(data, nrow)
n_cells
sum(n_cells)

# assay data: expression values
assay_data <- do.call("rbind", data)
stopifnot(nrow(assay_data) == sum(n_cells))

colnames(assay_data) <- unname(gsub("\\(.*", "", colnames(assay_data)))

# row data: sample information and cluster labels
stopifnot(length(n_cells) == length(sample_IDs), 
          length(n_cells) == length(group_IDs), 
          length(n_cells) == length(patient_IDs))

stopifnot(nrow(do.call("rbind", labels)) == sum(n_cells))

row_data <- data.frame(
  sample_IDs = rep(sample_IDs, n_cells), 
  group_IDs = rep(group_IDs, n_cells), 
  patient_IDs = rep(patient_IDs, n_cells), 
  population = do.call("rbind", labels)$population
)

# column data: marker information
stopifnot(ncol(assay_data) == nrow(marker_info))

col_data <- marker_info

rownames(col_data) <- marker_info$marker_names

# create SummarizedExperiment
d_se <- SummarizedExperiment(
  assay_data, 
  rowData = row_data, 
  colData = col_data
)



# -------------------------------------------------
# Check cluster / population labels with t-SNE plot
# -------------------------------------------------

# run t-SNE

# note: using cell type markers only
d_tsne <- assay(d_se)[, colData(d_se)$is_celltype_marker]
d_tsne <- as.matrix(d_tsne)

# remove any duplicate rows (required by Rtsne)
dups <- duplicated(d_tsne)

# subsample cells for t-SNE
n_sub <- 5000
set.seed(123)
# note: removing duplicated rows here
ix <- sample(seq_len(nrow(d_se))[!dups], n_sub)

d_tsne <- d_tsne[ix, ]
labels_tsne <- rowData(d_se)$population[ix]

stopifnot(nrow(d_tsne) == length(labels_tsne))

# run Rtsne
# (note: initial PCA step not required, since we do not have too many dimensions)
set.seed(123)
out_tsne <- Rtsne(d_tsne, pca = FALSE, verbose = TRUE)

tsne_coords <- as.data.frame(out_tsne$Y)
colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")


# create t-SNE plot

stopifnot(nrow(tsne_coords) == length(labels_tsne))

d_plot <- data.frame(
  tsne_coords, 
  population = labels_tsne
)

ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, color = population)) + 
  geom_point() + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("t-SNE plot: BCR-XL data set with manually merged populations") + 
  theme_bw() + 
  theme(aspect.ratio = 1)

# save plot
fn <- file.path(DIR_PLOTS, "BCR_XL_tSNE_plot.pdf")
ggsave(file = fn, width = 7, height = 7)



