# CyTOF-example-data


Example mass cytometry (CyTOF) data set: BCR-XL data set

For information about the data set, see comments at the top of the script [CyTOF_example_data_BCR_XL.R](CyTOF_example_data_BCR_XL.R).

Instructions to download and prepare the data set in SummarizedExperiment format:

- Set up directories listed at the top of both scripts.
- Run script [cell_population_labels_BCR_XL.R](cell_population_labels_BCR_XL.R) to reproduce the manually merged cell population labels from Nowicka et al. (2017). This script runs part of the code from the [CyTOF workflow](http://bioconductor.org/help/workflows/cytofWorkflow/) by Nowicka et al. (2017).
- Run script [CyTOF_example_data_BCR_XL.R](CyTOF_example_data_BCR_XL.R) to create a SummarizedExperiment and t-SNE plot.
