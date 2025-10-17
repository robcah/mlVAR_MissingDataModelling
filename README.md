# mlVAR Analysis of three datasets modelling missing data
## Code for replication
### Data will be accesible by request at discretion of the University of Manchester

Supporting materials for paper: Unsupervised identification of significant lineages of SARS-CoV-2 through scalable machine learning methods.

1. Install environment:\
`conda env create --file mlvar3datasets.yaml`

1. Install R and RStudio (recommended)

1. in RStudio run the file:
`BRC_mlVAR_Analysis_MultipleDatasets.R`

1. Activate environment and run jupyter notebook:\
`conda activate mlvar3datasets`\
`jupyter notebook`

1. Follow the instructions from jupyter notebook `mlVAR3Datasets_MissingDataModelling.ipynb` to produce the data with imputation and synthetic constructs modelling gaps.

1. Run R-code `BRC_mlVAR_Analysis_MultipleDatasets.R` either in R or RStudio. This produce two sets of mlVAR networks for the original and the imputed data, each set consist of networks: temporal, between-subjects and contemporaneous.

1. Follow the instructions from jupyter notebook `mlVAR3datasets_PsyNetworksPlot.ipynb` to plot results of both mlVAR set of networks.



