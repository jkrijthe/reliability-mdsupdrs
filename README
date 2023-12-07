Code to reproduce the analysis in Evers LJW, Krijthe JH, Meinders MJ, Bloem BR, Heskes TM. Measuring Parkinson's disease over time: The real-world within-subject reliability of the MDS-UPDRS. Mov Disord. 2019 Oct;34(10):1480-1487. doi: 10.1002/mds.27790.

# Data
The code relies on the Parkinson Progression Markers Initiative (PPMI) data, available at ppmi-info.org. Since this data cannot be freely redistributed, to run the code, the data needs to be downloaded and placed in the raw-data folder. For the paper, we used the version that was current on 27 June 2017.

# Running the analysis
To preprocess the data run:
`source("R/load_data.R")`
This generates a data.RData file.
This can then be used to generate a bootstrap-data.RData file by running knitr on measurement-error.Rmd
Th bootstrap estimates can then be calculated (this is a computationally expensive step, you can skip this step by using the provided estimates in the dlm_estimates_ppmi*.Rdata files) using:
`source("ppmi-ssm-bootstrap.R")`
`source("ppmi-ssm-bootstrap-compare.R")`
`source("ppmi-ssm-bootstrap-factors.R")`
Finally, the analysis and figures can be created by using knitr on the measurement-error.Rmd and supplementary.Rmd files.

