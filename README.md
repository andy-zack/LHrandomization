# Lawley-Hotelling Re-Randomization

This R package contains a function to rerandomize a dataset for conducting an experiment. It randomizes a dataset, n times, and selects the most balanced randomization using the Lawley-Hotelling statistic, checking balance on the specified covariates.

To install this package, run: 

```
# if devtools is not already installed, run:
# install.packages("devtools")

devtools::install_github("andy-zack/LHrandomization")
library(LHrandomization)
```

Then run `?balance_randomization` to view documentation.
