# SLICE
## Code Files:
### Functions.Rmd
Used to load packages and functions relevant for SLICE.

### SLICE - Data Creation.Rmd
Used to create the data for SLICE.

#### Requirments:
- Run 'Functions.Rmd' at the start of each session.

### SLICE - Data Loader.Rmd
Used to load the data created by the previous script into the R environment.

#### Requirments:
- Run 'Functions.Rmd' at the start of each session.
- Run 'SLICE - Data Creation.Rmd' at least once.

### SLICE - Analysing Results.Rmd
A script used to analyse the results of SLICE and provide further insight on the results.

#### Requirments:
- Run 'Functions.Rmd' and 'SLICE - Data Loader.Rmd' at the start of each session.
- Run 'SLICE - Data Creation.Rmd' at least once.

### shinySLICE.R
The shiny R application used for the exploration of the data and SLICE's results.

#### Requirments:
- Run 'Functions.Rmd' and 'SLICE - Data Loader.Rmd' at the start of each session.
- Run 'SLICE - Data Creation.Rmd' at least once.

## Things that need to be done:
### In the actual analysis
- Fix the ssGSEA function (ideally replace the implementation of the current functiion with a call to an existing function)
- Find some dataset to calculate driver mutation frequencies in pan-cancer

#### Potentially:
- Differential gene expression analysis between RB1-WT and RB1-MUT bladder cells.
- Survival analysis between RB1-WT and RB1-MUT bladder cells (currently using TCGA we got insignificant results).

### In the figures
- Changing the volcano visualization for the paper (currently using plotly which is interactive but ugly)
- Changing the network visualization for the paper (currently using plotly which is interactive but ugly)

### In the shiny application
- Fix the bug regarding the coloring of the network's edges (happens only when there are multiple "hubs")
- Improve on the BFS algorithm in the shiny app so it would use weights to find the best shortest path instead of just the first shortest path (take into consideration that it will make it run exponentially slower)
- Add a checkbox to filter all driver genes that are not mutated in over 15% of the patients (the frequencies are calculated in the analysis part, it probably would be best to first calculate them in the "Data Creation" script, save them to a file and then load them in the "Data Loader" because it takes a loooong time)

