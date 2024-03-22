# SLICE
## Code Files:
### SLICE - Data Creation.Rmd
Used to create the data for SLICE.

### shinySLICE - Data Loader.Rmd
Used to load the data created by the previous script into the R environment.

### shinySLICE.R
The shiny R application used for the exploration of the data and SLICE's results.

### SLICE - Analysing Results.Rmd
A script used to analyse the results of SLICE and provide further insight on the results.

## Things that need to be done:
### In the figures
- Changing the volcano visualiztion for the paper (currently using plotly which is interactive but ugly)
- Changing the network visualiztion for the paper (currently using plotly which is interactive but ugly)

### In the shiny application
- Fix the bug regarding the coloring of the network's edges (happens only when there are multiple "hubs")
- Improve on the BFS algorithm in the shiny app so it would use weights to find the best shortest path instead of just the first shortest path (take into consideration that it will make it run exponentially slower)
- Add a checkbox to filter all driver genes that are not mutated in over 15% of the patients (the frequencies are calculated in the analysis part, it probably would be best to first calculate them in the "Data Creation" script, save them to a file and then load them in the "Data Loader" because it takes a loooong time)

### In the actual analysis
- Fix the ssGSEA function (ideally replace the implementation of the current functiion with a call to an existing function)
- Find some dataset to calculate driver mutation frequencies in pan-cancer

#### Potentially:
- Differential gene expression analysis between RB1-WT and RB1-MUT bladder cells.
- Survival analysis between RB1-WT and RB1-MUT bladder cells.

