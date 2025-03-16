# SLAYER
## Code Files:
### Functions.Rmd
Used to load packages and functions relevant for SLAYER.

### SLAYER_Data_Creation.Rmd
Used to create the data for SLAYER.

#### Requirments:
- Run 'Functions.Rmd' at the start of each session.

### SLAYER_Data_Loader.Rmd
Used to load the data created by the previous script into the R environment.

#### Requirments:
- Run 'Functions.Rmd' at the start of each session.
- Run 'SLAYER_Data_Creation.Rmd' at least once.

### SLAYER_Analysing_Results.Rmd
A script used to analyse the results of SLAYER and provide further insight on the results.

#### Requirments:
- Run 'Functions.Rmd' and 'SLAYER_Data_Loader.Rmd' at the start of each session.
- Run 'SLAYER_Data_Creation.Rmd' at least once.

### Figures.Rmd
A script used to create figures from the paper.

### shinySLAYER.R
The shiny R application used for the exploration of the data and SLAYER's results.

#### Requirments:
- Run 'Functions.Rmd' and 'SLAYER_Data_Loader.Rmd' at the start of each session.
- Run 'Shiny_data_prep.R ' at least once.
