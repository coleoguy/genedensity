
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: combines the files in the results folder into a 
## single .csv file

library(data.table)
source("functions.R")
# directory containing individual result files
csv.dir.path <- "../results/vertebrates/individual_species_results"
csv.file.paths <- paste0(csv.dir.path, "/", list.files(csv.dir.path))
results.list <- lapply(csv.file.paths, fread)
results.table <- na.omit(rbindlist(results.list, fill = TRUE))
fwrite(results.table, "../results/vertebrates/combined_results.csv", row.names = FALSE)










