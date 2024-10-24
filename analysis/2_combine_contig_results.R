


# Zhaobo Hu
# zhaobohu2002@gmail.com

# description: combines the files in the results folder into a single .csv file

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

# load package
library(data.table)
# source function
source("functions.R")
# directory containing individual result files
csv.dir.path <- paste0("../results/", vert.invert, "/individual_species_results")
# paths of each individual result file
csv.file.paths <- paste0(csv.dir.path, "/", list.files(csv.dir.path))
# reads the results of every available species into a single list
results.list <- lapply(csv.file.paths, fread)
# converts the previous list into a table
results.table <- rbindlist(results.list, fill = TRUE)
# write csv
fwrite(results.table, paste0("../results/", vert.invert, "/combined_results.csv"), row.names = FALSE)










