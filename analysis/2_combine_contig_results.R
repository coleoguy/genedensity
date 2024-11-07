


# Zhaobo Hu
# zhaobohu2002@gmail.com

# description: combines the files in the results folder into a single .csv file

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"

# load package
library(data.table)
# directory containing individual result files
dir <- paste0("../results/", vert.invert, "/individual_species_results")
# paths of each individual result file
files <- paste0(dir, "/", list.files(dir))
# reads the results of every available species into a single list
results.list <- lapply(files, fread)
# converts the previous list into a table
results.table <- rbindlist(results.list, fill = TRUE)
# write csv
fwrite(results.table, paste0("../results/", vert.invert, "/combined_results.csv"), row.names = FALSE)










