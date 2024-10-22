
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: combines the files in the repeat landscape results folder
## into a single .csv file

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"

library(data.table)
source("functions.R")
# directory containing individual result files
csv.dir.path <- paste0("../results/", vert.invert, "/repeat_landscape_plotting")
csv.file.paths <- paste0(csv.dir.path, "/", list.files(csv.dir.path))
results.list <- lapply(csv.file.paths, fread)
results.table <- na.omit(rbindlist(results.list, fill = TRUE))
fwrite(results.table, paste0("../results/", vert.invert, "/all_repeat_landscapes.csv"), row.names = FALSE)










