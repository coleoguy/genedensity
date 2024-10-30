
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name, family, and assembly size

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"
assembly.sizes <- read.csv(
  paste0("../results/", vert.invert, "/assembly_sizes.csv"))
taxo.gnsz <- read.csv(paste0("../data/", vert.invert, "/taxo_gnsz.csv"))
taxo.gnsz <- taxo.gnsz[, c("species", "family", "est.gnsz_bp")]
results <- merge(taxo.gnsz, assembly.sizes, by = "species", all.x = TRUE)
results$species <- tolower(gsub(" ", "_", results$species))
results <- na.omit(results)
write.table(results, 
            paste0("../data/", vert.invert, "/species_family_asmblysz.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE,
            quote = FALSE)
