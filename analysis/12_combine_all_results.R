# zhaobo hu
# zhaobohu2002@gmail.com

# Description: combines all results into a single csv file for plotting

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"

chromnums <- read.csv(paste0("../data/", vert.invert, "/chromnums.csv"))
assembly.sizes <- read.csv(paste0("../results/", vert.invert, "/assembly_sizes.csv"))
contig.stats <- read.csv(paste0("../results/", vert.invert, "/contig_stats.csv"))
clades <- read.csv(paste0("../results/", vert.invert, "/clades.csv"))
taxo.gnsz <- read.csv(paste0("../data/", vert.invert, "/taxo_gnsz.csv"))
rep.landscape.stats <- read.csv(paste0("../results/", vert.invert, "/rep_landscape_stats.csv"))

a <- merge(chromnums, assembly.sizes, by = "species", all.x = TRUE)
b <- merge(a, contig.stats, by = "species", all.x = TRUE)
c <- merge(b, clades, by = "species", all.x = TRUE)
d <- merge(c, taxo.gnsz, by = "species", all.x = TRUE)
final.results <- merge(d, rep.landscape.stats, by = "species", all.x = TRUE)
write.csv(final.results, paste0("../results/", vert.invert, "/final_results.csv"), row.names = FALSE)
