

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

library(ggplot2)
parsed.results <- read.csv(paste0("../results/", vert.invert, "/parsed_results.csv"))
files <- list.files(paste0("../results/", vert.invert, "/repeat_landscape"))
for (file in files) {
  # read text file into lines
  divsum.vector <- readLines(paste0("../results/", vert.invert, "/repeat_landscape/", file))
  # look for the start of useful information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum.vector) + 1
  # condense the useful lines into a table
  divsum.vector <- divsum.vector[start.index:length(divsum.vector)]
  divsum.table <- read.table(textConnection(divsum.vector), sep = " ", header = TRUE)
  # drop NA columns
  divsum.table <- divsum.table[-c(which(sapply(divsum.table, 
                                              function(col) all(is.na(col)))))]
  # detect column names in the same repeat category
  repeat.cats <- sapply(colnames(divsum.table), 
                             function(names) strsplit(names, "\\.")[[1]][1])
  non.unique.cats <- names(table(repeat.cats)[table(repeat.cats) > 1])
  # get the indices of the non-unique colnames
  non.unique.indices <- sapply(non.unique.cats, 
                             function(str) which(as.character(repeat.cats) == str))
  # get the row sums of columns with non-unique names
  non.unique.sum <- lapply(non.unique.indices, 
                         function(index) rowSums(divsum.table[c(index)]))
  # remove non-unique columns from table
  results.table <- divsum.table[-c(as.numeric(unlist(non.unique.indices)))]
  # rename table
  colnames(results.table) <- repeat.cats[-c(as.numeric(unlist(non.unique.indices)))]
  # add summed columns to table
  results.table <- data.frame(results.table, non.unique.sum)
  # find indices of colnames of the largest repeat groups
  cats.sort.length <- sort(colSums(results.table), decreasing = TRUE)
  known.cats <- names(sort(colSums(results.table), decreasing = TRUE)) != "Unknown"
  top.4.cat.names <- names(head(cats.sort.length[known.cats], 4))
  top.4.cat.ind <- which(colnames(results.table) %in% top.4.cat.names)
  # find indices of columns with unknown and divergence
  unknown.ind <- which(colnames(results.table) %in% "Unknown")
  div.ind <- which(colnames(results.table) %in% "Div")
  # find and sum columns that don't belong to any aforementioned categories
  other <- rowSums(results.table[-c(top.4.cat.ind, div.ind, unknown.ind)])
  #remove "other" and replace with sums
  results.table <- results.table[c(top.4.cat.ind, div.ind, unknown.ind)]
  results.table <- data.frame(results.table, other)
  # grouped by columns
  divergence <- rep(results.table$Div, length(colnames(results.table))-1)
  sequence.lengths_bp <- as.numeric(unlist(results.table[-c(1)]))
  species.under <- gsub("_summary\\.divsum$", "", file)
  species.lower <- gsub("_", " ", species.under)
  species <- gsub("^(\\w)(.*)", "\\U\\1\\L\\2", species.lower, perl = TRUE)
  gnsz_bp <- unique(parsed.results$asmbly.size_bp[parsed.results$species == species])
  percent.of.genome <- sequence.lengths_bp / gnsz_bp * 100
  repeat.group <- as.character(
    sapply(colnames(results.table)[colnames(results.table) != "Div"], 
           function(rowname) rep(rowname, length(results.table$Div))))
  final <- data.frame(divergence, sequence.lengths_bp, percent.of.genome, repeat.group)
  final$repeat.group <- factor(final$repeat.group, 
                              levels = c(top.4.cat.names, "Other", "Unknown"))
  rsq <- signif(unique(parsed.results$species.rsquared[parsed.results$species == species]), 3)
  
  ggplot(final, aes(y = percent.of.genome, x = divergence, fill = repeat.group)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle(bquote("Repeat Landscape of"~italic(.(species))~~(italic(r)^2==.(rsq)))) +
    theme(plot.title = element_text(hjust = 0.45),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "#f2f2f2", 
                                           color = "black", 
                                           linewidth = 0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "black", 
                                          linetype = "dotted", 
                                          size = 0.25),
          legend.position = c(0.884, 0.557),
          legend.key.size = unit(23.5, "points")) +
    scale_fill_manual(
      breaks = c(top.4.cat.names, "other", "Unknown"), 
      values = c("#e41a1c", "#ffff33", "#4daf4a", "#377eb8", "#ff7f00", "#984ea3")) +
    labs(x = "K2P Distance (Substitutions per Site)", 
         y = "% of Genome Size Occupied by Repeats")
  ggsave(filename = paste0(gsub(" ", "_", species), ".jpg"), 
         plot = last_plot(), 
         path = paste0("../figures/", vert.invert, "/repeat_landscape"), 
         width = 7680, 
         height = 4320, 
         units = "px", 
         dpi = 1100)
}
