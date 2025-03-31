

library(MuMIn)

# columns in combined table
cols <- c("(Intercept)", "age.dna", "age.line", "age.ltr", "age.sine", 
          "age.unknown", "age.others", "prop.dna", "prop.line", "prop.ltr",
          "prop.sine", "prop.unknown", "prop.others", "age.dna:prop.dna", 
          "age.line:prop.line", "age.ltr:prop.ltr", "age.sine:prop.sine", 
          "age.unknown:prop.unknown", "age.others:prop.others", "shapirowilk.p", 
          "lambda.p", "df", "logLik", "AICc", "delta", "weight")

# factor for attributes(model)$column.types
ctypes <- c(rep("terms", 19), rep("extra", 2), "df", 
            "loglik", "ic", "delta", "weight")
levels <- c("terms", "varying", "extra", "df", 
            "loglik", "ic", "delta", "weight")
ctypes <- factor(ctypes, levels = levels)
names(ctypes) <- cols

# vector for attribute(model)$terms
terms <- cols[1:19]
attributes(terms)$interceptLabel <- "(Intercept)"

# add missing columns and sort by the order of cols
fix_columns <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  for(col in miss) df[[col]] <- NA
  return(df[, cols, drop = FALSE])
}

combined.models <- NULL
combined.coefTables <- list()

for(i in 1:896){
  if (!file.exists(paste0("../results/sauria.models/Sauria.", i, ".rds"))) {
    next
  }
  cur.models <- readRDS(paste0("../results/sauria.models/Sauria.", i, ".rds"))
  
  # fix attributes(cur.models)$column.types and $names
  cur.attr <- attributes(cur.models)
  cur.attr$column.types <- ctypes
  cur.attr$names <- cols
  
  # add missing columns and sort by the order of cols; 
  cur.models <- fix_columns(as.data.frame(cur.models), cols)
  
  # assign unique indices to each model
  cur.models$rn <- paste(i, rownames(cur.models), sep = ".")
  
  # save attributes(cur.models)$coefTables and assign unique indices to each
  cur.coefTables <- cur.attr$coefTables
  names(cur.coefTables) <- paste(i, names(cur.coefTables), sep = ".")
  
  # combine current and previous models
  if (is.null(combined.models)) {
    combined.models <- cur.models
  } else {
    combined.models <- rbind(combined.models, cur.models)
  }
  
  #combine current and previous attribute(cur.models)$coefTables
  combined.coefTables <- c(combined.coefTables, cur.coefTables)
}

# reassign rownames
rownames(combined.models) <- combined.models$rn
combined.models$rn <- NULL

# remove duplicate models based on parameter inclusion
# to accomodate rounding errors, actual values are not matched
# models with the same parameter inclusion patterns have the same parameter estimates
pattern <- apply(combined.models[, 1:19], 1, function(x) {
  paste(ifelse(is.na(x), "0", "1"), collapse = ".")
})
final.model <- combined.models[!duplicated(pattern), ]
final.coefTables <- combined.coefTables[!duplicated(pattern)]

# reassign rownames
rownames(final.model) <- as.character(seq_len(nrow(final.model)))
names(final.coefTables) <- as.character(seq_len(nrow(final.model)))

# set other attributes
attributes(final.model)$model.calls <- NULL
attributes(final.model)$coefTables <- final.coefTables
attributes(final.model)$column.types <- ctypes
attributes(final.model)$names <- cols
attributes(final.model)$class <- c("model.selection", "data.frame")
attributes(final.model)$terms <- terms

# reorder, recalculate weights, and save
final.model <- final.model[order(final.model$AICc), ]
final.model <- rbind(final.model, final.model[1,])
final.model <- final.model[-1, ]
saveRDS(final.model, "../results/sauria.09.rds")

