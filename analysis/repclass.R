


files <- list.files("../results/divsums")
sp <- gsub("_", " ", sub("\\..*", "", files))

dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$species %in% sp, ]
dat <- dat[order(dat$species == sp), ]

species <- clade <- LINE <- SINE <- LTR <- DNA <- RC <- satellite <- c()
for (i in 1:length(sp)) {
  s <- sp[i]
  species <- c(species, s)
  clade <- c(clade, dat[i, ]$clade)
  asmb.bp <- dat[i, ]$asmblysize.Mbp*1000000
  table <- read.delim(paste0("../results/divsums/", files[i]), 
                     sep = "\t", 
                     skip = 4, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) 
  drop1 <- which(table[[1]] == "Coverage for each repeat class and divergence (Kimura)")
  table <- table[4:drop1-1, ]
  drop2 <- which(table[[1]] == "Simple_repeat")
  table <- table[-drop2, ]
  table[, 3:5] <- sapply(table[, 3:5], as.numeric)
  table$score <- (table$absLen / asmb.bp) * (1 - (table$Kimura. / 100))
  # table$Subclass <- sub(".*\\/", "", table$Class)
  table$Class <- sub("/.*", "", table$Class)
  class <- c("LINE", "SINE", "LTR", "DNA", "RC", "Satellite")
  sums <- c()
  for (j in class) {
    sums <- c(sums, sum(table[table$Class == j, ]$score))
  }
  LINE <- c(LINE, sums[1])
  SINE <- c(SINE, sums[2])
  LTR <- c(LTR, sums[3])
  DNA <- c(DNA, sums[4])
  RC <- c(RC, sums[5])
  satellite <- c(satellite, sums[6])
}
df <- data.frame(species, clade, LINE, SINE, LTR, DNA, RC, satellite)

df$rsq <- dat$rsq
df$w <- 1 - (abs(dat$asmblysize.Mbp - dat$est.gnsz.Mbp) / dat$est.gnsz.Mbp)
df <- na.omit(df)
df <- df[df$w >= 0, ]

df <- df[df$clade == "Mammalia", ]

df$LINEn <- (df$LINE - min(df$LINE)) / (max(df$LINE) - min(df$LINE))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$LINEn, 
     df$rsq, 
     xlab = "LINE intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ LINEn, data = df, weight = w))
summary(glm(rsq ~ LINEn, data = df, weight = w))

df$SINEn <- (df$SINE - min(df$SINE)) / (max(df$SINE) - min(df$SINE))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$SINEn, 
     df$rsq, 
     xlab = "SINE intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ SINEn, data = df, weight = w))
summary(glm(rsq ~ SINEn, data = df, weight = w))

df$LTRn <- (df$LTR - min(df$LTR)) / (max(df$LTR) - min(df$LTR))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$LTRn, 
     df$rsq, 
     xlab = "LTR intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ LTRn, data = df, weight = w))
summary(glm(rsq ~ LTRn, data = df, weight = w))

df$DNAn <- (df$DNA - min(df$DNA)) / (max(df$DNA) - min(df$DNA))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$DNAn, 
     df$rsq, 
     xlab = "DNA intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ DNAn, data = df, weight = w))
summary(glm(rsq ~ DNAn, data = df, weight = w))

df$RCn <- (df$RC - min(df$RC)) / (max(df$RC) - min(df$RC))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$RCn, 
     df$rsq, 
     xlab = "RC intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ RCn, data = df, weight = w))
summary(glm(rsq ~ RCn, data = df, weight = w))

df$satelliten <- (df$satellite - min(df$satellite)) / (max(df$satellite) - min(df$satellite))
cols <- viridis(length(unique(df$w)), alpha = 0.45)[as.factor(df$w)]
plot(df$satelliten, 
     df$rsq, 
     xlab = "satellite intensity",
     ylab = "consistency", 
     col = cols, 
     pch = 16)
abline(glm(rsq ~ satelliten, data = df, weight = w))
summary(glm(rsq ~ satelliten, data = df, weight = w))
  