# Compress all .rda files in the data directory
setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70")

cat("Compressing data files with xz compression...\n")
tools::resaveRdaFiles("data", compress="xz")
cat("Data compression complete!\n")
