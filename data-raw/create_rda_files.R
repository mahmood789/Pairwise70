# Convert CSV files to RDA format
library(tools)

setwd('C:/Users/user/OneDrive - NHS/Documents/Pairwise70')

datasets <- list(
  list(name='CD000028_pub4_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD000028_pub4_CD000028-data-rows.csv'),
  list(name='CD002042_pub6_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD002042_pub6_CD002042-data-rows.csv'),
  list(name='CD002122_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD002122_pub3_CD002122-data-rows.csv'),
  list(name='CD006404_pub5_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD006404_pub5_CD006404-data-rows.csv'),
  list(name='CD006742_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD006742_pub3_CD006742-data-rows.csv'),
  list(name='CD006942_pub4_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD006942_pub4_CD006942-data-rows.csv'),
  list(name='CD008940_pub4_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD008940_pub4_CD008940-data-rows.csv'),
  list(name='CD009022_pub4_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD009022_pub4_CD009022-data-rows.csv'),
  list(name='CD011192_pub4_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD011192_pub4_CD011192-data-rows.csv'),
  list(name='CD012335_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD012335_pub2_CD012335-data-rows.csv'),
  list(name='CD012503_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD012503_pub3_CD012503-data-rows.csv'),
  list(name='CD012886_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD012886_pub3_CD012886-data-rows.csv'),
  list(name='CD013462_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD013462_pub3_CD013462-data-rows.csv'),
  list(name='CD013801_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD013801_pub2_CD013801-data-rows.csv'),
  list(name='CD014617_pub3_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD014617_pub3_CD014617-data-rows.csv'),
  list(name='CD014623_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD014623_pub2_CD014623-data-rows.csv'),
  list(name='CD015046_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD015046_pub2_CD015046-data-rows.csv'),
  list(name='CD015284_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD015284_pub2_CD015284-data-rows.csv'),
  list(name='CD015584_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD015584_pub2_CD015584-data-rows.csv'),
  list(name='CD015698_pub2_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD015698_pub2_CD015698-data-rows.csv'),
  list(name='CD016001_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD016001_CD016001-data-rows.csv'),
  list(name='CD016002_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD016002_CD016002-data-rows.csv'),
  list(name='CD016131_data', csv='C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise/10_1002_14651858_CD016131_CD016131-data-rows.csv')
)

for (ds in datasets) {
  cat(paste0("Converting ", ds$name, "...\n"))

  # Read CSV
  data <- read.csv(ds$csv, stringsAsFactors = FALSE)

  # Assign to variable name
  assign(ds$name, data)

  # Save as RDA
  save(list = ds$name, file = file.path('data', paste0(ds$name, '.rda')))

  cat(paste0("  Saved: ", nrow(data), " rows\n"))
}

cat("\nAll datasets converted to RDA format!\n")
