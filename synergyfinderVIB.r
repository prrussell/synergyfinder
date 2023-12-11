library(synergyfinder)
library(dplyr)
library(svDialogs)
library(writexl)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# choose file
filepath <- file.choose()

# enter assay round
user.input <- dlgInput("Enter assay round number", "")$res

# create folders for heat maps and synergy scores
if (!dir.exists('Heat Maps')){dir.create('Heat Maps')}
if (!dir.exists('Synergy Scores')){dir.create('Synergy Scores')}

# read selected file
readfile <- read.csv(filepath)

# file preprocessing
# combine cell line and kinase into target field
readfile$cell_kinase <- paste(readfile$CELL_LINE, readfile$KINASE)

# remove forward slashes from target
readfile$cell_kinase <- gsub("/","_",as.character(readfile$cell_kinase))

# assign block IDs and get unique list of kinases
df_pairs <- readfile %>% distinct(
  SAMPLE_ID_1, SAMPLE_ID_2, cell_kinase, .keep_all=FALSE)
targets <- unique(df_pairs$cell_kinase)
block_id <- seq_len(nrow(df_pairs))
df_pairs <- cbind(block_id, df_pairs)
df <- merge(readfile, df_pairs)

# rename columns to match synergyfinder expected input
df <- df %>%
  rename(
    drug_row = SAMPLE_ID_1,
    drug_col = SAMPLE_ID_2,
    conc_r = CONC_NM_1,
    conc_c = CONC_NM_2,
    inhibition = PCT_INH
    )
units = "nM"
df$conc_r_unit <- units
df$conc_c_unit <- units

# create empty dataframe for synergy scores
all_scores <- data.frame(block_id=integer(),
                  drug1=character(),
                  drug2=character(),
                  target=character(),
                  conc_unit1=character(),
                  conc_unit2=character(),
                  input_type=character(),
                  replicate=logical(),
                  response_p_value=character(),
                  response_origin_p_value=character(),
                  stringsAsFactors=FALSE)

# split file into sections by target
for (x in targets){

filter_df <- df %>% filter(cell_kinase == x)
filter_ids <- filter_df %>% distinct

# create folder for target
if (!dir.exists(paste('Heat Maps/', x, sep="")))
{dir.create(paste('Heat Maps/', x, sep=""))}

# synergyfinder data preprocessing function
res <- ReshapeData(
  data = filter_df,
  data_type = "inhibition",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  seed = 1)

# generate and save heat maps
PlotDoseResponse(
  data = res,
  block_ids = NULL,
  drugs = c(1,2),
  save_file = TRUE,
  file_type = "png"
)

# rename heat maps with target and round number
fnames <- list.files(pattern = "\\.png$")
new_names <- paste0(x, "_Rd", user.input, "_", fnames)
file.rename(paste0("./", fnames), paste0("./", new_names))
fnames_new <- list.files(pattern = "\\.png$")


# move heat maps to final folder
file.copy(from = paste0("./", fnames_new),
          to = paste0("./Heat Maps/", x, "/", fnames_new))
file.remove(from = paste0("./", fnames_new))

# calculated synergy scores for selected target
res <- CalculateSynergy(
  data = res,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  Emin = NA,
  Emax = NA,
  correct_baseline = "non")

# add target column to synergy scores
scores <- res$drug_pairs
scores$target <- x

# append synergy scores to complete list
all_scores <- rbind(all_scores, scores)
}

# save complete list of synergy scores
write_xlsx(all_scores, paste("./Synergy Scores/", "Rd", user.input, "_synscores.xlsx", sep=""))
