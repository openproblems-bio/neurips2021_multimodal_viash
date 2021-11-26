library(tidyverse)

dir <- "output/submissions"

dir.create(dir, showWarnings = FALSE)

valid_zip <- function(filepath) {
  tryCatch({
    unzip(filepath, list = TRUE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}
valid_json <- function(filepath) {
  tryCatch({
    jsonlite::validate(gsub("'", '"', readr::read_lines(filepath)))
  }, error = function(e) {
    FALSE
  })
}

all_submissions <- read_csv(paste0(dir, "/all_submissions.csv")) %>%
  mutate(
    url_zip = `Submitted File`,
    url_json = `Submission Result File`,
    dest_zip = paste0(dir, "/submission_", id, ".zip"),
    dest_json = ifelse(is.na(url_json), NA_character_, paste0(dir, "/submission_", id, ".json")),
    zip_exists = file.exists(dest_zip),
    zip_valid = map2_lgl(zip_exists, dest_zip, function(a, b) a && valid_zip(b)),
    json_exists = !is.na(url_json) & file.exists(dest_json),
    json_valid = map2_lgl(json_exists, dest_json, function(a, b) a && valid_json(b))
  )

write_rds(all_submissions, paste0(dir, "/all_submissions.rds"))

table(all_submissions$zip_exists)
table(all_submissions$zip_valid)
table(all_submissions$json_exists)
table(all_submissions$json_valid)

options(timeout = max(100000, getOption("timeout")))

# dl_zip <-
#   all_submissions %>%
#   filter(!zip_valid)
# zzz <- pbapply::pbmapply(
#   cl = 10,
#   FUN = download.file,
#   url = dl_zip$url_zip,
#   destfile = dl_zip$dest_zip,
#   quiet = TRUE
# )

dl_json <-
  all_submissions %>%
  filter(Status == "finished", !json_valid, !is.na(url_json))
zzz <- pbapply::pbmapply(
  cl = 10,
  FUN = download.file,
  url = dl_json$url_json,
  destfile = dl_json$dest_json,
  quiet = TRUE
)
