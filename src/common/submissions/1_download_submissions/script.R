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
    dest_zip = paste0(dir, "/", basename(url_zip)),
    dest_json = paste0(dir, "/", basename(url_json)),
    zip_exists = file.exists(dest_zip),
    zip_valid = map2_lgl(zip_exists, dest_zip, function(a, b) a && valid_zip(b)),
    json_exists = file.exists(dest_json),
    json_valid = map2_lgl(json_exists, dest_json, function(a, b) a && valid_json(b))
  )

write_rds(all_submissions, paste0(dir, "/all_submissions.rds"))

table(all_submissions$zip_exists)
table(all_submissions$zip_valid)
table(all_submissions$json_exists)
table(all_submissions$json_valid)

options(timeout = max(10000, getOption("timeout")))

dl_zip <-
  all_submissions %>%
  filter(!zip_valid)
zzz <- pbapply::pbmapply(
  cl = 10,
  FUN = download.file,
  url = dl_zip$url_zip,
  destfile = dl_zip$dest_zip,
  quiet = TRUE
)

dl_json <-
  all_submissions %>%
  filter(!json_exists)
zzz <- pbapply::pbmapply(
  cl = 10,
  FUN = download.file,
  url = dl_zip$url_zip,
  destfile = dl_zip$dest_zip,
  quiet = TRUE
)

zzz <- pbapply::pblapply(seq_len(nrow(all_submissions)), cl = 10, function(i) {
  id <- all_submissions$id[[i]]
  dest_zip <- all_submissions$dest_zip[[i]]
  dest_json <- all_submissions$dest_json[[i]]
  status <- all_submissions$Status[[i]]

  tryCatch({
    if (!file.exists(dest_zip) || !valid_zip(dest_zip)) {
      download.file(all_submissions$`Submitted File`[[i]], dest_zip, quiet = TRUE)
    }
  }, error = function(e) {
    cat("Error: ", e$message, "\n", sep = "")
  })

  tryCatch({
    # check if already downloaded
    if (!file.exists(dest_json)) {
      download.file(all_submissions$`Submission Result File`[[i]], dest_json, quiet = TRUE)
    }
  }, error = function(e) {
    cat("Error: ", e$message, "\n", sep = "")
  })
})
