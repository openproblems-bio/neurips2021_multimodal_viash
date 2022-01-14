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
    # zip
    url_zip = `Submitted File`,
    dest_zip = paste0(dir, "/submission_", id, ".zip"),
    zip_exists = file.exists(dest_zip),
    zip_valid = map2_lgl(zip_exists, dest_zip, function(a, b) a && valid_zip(b)),
    # result
    url_json = `Submission Result File`,
    dest_json = ifelse(is.na(url_json), NA_character_, paste0(dir, "/submission_", id, ".json")),
    json_exists = !is.na(url_json) & file.exists(dest_json),
    json_valid = map2_lgl(json_exists, dest_json, function(a, b) a && valid_json(b))#,
    # stdout
    # url_stdout = `Stdout File`,
    # dest_stdout = ifelse(is.na(url_stdout), NA_character_, paste0(dir, "/submission_", id, ".stdout.txt")),
    # stdout_exists = file.exists(dest_stdout),
    # # stderr
    # url_stderr = `Stderr File`,
    # dest_stderr = ifelse(is.na(url_stderr), NA_character_, paste0(dir, "/submission_", id, ".stderr.txt")),
    # stderr_exists = file.exists(dest_stderr),
  )

write_rds(all_submissions, paste0(dir, "/all_submissions.rds"))

table(all_submissions$zip_exists)
table(all_submissions$zip_valid)
table(all_submissions$json_exists)
table(all_submissions$json_valid)

options(timeout = max(100000, getOption("timeout")))

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

dl_zip <-
  all_submissions %>%
  filter(!zip_valid, grepl("Phase 2", `Challenge Phase`))
zzz <- pbapply::pbmapply(
  cl = 10,
  FUN = download.file,
  url = dl_zip$url_zip,
  destfile = dl_zip$dest_zip,
  quiet = TRUE
)

# dl_stdout <-
#   all_submissions %>%
#   filter(Status == "failed", !is.na(url_stdout), !stdout_exists)
# zzz <- pbapply::pbmapply(
#   cl = 10,
#   FUN = download.file,
#   url = dl_stdout$url_stdout,
#   destfile = dl_stdout$dest_stdout,
#   quiet = TRUE
# )
#
# read_stdout <- all_submissions %>%
#   filter(Status == "failed", !is.na(url_stdout), stdout_exists) %>%
#   mutate(text = map_chr(dest_stdout, function(fil) readr::read_lines(fil) %>% paste(collapse = "\n")))
#
# read_stdout %>% filter(grepl("No space left on device", text))
#
#
#
#
#
# dl_stderr <-
#   all_submissions %>%
#   filter(Status == "failed", !is.na(url_stderr), !stderr_exists)
# zzz <- pbapply::pbmapply(
#   cl = 10,
#   FUN = download.file,
#   url = dl_stderr$url_stderr,
#   destfile = dl_stderr$dest_stderr,
#   quiet = TRUE
# )
#
# read_stderr <- all_submissions %>%
#   filter(Status == "failed", !is.na(url_stderr), stderr_exists) %>%
#   mutate(text = map_chr(dest_stderr, function(fil) readr::read_lines(fil) %>% paste(collapse = "\n")))
#
# read_stderr %>% filter(grepl("No space left on device", text))
