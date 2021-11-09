library(tidyverse)

dir <- "output/submissions"

all_submissions <- read_rds(paste0(dir, "/all_submissions.rds"))

out <- pbapply::pblapply(seq_len(nrow(all_submissions)), cl = 30, function(i) {
  id <- all_submissions$id[[i]]
  dest_zip <- all_submissions$dest_zip[[i]]
  zip_valid <- all_submissions$zip_valid[[i]]
  dest_json <- all_submissions$dest_json[[i]]
  json_valid <- all_submissions$json_valid[[i]]
  status <- all_submissions$Status[[i]]

  df <- tibble(id = id)

  df2 <-
    tryCatch({
      # check if yaml exists
      if (!zip_valid) stop("Zip file was not valid")

      conn <- unz(dest_zip, "config.vsh.yaml")
      yaml_obj <-
        tryCatch({
          yaml::yaml.load(readLines(conn))
        }, finally = {
          close(conn)
        }, error = function(e) {
          stop("Config file was not valid")
        })
      if (is.null(yaml_obj)) stop("Config is empty?")

      authors_df <- map_df(yaml_obj$functionality$authors, function(aut) {
        aut$roles <- paste(aut$roles, collapse = ", ")
        props <- aut$props
        aut[names(props)] <- unlist(props)
        aut$props <- NULL
        as.data.frame(aut)
      })
      resources <- yaml_obj$functionality$resources
      res_types <- map_chr(resources, function(res) {
        if ("type" %in% names(res)) {
          res[["type"]]
        } else {
          "file"
        }
      })
      script <- resources[[first(grep("script", res_types))]]

      len1 <- function(x, err = TRUE) {
        if (length(x) == 0) {
          z <- NA
          class(z) <- class(x)
          z
        } else if (length(x) > 1) {
          comb <- paste0(x, collapse = ", ")
          if (err) {
            stop("Length of ", comb, " is >1")
          } else {
            comb
          }
        } else {
          x
        }
      }

      df. <- tibble(
        task = gsub("_methods", "", yaml_obj$functionality$namespace) %>% len1,
        method_id = yaml_obj$functionality$name %>% len1,
        description = yaml_obj$functionality$description %>% len1,
        maintainer = authors_df %>% filter(grepl("maintainer", roles)) %>% pull(name) %>% len1(err = FALSE),
        authors = list(authors_df),
        num_authors = nrow(authors_df) %>% len1,
        language = gsub("_script", "", script$type) %>% len1
      )

      info <- yaml_obj$functionality$info
      if (length(info) > 0) {
        df.[names(info)] <- info
      }
      df.
    }, error = function(e) {
      tibble(zip_error = e$message)
    })

  df3 <-
    tryCatch({
      if (!json_valid) stop("Json file was not valid")
      scores <-
        if (status == "finished") {
          jsonlite::parse_json(gsub("'", '"', readLines(dest_json, warn = FALSE)), simplifyVector = TRUE)
        } else {
          tibble(a = 1)[, -1]
        }
      tibble(scores = list(scores))
    }, error = function(e) {
      tibble(json_error = e$message)
    })

  bind_cols(df, df2, df3)
})

df <- all_submissions %>% left_join(bind_rows(out), by = "id")

expected_zip_errs <- c("Zip file was not valid", "Config file was not valid", "Config is empty?")
df %>% filter(!is.na(json_error) & json_error != "Json file was not valid") %>% pull(json_error) %>% table
df %>% filter(!is.na(zip_error) & !zip_error %in% expected_zip_errs) %>% pull(zip_error) %>% table
df %>% mutate(i = row_number()) %>% filter(!is.na(zip_error) & !zip_error %in% expected_zip_errs) %>% pull(i)

write_rds(df, paste0(dir, "/all_submissions_read.rds"))
