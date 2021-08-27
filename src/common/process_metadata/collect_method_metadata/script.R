cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("yaml", quietly = TRUE)

## VIASH START
par <- list(
  input = list.files("src/", recursive = TRUE, full.names = TRUE, pattern = "*.vsh.yaml"),
  output = "results/meta_methods.tsv"
)
## VIASH END

method_paths <- par$input[grepl("/methods/", par$input)]

meta_method <- map_df(method_paths, function(x) {
  yaml_str <- system(paste0("viash config view ", x), intern = TRUE, ignore.stderr = TRUE)
  yaml_obj <- yaml::yaml.load(string = yaml_str)
  authors_df <- map_df(yaml_obj$functionality$authors, function(aut) {
    aut$roles <- paste(aut$roles, collapse = ", ")
    props <- aut$props
    aut[names(props)] <- unlist(props)
    aut$props <- NULL
    as.data.frame(aut)
  })
  resources <- yaml_obj$functionality$resources
  script <- resources[[grep("script", map_chr(resources, "type"))]]
  
  df <- tibble(
    task = gsub("_methods", "", yaml_obj$functionality$namespace),
    name = yaml_obj$functionality$name,
    description = yaml_obj$functionality$description,
    maintainer = authors_df %>% filter(grepl("maintainer", roles)) %>% pull(name),
    authors = list(authors_df),
    num_authors = nrow(authors_df),
    language = gsub("_script", "", script$type)
  )

  info <- yaml_obj$functionality$info
  if (length(info) > 0) {
    df[names(info)] <- unlist(info)
  }

  df
})

cat("Writing output file\n")
write_tsv(meta_method %>% select(-authors), par$output)
