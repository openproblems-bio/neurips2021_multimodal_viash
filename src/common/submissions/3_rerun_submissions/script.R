library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds")) %>%
  mutate(
    task_id = `Challenge Phase` %>% gsub(" - .*", "", .) %>% gsub(" ", "_", .) %>% tolower(),
    dest_scores = gsub("\\.zip$", "_scores_1.2.0.rds", dest_zip)
  )

df <- df %>% filter(task_id == "joint_embedding", zip_valid)

# fixed params
pipeline_repo <- "openproblems-bio/neurips2021_multimodal_viash"
pipeline_version <- "1.2.0"
pipeline_tmp_dir <- paste0(Sys.getenv("VIASH_TEMP"), "/neurips2021_work")
data_loc <- Sys.getenv("NEURIPS_DATA_DIR")

zzz <- pbapply::pblapply(seq_len(nrow(df)), function(i) {
  id <- df$id[[i]]
  task_id <- df$task_id[[i]]
  dest_zip <- df$dest_zip[[i]]
  zip_valid <- df$zip_valid[[i]]
  dest_scores <- df$dest_scores[[i]]

  if (!file.exists(dest_scores)) {
    out <- tryCatch({
      if (!zip_valid) {
        stop("Invalid zip file")
      } else {
        dir_path <- tempfile("process_submission")
        on.exit(unlink(dir_path, recursive = TRUE, force = TRUE))
        dir.create(dir_path)
        
        unzip(dest_zip, exdir = dir_path)

        config_path <- paste0(dir_path, "/config.vsh.yaml")

        # build docker
        target_docker_path <- paste0(dir_path, "/target/docker")
        system2(
          "viash",
          args = c(
            "build", conf_path,
            "--config_mod", ".functionality.name := 'method'",
            "--platform", "docker",
            "--output", target_docker_path,
            "--setup", "cachedbuild"
          )
        )

        # build nextflow
        target_nextflow_path <- paste0(dir_path, "/target/nextflow")
        system2(
          "viash",
          args = c(
            "build", conf_path,
            "--config_mod", ".functionality.name := 'method'",
            "--platform", "nextflow",
            "--output", target_nextflow_path
          )
        )

        # run nextflow
        solution_path <- paste0(data_loc, task_id)
        output_dir <- paste0(dir_path, "/output/evaluation/", task_id, "/")
        predictions <- paste0(dir_path, "/output/predictions/", task_id, "/**.h5ad")
        pipeline_script <- paste0("src/", task_id, "/workflows/evaluate_submission/main.nf")

        num_preds <- list.files(dirname(predictions), pattern = ".*\\.h5ad$", recursive = TRUE) %>% length
        if (num_preds == 0) {
          stop("No prediction files found")
        }

        system2(
          "nextflow",
          args = c(
                "run", pipeline_repo,
                "-r", pipeline_version,
                "-main-script", pipeline_script,
                "-work-dir", pipeline_tmp_dir,
                "--solutionDir", solution_path,
                "--predictions", predictions,
                "--publishDir", output_dir,
                "-latest",
                "-resume"
          )
        )

        # read output
        summary_tsv <- paste0(output_dir, "output.final_scores.output_summary.tsv")
        if (!file.exists(summary_tsv)) {
          stop("No output was found")
        }
        
        summary <- read_tsv(summary_tsv)

        summary %>% mutate(id = id) %>% select(-method_id, -var) %>% select(id, everything())
        
      }
    }, error = function(e) {
      tibble(id, run_error = e$message)
    })

    write_rds(out, dest_scores)
  }

  NULL
})


out <- bind_rows(pbapply::pblapply(seq_len(nrow(df)), function(i) {
  id <- df$id[[i]]
  dest_scores <- df$dest_scores[[i]]

  if (file.exists(dest_scores)) {
    tryCatch({ read_rds(dest_scores) }, error = function(e) NULL)
  } else {
    NULL
  }
}))
