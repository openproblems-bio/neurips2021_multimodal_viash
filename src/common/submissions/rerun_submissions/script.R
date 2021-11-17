library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds")) %>%
  filter(id >= 157995) %>%
  mutate(
    task_id = `Challenge Phase` %>% gsub(" - .*", "", .) %>% gsub(" ", "_", .) %>% tolower(),
    # dest_scores = gsub("\\.zip$", "_scores_1.2.0.rds", dest_zip)
    dest_scores = gsub("\\.zip$", "_scores_mainbuild.rds", dest_zip),
  )

# # remove all errorred
# walk(seq_len(nrow(df)), function(i) {
#   id <- df$id[[i]]
#   dest_scores <- df$dest_scores[[i]]

#   if (file.exists(dest_scores)) {
#     try({
#       scor <- read_rds(dest_scores)
#       # if (!is.null(scor$run_error)) {
#         cat("Removing ", i, ", ", id, "\n", sep = "")
#         file.remove(dest_scores)
#       # }
#     })
#   }
# })

# df <- df %>% filter(task_id == "joint_embedding", zip_valid)
# df <- df %>% filter(task_id == "predict_modality", zip_valid)
# df <- df %>% filter(task_id == "match_modality", zip_valid)

# fixed params
pipeline_repo <- "openproblems-bio/neurips2021_multimodal_viash"
# pipeline_version <- "1.2.0"
pipeline_version <- "main_build"
# pipeline_tmp_dir <- paste0(Sys.getenv("VIASH_TEMP"), "/neurips2021_work")
pipeline_tmp_dir <- "work"
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
            "build", config_path,
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
            "build", config_path,
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

        any_dup <- summary %>% select(metric_id, dataset_subtask) %>% { any(duplicated(paste0(.$metric_id, "_", .$dataset_subtask)))}
        if (any_dup) stop("More than 1 method detected")

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
    tryCatch(
      {
        df <- read_rds(dest_scores)

        any_dup <- df %>% select(metric_id, dataset_subtask) %>% { any(duplicated(paste0(.$metric_id, "_", .$dataset_subtask)))}

        if (any_dup) stop("More than 1 method detected")

        df
      }, error = function(e) {
        NULL
      })
  } else {
    NULL
  }
}))

ids <- df %>% filter(task_id == "predict_modality") %>% pull(id)
yy <- out %>% filter(id %in% ids) %>% spread(metric_id, mean) %>% filter(correct_format == 1) %>% mutate(rmse = ifelse(rmse > 1000, Inf, rmse))
ggplot(yy) + geom_point(aes(rmse, mean_pearson_per_cell)) + facet_wrap(~dataset_subtask) + scale_x_log10()
ggplot(yy) + geom_point(aes(rmse, mean_spearman_per_cell)) + facet_wrap(~dataset_subtask) + scale_x_log10() +
ggplot(yy) + geom_point(aes(rmse, mean_pearson_per_gene)) + facet_wrap(~dataset_subtask) + scale_x_log10()
ggplot(yy) + geom_point(aes(rmse, mean_spearman_per_gene)) + facet_wrap(~dataset_subtask) + scale_x_log10()

orig_score <- df %>%
  filter(zip_valid, Status != "failed") %>%
  select(id, Status, task_id, scores) %>%
  mutate(scores = map(scores, function(sc) {
    if (!is.null(sc)) {
      sc %>% gather(metric, orig_score)
    } else {
      NULL
    }
  })) %>%
  unnest(scores) %>%
  mutate(
    metric_prefix = case_when(
      task_id == "predict_modality" ~ "rmse_",
      task_id == "match_modality" ~ "match_probability_",
      TRUE ~ ""
    ),
    metric = paste0(metric_prefix, metric)
  ) %>%
  select(-Status, -metric_prefix, -task_id)
new_score <- out %>%
  mutate(metric = paste0(metric_id, "_", dataset_subtask)) %>%
  filter(metric %in% unique(orig_score$metric) | metric_id %in% c("ti_cons_batch_mean", "pairing_aupr")) %>%
  transmute(id, metric, new_score = mean, error = run_error)
scores <- full_join(orig_score, new_score, by = c("id", "metric")) %>%
  filter(!is.na(orig_score) | !is.na(new_score)) %>%
  left_join(df %>% select(1:10, task_id), by = "id") %>%
  mutate(
    orig_score = ifelse(!is.na(orig_score), ifelse(orig_score > 1000, Inf, orig_score), ifelse(grepl("rmse", metric), Inf, 0)),
    new_score = ifelse(!is.na(new_score), ifelse(new_score > 1000, Inf, new_score), ifelse(grepl("rmse", metric), Inf, 0)),
    status = case_when(
      orig_score < 1e-4 & new_score > 1e-4 ~ "better",
      new_score < 1e-4 & orig_score > 1e-4 ~ "worse",
      abs(orig_score - new_score) / (new_score + orig_score) > .2 ~ "different",
      TRUE ~ "good"
    )
  )
nr <- c(predict_modality = 1, joint_embedding = 3, match_modality = 2)
patchwork::wrap_plots(
  list = map(unique(scores$task_id), function(tid) {
    ggplot(scores %>% filter(task_id == tid)) +
      geom_point(aes(orig_score, new_score, colour = status)) +
      facet_wrap(~metric, scales = "free", nrow = nr[[tid]]) +
      theme_bw() +
      labs(title = tid)
  }),
  ncol = 1,
  heights = nr
)
ggplot(scores) +
  geom_point(aes(orig_score, new_score, colour = status)) +
  facet_wrap(task_id~metric, scales = "free") +
  theme_bw()

ids <- scores %>% filter((new_score == 0 & orig_score > 0) | (is.finite(orig_score) & !is.finite(new_score))) %>% pull(id) %>% unique
df %>% filter(id %in% !!ids)








je_baseline <- read_tsv("results/inhouse_joint_embedding_scores.tsv") %>%
  mutate(metric = paste0(metric_id, "_", dataset_subtask)) %>%
  filter(metric %in% unique(scores$metric))


weights <-
  je_scores_g %>%
  filter(value > 0) %>%
  group_by(metric) %>%
  summarise(
    lower = floor(min(value) * 20) / 20,
    upper = ceiling(max(value) * 20) / 20
  )

je_comb <- bind_rows(
  je_baseline %>% transmute(method_id, team_id = "baseline", method_type, metric, value),
  je_scores_g %>% transmute(method_id = paste0("submission_", id), team_id = `Team Name`, method_type = "submission", metric, value)
) %>%
  mutate(
    value = ifelse(is.na(value), 0, value)
  )

weight_values <- c(
  asw_batch_ADT = 2,
  asw_batch_ATAC = 2,
  asw_label_ADT = 1,
  asw_label_ATAC = 1,
  cc_cons_ADT = 1,
  cc_cons_ATAC = 1,
  graph_conn_ADT = 2,
  graph_conn_ATAC = 2,
  nmi_ADT = 1,
  nmi_ATAC = 1,
  ti_cons_mean_ADT = 1,
  ti_cons_mean_ATAC = 1
)
weights <-
  je_comb %>%
  filter(value > 0) %>%
  group_by(metric) %>%
  summarise(
    lower = floor(min(value) * 20) / 20,
    upper = ceiling(max(value) * 20) / 20,
    weight = weight_values[metric[[1]]],
    .groups = "drop"
  )

je_comb2 <- je_comb %>%
  left_join(weights, by = "metric") %>%
  mutate(
    scaled_value = value %>% pmin(upper) %>% pmax(lower) %>% {(. - lower) / (upper - lower)}
  )

je_out <- je_comb2 %>%
  group_by(method_id, team_id, method_type) %>%
  summarise(
    arith_mean = mean(value),
    scaled_arith_mean = mean(scaled_value),
    weighted_scaled_arith_mean = sum(scaled_value * weight) / sum(weight),
    .groups = "drop"
  ) %>%
  left_join(je_comb2 %>% select(method_id, metric, value) %>% spread(metric, value), by = "method_id")

je_out %>% arrange(desc(arith_mean))
je_out %>% arrange(desc(scaled_arith_mean))
je_out %>% arrange(desc(weighted_scaled_arith_mean))

je_out %>%
  select(1:6) %>%
  gather(agg, score, 4:6) %>%
  ggplot() +
  geom_point(aes(score, agg, colour = method_type)) +
  geom_path(aes(score, agg, group = method_id, colour = method_type)) +
  theme_bw()

