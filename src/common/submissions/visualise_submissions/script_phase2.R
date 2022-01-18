library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds"))



# Joint embedding ---------------------------------------------------------
je_metrics <- c("asw_batch", "asw_label", "cc_cons", "graph_conn", "nmi", "ti_cons_batch_mean", "arithmetic_mean")
je_metrics2 <- c(paste0(je_metrics, "_ADT"), paste0(je_metrics, "_ATAC"), "arithmetic_mean")
je_scores <- df %>%
  filter(`Challenge Phase` == "Joint Embedding - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

je_scores %>% select(id, language, Status, ends_with("_ADT"), ends_with("_ATAC"))

je_scores_g <- je_scores %>% gather(metric, value, arithmetic_mean_ADT, arithmetic_mean_ATAC, arithmetic_mean)

je_baseline <- read_tsv("results/phase2_private/inhouse_joint_embedding_scores.tsv") %>%
  mutate(metric = paste0(metric_id, "_", dataset_subtask)) %>%
  filter(metric %in% c("arithmetic_mean_ADT", "arithmetic_mean_ATAC"))
je_baseline <- bind_rows(je_baseline, je_baseline %>% group_by(method_id, method_type) %>% summarise(value = mean(value), metric = "arithmetic_mean"))



gje <- ggplot() +
  ggbeeswarm::geom_quasirandom(aes(1, value), je_scores_g %>% filter(value > 0)) +
  geom_hline(aes(yintercept = value, colour = method_type), je_baseline %>% filter(value > 0)) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_bw()
gje












# Match modality ----------------------------------------------------------

mm_scores <- df %>%
  filter(`Challenge Phase` == "Match Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

mm_scores %>% select(id, language, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

mm_scores_g <- mm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

mm_baseline <- read_tsv("results/phase2_private/inhouse_match_modality_scores.tsv") %>%
  filter(metric_id == "match_probability") %>%
  mutate(metric = dataset_subtask)
mm_baseline <- bind_rows(
  mm_baseline,
  mm_baseline %>% group_by(method_id, method_type) %>% summarise(metric = "Overall", value = mean(value))
)


gmm <- ggplot() +
  ggbeeswarm::geom_quasirandom(aes(1, value), mm_scores_g %>% filter(value > 0)) +
  geom_hline(aes(yintercept = value, colour = method_type), mm_baseline %>% filter(value > 0)) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_bw() +
  scale_y_log10()
gmm












# Predict modality --------------------------------------------------------


pm_scores <- df %>%
  filter(`Challenge Phase` == "Predict Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

pm_scores %>% select(id, language, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pm_scores_g <- pm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pm_baseline <- read_tsv("results/phase2_private/inhouse_predict_modality_scores.tsv") %>%
  filter(metric_id == "rmse") %>%
  mutate(metric = dataset_subtask)
pm_baseline <- bind_rows(
  pm_baseline,
  pm_baseline %>% group_by(method_id, method_type) %>% summarise(metric = "Overall", value = mean(value))
)




gpm <- ggplot() +
  ggbeeswarm::geom_quasirandom(aes(1, value), pm_scores_g %>% filter(value < 1e4)) +
  geom_hline(aes(yintercept = value, colour = method_type), pm_baseline %>% filter(value < 1e4)) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "RMSE (lower is better)")
  # scale_y_log10()
gpm




patchwork::wrap_plots(
  gpm,
  gmm,
  gje,
  ncol = 1
)


patchwork::wrap_plots(
  gpm + coord_flip(),
  gmm + coord_flip(),
  gje + coord_flip(),
  ncol = 1
)
