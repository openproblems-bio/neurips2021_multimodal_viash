library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds"))
pal <- setNames(RColorBrewer::brewer.pal(3, "Set1")[c(3, 2, 1)], c("positive_control", "baseline", "negative_control"))

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

cutoff <- c(ADT2GEX=.5, ATAC2GEX=1, GEX2ADT=1, GEX2ATAC=.5, Overall=1)
gpm <- ggplot() +
  ggridges::geom_density_ridges(aes(value, 1), pm_scores_g %>% filter(value < cutoff[metric]), jittered_points = TRUE, point_size = .5) +
  geom_vline(aes(xintercept = value, colour = method_type), data = pm_baseline %>% filter(value < 1) %>% mutate(value = ifelse(method_type == "positive_control", -Inf, value)), linetype = "dashed", size = .5) +
  theme_bw() +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  labs(x = "RMSE (lower is better)", y = "Density", colour = "Pilot method") +
  scale_colour_manual(values = pal)

gpm


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
  ggridges::geom_density_ridges(aes(value, 1), mm_scores_g %>% filter(value > 0), jittered_points = TRUE, point_size = .5) +
  geom_vline(aes(xintercept = value, colour = method_type), data = mm_baseline %>% filter(value > 0, method_id != "dummy_semisolution"), linetype = "dashed", size = .5) +
  theme_bw() +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  labs(x = "Match probability score (higher is better)", y = "Density", colour = "Pilot method") +
  scale_x_log10() +
  scale_colour_manual(values = pal)


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
  ggridges::geom_density_ridges(aes(value, 1), je_scores_g %>% filter(value > 0), jittered_points = TRUE, point_size = .5) +
  geom_vline(aes(xintercept = value, colour = method_type), data = je_baseline %>% filter(value > 0), linetype = "dashed", size = .5) +
  theme_bw() +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  labs(x = "Mean of metrics (higher is better)", y = "Density", colour = "Pilot method") +
  scale_colour_manual(values = pal)


# Combined plot ---------------------------------------------------------
g <- patchwork::wrap_plots(
  gpm + labs(title = "Task 1: Predict Modality"),
  gmm + labs(title = "Task 2: Match Modality"),
  gje + labs(title = "Task 3: Joint Embedding"),
  ncol = 1,
  guides = "collect"
) & theme(legend.position = "bottom")

ggsave("results/phase2_private/compare_pilot_submissions.pdf", g, width = 12, height = 8)
ggsave("results/phase2_private/compare_pilot_submissions.png", g, width = 12, height = 8)

ggsave("results/phase2_private/compare_pilot_submissions_task1.png", gpm + facet_wrap(~metric, scales = "free"), width = 8, height = 4)
ggsave("results/phase2_private/compare_pilot_submissions_task2.png", gmm + facet_wrap(~metric, scales = "free"), width = 8, height = 4)
ggsave("results/phase2_private/compare_pilot_submissions_task3.png", gje + facet_wrap(~metric, scales = "free", nrow = 2), width = 8, height = 4)
