library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds"))







# Joint embedding ---------------------------------------------------------

je_scores <- df %>%
  filter(`Challenge Phase` == "Joint Embedding - Phase 1", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

je_scores %>% select(id, language, Status, ends_with("_ADT"), ends_with("_ATAC"))

je_scores_g <- je_scores %>% gather(metric, value, ends_with("_ADT"), ends_with("_ATAC"))

je_baseline <- read_tsv("results/inhouse_joint_embedding_scores.tsv") %>%
  mutate(metric = paste0(metric_id, "_", dataset_subtask)) %>%
  filter(metric %in% unique(je_scores_g$metric))

ggplot() +
  geom_histogram(aes(value), je_scores_g, binwidth = .1) +
  geom_vline(aes(xintercept = value, colour = method_type), je_baseline, binwidth = .1) +
  facet_wrap(~metric, scales = "free") +
  theme_bw()

ggplot() +
  geom_density(aes(value), je_scores_g) +
  geom_vline(aes(xintercept = value, colour = method_type), je_baseline, binwidth = .1) +
  facet_wrap(~metric, scales = "free")

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











# Match modality ----------------------------------------------------------

mm_scores <- df %>%
  filter(`Challenge Phase` == "Match Modality - Phase 1", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

mm_scores %>% select(id, language, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

mm_scores_g <- mm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

mm_baseline <- read_tsv("results/inhouse_match_modality_scores.tsv") %>%
  filter(metric_id == "match_probability") %>%
  mutate(metric = dataset_subtask)
mm_baseline <- bind_rows(
  mm_baseline,
  mm_baseline %>% group_by(method_id, method_type) %>% summarise(metric = "Overall", value = mean(value))
)


mm_comb <- bind_rows(
  mm_baseline %>% transmute(method_id, team_id = "baseline", method_type, metric, value),
  mm_scores_g %>% transmute(method_id = paste0("submission_", id), team_id = `Team Name`, method_type = "submission", metric, value = value / 1000)
)




mm_out <- mm_comb %>%
  spread(metric, value)

mm_out %>% arrange(desc(Overall))
mm_out %>% arrange(desc(ADT2GEX))
mm_out %>% arrange(desc(ATAC2GEX))
mm_out %>% arrange(desc(GEX2ADT))
mm_out %>% arrange(desc(GEX2ATAC))

mm_out %>%
  gather(agg, score, -1:-3) %>%
  mutate(agg = factor(agg, levels = c("Overall", "GEX2ATAC", "ATAC2GEX", "GEX2ADT", "ADT2GEX"))) %>%
  arrange(agg, method_id) %>%
  ggplot() +
  geom_point(aes(score, agg, colour = method_type)) +
  geom_path(aes(score, agg, group = method_id, colour = method_type)) +
  theme_bw() +
  scale_x_log10()

















# Predict modality --------------------------------------------------------


pm_scores <- df %>%
  filter(`Challenge Phase` == "Predict Modality - Phase 1", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

pm_scores %>% select(id, language, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pm_scores_g <- pm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pm_baseline <- read_tsv("results/inhouse_predict_modality_scores.tsv") %>%
  filter(metric_id == "rmse") %>%
  mutate(metric = dataset_subtask)
pm_baseline <- bind_rows(
  pm_baseline,
  pm_baseline %>% group_by(method_id, method_type) %>% summarise(metric = "Overall", value = mean(value))
)


pm_comb <- bind_rows(
  pm_baseline %>% transmute(method_id, team_id = "baseline", method_type, metric, value),
  pm_scores_g %>% transmute(method_id = paste0("submission_", id), team_id = `Team Name`, method_type = "submission", metric, value)
) %>%
  mutate(
    value = ifelse(is.na(value), 99999, value)
  )




pm_out <- pm_comb %>%
  spread(metric, value)

pm_out %>% arrange(Overall)
pm_out %>% arrange(ADT2GEX)
pm_out %>% arrange(ATAC2GEX)
pm_out %>% arrange(GEX2ADT)
pm_out %>% arrange(GEX2ATAC)

pm_out %>%
  gather(agg, score, -1:-3) %>%
  arrange(agg, method_id) %>%
  ggplot() +
  geom_point(aes(score, agg, colour = method_type)) +
  geom_path(aes(score, agg, group = method_id, colour = method_type), alpha = .1) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1))

