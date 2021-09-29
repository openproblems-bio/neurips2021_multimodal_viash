library(tidyverse)

mt_map <- c(baseline = "B", negative_control = "NC", positive_control = "PC")
mt_map2 <- c(baseline = "Baseline", negative_control = "Positive control", positive_control = "Negative control")

summ_pm <- read_tsv("results/inhouse_predict_modality_scores.tsv") %>%
  filter(metric_id == "rmse") %>%
  mutate(
    value = ifelse(value > 100, Inf, value),
    mt = mt_map[method_type],
    mt2 = mt_map2[method_type],
    ml = paste0(mt, " - ", method_label),
    dataset_subtask = gsub("2", "->", dataset_subtask)
  )

summ_mm <- read_tsv("results/inhouse_match_modality_scores.tsv") %>%
  filter(method_id != "dummy_semisolution") %>%
  filter(metric_id == "match_probability") %>%
  mutate(
    value = ifelse(value > .0025, Inf, value),
    mt = mt_map[method_type],
    mt2 = mt_map2[method_type],
    ml = paste0(mt, " - ", method_label),
    dataset_subtask = gsub("2", "->", dataset_subtask)
  )
summ_je <- read_tsv("results/inhouse_joint_embedding_scores.tsv") %>%
  filter(metric_id == "geometric_mean") %>%
  mutate(
    mt = mt_map[method_type],
    mt2 = mt_map2[method_type],
    ml = paste0(mt, " - ", method_label),
    dataset_subtask = gsub("2", "→", dataset_subtask),
    value = ifelse(is.na(value), 0, value)
  )

gpm <-
  ggplot(summ_pm) +
  geom_point(aes(value, ml, colour = dataset_subtask)) +
  geom_point(aes(value, ml, colour = "Mean"), summ_pm %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean"), size = 4) +
  theme_bw() +
  labs(x = "RMSE", y = NULL, colour = "Subtask", tag = "A") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2"), "#222222")) +
  xlim(0, 1.5)

gmm <-
  ggplot(summ_mm) +
  geom_point(aes(value, ml, colour = dataset_subtask)) +
  geom_point(aes(value, ml, colour = "Mean"), summ_mm %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean"), size = 4) +
  # ggforce::facet_zoom(xlim = c(0, .0025)) +
  xlim(0, .002) +
  theme_bw() +
  labs(x = "Match Probability", y = NULL, colour = "Subtask", tag = "B") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2"), "#222222"))

gje <-
  ggplot(summ_je) +
  geom_point(aes(value, ml, colour = dataset_subtask)) +
  geom_point(aes(value, ml, colour = "Mean"), summ_je %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean"), size = 4) +
  theme_bw() +
  labs(x = "Geometric Mean", y = NULL, colour = "Subtask", tag = "C") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2")[1:2], "#222222"))

g <- patchwork::wrap_plots(gpm, gmm, gje, ncol = 1)

ggsave("results/baseline/figure_baseline.pdf", height = 8, width = 8)
