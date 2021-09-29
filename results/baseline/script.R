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
  ggplot(summ_pm, aes(value, ml)) +
  geom_point(aes(colour = dataset_subtask)) +
  geom_point(aes(colour = dataset_subtask), function(df) df %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean") %>% filter(is.finite(value)), size = 4) +
  ggrepel::geom_label_repel(aes(label = "DNF"), function(df) df %>% filter(!is.finite(value)), nudge_x = -1, box.padding = 1, direction = "x") +
  theme_bw() +
  labs(x = "RMSE (lower is better)", y = NULL, colour = "Subtask", tag = "a") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2"), "#222222")) +
  xlim(0, 1.5)

gmm <-
  ggplot(summ_mm, aes(value, ml)) +
  geom_point(aes(colour = dataset_subtask)) +
  geom_point(aes(colour = dataset_subtask), function(df) df %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean"), size = 4) +
  ggrepel::geom_label_repel(aes(label = "Score = 1"), function(df) df %>% filter(!is.finite(value)) %>% select(value, ml) %>% unique(), box.padding = 1, direction = "x") +
  xlim(0, .002) +
  theme_bw() +
  labs(x = "Match Probability (higher is better)", y = NULL, colour = "Subtask", tag = "b") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2"), "#222222"))

gje <-
  ggplot(summ_je) +
  geom_point(aes(value, ml, colour = dataset_subtask)) +
  geom_point(aes(value, ml, colour = "Mean"), function(df) df %>% group_by(ml) %>% summarise(value = mean(value), dataset_subtask = "Mean"), size = 4) +
  theme_bw() +
  labs(x = "Geometric Mean (higher is better)", y = NULL, colour = "Subtask", tag = "c") +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(4, "Set2")[1:2], "#222222"))

g <- patchwork::wrap_plots(gpm, gmm, gje, ncol = 1, heights = c(8, 7, 6))

ggsave("results/baseline/figure_baseline.pdf", height = 8, width = 8)

system2(command = "pdfcrop", args = c("results/baseline/figure_baseline.pdf", "results/baseline/figure_baseline.pdf"))
