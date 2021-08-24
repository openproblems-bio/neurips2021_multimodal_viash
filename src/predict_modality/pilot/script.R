library(tidyverse)

df <- readr::read_tsv("output/pilot/predict_modality/output.extract_scores.output.tsv", col_types = c(.default = "c")) %>%
  mutate_at(vars(-one_of(c("method_id", "dataset_id"))), as.numeric) %>%
  mutate(
    dataset_loader = gsub("_.*", "", dataset_id),
    dataset_group = gsub(".*_", "", dataset_id),
    method_group = gsub("_.*", "", method_id),
    format_check = (correct_format + finished) / 2
  )


df %>% group_by(method_id) %>% filter(correct_format > 0) %>% summarise_if(is.numeric, mean)
df %>% group_by(method_id) %>% summarise_if(is.numeric, mean)

ggplot(df) +
  geom_point(aes(mean_spearman, msle, colour = method_id)) +
  facet_wrap(~dataset_group) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1")

ggplot(df) +
  geom_point(aes(mean_spearman, msle, colour = dataset_group)) +
  facet_wrap(~method_id) +
  theme_bw() +
  scale_colour_brewer(palette = "Set2")


bind_rows(df, df %>% mutate(dataset_group = "overall")) %>%
  group_by(method_id, dataset_group) %>%
  filter(correct_format > 0) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup()


df2 <- df %>%
  mutate(filter = finished > 0) %>%
  select(-correct_format, -finished) %>%
  filter(method_id != "dummy_solution") %>%
  gather(metric_id, value, -starts_with("dataset_"), -starts_with("method"), -filter) %>%
  filter(filter | metric_id %in% c("finished", "correct_format", "format_check"))
ggplot(df2) +
  geom_violin(aes(method_id, value, fill = method_id), alpha = .25) +
  #geom_path(aes(method_id, value, group = dataset_id), colour = "gray", size = .5) +
  geom_point(aes(method_id, value, colour = method_id)) +
  facet_wrap(~metric_id, scales = "free_y") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Metric value", title = "Predict Modality scores overview")



ggplot(df2) +
  ggbeeswarm::geom_quasirandom(aes(method_id, value, colour = dataset_loader, shape = dataset_group), groupOnX = TRUE) +
  facet_wrap(~metric_id, scales = "free_x") +
  theme_bw() +
  coord_flip() +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "Metric Value", title = "Predict Modality scores overview")
