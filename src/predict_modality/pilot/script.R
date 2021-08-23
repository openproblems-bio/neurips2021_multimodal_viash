library(tidyverse)

df <- readr::read_tsv("output/pilot/predict_modality/output.extract_scores.output.tsv", col_types = c(.default = "c")) %>%
  mutate_at(vars(-one_of(c("method_id", "dataset_id"))), as.numeric) %>%
  mutate(
    dataset_group = gsub(".*_", "", dataset_id),
    method_group = gsub("_.*", "", method_id)
  )

df %>% group_by(method_id) %>% filter(correct_format > 0) %>% summarise_if(is.numeric, mean)
df %>% group_by(method_id) %>% summarise_if(is.numeric, mean)

ggplot(df) +
  geom_point(aes(mean_spearman, msle, colour = method_id)) +
  facet_wrap(~dataset_group) +
  theme_classic() +
  scale_colour_brewer(palette = "Set2")

ggplot(df) +
  geom_point(aes(mean_spearman, msle, colour = dataset_group)) +
  facet_wrap(~method_id) +
  theme_bw() +
  scale_colour_brewer(palette = "Set2")
