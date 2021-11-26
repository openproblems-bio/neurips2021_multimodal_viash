library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds")) %>%
  mutate(
    execution_time = ifelse(`Execution Time(sec.)` == "None", NA_real_, as.numeric(`Execution Time(sec.)`)),
    phase = factor(`Challenge Phase`, levels = c("Predict Modality - Phase 1", "Match Modality - Phase 1", "Joint Embedding - Phase 1", "Predict Modality - Phase 2", "Match Modality - Phase 2", "Joint Embedding - Phase 2"))
  )




g1 <- ggplot(df %>% filter(execution_time > 0)) +
  geom_histogram(aes(execution_time / 60)) +
  facet_wrap(~phase, scales = "free") +
  theme_bw() +
  labs(x = "Execution time (mins)", title = "Execution time of submissions")


g2 <- ggplot(df %>% filter(Status %in% c("submitted", "running"), `Submitted At` > lubridate::ymd_hms("2021-11-22 00:01:00"))) +
  geom_bar(aes(forcats::fct_rev(phase))) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, title = "Number of submissions in queue", y = "# submissions")

g <- patchwork::wrap_plots(
  g1,
  g2,
  nrow = 1,
  widths = c(2,1)
)

filename <- format(Sys.time(), "output/queue_%Y-%m-%d_%H:%M.png")
ggsave(filename, g, width = 12, height = 5)

