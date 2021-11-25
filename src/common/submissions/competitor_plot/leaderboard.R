library(tidyverse)

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds"))




pal <- RColorBrewer::brewer.pal(3, "Dark2")



# Predict modality --------------------------------------------------------


pm_scores <- df %>%
  filter(`Challenge Phase` == "Predict Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

pm_scores %>% select(id, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pm_scores_g <- pm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

pmzz <- pm_scores_g %>%
  group_by(metric, `Team Name`) %>%
  summarise(value = min(value, na.rm = TRUE)) %>%
  arrange((value)) %>%
  mutate(position = row_number()) %>%
  slice(1:9) %>%
  ungroup() %>%
  mutate(str = paste0("#", position, " - ", `Team Name`))

pm_met_lab <- c(
  Overall = "Mean RMSE",
  ADT2GEX = "RMSE for ADT2GEX",
  GEX2ADT = "RMSE for GEX2ADT",
  ATAC2GEX = "RMSE for ATAC2GEX",
  GEX2ATAC = "RMSE for GEX2ATAC"
)

pms <- map(names(pm_met_lab), function(met) {
  pmzz_ <- pmzz %>% filter(metric == met) %>% mutate(str = forcats::fct_inorder(str))
  min <- min(pmzz_$value)
  max <- max(pmzz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .5
  ggplot(pmzz %>% filter(metric == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[1]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[1]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(y = pm_met_lab[[met]], x = NULL)
})






# Match modality ----------------------------------------------------------

mm_scores <- df %>%
  filter(`Challenge Phase` == "Match Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

mm_scores %>% select(id, Status, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)

mm_scores_g <- mm_scores %>% gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall)


mmzz <- mm_scores_g %>%
  group_by(metric, `Team Name`) %>%
  summarise(value = max(value, na.rm = TRUE)) %>%
  arrange(desc(value)) %>%
  mutate(position = row_number()) %>%
  slice(1:9) %>%
  ungroup() %>%
  mutate(str = paste0("#", position, " - ", `Team Name`))

mm_met_lab <- c(
  Overall = "Mean match probability",
  ADT2GEX = "Match probability for ADT2GEX",
  GEX2ADT = "Match probability for GEX2ADT",
  ATAC2GEX = "Match probability for ATAC2GEX",
  GEX2ATAC = "Match probability for GEX2ATAC"
)

mms <- map(names(mm_met_lab), function(met) {
  mmzz_ <- mmzz %>% filter(metric == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))
  min <- min(mmzz_$value)
  max <- max(mmzz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .5
  # nud <- max(mmzz_$value) * .5
  ggplot(mmzz_) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[2]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[2]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    # expand_limits(y = c(.7, .85)) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(y = mm_met_lab[[met]], x = NULL)
})







# Joint embedding ---------------------------------------------------------

je_scores <- df %>%
  filter(`Challenge Phase` == "Joint Embedding - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  )

je_scores %>% select(id, Status, ends_with("_ADT"), ends_with("_ATAC"), arithmetic_mean)

je_scores_g <- je_scores %>% gather(metric, value, ends_with("_ADT"), ends_with("_ATAC"), arithmetic_mean)

je_scores_g %>%
  filter(`Team Name` == "LiuZ_Lab_BCM", metric == "arithmetic_mean_ATAC") %>%
  select(id, Status, metric, value)

jezz <- je_scores_g %>%
  group_by(metric, `Team Name`) %>%
  summarise(value = max(value, na.rm = TRUE)) %>%
  arrange(desc(value)) %>%
  mutate(position = row_number()) %>%
  slice(1:9) %>%
  ungroup() %>%
  filter(grepl("arithmetic_mean", metric)) %>%
  mutate(str = paste0("#", position, " - ", `Team Name`))

je_met_lab <- c(
  arithmetic_mean = "Arithmetic mean",
  arithmetic_mean_ADT = "Arithmetic mean for ADT",
  arithmetic_mean_ATAC = "Arithmetic mean for ATAC"
)

jes <- map(names(je_met_lab), function(met) {
  jezz_ <- jezz %>% filter(metric == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))
  min <- min(jezz_$value)
  max <- max(jezz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .5
  ggplot(jezz_) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[3]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[3]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    # expand_limits(y = c(.7, .85)) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(y = je_met_lab[[met]], x = NULL)
})



pm_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Predict Modality"), fontface = "bold") + ggplot2::theme_void()
mm_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Match Modality"), fontface = "bold") + ggplot2::theme_void()
je_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Joint Embedding"), fontface = "bold") + ggplot2::theme_void()
date <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "OpenProblems-NeurIPS2021 Leaderboard on 22 Nov 2021"), fontface = "bold") + ggplot2::theme_void()
em <- patchwork::plot_spacer()

g <- patchwork::wrap_plots(
  c(
    list(em, em, pm_title, em, em),
    pms,
    list(em, em, mm_title, em, em),
    mms,
    list(em, em, je_title, em, em),
    list(em), jes, list(em),
    list(date, em, em, em, em)
  ),
  ncol = 5,
  byrow = TRUE,
  heights = c(.1, 1, .1, 1, .1, 1, .1)
)

layout <- "
QQQQQ
AAAAA
BCDEF
GGGGG
HIJKL
MMMMM
#NOP#
"

g <-
  pm_title + pms[[1]] + pms[[2]] + pms[[3]] + pms[[4]] + pms[[5]] +
  mm_title + mms[[1]] + mms[[2]] + mms[[3]] + mms[[4]] + mms[[5]] +
  je_title + jes[[1]] + jes[[2]] + jes[[3]] +
  date +
  patchwork::plot_layout(design = layout, heights = c(.1, .1, 1, .1, 1, .1, 1))
ggsave("output/leaderboard_2021-11-25.pdf", g, width = 20, height = 10)
