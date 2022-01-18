library(tidyverse)

googlesheets4::gs4_auth(email = "brechtus@gmail.com")

dir <- "output/submissions"
df <- read_rds(paste0(dir, "/all_submissions_read.rds")) %>%
  filter(`Team Name` != "OpenProblems") %>%
  mutate(
    submission_time = `Submitted At`,
    team_name = `Team Name`
  )

paste0(sl$`Team Name`, collapse = "\n") %>% cat
paste0(sl$`Team Members Email Id`, collapse = "\n") %>% cat
paste0(sl$`Submitted File`, collapse = "\n") %>% cat
paste0(sl$submission_time, collapse = "\n") %>% cat
paste0(sl$`Execution Time(sec.)`, collapse = "\n") %>% cat
sl %>% select(1:10)


pto <- googlesheets4::read_sheet(
  "https://docs.google.com/spreadsheets/d/1SPTQvgmw2-1otZuwpGMB1R57qnmVMOzKm2IGh63P7UY/edit#gid=277201005",
  sheet = "JE_submissions"
) %>% select(id, type)
pto %>% filter(duplicated(id))


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
pm_scores_g <- bind_rows(
  pm_scores_g,
  pm_scores_g %>%
    filter(metric != "Overall") %>%
    group_by(`Team Name`, metric) %>%
    summarise(value = min(value, na.rm = TRUE)) %>%
    summarise(value = sum(value)/4, metric = "Best")
)

pmzz <- pm_scores_g %>%
  group_by(metric, `Team Name`) %>%
  summarise(i = which.min(value), value = value[[i]], id = id[[i]]) %>%
  arrange((value)) %>%
  mutate(position = row_number()) %>%
  slice(1:10) %>%
  ungroup() %>%
  mutate(str = paste0("#", position, " - ", `Team Name`, ifelse(is.na(id), "", paste0("\\nid: ", id))))
unique(pmzz$id)

pm_met_lab <- c(
  Overall = "Mean RMSE",
  ADT2GEX = "RMSE for ADT2GEX",
  GEX2ADT = "RMSE for GEX2ADT",
  ATAC2GEX = "RMSE for ATAC2GEX",
  GEX2ATAC = "RMSE for GEX2ATAC",
  Best = "Mean of best RMSE"
)

pms <- map(names(pm_met_lab), function(met) {
  pmzz_ <- pmzz %>% filter(metric == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))
  min <- min(pmzz_$value)
  max <- max(pmzz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .25
  ggplot(pmzz_) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[1]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[1]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks.x = element_line()) +
    labs(y = pm_met_lab[[met]], x = NULL) +
    scale_x_discrete(labels = function(x){sub("\\\\n", "\n", x)})
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

mm_scores_g <- mm_scores %>%
  gather(metric, value, ADT2GEX, GEX2ADT, ATAC2GEX, GEX2ATAC, Overall) %>%
  filter(is.finite(value))

mm_scores_g <- bind_rows(
  mm_scores_g,
  mm_scores_g %>%
    filter(metric != "Overall") %>%
    group_by(`Team Name`, metric) %>%
    summarise(value = max(value, na.rm = TRUE)) %>%
    summarise(value = sum(value)/4, metric = "Best")
)

mmzz <- mm_scores_g %>%
  group_by(metric, `Team Name`) %>%
  summarise(i = which.max(value), value = value[[i]], id = id[[i]]) %>%
  arrange(desc(value)) %>%
  mutate(position = row_number()) %>%
  slice(1:10) %>%
  ungroup() %>%
  mutate(str = paste0("#", position, " - ", `Team Name`, ifelse(is.na(id), "", paste0("\\nid: ", id))))
unique(mmzz$id)

mm_met_lab <- c(
  Overall = "Mean match probability",
  ADT2GEX = "Match probability for ADT2GEX",
  GEX2ADT = "Match probability for GEX2ADT",
  ATAC2GEX = "Match probability for ATAC2GEX",
  GEX2ATAC = "Match probability for GEX2ATAC",
  Best = "Mean best match probability"
)

mms <- map(names(mm_met_lab), function(met) {
  mmzz_ <- mmzz %>% filter(metric == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))
  min <- min(mmzz_$value)
  max <- max(mmzz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .25
  # nud <- max(mmzz_$value) * .5
  ggplot(mmzz_) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[2]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[2]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    # expand_limits(y = c(.7, .85)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks.x = element_line()) +
    labs(y = mm_met_lab[[met]], x = NULL) +
    scale_x_discrete(labels = function(x){sub("\\\\n", "\n", x)})
})







# Joint embedding ---------------------------------------------------------

je_scores <- df %>%
  filter(`Challenge Phase` == "Joint Embedding - Phase 2", Status == "finished") %>%
  unnest(scores) %>%
  select(
    -`Team Members Email Id`, -`Team Members Affiliaton`, -`Challenge Phase`, -`Submission Number`,
    -`Submitted File`, -`Stdout File`, -`Stderr File`, -dest_zip, -dest_json
  ) %>%
  left_join(pto, by = "id") %>%
  mutate(type = ifelse(is.na(type), "", type))

je_scores %>%
  arrange(desc(arithmetic_mean_ATAC)) %>%
  select(1:6, type, contains("arithmetic_mean"))

je_scores %>% select(id, Status, ends_with("_ADT"), ends_with("_ATAC"), arithmetic_mean)

je_scores_g <- je_scores %>% gather(metric, value, ends_with("_ADT"), ends_with("_ATAC"), arithmetic_mean)

je_scores_g <- bind_rows(
  je_scores_g,
  je_scores_g %>%
    filter(grepl("arithmetic_mean_", metric)) %>%
    group_by(`Team Name`, metric) %>%
    summarise(value = max(value, na.rm = TRUE)) %>%
    summarise(value = sum(value)/2, metric = "arithmetic_mean_best")
)

jezz <- je_scores_g %>%
  filter(grepl("arithmetic_mean", metric)) %>%
  mutate(metric2 = ifelse(metric == "arithmetic_mean", metric, paste0(metric, "_", type))) %>%
  group_by(metric2, `Team Name`) %>%
  summarise(i = which.max(value), value = value[[i]], id = id[[i]], type = type[[i]], metric = metric[[i]]) %>%
  arrange(desc(value)) %>%
  mutate(position = row_number()) %>%
  slice(1:10) %>%
  ungroup() %>%
  mutate(str = paste0("#", position, " - ", `Team Name`, ifelse(is.na(id), "", paste0("\\nid: ", id))))
unique(jezz$id)

table(jezz$metric2)

je_met_lab <- c(
  arithmetic_mean = "Mean JE mean",
  arithmetic_mean_ADT_online = "JE mean for ADT (online only)",
  arithmetic_mean_ADT_pretrained = "JE mean for ADT (pre-trained only)",
  arithmetic_mean_ATAC_online = "JE mean for ATAC (online only)",
  arithmetic_mean_ATAC_pretrained = "JE mean for ATAC (pre-trained only)",
  arithmetic_mean_ADT_ = "JE mean for ADT",
  arithmetic_mean_ATAC_ = "JE mean for ATAC",
  arithmetic_mean_best_NA = "Mean of best JE mean"
)

jes <- map(names(je_met_lab), function(met) {
  jezz_ <- jezz %>% filter(metric2 == met) %>% mutate(str = forcats::fct_rev(forcats::fct_inorder(str)))
  min <- min(jezz_$value)
  max <- max(jezz_$value)
  ori <- max(min-(max-min)*2, 0)
  nud <- (max - ori) * .25
  ggplot(jezz_) +
    geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[3]]) +
    geom_point(aes(str, value), stat = "identity", colour = pal[[3]]) +
    geom_text(aes(str, value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = nud, size = 3) +
    coord_flip() +
    # expand_limits(y = c(.7, .85)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks.x = element_line()) +
    labs(y = je_met_lab[[met]], x = NULL) +
    scale_x_discrete(labels = function(x) gsub("\\\\n", "\n", x))
})

jes[[6]]
jes[[7]]

pm_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Predict Modality"), fontface = "bold") + ggplot2::theme_void()
mm_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Match Modality"), fontface = "bold") + ggplot2::theme_void()
je_title <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = "Joint Embedding"), fontface = "bold") + ggplot2::theme_void()
date <- ggplot(tibble(a = 1)) + geom_text(aes(a, a, label = format(Sys.time(), "OpenProblems-NeurIPS2021 Leaderboard on %d %h %Y at %H:%M")), fontface = "bold") + ggplot2::theme_void()
em <- patchwork::plot_spacer()

layout <- "
AAAAA
BBBBB
CDEFG
HHHHH
IJKLM
NNNNN
OPQRS
"

# layout <- "
# AAAAAAA
# BBBBB##
# CDEFG#T
# HHHHH##
# IJKLM#U
# NNNNN##
# OPQRSWV
# "




g <-
  date +
  pm_title + pms[[6]] + pms[[2]] + pms[[3]] + pms[[4]] + pms[[5]] +
  mm_title + mms[[6]] + mms[[2]] + mms[[3]] + mms[[4]] + mms[[5]] +
  je_title + jes[[8]] + jes[[2]] + jes[[3]] + jes[[4]] + jes[[5]] +
  patchwork::plot_layout(
    design = layout,
    heights = c(.1, .1, 1, .1, 1, .1, 1),
    widths = c(1, 1, 1, 1, 1, .5, 1)
  )




filename <- format(Sys.Date(), "output/leaderboard_%Y-%m-%d.pdf")
ggsave(filename, g, width = 4*length(pms), height = 14)
filename <- format(Sys.Date(), "output/leaderboard_%Y-%m-%d.png")
ggsave(filename, g, width = 4*length(pms), height = 14)






layout <- "
AAAAAAA
BBBBB##
CDEFG#T
HHHHH##
IJKLM#U
NNNNN##
OPQRSWV
"

g <-
  date +
  pm_title + pms[[1]] + pms[[2]] + pms[[3]] + pms[[4]] + pms[[5]] +
  mm_title + mms[[1]] + mms[[2]] + mms[[3]] + mms[[4]] + mms[[5]] +
  je_title + jes[[1]] + jes[[2]] + jes[[3]] + jes[[4]] + jes[[5]] +
  pms[[6]] + mms[[6]] + jes[[8]] +
  patchwork::plot_layout(
    design = layout,
    heights = c(.1, .1, 1, .1, 1, .1, 1),
    widths = c(1, 1, 1, 1, 1, .5, 1)
  )




filename <- format(Sys.Date(), "output/leaderboard_%Y-%m-%d_all.pdf")
ggsave(filename, g, width = 4*length(pms), height = 14)








#
# library(lubridate)
#
# met <- "Overall"
# start <- lubridate::floor_date(min(pm_scores_g$submission_time), "hour")
# end <- lubridate::ceiling_date(max(pm_scores_g$submission_time), "hour")
#
# timeseq <- seq(start, end, by = "1 hour")
# timeline <- map_df(seq_along(timeseq), function(i) {
#   time <- timeseq[[i]]
#   pm_scores_g %>%
#     filter(submission_time <= time, value < 1000) %>%
#     group_by(team_name, metric) %>%
#     arrange(value) %>%
#     slice(1) %>%
#     ungroup() %>%
#     group_by(metric) %>%
#     arrange(value) %>%
#     mutate(rank = row_number(), time = time, ti = i) %>%
#     ungroup()
# }) %>%
#   arrange(metric, time, rank) %>%
#   filter(rank <= 10)
#
#
# timeline %>% select(time, rank, team_name, metric)
#
# min <- min(timeline$value)
# max <- max(timeline$value)
# ori <- max(min-(max-min)*2, 0)
# nud <- (max - ori) * .02
# g <- ggplot(timeline, aes(rank, group = team_name, fill = team_name, colour = team_name)) +
#   geom_bar(aes(y = value), colour = NA, stat = "identity") +
#   # geom_segment(aes(x = str, xend = str, y = ori, yend = value), colour = pal[[1]]) +
#   # geom_point(aes(str, value), stat = "identity", colour = pal[[1]]) +
#   geom_text(aes(y = value, label = team_name), hjust = 0, nudge_y = nud, size = 3) +
#   geom_text(aes(y = value, label = sprintf("%.4f", value)), hjust = 1, nudge_y = -nud, size = 3, colour = "white") +
#   facet_wrap(~metric) +
#   coord_flip() +
#   theme_minimal() +
#   theme(panel.grid = element_blank(), legend.position = "none") +
#   labs(y = pm_met_lab[[met]], x = NULL) +
#   scale_x_reverse()
# g
#
# library(gganimate)
# anim <- g + transition_states(ti, transition_length = 4, state_length = 1) +
#   view_follow(fixed_x = TRUE)# +
#   # labs(title = 'GDP per Year : {time}',
#   #      subtitle  =  "Top 10 Countries",
#   #      caption  = "GDP in Billions USD | Data Source: World Bank Data")
# animate(anim, 200, fps = 20,  width = 1200, height = 1000,
#         renderer = gifski_renderer("gganim.gif"))
#
# animate(anim, 200, fps = 20,  width = 1200, height = 1000,
#         renderer = ffmpeg_renderer()) -> for_mp4anim_save("animation.mp4", animation = for_mp4 )

df %>%
  filter(`Challenge Phase` == "Match Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>% arrange(desc(Overall)) %>% write_csv("output/submissions_mm.csv")
df %>%
  filter(`Challenge Phase` == "Joint Embedding - Phase 2", Status == "finished") %>%
  unnest(scores) %>% arrange(desc(arithmetic_mean)) %>% write_csv("output/submissions_je.csv")
df %>%
  filter(`Challenge Phase` == "Predict Modality - Phase 2", Status == "finished") %>%
  unnest(scores) %>% arrange(Overall) %>% write_csv("output/submissions_pm.csv")


# jeyy <- df %>%
#   filter(`Challenge Phase` == "Joint Embedding - Phase 2", Status == "finished") %>%
#   select(id, team = `Team Name`, scores) %>%
#   unnest(scores) %>%
#   filter(arithmetic_mean > 0)
# jeyy1 <- jeyy %>% dplyr::select(id, team) %>% left_join(jezz %>% filter(metric == "arithmetic_mean") %>% select(team = `Team Name`, position), by = "team")
# jeyy2 <- jeyy %>% dplyr::select(-id, -team) %>% mutate_all(function(x) ifelse(is.na(x), 0, x)) %>% as.matrix %>% dyndimred::dimred_mds(distance_method = "pearson")
# jeyy3 <- data.frame(jeyy1, jeyy2) %>% as_tibble
# ggplot(jeyy3) +
#   geom_point(aes(comp_1, comp_2, colour = team)) +
#   geom_text(aes(comp_1, comp_2, label = stringr::str_sub(team, 0, 4), colour = team), nudge_y = .005) +
#   # ggrepel::geom_label_repel(aes(comp_1, comp_2, label = stringr::str_sub(team, 0, 4), colour = team), max.overlaps = 50, size = 2) +
#   theme_bw()

