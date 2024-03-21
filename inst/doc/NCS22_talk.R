## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(designit)
library(tidyverse)

## ----get_data, include = TRUE-------------------------------------------------
data("longitudinal_subject_samples")

dat <- longitudinal_subject_samples |> 
  filter(Group %in% 1:5, Week %in% c(1, 4)) |>
  select(SampleID, SubjectID, Group, Sex, Week)

# for simplicity: remove two subjects that don't have both visits
dat <- dat |>
  filter(SubjectID %in%
    (dat |> count(SubjectID) |> filter(n == 2) |> pull(SubjectID)))


subject_data <- dat |>
  select(SubjectID, Group, Sex) |>
  unique()

## ----fig.width= 4, fig.height=3, echo = FALSE---------------------------------
data("plate_effect_example")
plate_effect_example |>
  ggplot() +
  aes(x = column, y = row, fill = treatment, alpha = log_conc) +
  geom_tile() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_discrete(limits = rev) +
  scale_fill_brewer(palette = "Set1") +
  # make transparency more visible
  scale_alpha_continuous(range = c(0.2, 1)) +
  ggtitle("Design")

## ----fig.width= 4, fig.height=5, echo = FALSE---------------------------------
p1 <- plate_effect_example |>
  ggplot() +
  aes(x = column, y = row, fill = readout) +
  geom_tile() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis_c() +
  ggtitle("Readout")

p2 <- plate_effect_example |>
  filter(treatment == "control") |>
  mutate(column = as.numeric(column)) |>
  ggplot() +
  aes(x = column, y = readout, color = row) +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  ggtitle("Control")

cowplot::plot_grid(p1, p2, nrow = 2)

## ----echo = FALSE-------------------------------------------------------------
set.seed(17) # gives `bad` random assignment

bc <- BatchContainer$new(
  dimensions = list("batch" = 3, "location" = 11)
) |>
  assign_random(subject_data)

## ----fig.width= 3, fig.height=3, echo = FALSE---------------------------------
bc$get_samples() |>
  ggplot(aes(x = batch, fill = Group)) +
  geom_bar() +
  labs(y = "subject count")

## ----include=FALSE------------------------------------------------------------
set.seed(17) # gives `bad` random assignment

## -----------------------------------------------------------------------------
bc <- BatchContainer$new(
  dimensions = list("batch" = 3, "location" = 11)
) |>
  assign_random(subject_data)

## ----fig.width= 5.5, fig.height=3, echo = FALSE-------------------------------
cowplot::plot_grid(
  plotlist = list(
    bc$get_samples() |>
      ggplot(aes(x = batch, fill = Group)) +
      geom_bar() +
      labs(y = "subject count"),
    bc$get_samples() |>
      ggplot(aes(x = batch, fill = Sex)) +
      geom_bar() +
      labs(y = "subject count")
  ),
  nrow = 1
)

## ----eval = FALSE-------------------------------------------------------------
#  bc$get_samples()

## ----echo=FALSE---------------------------------------------------------------
bind_rows(
  head(bc$get_samples(), 3) |>
    mutate(across(everything(), as.character)),
  tibble(
    batch = "...",
    location = " ...",
    SubjectID = "...",
    Group = "...", Sex = "..."
  ),
  tail(bc$get_samples(), 3) |>
    mutate(across(everything(), as.character))
) |>
  gt::gt() |>
  gt::tab_options(
    table.font.size = 11,
    data_row.padding = 0.1
  )

## ----warning=FALSE------------------------------------------------------------
bc <- optimize_design(
  bc,
  scoring = list(
    group = osat_score_generator(
      batch_vars = "batch",
      feature_vars = "Group"
    ),
    sex = osat_score_generator(
      batch_vars = "batch",
      feature_vars = "Sex"
    )
  ),
  n_shuffle = 1,
  acceptance_func =
    ~ accept_leftmost_improvement(..., tolerance = 0.01),
  max_iter = 150,
  quiet = TRUE
)

## ----fig.width= 8, fig.height=3, echo = FALSE---------------------------------
cowplot::plot_grid(
  plotlist = list(
    bc$get_samples() |>
      ggplot(aes(x = batch, fill = Group)) +
      geom_bar() +
      labs(y = "subject count"),
    bc$get_samples() |>
      ggplot(aes(x = batch, fill = Sex)) +
      geom_bar() +
      labs(y = "subject count"),
    bc$plot_trace(include_aggregated = TRUE)
  ),
  ncol = 3
)

## ----echo=FALSE---------------------------------------------------------------
bind_rows(
  head(bc$get_samples(), 3) |>
    mutate(across(everything(), as.character)),
  tibble(
    batch = "...",
    location = " ...",
    SubjectID = "...",
    Group = "...", Sex = "..."
  ),
  tail(bc$get_samples(), 3) |>
    mutate(across(everything(), as.character))
) |>
  gt::gt() |>
  gt::tab_options(
    table.font.size = 11,
    data_row.padding = 0.1
  )

## -----------------------------------------------------------------------------
set.seed(4)

bc <- BatchContainer$new(
  dimensions = list("plate" = 3, "row" = 4, "col" = 6)
) |>
  assign_in_order(dat)

## ----fig.width= 5, fig.height=4.5, eval=FALSE---------------------------------
#  plot_plate(bc,
#    plate = plate, row = row, column = col,
#    .color = Group, title = "Initial layout by Group"
#  )
#  plot_plate(bc,
#    plate = plate, row = row, column = col,
#    .color = Sex, title = "Initial layout by Sex"
#  )

## ----fig.width= 5, fig.height=4.5, echo=FALSE---------------------------------
cowplot::plot_grid(
  plotlist = list(
    plot_plate(bc,
      plate = plate, row = row, column = col,
      .color = Group, title = "Initial layout by Group"
    ),
    plot_plate(bc,
      plate = plate, row = row, column = col,
      .color = Sex, title = "Initial layout by Sex"
    )
  ),
  nrow = 2
)

## ----warning=FALSE------------------------------------------------------------
bc1 <- optimize_design(
  bc,
  scoring = list(
    group = osat_score_generator(
      batch_vars = "plate",
      feature_vars = "Group"
    ),
    sex = osat_score_generator(
      batch_vars = "plate",
      feature_vars = "Sex"
    )
  ),
  n_shuffle = 1,
  acceptance_func =
    ~ accept_leftmost_improvement(..., tolerance = 0.01),
  max_iter = 150,
  quiet = TRUE
)

## ----fig.width= 5, fig.height=4.5, echo=FALSE---------------------------------
cowplot::plot_grid(
  plotlist = list(
    plot_plate(bc1,
      plate = plate, row = row, column = col,
      .color = Group, title = "Layout after the first step, Group"
    ),
    plot_plate(bc1,
      plate = plate, row = row, column = col,
      .color = Sex, title = "Layout after the first step, Sex"
    )
  ),
  nrow = 2
)

## ----warning=FALSE------------------------------------------------------------
bc2 <- optimize_design(
  bc1,
  scoring = mk_plate_scoring_functions(
    bc1,
    plate = "plate", row = "row", column = "col",
    group = "Group"
  ),
  shuffle_proposal_func = shuffle_with_constraints(dst = plate == .src$plate),
  max_iter = 150,
  quiet = TRUE
)

## ----fig.width= 5, fig.height=4.5, echo=FALSE---------------------------------
cowplot::plot_grid(
  plotlist = list(
    plot_plate(bc2,
      plate = plate, row = row, column = col,
      .color = Group, title = "Layout after the second step, Group"
    ),
    plot_plate(bc2,
      plate = plate, row = row, column = col,
      .color = Sex, title = "Layout after the second step, Sex"
    )
  ),
  nrow = 2
)

## ----warning=FALSE, message=FALSE---------------------------------------------
bc <- optimize_multi_plate_design(
  bc,
  across_plates_variables = c("Group", "Sex"),
  within_plate_variables = c("Group"),
  plate = "plate", row = "row", column = "col",
  n_shuffle = 2,
  max_iter = 500 # 2000
)

## ----fig.width= 5, fig.height=4.5, echo=FALSE---------------------------------
cowplot::plot_grid(
  plotlist = list(
    plot_plate(bc,
      plate = plate, row = row, column = col,
      .color = Group, title = "After optimization, Group"
    ),
    plot_plate(bc,
      plate = plate, row = row, column = col,
      .color = Sex, title = "After optimization, Sex"
    )
  ),
  nrow = 2
)

## ----fig.width=5, fig.height=4, echo=FALSE------------------------------------
bc$plot_trace()

## ----fig.width=4.0, fig.hight = 5.0, echo = FALSE-----------------------------
layout <- crossing(row = 1:9, column = 1:12) |>
  mutate(Questions = "no")
layout$Questions[c(
  16, 17, 18, 19, 20, 21,
  27, 28, 33, 34,
  45, 46,
  55, 56, 66, 67, 90, 91
)] <- "yes"

plot_plate(layout, .color = Questions, title = "Thank you")

