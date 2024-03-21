## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(designit)
library(ggplot2)
library(dplyr)
library(tidyr)

## -----------------------------------------------------------------------------
# conditions to use
conditions <- data.frame(
  group = c(1, 2, 3, 4, 5),
  treatment = c(
    "vehicle", "TRT1", "TRT2",
    "TRT1", "TRT2"
  ),
  dose = c(0, 25, 25, 50, 50)
)

gt::gt(conditions)

## -----------------------------------------------------------------------------
# sample table (2 animals per group with 3 replicates)
n_reps <- 4
n_animals <- 3
animals <- bind_rows(replicate(n_animals, conditions, simplify = FALSE),
  .id = "animal"
)
samples <- bind_rows(replicate(n_reps, animals, simplify = FALSE),
  .id = "replicate"
) |>
  mutate(
    SampleID = paste0(treatment, "_", animal, "_", replicate),
    AnimalID = paste0(treatment, "_", animal)
  ) |>
  mutate(dose = factor(dose))

samples |>
  head(10) |>
  gt::gt()

## -----------------------------------------------------------------------------
n_samp <- nrow(samples)
n_loc_per_plate <- 48 - 4
n_plates <- ceiling(n_samp / n_loc_per_plate)

exclude_wells <- expand.grid(plate = seq(n_plates), column = c(1, 8), row = c(1, 6))

## -----------------------------------------------------------------------------
bc <- BatchContainer$new(
  dimensions = c("plate" = n_plates, "column" = 8, "row" = 6),
  exclude = exclude_wells
)
bc

bc$n_locations
bc$exclude
bc$get_locations() |> head()

## -----------------------------------------------------------------------------
bc <- assign_random(bc, samples)

bc$get_samples()
bc$get_samples(remove_empty_locations = TRUE)

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_plate(bc,
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose
)

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_plate(bc$get_samples(remove_empty_locations = TRUE),
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose
)

## ----fig.width=6, fig.height=3.5----------------------------------------------
bc$move_samples(src = c(1L, 2L), dst = c(2L, 1L))

plot_plate(bc$get_samples(remove_empty_locations = TRUE),
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose
)

## ----fig.width=6, fig.height=3.5----------------------------------------------
bc$move_samples(
  location_assignment = c(
    1:nrow(samples),
    rep(NA, (bc$n_locations - nrow(samples)))
  )
)

plot_plate(bc$get_samples(remove_empty_locations = TRUE, include_id = TRUE),
  plate = plate, column = column, row = row,
  .color = .sample_id
)

## -----------------------------------------------------------------------------
bc <- optimize_design(bc,
  scoring = osat_score_generator(
    batch_vars = "plate",
    feature_vars = c("treatment", "dose")
  ),
  # shuffling schedule
  n_shuffle = c(rep(10, 200), rep(2, 400))
)

## ----fig.width=3.5, fig.height=3----------------------------------------------
bc$plot_trace()

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_plate(bc$get_samples(remove_empty_locations = TRUE),
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose
)

## ----fig.width=4, fig.height=3.5----------------------------------------------
ggplot(
  bc$get_samples(remove_empty_locations = TRUE),
  aes(x = treatment, fill = treatment)
) +
  geom_bar() +
  facet_wrap(~plate)

## ----fig.width=6, fig.height=3.5----------------------------------------------
color_palette <- c(
  TRT1 = "blue", TRT2 = "purple",
  vehicle = "orange", empty = "white"
)

plot_plate(bc,
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose,
  add_excluded = TRUE, rename_empty = TRUE
) +
  scale_fill_manual(values = color_palette, na.value = "darkgray")

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_plate(bc$get_samples(remove_empty_locations = TRUE),
  plate = plate, column = column, row = row,
  .color = treatment, .alpha = dose
) +
  scale_fill_viridis_d()

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_plate(bc$get_samples(remove_empty_locations = TRUE) |>
  filter(column != 2),
plate = plate, column = column, row = row,
.color = treatment, .alpha = dose
) +
  scale_fill_viridis_d()

