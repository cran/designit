params <-
list(iterations = 100L, n_samples = 50L)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----setup--------------------------------------------------------------------
library(designit)
library(tidyverse)

## -----------------------------------------------------------------------------
set.seed(43)
samples <- tibble(
  id = 1:params$n_samples,
  sex = sample(c("F", "M"), params$n_samples, replace = TRUE, prob = c(0.8, 0.2)),
  group = sample(c("treatment", "control"), params$n_samples, replace = TRUE),
  weight = runif(params$n_samples, 20, 30)
)
samples |>
  head()

## -----------------------------------------------------------------------------
samples |>
  count(sex)

## -----------------------------------------------------------------------------
set.seed(42)
bc <- BatchContainer$new(
  dimensions = c("cage" = 11, "position" = 5)
) |>
  assign_random(samples)
bc

## -----------------------------------------------------------------------------
males_per_cage <- function(bc) {
  bc$get_samples() |>
    filter(sex == "M") |>
    count(cage) |>
    ggplot(aes(cage, n)) +
    geom_col()
}

weight_d <- function(bc) {
  bc$get_samples() |>
    ggplot(aes(factor(cage), weight)) +
    geom_violin() +
    geom_point() +
    stat_summary(fun = mean, geom = "point", size = 2, shape = 23, color = "red")
}

group_d <- function(bc) {
  bc$get_samples(remove_empty_locations = TRUE) |>
    ggplot(aes(factor(cage), fill = group)) +
    geom_bar(position = "fill")
}

## -----------------------------------------------------------------------------
males_per_cage(bc)
weight_d(bc)
group_d(bc)

## -----------------------------------------------------------------------------
set.seed(10)

bc <- optimize_design(
  bc,
  scoring = osat_score_generator(
    "cage",
    "sex"
  ),
  shuffle_proposal_func = shuffle_with_constraints(
    sex == "M",
    cage != .src$cage
  ),
  max_iter = 10
)

bc$plot_trace()

## -----------------------------------------------------------------------------
males_per_cage(bc)
weight_d(bc)
group_d(bc)

## -----------------------------------------------------------------------------
scoring_f <- function(bc) {
  samples <- bc$get_samples(include_id = TRUE, as_tibble = FALSE)
  avg_w <- samples[, mean(weight, na.rm = TRUE)]
  avg_w_per_cage <- samples[!is.na(weight), mean(weight), by = cage]$V1
  trt_per_cage <- samples[!is.na(group), sum(group == "treatment") / .N, by = cage]$V1

  w_score <- mean((avg_w - avg_w_per_cage)**2)
  trt_score <- mean((trt_per_cage - 0.5)**2)
  w_score + 10 * trt_score
}

set.seed(12)
bc <- optimize_design(bc,
  scoring = scoring_f,
  shuffle_proposal = shuffle_with_constraints(
    sex == "F",
    cage != .src$cage & (is.na(sex) | sex != "M")
  ),
  n_shuffle = c(rep(10, 20), rep(5, 20), rep(3, 20), rep(1, 140)),
  max_iter = 200
)
bc$plot_trace()
scoring_f(bc)

## -----------------------------------------------------------------------------
males_per_cage(bc)
weight_d(bc)
group_d(bc)

