## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(designit)
library(tidyverse)

## -----------------------------------------------------------------------------
data("multi_trt_day_samples")

## -----------------------------------------------------------------------------
multi_trt_day_samples |>
  count(Time, Treatment) |>
  gt::gt()

## -----------------------------------------------------------------------------
# Setting up the batch container
bc <- BatchContainer$new(
  dimensions = c(
    batch = ceiling(nrow(multi_trt_day_samples) / 8),
    run = 2, position = 4
  )
)

# Initial random assignment
bc <- assign_in_order(bc, multi_trt_day_samples)
bc

## -----------------------------------------------------------------------------
n_shuffle <- rep(c(32, 10, 2), c(100, 80, 20))
n_iterations <- length(n_shuffle)

set.seed(42) # should we have conventions for this?

scoring_f <- osat_score_generator(c("batch"), c("Treatment", "Time"))
bc <- optimize_design(
  bc,
  scoring = scoring_f,
  n_shuffle = n_shuffle,
  max_iter = n_iterations
) # default is 10000

## ----fig.width=5, fig.height= 4-----------------------------------------------
qplot(
  x = bc$trace$scores[[1]]$step,
  y = bc$trace$scores[[1]]$score_1,
  color = factor(c(32, n_shuffle)),
  main = str_glue("Final score={bc$score(scoring_f)}"), geom = "point"
)

## ----fig.width=6, fig.height=5------------------------------------------------
bc$get_samples(assignment = TRUE) |>
  mutate(batch = factor(batch)) |>
  ggplot(aes(x = batch, fill = Treatment, alpha = factor(Time))) +
  geom_bar()

## ----fig.width=6, fig.height= 4-----------------------------------------------
# copy batch container for second optimization
bc2 <- assign_in_order(bc)

n_iterations <- 200

set.seed(42) # should we have conventions for this?

bc2 <- optimize_design(
  bc2,
  scoring = scoring_f,
  shuffle_proposal = shuffle_with_constraints(
    src = TRUE,
    # batch needs to change for shuffle to be accepted
    dst = .src$batch != batch
  ),
  max_iter = n_iterations
)

qplot(
  x = bc2$trace$scores[[1]]$step,
  y = bc2$trace$scores[[1]]$score_1,
  main = str_glue("Final score={bc2$score(scoring_f)}"), geom = "point"
)

bc2$get_samples(assignment = TRUE) |>
  mutate(batch = factor(batch)) |>
  ggplot(aes(x = batch, fill = Treatment, alpha = factor(Time))) +
  geom_bar()

## -----------------------------------------------------------------------------
n_iterations <- 100

# new optimization function
scoring_f <- osat_score_generator(c("run"), c("Treatment", "Time"))
# like this the optimization score is wrong because it tries to optimize across Batches.
# Possible ways to go:
# - we'd need something like c("batch", batch/run") for optimize by batch and run within batch.
# - or we add "batch/run" to the constraints somehow.
bc$score(scoring_f)

bc <- optimize_design(
  bc,
  scoring = scoring_f,
  shuffle_proposal = shuffle_with_constraints(
    src = TRUE,
    # batch remains the same and run needs to change
    dst = batch == .src$batch & run != .src$run
  ),
  max_iter = n_iterations
)

## ----fig.width=6, fig.height= 4-----------------------------------------------
qplot(
  x = bc$trace$scores[[1]]$step,
  y = bc$trace$scores[[1]]$score_1,
  color = factor(n_iterations),
  main = str_glue("Final score={bc$score(scoring_f)}"), geom = "point"
)

## ----fig.width=6, fig.height= 4-----------------------------------------------
bc$get_samples() |>
  mutate(run = factor(run)) |>
  ggplot(aes(x = run, fill = Treatment, alpha = factor(Time))) +
  geom_bar() +
  facet_wrap(~batch)

