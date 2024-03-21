## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(designit)

## -----------------------------------------------------------------------------
data("multi_trt_day_samples")

## -----------------------------------------------------------------------------
multi_trt_day_samples |>
  dplyr::count(Time, Treatment) |>
  gt::gt()

## -----------------------------------------------------------------------------
# Setting up the batch container
bc <- BatchContainer$new(
  dimensions = c(
    batch = ceiling(nrow(multi_trt_day_samples) / 8),
    run = 2, position = 5
  ),
  exclude = tibble::tibble(batch = 4, run = c(1, 2), position = c(5, 5))
) |>
  # Add samples to container
  assign_in_order(samples = multi_trt_day_samples)

bc

## -----------------------------------------------------------------------------
n_shuffle <- rep(c(32, 10, 5, 2, 1), c(20, 40, 40, 50, 50))

scoring_f <- osat_score_generator(c("batch"), c("Treatment", "Time"))

bc1 <- optimize_design(
  bc,
  scoring = scoring_f,
  n_shuffle = n_shuffle # will implicitly generate a shuffling function according to the provided schedule
)

bc1$trace$elapsed

## ----fig.width=5, fig.height= 4-----------------------------------------------
bc1$scores_table() |>
  dplyr::mutate(
    n_shuffle = c(NA, n_shuffle)
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(step, value, color = factor(n_shuffle))
  ) +
  ggplot2::geom_point() +
  ggplot2::labs(
    title = "Score 1 tracing",
    subtitle = stringr::str_glue("Final score = {bc1$score(scoring_f)}"),
    x = "Iteration",
    y = "Score",
    color = "n_shuffle"
  )

## ----fig.width=5, fig.height= 4-----------------------------------------------
bc1$plot_trace()

## ----fig.width=6, fig.height=5------------------------------------------------
bc1$score(scoring_f)

bc1$get_samples(assignment = TRUE) |>
  dplyr::filter(!is.na(Treatment)) |>
  dplyr::mutate(anno = stringr::str_c(Time, " hr")) |>
  ggplot2::ggplot(ggplot2::aes(x = batch, y = interaction(position, run), fill = Treatment)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_hline(yintercept = 5.5, size = 1) +
  ggplot2::geom_text(ggplot2::aes(label = anno)) +
  ggplot2::labs(x = "Batch", y = "Position . Run")

## -----------------------------------------------------------------------------
n_shuffle <- rep(c(5, 2, 1), c(30, 30, 30))

bc1 <- optimize_design(
  bc1,
  scoring = scoring_f,
  n_shuffle = n_shuffle
)

## -----------------------------------------------------------------------------
bc2 <- optimize_design(
  bc,
  scoring = scoring_f,
  n_shuffle = 3, # will implicitly generate a shuffling function that will do 3 swaps at each iteration
  max_iter = 2000,
  min_delta = 0.1
)

## -----------------------------------------------------------------------------
multi_scoring_f <- list(
  osat_score_generator(c("batch"), c("Treatment", "Time")),
  osat_score_generator(c("batch"), c("Treatment"))
)


bc3 <- optimize_design(
  bc,
  scoring = multi_scoring_f,
  n_shuffle = 3,
  max_iter = 200,
  min_delta = 0.1
)

## -----------------------------------------------------------------------------
bc3_as <- optimize_design(
  bc,
  scoring = multi_scoring_f,
  n_shuffle = 3,
  max_iter = 200,
  min_delta = 0.01,
  autoscale_scores = T,
  autoscaling_permutations = 200
)

## -----------------------------------------------------------------------------
bc4 <- optimize_design(
  bc,
  scoring = multi_scoring_f,
  n_shuffle = 3,
  aggregate_scores_func = worst_score,
  max_iter = 200,
  autoscale_scores = TRUE,
  autoscaling_permutations = 200
)

## ----eval = FALSE-------------------------------------------------------------
#  bc5 <- optimize_design(
#    bc,
#    scoring = multi_scoring_f,
#    aggregate_scores_func = sum_scores,
#    max_iter = 200,
#    autoscale_scores = TRUE,
#    autoscaling_permutations = 200
#  )

## ----eval = FALSE-------------------------------------------------------------
#  bc5_2 <- optimize_design(
#    bc,
#    scoring = multi_scoring_f,
#    aggregate_scores_func = L2s_norm,
#    max_iter = 200,
#  )

## -----------------------------------------------------------------------------
bc6 <- optimize_design(
  bc,
  scoring = scoring_f,
  shuffle_proposal_func = complete_random_shuffling,
  max_iter = 200
)

## -----------------------------------------------------------------------------
bc7 <- optimize_design(
  bc,
  scoring = scoring_f,
  n_shuffle = 1,
  acceptance_func = mk_simanneal_acceptance_func(),
  max_iter = 200
)

## ----fig.width=5, fig.height= 4-----------------------------------------------
bc7$plot_trace()

## -----------------------------------------------------------------------------
bc8 <- optimize_design(
  bc,
  scoring = scoring_f,
  n_shuffle = 1,
  acceptance_func = mk_simanneal_acceptance_func(mk_simanneal_temp_func(T0 = 100, alpha = 2)),
  max_iter = 150
)

bc8$plot_trace()

## -----------------------------------------------------------------------------
n_shuffle <- rep(c(3, 2, 1), c(20, 20, 200))

bc9 <- optimize_design(
  bc,
  scoring = list(
    osat_score_generator(c("batch"), c("Treatment", "Time")),
    osat_score_generator(c("batch"), c("Treatment")),
    osat_score_generator(c("batch"), c("Time"))
  ),
  n_shuffle = n_shuffle,
  aggregate_scores_func = sum_scores,
  acceptance_func = mk_simanneal_acceptance_func(mk_simanneal_temp_func(T0 = 500, alpha = 1)),
  max_iter = 200,
  min_delta = 1e-8,
  autoscale_scores = T
)

bc9$plot_trace()

bc9$get_samples(assignment = TRUE) |>
  dplyr::mutate(batch = factor(batch)) |>
  ggplot2::ggplot(ggplot2::aes(x = batch, fill = Treatment, alpha = factor(Time))) +
  ggplot2::geom_bar()

