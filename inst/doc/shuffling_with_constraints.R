## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE,
  fig.width = 6,
  fig.height = 6
)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(designit)
library(tidyverse)

## ----echo=TRUE----------------------------------------------------------------
data("invivo_study_samples")

invivo_study_samples <- dplyr::mutate(invivo_study_samples,
  Litter_combine_females = ifelse(Sex == "F", "female_all", Litter)
)
str(invivo_study_samples)

invivo_study_samples |>
  dplyr::count(Strain, Sex, Litter_combine_females) |>
  gt::gt()

## ----echo=TRUE----------------------------------------------------------------
treatments <- factor(rep(c("Treatment A", "Treatment B"), c(30, 29)))
table(treatments)

bc <- BatchContainer$new(locations_table = data.frame(Treatment = treatments, Position = seq_along(treatments)))

bc <- assign_in_order(bc, invivo_study_samples)

scoring_f <- osat_score_generator(batch_vars = "Treatment", feature_vars = c("Strain", "Sex"))

bc

## -----------------------------------------------------------------------------
bc2 <- optimize_design(
  bc,
  scoring = scoring_f,
  shuffle_proposal_func = shuffle_grouped_data(bc,
    allocate_var = "Treatment",
    keep_together_vars = c("Strain", "Sex"),
    keep_separate_vars = c("Earmark"),
    subgroup_var_name = "Cage",
    n_min = 2, n_ideal = 3, n_max = 5,
    strict = TRUE,
    report_grouping_as_attribute = TRUE
  ),
  max_iter = 600
)

design <- bc2$get_samples()

## ----echo=TRUE----------------------------------------------------------------
design |>
  dplyr::count(Cage, Strain, Sex, Treatment) |>
  gt::gt()

## ----echo=TRUE----------------------------------------------------------------
subg <- form_homogeneous_subgroups(
  batch_container = bc, allocate_var = "Treatment",
  keep_together_vars = c("Strain", "Sex", "Litter_combine_females"),
  subgroup_var_name = "Cage",
  n_min = 2, n_ideal = 3, n_max = 5
)

## -----------------------------------------------------------------------------
subg$Subgroup_Sizes

## -----------------------------------------------------------------------------
possible <- compile_possible_subgroup_allocation(subg)

## -----------------------------------------------------------------------------
shuffle_proposal <- shuffle_with_subgroup_formation(subg, possible, report_grouping_as_attribute = TRUE)

shuffle_proposal()

## ----echo=TRUE----------------------------------------------------------------
bc3 <- optimize_design(
  bc,
  scoring = scoring_f,
  shuffle_proposal_func = shuffle_proposal,
  max_iter = 300
)

design <- bc3$get_samples()

# Obeying all constraints does not lead to a very balanced sample allocation:
dplyr::count(design, Treatment, Strain) |> gt::gt()

dplyr::count(design, Treatment, Sex) |> gt::gt()

