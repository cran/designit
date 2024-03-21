## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE,
  fig.width = 6,
  fig.height = 6
)

## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
library(designit)
library(dplyr)
library(assertthat)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)

## ----helper_functions, echo=FALSE---------------------------------------------
# -----------------------------------------------------------------------------
# Helper functions used by the in vivo wrappers

DoE_multipleLevelCols <- function(df, na.rm = FALSE) {
  tmp <- dplyr::summarise(df, across(everything(), ~ dplyr::n_distinct(.x, na.rm = na.rm)))

  colnames(tmp)[tmp > 1]
}

DoE_numericalCols <- function(df, na.rm = FALSE) {
  dplyr::select(df, where(is.numeric)) |>
    colnames()
}

## ----scoring_functions, echo=FALSE--------------------------------------------
# -----------------------------------------------------------------------------
# Dedicated scoring functions used by the in vivo wrappers

# Scoring function to assess more specific constraints of the type Column A == Column B
bc_misaligned <- function(batch_vars, feature_vars, weight = 1, test_data = NULL) {
  force(batch_vars)
  force(feature_vars)
  force(weight)

  assertthat::assert_that(length(batch_vars) > 0, length(batch_vars) == length(feature_vars),
    msg = "Batch and feature variables must be of equal length (>0)"
  )
  if (!is.null(test_data)) {
    batch_match <- match(batch_vars, colnames(test_data))
    feature_match <- match(feature_vars, colnames(test_data))
    assertthat::assert_that(!any(is.na(batch_match)), !any(is.na(feature_match)),
      msg = "All columns must exist in passed sample data"
    )
    assertthat::assert_that(length(intersect(batch_vars, feature_vars)) == 0,
      msg = "There should be no overlap between batch and feature variables"
    )
    # Could do more checks, like matching number of levels
  }


  function(bc) {
    samples <- bc$get_samples(include_id = TRUE, as_tibble = FALSE)
    misaligns <- 0
    for (i in seq_along(batch_vars)) {
      misaligns <- misaligns + sum(samples[[batch_vars[i]]] != samples[[feature_vars[i]]])
    }
    weight * misaligns
  }
}

# Scoring function to assess more specific constraints of the type 'Feature levels unique in Container'
bc_homogenicity_violations <- function(container_var, feature_var, weight = 1, test_data = NULL) {
  force(container_var)
  force(feature_var)
  force(weight)

  assertthat::assert_that(length(container_var) == 1, length(feature_var) == 1,
    msg = "Container and feature variables must be of length 1"
  )
  if (!is.null(test_data)) {
    container_match <- match(container_var, colnames(test_data))
    feature_match <- match(feature_var, colnames(test_data))
    assertthat::assert_that(!any(is.na(container_match)), !any(is.na(feature_match)),
      msg = "All columns must exist in passed sample data"
    )
    assertthat::assert_that(container_var != feature_var,
      msg = "Container and feature variable must be different"
    )
  }

  function(bc) {
    samples <- bc$get_samples(include_id = TRUE, as_tibble = FALSE)
    viol <- table(samples[[container_var]], samples[[feature_var]]) |>
      apply(1, function(a) {
        sum(a > 0) - 1
      }) |>
      sum()
    weight * viol
  }
}

# Scoring function with one-way-ANOVA; Only model with one descriptor variable possible yet
anova_logP <- function(batch_var, feature_var, weight = 1, test_data = NULL) {
  force(batch_var)
  force(feature_var)
  force(weight)

  assertthat::assert_that(length(batch_var) == 1, length(feature_var) == 1,
    msg = "Batch and feature variable must be of length 1"
  )
  if (!is.null(test_data)) {
    assertthat::assert_that(batch_var != feature_var, !any(is.na(match(c(batch_var, feature_var), colnames(test_data)))),
      msg = "All columns must exist in passed sample data and be different from each other"
    )
    assertthat::assert_that(is.numeric(test_data[[feature_var]]), msg = "Feature variable must be numeric")
    # assertthat::assert_that(is.factor(test_data[[batch_var]]) || is.character(test_data[[batch_var]]), msg="Batch variable must be categorical")
  }

  function(bc) {
    samples <- bc$get_samples(include_id = TRUE, as_tibble = TRUE)
    mod <- lm(samples[[feature_var]] ~ factor(samples[[batch_var]]), na.action = na.omit)

    f <- summary(mod)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail = F)
    weight * unname(-log10(p))
  }
}

## ----invivo_functions, echo=FALSE---------------------------------------------
# -----------------------------------------------------------------------------
# User callable functions for helping with the design of in vivo studies

# Assignment of treatments to an animal list
InVivo_assignTreatments <- function(animal_list, treatments,
                                    balance_treatment_vars = c(),
                                    form_cages_by = c(),
                                    n_shuffle = c(rep(5, 100), rep(3, 200), rep(2, 500), rep(1, 20000)),
                                    quiet_process = FALSE, quiet_optimize = TRUE) {
  if (is.vector(treatments)) treatments <- data.frame(Treatment = treatments)

  if (!quiet_process) message("Performing treatment assignment with constrained animal selection.")

  assertthat::assert_that("Treatment" %in% colnames(treatments), msg = "Column 'Treatment' missing from treatment object")

  # Make sure we have defined factor levels for categorical variables
  treatments_sortedfac <- dplyr::mutate(treatments, across(where(is_character), as.factor)) # here we want alphabetic sort order!
  treatments <- dplyr::mutate(treatments, across(where(is_character), function(v) factor(v, levels = unique(v)))) # Keep orig sort order in treatment object!

  # Create a CageGroup variable if not yet provided in animal list or overridden by specific argument
  if (!is.null(form_cages_by) && length(form_cages_by) > 0) {
    animal_list <- tidyr::unite(animal_list, "CageGroup", any_of(form_cages_by), sep = "_", remove = F)
  }

  # Do we have batch container columns that are multi-value and have to align with identically named sample columns?
  # In this case, make those columns factors, check that numbers match and thus a solution is possible at all
  # Append '_bc' to the common variables in the batch container to allow both sets of variables in the same container
  bc_columns <- intersect(colnames(treatments), colnames(animal_list)) |> intersect(DoE_multipleLevelCols(treatments))

  ani_fac <- dplyr::mutate(animal_list, across(all_of(bc_columns), as.factor)) # here we want alphabetic sort order!

  if (length(bc_columns) > 0) {
    if (!quiet_process) {
      message(
        "Using constraints in variables: ", stringr::str_c(bc_columns, collapse = ", "),
        "\nChecking if solution is possible:"
      )
    }

    ani_counts <- dplyr::count(as_tibble(ani_fac), across(all_of(bc_columns))) |> dplyr::arrange(across(all_of(bc_columns)))
    treat_counts <- dplyr::count(as_tibble(treatments_sortedfac), across(all_of(bc_columns))) |> dplyr::arrange(across(all_of(bc_columns)))
    assertthat::assert_that(identical(ani_counts, treat_counts), msg = "Levels of common variables don't match between animal list and treatment list")
    if (!quiet_process) message("   ... Yes!")
    for (i in seq_along(bc_columns)) { # add bc_ prefix to constraint variables to allow both variable sets in batch container
      treatments_sortedfac <- dplyr::rename(treatments_sortedfac, !!stringr::str_c(bc_columns[i], "_bc") := bc_columns[i])
    }
  }

  if (!quiet_process) message("Setting up batch container.")

  # Set up batch container. Now it would be VERY handy to be able to just pass the treatment object directly. ;)
  # Instead, have to set up dimension vector and exclude table to reflect the valid container positions.
  n_lev <- dplyr::summarize(treatments_sortedfac, across(everything(), nlevels))
  all_fac_combs <- dplyr::count(treatments_sortedfac, across(everything()), .drop = F)
  positions <- max(all_fac_combs$n)

  dimension_vector <- c(setNames(as.integer(n_lev), nm = colnames(n_lev)), Position = positions)

  exclude_table <- tidyr::expand_grid(dplyr::select(all_fac_combs, -c("n")), Position = 1:positions) |>
    dplyr::left_join(all_fac_combs, by = setdiff(colnames(all_fac_combs), "n")) |>
    dplyr::filter(Position > n) |>
    dplyr::select(-c("n")) |>
    dplyr::mutate(across(everything(), as.integer))


  # Prepare all constraint columns to hold integer level numbers --> comparable to batch container
  ani_bclevels <- dplyr::mutate(ani_fac, across(all_of(bc_columns), as.integer))

  # Initial batch container, treatments will be assigned in first step
  bc_treatment <- BatchContainer$new(
    dimensions = dimension_vector,
    exclude = exclude_table
  )

  bc_treatment <- assign_in_order(bc_treatment, ani_bclevels)

  if (!quiet_process) message("Constructing scoring functions:")

  # Set up required scoring functions; most relevant ones come first to allow later application of relevance ranking
  scoring_functions <- list()
  bc_data <- bc_treatment$get_samples()

  if (length(bc_columns) > 0) {
    if (!quiet_process) {
      message(
        "     ... user specified treatment allocation constraint (Treatment-",
        stringr::str_c(bc_columns, collapse = "-"), ")"
      )
    }
    scoring_functions <- c(scoring_functions, trt_constraints = bc_misaligned(
      batch_vars = str_c(bc_columns, "_bc"),
      feature_vars = bc_columns,
      weight = 1,
      test_data = bc_data
    ))
  }

  # This constraint is there to restrict number of distinct treatments within the cage groups, in order to find
  # viable solutions if treatment should be homogeneous in every cage
  if ("Treatment" %in% colnames(bc_data) && "CageGroup" %in% colnames(bc_data) &&
    "CageGroup" %in% DoE_multipleLevelCols(bc_data)) {
    if (!quiet_process) message("     ... facilitating homogeneity of treatment in cages (CageGroup)")
    scoring_functions <- c(scoring_functions, cage_trt_homogeneity = bc_homogenicity_violations(
      container_var = "CageGroup",
      feature_var = "Treatment",
      weight = 1,
      test_data = bc_data
    ))
  }

  trt_vars_cat <- intersect(balance_treatment_vars, colnames(bc_data)) |>
    setdiff(bc_columns) |>
    intersect(DoE_multipleLevelCols(bc_data)) |>
    setdiff(DoE_numericalCols(bc_data))
  if (length(trt_vars_cat) > 0) {
    if (!quiet_process) message("     ... OSAT for categorical variables balanced across treatment (", stringr::str_c(trt_vars_cat, collapse = ", "), ")")
    scoring_functions <- c(scoring_functions, balanced_categorical = osat_score_generator(batch_vars = "Treatment", feature_vars = trt_vars_cat))
  }

  trt_vars_num <- intersect(balance_treatment_vars, colnames(bc_data)) |>
    setdiff(bc_columns) |>
    intersect(DoE_multipleLevelCols(bc_data)) |>
    intersect(DoE_numericalCols(bc_data))
  if (length(trt_vars_num) > 0) {
    if (!quiet_process) message("     ... ANOVA -logP for numerical variables balanced across treatment (", stringr::str_c(trt_vars_num, collapse = ", "), ")")
    sf <- list()
    for (i in seq_along(trt_vars_num)) {
      sf <- c(sf, anova_logP(batch_var = "Treatment", feature_var = trt_vars_num[i], test_data = bc_data))
    }
    names(sf) <- stringr::str_c("balanced_", trt_vars_num)
    scoring_functions <- c(scoring_functions, sf)
  }

  assertthat::assert_that(length(scoring_functions) > 0, msg = "No variables for scoring found or all have only one level. Nothing to do.")
  bc_treatment$score(scoring_functions)

  bc_treatment <- optimize_design(
    bc_treatment,
    scoring = scoring_functions,
    n_shuffle = n_shuffle,
    acceptance_func = ~ accept_leftmost_improvement(..., tolerance = 0.1),
    quiet = quiet_optimize
  )

  # Check if user given constraints (if provided) could be satisfied
  if ("trt_constraints" %in% names(bc_treatment$score(scoring_functions))) {
    if (bc_treatment$score(scoring_functions)[["trt_constraints"]] > 0) {
      message("CAUTION: User defined constraints could not be fully met (remaining score ", bc_treatment$score(scoring_functions)[["trt_constraints"]], ")")
    } else {
      if (!quiet_process) message("Success. User provided constraints could be fully met.")
    }
  }

  design_trt <- bc_treatment$get_samples()

  # Translate back integer codes for 'used' factors into meaningful labels and remove intermediate helper columns
  design_trt[["Treatment"]] <- levels(treatments_sortedfac[["Treatment"]])[design_trt[["Treatment"]]]
  for (i in seq_along(bc_columns)) {
    design_trt[[bc_columns[i]]] <- levels(treatments_sortedfac[[stringr::str_c(bc_columns[i], "_bc")]])[design_trt[[bc_columns[i]]]]
  }

  dplyr::select(design_trt, -all_of(ends_with("_bc"))) |>
    dplyr::select(-any_of(c("Position")))
}

# Form cages from animal list with treatment groups assigned

Invivo_assignCages <- function(design_trt,
                               cagegroup_vars,
                               unique_vars = c(),
                               balance_cage_vars = c(),
                               n_min = 2, n_max = 5, n_ideal = 2, prefer_big_groups = TRUE, strict = TRUE,
                               maxiter = 5e3,
                               quiet_process = FALSE, quiet_optimize = TRUE) {

  # Set up batch container
  # We just need the subgrouping functionality of shuffle_grouped_data(), so the obligatory assignment variable
  # is just a dummy with one factor level. (Initially thought that Treatment would be assigned at this step, too.)

  if (!quiet_process) message("Setting up batch container.")

  bc_cage <- BatchContainer$new(
    dimensions = c("Dummy" = 1, ID = nrow(design_trt))
  )
  bc_cage <- assign_in_order(bc_cage, design_trt)

  shuffle_proposal <- shuffle_grouped_data(bc_cage,
    allocate_var = "Dummy",
    keep_together_vars = cagegroup_vars,
    keep_separate_vars = unique_vars,
    subgroup_var_name = "Cage",
    report_grouping_as_attribute = TRUE,
    n_min = n_min, n_max = n_max, n_ideal = n_ideal,
    prefer_big_groups = prefer_big_groups, strict = strict
  )

  if (!quiet_process) {
    message(
      "\nExpecting ", dplyr::n_distinct(shuffle_proposal()$samples_attr$Cage),
      " cages to be created and ", sum(table(shuffle_proposal()$samples_attr$Cage) == 1),
      " single-housed animals.\n"
    )
  }

  if (!quiet_process) message("Constructing scoring functions:")

  # Set up required scoring functions; most relevant ones come first to allow later application of relevance ranking
  scoring_functions <- list()
  bc_data <- bc_cage$get_samples()

  trt_vars_cat <- intersect(balance_cage_vars, colnames(bc_data)) |>
    intersect(DoE_multipleLevelCols(bc_data)) |>
    setdiff(DoE_numericalCols(bc_data))
  if (length(trt_vars_cat) > 0) {
    if (!quiet_process) message("     ... OSAT for categorical variables balanced across cages (", stringr::str_c(trt_vars_cat, collapse = ", "), ")")
    scoring_functions <- c(scoring_functions, balanced_categorical = osat_score_generator(batch_vars = "Cage", feature_vars = trt_vars_cat))
  }

  trt_vars_num <- intersect(balance_cage_vars, colnames(bc_data)) |>
    intersect(DoE_multipleLevelCols(bc_data)) |>
    intersect(DoE_numericalCols(bc_data))
  if (length(trt_vars_num) > 0) {
    if (!quiet_process) message("     ... ANOVA -logP for numerical variables balanced across cages (", stringr::str_c(trt_vars_num, collapse = ", "), ")")
    sf <- list()
    for (i in seq_along(trt_vars_num)) {
      sf <- c(sf, anova_logP(batch_var = "Cage", feature_var = trt_vars_num[i]))
    }
    names(sf) <- stringr::str_c("balanced_", trt_vars_num)
    scoring_functions <- c(scoring_functions, sf)
  }

  if (length(scoring_functions) == 0) {
    if (!quiet_process) message("     ... just a dummy score as there are no user provided balancing variables")
    scoring_functions <- osat_score_generator(batch_vars = "Dummy", feature_vars = c("Treatment"))
  }

  bc_cage <- optimize_design(
    bc_cage,
    scoring = scoring_functions,
    shuffle_proposal_func = shuffle_proposal,
    acceptance_func = accept_leftmost_improvement,
    max_iter = maxiter,
    sample_attributes_fixed = F,
    quiet = quiet_optimize
  )

  bc_cage$get_samples() |>
    dplyr::select(-any_of(c("Dummy", "ID", "alloc_var_level", "group", "subgroup")))
}

# Arrange cages in a rack

Invivo_arrangeCages <- function(design_cage,
                                distribute_cagerack_vars = "Treatment",
                                rack_size_x = 4,
                                rack_size_y = 4,
                                n_shuffle = c(rep(5, 100), rep(3, 400), rep(2, 500), rep(1, 4000)),
                                quiet_process = FALSE, quiet_optimize = TRUE) {

  # How many racks do we need?
  nr_cages <- dplyr::n_distinct(design_cage[["Cage"]])
  nr_racks <- ceiling(nr_cages / (rack_size_x * rack_size_y))

  # Identify a minimal n x m rectangle to use from every rack
  approx_cages_per_rack <- ceiling(nr_cages / nr_racks)
  n_cage_x <- rack_size_x
  n_cage_y <- rack_size_y
  may_shrink <- TRUE
  while (may_shrink) {
    may_shrink <- FALSE
    if ((n_cage_y - 1) * n_cage_x >= approx_cages_per_rack) {
      n_cage_y <- n_cage_y - 1
      may_shrink <- TRUE
    }
    if ((n_cage_x - 1) * n_cage_y >= approx_cages_per_rack) {
      n_cage_x <- n_cage_x - 1
      may_shrink <- TRUE
    }
  }

  if (!quiet_process) message("Needing ", nr_racks, " rack", ifelse(nr_racks == 1, "", "s"), " with a grid of ", n_cage_x, " x ", n_cage_y, " cages.")

  empty <- nr_racks * (n_cage_x * n_cage_y) - nr_cages
  if (!quiet_process && empty > 0) message("There will be ", empty, " empty position", ifelse(empty == 1, "", "s"), " overall.")

  if (!quiet_process) message("Setting up batch container.")

  design_rack <- dplyr::select(design_cage, all_of(c("Cage", distribute_cagerack_vars))) |> dplyr::distinct()
  assertthat::assert_that(!any(duplicated(design_rack[["Cage"]])), msg = "'Distribute rack' variables are not homogeneous within cages")

  bc_rack <- BatchContainer$new(
    dimensions = c(Rack = nr_racks, CageRow = n_cage_x, CageCol = n_cage_y)
  )

  bc_rack <- assign_random(bc_rack, design_rack)

  # Firstly, distribute variables across racks if necessary
  if (nr_racks > 1) {
    if (!quiet_process) {
      message(
        " \nDistributing target variables evenly across racks (OSAT score for ",
        stringr::str_c(distribute_cagerack_vars, collapse = ", "), ")"
      )
    }

    scoring_functions <- list(across_rack = osat_score_generator(batch_vars = c("Rack"), feature_vars = distribute_cagerack_vars))

    bc_rack <- optimize_design(
      bc_rack,
      scoring = scoring_functions,
      quiet = quiet_optimize,
      min_score = 0, max_iter = 1e3,
      n_shuffle = 2,
      acceptance_func = mk_simanneal_acceptance_func(mk_simanneal_temp_func(T0 = 10000, alpha = 0.5))
    )

    if (!quiet_process) message("   ... final score: ", bc_rack$score(scoring_f))
  }

  if (!quiet_process) {
    message(
      "\nDistributing target variables (", stringr::str_c(distribute_cagerack_vars, collapse = ", "),
      ") within rack", ifelse(nr_racks == 1, "", "s")
    )
  }

  scoring_functions <- list()
  for (tv in distribute_cagerack_vars) {
    scoring_functions <- c(scoring_functions, mk_plate_scoring_functions(bc_rack,
      plate = "Rack", row = "CageRow", column = "CageCol",
      group = tv, penalize_lines = "none"
    ))
  }
  names(scoring_functions) <- stringr::str_c(names(scoring_functions), rep(distribute_cagerack_vars, each = nr_racks), sep = "_")

  for (i in 1:nr_racks) {
    if (!quiet_process) message("   ... Rack ", i)

    bc_rack <- optimize_design(
      bc_rack,
      scoring = scoring_functions,
      shuffle_proposal_func = mk_subgroup_shuffling_function(
        subgroup_vars = "Rack",
        restrain_on_subgroup_levels = c(i),
        n_swaps = n_shuffle
      ),
      # aggregate_scores_func = sum_scores,  # L2s aggregation doesn't make sense for autoscaled scores!
      # acceptance_func = mk_simanneal_acceptance_func(mk_simanneal_temp_func(T0 = 10000, alpha = 0.5)),
      autoscale_scores = TRUE,
      acceptance_func = ~ accept_leftmost_improvement(..., tolerance = 0.1),
      quiet = quiet_optimize
    )
  }

  if (!quiet_process) message("   ... final scores: ", paste(names(bc_rack$score(scoring_functions)), round(bc_rack$score(scoring_functions), 2), sep = ": ", collapse = ", "))

  # Translate Rack numbers to some text output and assign CageNr
  design_rack <- bc_rack$get_samples() |>
    dplyr::filter(!is.na(Cage)) |> # strip empty locations in rack
    dplyr::arrange(Rack, CageRow, CageCol) |>
    dplyr::mutate(
      CageNr = 1:nr_cages,
      Rack = stringr::str_c("Rack ", Rack)
    )

  # Join animal list with cage information and add Group variable needed for the ScoreSheet generator
  dplyr::inner_join(design_cage, design_rack, by = intersect(colnames(design_cage), colnames(design_rack)))
}

## ----plotting_functions, echo=FALSE-------------------------------------------
Invivo_plotByRack <- function(design, colorBy = "Treatment", showAnimals = FALSE, animalLabel = "AnimalID", showLegend = FALSE) {
  ps <- list()

  paste_animals <- function(animalids, earmarks) {
    if (any(earmarks != "")) {
      return(stringr::str_c(animalids, " (", earmarks, ")", collapse = "\n"))
    }
    stringr::str_c(animalids, collapse = "\n")
  }

  if (!"Earmark" %in% colnames(design)) {
    design$Earmark <- "-"
  }

  if (!"Rack" %in% colnames(design)) {
    design$Rack <- "Rack 1"
  }

  for (rack in unique(design$Rack)) {
    design_f <- dplyr::filter(design, Rack == rack)
    design_u <- dplyr::select(design_f, all_of(c("CageNr", "CageRow", "CageCol", colorBy))) |> dplyr::distinct()

    p <- ggplot2::ggplot(design_f, aes(x = NA, y = NA)) +
      facet_grid(CageRow ~ CageCol) +
      geom_raster(data = design_u, aes(fill = !!as.symbol(colorBy)), alpha = 0.9) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )

    if (showAnimals) {
      cage_list <- dplyr::group_by(design_f, across(all_of(c("CageNr", "CageRow", "CageCol")))) |>
        dplyr::summarize(Animals = paste_animals(!!as.symbol(animalLabel), Earmark), .groups = "drop")
      p <- p + geom_text(data = cage_list, aes(label = Animals), size = 3) +
        theme(legend.position = "none") +
        labs(title = rack, subtitle = paste("Animals labeled by", animalLabel, "(Earmark)"))
    } else {
      p <- p + geom_text(data = design_u, aes(label = CageNr), size = 3) +
        labs(title = rack, fill = "", subtitle = paste("Cage by", colorBy))
      if (showLegend) {
        p <- p + theme(legend.position = "bottom")
      } else {
        p <- p + theme(legend.position = "none")
      }
    }

    ps[[rack]] <- p
  }

  ps
}

## -----------------------------------------------------------------------------
data("invivo_study_samples")
data("invivo_study_treatments")

## ----echo=TRUE----------------------------------------------------------------
str(invivo_study_samples)

invivo_study_samples |>
  dplyr::count(Strain, Sex, BirthDate) |>
  gt::gt()

## ----echo=TRUE----------------------------------------------------------------
invivo_study_samples |>
  dplyr::count(Strain, Litter, BirthDate) |>
  gt::gt()

## ----echo=TRUE----------------------------------------------------------------
str(invivo_study_treatments)

invivo_study_treatments |>
  dplyr::count(Treatment, Strain, Sex) |>
  gt::gt()

## ----echo=TRUE----------------------------------------------------------------
invivo_study_samples <- dplyr::mutate(invivo_study_samples,
  AgeGroup = as.integer(factor(BirthDate, exclude = NULL)),
  Litter_combine_females = ifelse(Sex == "F", "female_all", Litter)
)

invivo_study_samples |>
  dplyr::count(Strain, Litter_combine_females, BirthDate, AgeGroup) |>
  gt::gt()

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  InVivo_assignTreatments <- function(animal_list, treatments,
#                                      balance_treatment_vars = c(),
#                                      form_cages_by = c(),
#                                      n_shuffle = c(rep(5, 100), rep(3, 200), rep(2, 500), rep(1, 20000)),
#                                      quiet_process = FALSE, quiet_optimize = TRUE) {
#    (...)
#  }

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  Invivo_assignCages <- function(design_trt,
#                                 cagegroup_vars,
#                                 unique_vars = c(),
#                                 balance_cage_vars = c(),
#                                 n_min = 2, n_max = 5, n_ideal = 2, prefer_big_groups = TRUE, strict = TRUE,
#                                 maxiter = 5e3,
#                                 quiet_process = FALSE, quiet_optimize = TRUE) {
#    (...)
#  }

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  Invivo_arrangeCages <- function(design_cage,
#                                  distribute_cagerack_vars = "Treatment",
#                                  rack_size_x = 4,
#                                  rack_size_y = 4,
#                                  n_shuffle = c(rep(5, 100), rep(3, 400), rep(2, 500), rep(1, 4000)),
#                                  quiet_process = FALSE, quiet_optimize = TRUE) {
#    (...)
#  }

## ----include=FALSE------------------------------------------------------------
# stop if called from precompiling the child
if (exists(".precompile")) knitr::knit_exit()

