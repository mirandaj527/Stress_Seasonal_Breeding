#!/usr/bin/env Rscript

# ==============================================================================
# PURPOSE
# ==============================================================================
#
# This script reads one or more dynamic-programming strategy files and
# summarises the optimal hormone decisions across:
#
#   t = time since the most recent predator attack
#   d = damage state
#   s = breeding phase
#
# The script produces four linked, detailed tables:
#
#   1. strategy_summary_overall.csv
#   2. strategy_summary_by_phase.csv
#   3. strategy_hormone_frequencies.csv
#   4. strategy_switching_damage_by_time.csv
#
# These tables deliberately contain overlapping and exploratory measures.
# After examining a wider range of model outputs, some columns may be removed
# or moved into condensed or supplementary tables.
#
# IMPORTANT TERMINOLOGY
#
# h0:
#   The built-in baseline hormone value supplied to the model.
#
# h_max:
#   The highest hormone value permitted by the model.
#   This is calculated as maxH - 1.
#
# predator_presence_risk:
#   The probability that the predator is present, calculated from lambdaA and
#   lambdaL. This is different from pAtt, which is the probability that a
#   present predator attacks.


# ==============================================================================
# BASIC SCRIPT SET-UP
# ==============================================================================

# Remove existing objects so they cannot interfere with the script.

rm(list = ls())


# Load the required packages.

suppressPackageStartupMessages(
  library("ggplot2")
)

suppressPackageStartupMessages(
  library("patchwork")
)

suppressPackageStartupMessages(
  library("readr")
)

suppressPackageStartupMessages(
  library("filesstrings")
)

suppressPackageStartupMessages(
  library("tidyverse")
)

suppressPackageStartupMessages(
  library("viridis")
)

suppressPackageStartupMessages(
  library("here")
)

suppressPackageStartupMessages(
  library("doParallel")
)


# Number of cores available for parallel processing.
#
# From Bram: please adjust the number of cores to what you have on your
# computer. If you run this in RStudio, you could go quite high, assuming that
# it will not take mega long. Try it out!
#
# NOTE:
# The current version uses lapply(), so it processes files sequentially.
# The registered cores are retained in case the file-processing step is changed
# back to foreach() and %dopar% later.

registerDoParallel(
  cores = 5
)


# ==============================================================================
# SCRIPT INPUTS AND DIRECTORIES
# ==============================================================================

# Folder containing this script and stress_file_functions.r.

code_dir <- path.expand(
  "~/Documents/StressModel/code/Stress_Seasonal_Breeding"
)

# The script requires three command-line arguments:
#
# Rscript summarise_strategies.R DATA_DIR SWEEP_NAME VARIED_PARAMETERS
#
# VARIED_PARAMETERS can contain one name, such as "mu", or a comma-separated
# list, such as "pAtt,alpha".

script_arguments <- commandArgs(trailingOnly = TRUE)

if (length(script_arguments) == 0) {
  stop(
    paste0(
      "No arguments specified. Please run the script using:\n\n",
      "Rscript summarise_strategies.R ",
      "\"DATA_DIR\" \"SWEEP_NAME\" \"VARIED_PARAMETERS\"\n\n",
      "Example for one varied parameter:\n",
      "Rscript summarise_strategies.R ",
      "\"~/Documents/StressModel/data/mu_sweep_01\" ",
      "\"mu_sweep_01\" \"mu\"\n\n",
      "For multiple varied parameters, use a comma-separated list, ",
      "for example: \"mu,pAtt\"."
    )
  )
}

if (length(script_arguments) != 3) {
  stop(
    paste0(
      "Incorrect number of arguments supplied. The script requires exactly ",
      "three arguments:\n\n",
      "Rscript summarise_strategies.R ",
      "\"DATA_DIR\" \"SWEEP_NAME\" \"VARIED_PARAMETERS\""
    )
  )
}

data_dir <- path.expand(script_arguments[1])
sweep_name <- script_arguments[2]
varied_parameters_text <- script_arguments[3]

varied_parameters <- trimws(
  unlist(
    strsplit(
      varied_parameters_text,
      split = ",",
      fixed = TRUE
    )
  )
)

varied_parameters <- varied_parameters[varied_parameters != ""]

if (!dir.exists(data_dir)) {
  stop("The data directory does not exist: ", data_dir)
}

if (sweep_name == "") {
  stop("sweep_name cannot be empty.")
}

if (length(varied_parameters) == 0) {
  stop("At least one varied parameter must be supplied.")
}

message("Data directory: ", data_dir)
message("Sweep name: ", sweep_name)
message("Varied parameters: ", paste(varied_parameters, collapse = ", "))


# Load the existing functions used to read the DP output files.

source(
  file.path(
    code_dir,
    "stress_file_functions.r"
  )
)


# ==============================================================================
# FIND THE STRATEGY FILES
# ==============================================================================

strategy_files <- data.frame(
  file = list.files(
    path = data_dir,
    pattern = "^sim_stress_life_.*\\.csv$",
    full.names = TRUE
  )
)


# Stop with a clear message if no strategy files were found.

if (nrow(strategy_files) == 0) {
  stop(
    "No strategy files were found in: ",
    data_dir
  )
}


# Show which files will be processed.

print(strategy_files)


# ==============================================================================
# SMALL HELPER FUNCTIONS
# ==============================================================================

# SMALL HELPER FUNCTION 1: FIND THE LOWEST DAMAGE USING A HORMONE

# This function finds the lowest damage state at which a particular hormone
# value is used.
#
# It returns NA when the requested hormone does not occur. This prevents
# min() from returning an invalid result when, for example, h_max is absent.

find_lowest_damage_using_hormone <- function(
    data_strategy_subset,
    target_hormone
) {
  
  matching_states <-
    data_strategy_subset %>%
    filter(
      hormone == target_hormone
    )
  
  if (nrow(matching_states) == 0) {
    
    lowest_damage <- NA_real_
    
  } else {
    
    lowest_damage <- min(
      matching_states$d
    )
  }
  
  return(lowest_damage)
}

# SMALL HELPER FUNCTION 2: FIND THE HIGHEST HORMONE BELOW H_MAX

# This finds the highest observed hormone below h_max.
#
# If every state uses h_max, there is no observed hormone below h_max, so the
# function returns NA.

find_maximum_hormone_below_h_max <- function(
    data_strategy_subset,
    h_max
) {
  
  states_below_h_max <-
    data_strategy_subset %>%
    filter(
      hormone < h_max
    )
  
  if (nrow(states_below_h_max) == 0) {
    
    maximum_below_h_max <- NA_real_
    
  } else {
    
    maximum_below_h_max <- max(
      states_below_h_max$hormone
    )
  }
  
  return(maximum_below_h_max)
}


# ==============================================================================
# EXTRACT THE PARAMETERS AND STRATEGY FROM ONE FILE
# ==============================================================================

extract_strategy_and_params_single_file <- function(filename) {
  
  # Get the parameters for this model run.
  
  params <- read.parameters(
    filename
  )
  
  
  # Autocorrelation in predator presence.
  
  params$autocorrelation <-
    1 -
    params$lambdaA -
    params$lambdaL
  
  
  # Probability that the predator is present.
  #
  # This used to be called "risk", but predator_presence_risk is more precise.
  # It should not be confused with pAtt, which describes the probability that
  # a present predator attacks.
  
  params$predator_presence_risk <-
    params$lambdaA /
    (
      params$lambdaL +
      params$lambdaA
    )

  # Check that every named varied parameter is present in the model output.

  missing_varied_parameters <- setdiff(
    varied_parameters,
    names(params)
  )

  if (length(missing_varied_parameters) > 0) {
    stop(
      "These varied parameters are missing from ",
      basename(filename),
      ": ",
      paste(missing_varied_parameters, collapse = ", ")
    )
  }

  # Record the varied value or combination of values for this model run.

  varied_parameter_values <- paste(
    paste0(
      varied_parameters,
      "=",
      unlist(
        params[1, varied_parameters, drop = FALSE],
        use.names = FALSE
      )
    ),
    collapse = ", "
  )

  params$sweep_name <- sweep_name
  params$varied_parameters <- paste(varied_parameters, collapse = ",")
  params$varied_parameter_values <- varied_parameter_values
  
  
  # Get the raw strategy data: the optimal hormone for every s, d, and t.
  
  data_strategy <- as_tibble(
    read.stress.file(
      filename
    )
  )
  
  
  # Put the filename, parameters, and strategy together in a list so that the
  # table-making functions can all use the same extracted information.
  
  result <- list(
    filename = filename,
    params = params,
    data_strategy = data_strategy
  )
  
  return(result)
}


# ==============================================================================
# TABLE 1 FUNCTION: summarise_overall_strategy()
# OUTPUT: strategy_summary_overall.csv
# ==============================================================================

# One row per complete model run.
#
# This table provides a whole-strategy summary across all times, damage states,
# and breeding phases.

summarise_overall_strategy <- function(
    data_strategy,
    params,
    filename
) {
  
  # --------------------------------------------------------------------------
  # Important reference values
  
  t_stress <- 0
  
  t_far <-
    as.numeric(params$maxT) - 1
  
  h0 <-
    as.numeric(params$h0)
  
  h_max <-
    as.numeric(params$maxH) - 1
  
  
  # --------------------------------------------------------------------------
  # Basic strategy-size and hormone-range measures
  
  number_of_states <-
    nrow(data_strategy)
  
  
  number_of_hormone_levels_used <-
    data_strategy %>%
    distinct(
      hormone
    ) %>%
    nrow()
  
  # This is the number of different hormone levels actually used by the
  # complete strategy.
  
  
  min_hormone <-
    min(
      data_strategy$hormone
    )
  
  
  max_hormone <-
    max(
      data_strategy$hormone
    )
  
  
  max_hormone_below_h_max <-
    find_maximum_hormone_below_h_max(
      data_strategy_subset = data_strategy,
      h_max = h_max
    )
  
  # This is the highest hormone used below the maximum permitted by the model.
  
  
  # --------------------------------------------------------------------------
  # Mean post-stressor and far-from-stressor hormone
  
  post_stressor_mean_h <-
    data_strategy %>%
    filter(
      t == t_stress
    ) %>%
    summarise(
      value = mean(hormone)
    ) %>%
    pull(value)
  
  # Mean hormone across the complete strategy at t = 0.
  
  
  far_from_stressor_mean_h <-
    data_strategy %>%
    filter(
      t == t_far
    ) %>%
    summarise(
      value = mean(hormone)
    ) %>%
    pull(value)
  
  # Mean hormone across the complete strategy at t = maxT - 1.
  
  
  response_amplitude <-
    post_stressor_mean_h -
    far_from_stressor_mean_h
  
  
  # --------------------------------------------------------------------------
  # Frequency measures
  
  percentage_at_h0 <-
    (
      data_strategy %>%
        filter(
          hormone == h0
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  percentage_at_min_h <-
    (
      data_strategy %>%
        filter(
          hormone == min_hormone
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  percentage_at_max_h <-
    (
      data_strategy %>%
        filter(
          hormone == max_hormone
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  percentage_below_max_h <-
    (
      data_strategy %>%
        filter(
          hormone < max_hormone
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  percentage_at_h_max <-
    (
      data_strategy %>%
        filter(
          hormone == h_max
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  percentage_below_h_max <-
    (
      data_strategy %>%
        filter(
          hormone < h_max
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  # --------------------------------------------------------------------------
  # Undamaged far-from-stressor reference
  
  undamaged_far_from_stressor_hormone <-
    data_strategy %>%
    filter(
      t == t_far &
        d == 0
    ) %>%
    summarise(
      value = mean(hormone)
    ) %>%
    pull(value)
  
  # This is the mean undamaged far-from-stressor hormone across all phases.
  
  
  percentage_at_or_below_undamaged_far_reference <-
    (
      data_strategy %>%
        filter(
          hormone <=
            undamaged_far_from_stressor_hormone
        ) %>%
        nrow()
    ) /
    nrow(data_strategy) *
    100
  
  
  # --------------------------------------------------------------------------
  # Here's where I'd put my K means clustering & largest gap analysis
  # IF I HAD ANY.
  
  
  # --------------------------------------------------------------------------
  # Assemble the one-row result
  
  overall_summary <- tibble(
    
    number_of_strategy_states =
      number_of_states,
    
    number_of_distinct_hormone_levels =
      number_of_hormone_levels_used,
    
    overall_mean_hormone =
      mean(data_strategy$hormone),
    
    overall_post_stressor_mean_hormone =
      post_stressor_mean_h,
    
    overall_far_from_stressor_mean_hormone =
      far_from_stressor_mean_h,
    
    overall_response_amplitude =
      response_amplitude,
    
    minimum_observed_hormone =
      min_hormone,
    
    maximum_observed_hormone =
      max_hormone,
    
    h_max =
      h_max,
    
    maximum_observed_below_h_max =
      max_hormone_below_h_max,
    
    percentage_at_h0 =
      percentage_at_h0,
    
    percentage_at_observed_minimum =
      percentage_at_min_h,
    
    percentage_at_observed_maximum =
      percentage_at_max_h,
    
    percentage_below_observed_maximum =
      percentage_below_max_h,
    
    percentage_at_h_max =
      percentage_at_h_max,
    
    percentage_below_h_max =
      percentage_below_h_max,
    
    undamaged_far_from_stressor_hormone =
      undamaged_far_from_stressor_hormone,
    
    percentage_at_or_below_undamaged_far_reference =
      percentage_at_or_below_undamaged_far_reference
  )
  
  
  # Add the filename and all model parameters.
  
  overall_summary <- bind_cols(
    tibble(
      filename = basename(filename)
    ),
    params,
    overall_summary
  )
  
  
  return(overall_summary)
  
} # End of summarise_overall_strategy()


# ==============================================================================
# TABLE 2 FUNCTION: summarise_strategy_by_breeding_phase()
# OUTPUT: strategy_summary_by_phase.csv
# ==============================================================================

# One row per:
#
#   model run * breeding phase
#
# This is the main large exploratory summary table.

summarise_strategy_by_breeding_phase <- function(
    data_strategy,
    params,
    filename
) {
  
  # --------------------------------------------------------------------------
  # Important reference values
  
  t_stress <- 0
  
  t_far <-
    as.numeric(params$maxT) - 1
  
  t_middle <-
    floor(
      t_far / 2
    )
  
  h0 <-
    as.numeric(params$h0)
  
  h_max <-
    as.numeric(params$maxH) - 1
  
  
  # Find all breeding phases present in the strategy.
  
  breeding_phases <- sort(
    unique(data_strategy$s)
  )
  
  
  # Create an empty tibble.
  # Each breeding phase will produce one row.
  
  phase_summary <- tibble()
  
  
  # Loop through the breeding phases one at a time.
  
  for (focal_phase in breeding_phases) {
    
    # ------------------------------------------------------------------------
    # Select only states in the focal breeding phase
    
    data_strategy_for_phase <-
      data_strategy %>%
      filter(
        s == focal_phase
      )
    
    
    # ------------------------------------------------------------------------
    # Basic strategy-size and hormone-range measures
    
    number_of_states_in_phase <-
      nrow(data_strategy_for_phase)
    
    
    number_of_hormone_levels_used <-
      data_strategy_for_phase %>%
      distinct(
        hormone
      ) %>%
      nrow()
    
    
    overall_mean_hormone_in_phase <-
      mean(
        data_strategy_for_phase$hormone
      )
    
    
    min_hormone <-
      min(
        data_strategy_for_phase$hormone
      )
    
    
    max_hormone <-
      max(
        data_strategy_for_phase$hormone
      )
    
    
    max_hormone_below_h_max <-
      find_maximum_hormone_below_h_max(
        data_strategy_subset =
          data_strategy_for_phase,
        
        h_max =
          h_max
      )
    
    
    # ------------------------------------------------------------------------
    # Mean post-stressor and far-from-stressor hormone
    
    post_stressor_mean_h <-
      data_strategy_for_phase %>%
      filter(
        t == t_stress
      ) %>%
      summarise(
        value = mean(hormone)
      ) %>%
      pull(value)
    
    
    far_from_stressor_mean_h <-
      data_strategy_for_phase %>%
      filter(
        t == t_far
      ) %>%
      summarise(
        value = mean(hormone)
      ) %>%
      pull(value)
    
    
    response_amplitude <-
      post_stressor_mean_h -
      far_from_stressor_mean_h
    
    
    # ------------------------------------------------------------------------
    # Frequency measures
    
    percentage_at_h0 <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone == h0
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    percentage_at_min_h <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone == min_hormone
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    percentage_at_max_h <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone == max_hormone
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    percentage_below_max_h <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone < max_hormone
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    percentage_at_h_max <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone == h_max
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    percentage_below_h_max <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone < h_max
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    # ------------------------------------------------------------------------
    # Undamaged far-from-stressor reference
    
    undamaged_far_from_stressor_hormone <-
      data_strategy_for_phase %>%
      filter(
        t == t_far &
          d == 0
      ) %>%
      summarise(
        value = mean(hormone)
      ) %>%
      pull(value)
    
    
    percentage_at_or_below_undamaged_far_reference <-
      (
        data_strategy_for_phase %>%
          filter(
            hormone <=
              undamaged_far_from_stressor_hormone
          ) %>%
          nrow()
      ) /
      nrow(data_strategy_for_phase) *
      100
    
    
    # ------------------------------------------------------------------------
    # Damage switching shenanigans
    #
    # Most strategies appear to jump between h_max and lower hormone values,
    # depending on t, d, and s.
    #
    # This section examines the damage at which the strategy first switches to
    # h_max within each breeding phase.
    
    
    # ------------------------------------------------------------------------
    # H_MAX SWITCHING DAMAGE AT SELECTED TIMES
    # ------------------------------------------------------------------------
    
    # Lowest damage using h_max at t = 0.
    
    h_max_switch_damage_at_t0 <-
      data_strategy_for_phase %>%
      filter(
        t == t_stress
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = h_max
      )
    
    
    # Lowest damage using h_max at the middle time.
    
    h_max_switch_damage_at_middle_time <-
      data_strategy_for_phase %>%
      filter(
        t == t_middle
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = h_max
      )
    
    
    # Lowest damage using h_max at t = maxT - 1.
    
    h_max_switch_damage_at_far_time <-
      data_strategy_for_phase %>%
      filter(
        t == t_far
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = h_max
      )
    
    
    # ------------------------------------------------------------------------
    # H_MAX SWITCHING DAMAGE AT EVERY TIME
    # ------------------------------------------------------------------------
    
    h_max_switch_damage_by_time <-
      data_strategy_for_phase %>%
      filter(
        hormone == h_max
      ) %>%
      group_by(
        t
      ) %>%
      summarise(
        switch_damage = min(d),
        .groups = "drop"
      )
    
    
    # If h_max is never used in this phase, the table above contains no rows.
    # In that case, store NA for the summary measurements.
    
    if (nrow(h_max_switch_damage_by_time) == 0) {
      
      mean_h_max_switch_damage <-
        NA_real_
      
      median_h_max_switch_damage <-
        NA_real_
      
      minimum_h_max_switch_damage <-
        NA_real_
      
      maximum_h_max_switch_damage <-
        NA_real_
      
      h_max_switch_damage_range <-
        NA_real_
      
    } else {
      
      mean_h_max_switch_damage <-
        mean(
          h_max_switch_damage_by_time$switch_damage
        )
      
      median_h_max_switch_damage <-
        median(
          h_max_switch_damage_by_time$switch_damage
        )
      
      minimum_h_max_switch_damage <-
        min(
          h_max_switch_damage_by_time$switch_damage
        )
      
      maximum_h_max_switch_damage <-
        max(
          h_max_switch_damage_by_time$switch_damage
        )
      
      h_max_switch_damage_range <-
        maximum_h_max_switch_damage -
        minimum_h_max_switch_damage
    }
    
    
    h_max_switch_change_t0_to_far <-
      h_max_switch_damage_at_far_time -
      h_max_switch_damage_at_t0
    
    
    # ------------------------------------------------------------------------
    # OBSERVED-MAXIMUM SWITCHING DAMAGE
    # ------------------------------------------------------------------------
    
    # The observed maximum is the highest hormone actually used in this phase.
    # It may or may not equal h_max.
    
    
    observed_switch_damage_at_t0 <-
      data_strategy_for_phase %>%
      filter(
        t == t_stress
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = max_hormone
      )
    
    
    observed_switch_damage_at_middle_time <-
      data_strategy_for_phase %>%
      filter(
        t == t_middle
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = max_hormone
      )
    
    
    observed_switch_damage_at_far_time <-
      data_strategy_for_phase %>%
      filter(
        t == t_far
      ) %>%
      find_lowest_damage_using_hormone(
        target_hormone = max_hormone
      )
    
    
    # Find observed-maximum switching damage separately at every time.
    
    observed_switch_damage_by_time <-
      data_strategy_for_phase %>%
      filter(
        hormone == max_hormone
      ) %>%
      group_by(
        t
      ) %>%
      summarise(
        switch_damage = min(d),
        .groups = "drop"
      )
    
    
    # The phase maximum must occur somewhere in the phase, so this table should
    # normally contain at least one row.
    
    mean_observed_switch_damage <-
      mean(
        observed_switch_damage_by_time$switch_damage
      )
    
    
    median_observed_switch_damage <-
      median(
        observed_switch_damage_by_time$switch_damage
      )
    
    
    minimum_observed_switch_damage <-
      min(
        observed_switch_damage_by_time$switch_damage
      )
    
    
    maximum_observed_switch_damage <-
      max(
        observed_switch_damage_by_time$switch_damage
      )
    
    
    observed_switch_damage_range <-
      maximum_observed_switch_damage -
      minimum_observed_switch_damage
    
    
    observed_switch_change_t0_to_far <-
      observed_switch_damage_at_far_time -
      observed_switch_damage_at_t0
    
    
    # ------------------------------------------------------------------------
    # Assemble the one-row result for this breeding phase
    
    result_for_focal_phase <- tibble(
      
      breeding_phase =
        focal_phase,
      
      number_of_strategy_states_in_phase =
        number_of_states_in_phase,
      
      number_of_distinct_hormone_levels =
        number_of_hormone_levels_used,
      
      overall_mean_hormone_in_phase =
        overall_mean_hormone_in_phase,
      
      post_stressor_mean_hormone =
        post_stressor_mean_h,
      
      far_from_stressor_mean_hormone =
        far_from_stressor_mean_h,
      
      response_amplitude =
        response_amplitude,
      
      minimum_observed_hormone =
        min_hormone,
      
      maximum_observed_hormone =
        max_hormone,
      
      h_max =
        h_max,
      
      maximum_observed_below_h_max =
        max_hormone_below_h_max,
      
      percentage_at_h0 =
        percentage_at_h0,
      
      percentage_at_observed_minimum =
        percentage_at_min_h,
      
      percentage_at_observed_maximum =
        percentage_at_max_h,
      
      percentage_below_observed_maximum =
        percentage_below_max_h,
      
      percentage_at_h_max =
        percentage_at_h_max,
      
      percentage_below_h_max =
        percentage_below_h_max,
      
      undamaged_far_from_stressor_hormone =
        undamaged_far_from_stressor_hormone,
      
      percentage_at_or_below_undamaged_far_reference =
        percentage_at_or_below_undamaged_far_reference,
      
      h_max_switch_damage_at_t0 =
        h_max_switch_damage_at_t0,
      
      h_max_switch_damage_at_middle_time =
        h_max_switch_damage_at_middle_time,
      
      h_max_switch_damage_at_far_time =
        h_max_switch_damage_at_far_time,
      
      mean_h_max_switch_damage =
        mean_h_max_switch_damage,
      
      median_h_max_switch_damage =
        median_h_max_switch_damage,
      
      minimum_h_max_switch_damage =
        minimum_h_max_switch_damage,
      
      maximum_h_max_switch_damage =
        maximum_h_max_switch_damage,
      
      h_max_switch_damage_range =
        h_max_switch_damage_range,
      
      h_max_switch_change_t0_to_far =
        h_max_switch_change_t0_to_far,
      
      observed_switch_damage_at_t0 =
        observed_switch_damage_at_t0,
      
      observed_switch_damage_at_middle_time =
        observed_switch_damage_at_middle_time,
      
      observed_switch_damage_at_far_time =
        observed_switch_damage_at_far_time,
      
      mean_observed_switch_damage =
        mean_observed_switch_damage,
      
      median_observed_switch_damage =
        median_observed_switch_damage,
      
      minimum_observed_switch_damage =
        minimum_observed_switch_damage,
      
      maximum_observed_switch_damage =
        maximum_observed_switch_damage,
      
      observed_switch_damage_range =
        observed_switch_damage_range,
      
      observed_switch_change_t0_to_far =
        observed_switch_change_t0_to_far
    )
    
    
    # Add identifying information for the model run.
    
    result_for_focal_phase <- bind_cols(
      tibble(
        filename = basename(filename)
      ),
      params,
      result_for_focal_phase
    )
    
    
    # Add this breeding phase to the accumulating results.
    
    phase_summary <- bind_rows(
      phase_summary,
      result_for_focal_phase
    )
    
  } # End of breeding-phase loop
  
  
  # Return one row per breeding phase.
  
  return(phase_summary)
  
} # End of summarise_strategy_by_breeding_phase()


# ==============================================================================
# TABLE 3 FUNCTION: summarise_hormone_frequencies()
# OUTPUT: strategy_hormone_frequencies.csv
# ==============================================================================

# One row per:
#
#   model run * breeding-phase category * observed hormone level
#
# The breeding_phase column contains either:
#
#   - the number of a particular breeding phase; or
#   - "all" for a calculation across the complete strategy.
#
# This table preserves the complete distribution of optimal hormone decisions.

summarise_hormone_frequencies <- function(
    data_strategy,
    params,
    filename
) {
  
  # --------------------------------------------------------------------------
  # Phase-specific hormone frequencies
  
  # Group the strategy by breeding phase and hormone.
  #
  # This creates one group for every combination of:
  #
  #   breeding phase * observed hormone level
  
  phase_frequencies <-
    data_strategy %>%
    group_by(
      s,
      hormone
    ) %>%
    summarise(
      
      # Number of states in this phase using the focal hormone.
      
      number_of_states =
        n(),
      
      # ----------------------------------------------------------------------
      # Damage distribution
      
      minimum_damage_where_used =
        min(d),
      
      maximum_damage_where_used =
        max(d),
      
      damage_range_where_used =
        max(d) - min(d),
      
      number_of_distinct_damage_levels_where_used =
        n_distinct(d),
      
      mean_damage_where_used =
        mean(d),
      
      median_damage_where_used =
        median(d),
      
      # ----------------------------------------------------------------------
      # Time distribution
      
      minimum_time_where_used =
        min(t),
      
      maximum_time_where_used =
        max(t),
      
      number_of_distinct_times_where_used =
        n_distinct(t),
      
      .groups = "drop"
    )
  
  
  # --------------------------------------------------------------------------
  # Calculate percentages and ranks within each breeding phase
  
  phase_frequencies <-
    phase_frequencies %>%
    group_by(
      s
    ) %>%
    mutate(
      
      # The denominator is the number of states in the focal phase.
      
      percentage_of_states =
        100 *
        number_of_states /
        sum(number_of_states),
      
      # Rank hormones from most frequently used to least frequently used.
      # Tied hormones receive the same rank.
      
      frequency_rank =
        min_rank(
          desc(number_of_states)
        )
    ) %>%
    ungroup() %>%
    rename(
      breeding_phase = s
    ) %>%
    mutate(
      
      # Convert the phase numbers to text so they can later be combined with
      # the overall category called "all".
      
      breeding_phase =
        as.character(breeding_phase)
    )
  
  
  # --------------------------------------------------------------------------
  # Overall hormone frequencies
  
  # Repeat the same calculations across the complete strategy.
  # The data are grouped only by hormone, not by breeding phase.
  
  overall_frequencies <-
    data_strategy %>%
    group_by(
      hormone
    ) %>%
    summarise(
      
      number_of_states =
        n(),
      
      # ----------------------------------------------------------------------
      # Damage distribution
      
      minimum_damage_where_used =
        min(d),
      
      maximum_damage_where_used =
        max(d),
      
      damage_range_where_used =
        max(d) - min(d),
      
      number_of_distinct_damage_levels_where_used =
        n_distinct(d),
      
      mean_damage_where_used =
        mean(d),
      
      median_damage_where_used =
        median(d),
      
      # ----------------------------------------------------------------------
      # Time distribution
      
      minimum_time_where_used =
        min(t),
      
      maximum_time_where_used =
        max(t),
      
      number_of_distinct_times_where_used =
        n_distinct(t),
      
      .groups = "drop"
    ) %>%
    mutate(
      
      # Label these rows as whole-strategy summaries.
      
      breeding_phase =
        "all",
      
      percentage_of_states =
        100 *
        number_of_states /
        sum(number_of_states),
      
      frequency_rank =
        min_rank(
          desc(number_of_states)
        )
    )
  
  
  # --------------------------------------------------------------------------
  # Put the two tables into the same column order
  
  phase_frequencies <-
    phase_frequencies %>%
    select(
      breeding_phase,
      hormone,
      number_of_states,
      percentage_of_states,
      frequency_rank,
      
      minimum_damage_where_used,
      maximum_damage_where_used,
      damage_range_where_used,
      number_of_distinct_damage_levels_where_used,
      mean_damage_where_used,
      median_damage_where_used,
      
      minimum_time_where_used,
      maximum_time_where_used,
      number_of_distinct_times_where_used
    )
  
  
  overall_frequencies <-
    overall_frequencies %>%
    select(
      breeding_phase,
      hormone,
      number_of_states,
      percentage_of_states,
      frequency_rank,
      
      minimum_damage_where_used,
      maximum_damage_where_used,
      damage_range_where_used,
      number_of_distinct_damage_levels_where_used,
      mean_damage_where_used,
      median_damage_where_used,
      
      minimum_time_where_used,
      maximum_time_where_used,
      number_of_distinct_times_where_used
    )
  
  
  # --------------------------------------------------------------------------
  # Combine the phase-specific and overall results
  
  frequency_summary <-
    bind_rows(
      phase_frequencies,
      overall_frequencies
    )
  
  
  # Within each phase category, put the most common hormones first.
  # Hormone is used as the second ordering variable when frequencies are tied.
  
  frequency_summary <-
    frequency_summary %>%
    arrange(
      breeding_phase,
      frequency_rank,
      hormone
    )
  
  
  # --------------------------------------------------------------------------
  # Add the filename and model parameters
  
  # frequency_summary contains multiple rows, but params contains one row.
  # Repeat the parameter row once for every frequency row.
  
  number_of_frequency_rows <-
    nrow(frequency_summary)
  
  
  repeated_params <-
    params[
      rep(
        1,
        number_of_frequency_rows
      ),
      ,
      drop = FALSE
    ]
  
  
  repeated_filename <- tibble(
    filename = rep(
      basename(filename),
      number_of_frequency_rows
    )
  )
  
  
  frequency_summary <- bind_cols(
    repeated_filename,
    repeated_params,
    frequency_summary
  )
  
  
  return(frequency_summary)
  
} # End of summarise_hormone_frequencies()


# ==============================================================================
# TABLE 4 FUNCTION: summarise_switching_damage_by_time()
# OUTPUT: strategy_switching_damage_by_time.csv
# ==============================================================================

# One row per:
#
#   model run * breeding phase * time since stressor
#
# This table retains the complete time series of h_max switching damage instead
# of immediately reducing it to a phase mean, median, minimum, or maximum.
#
# The condensed-summary script will later use this table to calculate:
#
#   - switching damage at t = 0;
#   - switching damage at t = 1;
#   - switching damage at the far time;
#   - minimum and maximum switching damage;
#   - switching-damage range;
#   - number of distinct switching-damage values;
#   - switching-damage stabilisation time.

summarise_switching_damage_by_time <- function(
    data_strategy,
    params,
    filename
) {
  
  # --------------------------------------------------------------------------
  # Important reference value
  
  h_max <-
    as.numeric(params$maxH) - 1
  
  
  # Find all breeding phases and time steps represented in the strategy.
  
  breeding_phases <- sort(
    unique(data_strategy$s)
  )
  
  time_steps <- sort(
    unique(data_strategy$t)
  )
  
  
  # Create an empty tibble.
  # Each phase and time combination will produce one row.
  
  switching_summary <- tibble()
  
  
  # --------------------------------------------------------------------------
  # Loop through breeding phases
  
  for (focal_phase in breeding_phases) {
    
    # Select only the focal breeding phase.
    
    data_strategy_for_phase <-
      data_strategy %>%
      filter(
        s == focal_phase
      )
    
    
    # ------------------------------------------------------------------------
    # Loop through time steps within the focal phase
    
    for (focal_timestep in time_steps) {
      
      # Select one complete damage profile:
      # one phase at one time since the stressor.
      
      data_strategy_for_phase_and_time <-
        data_strategy_for_phase %>%
        filter(
          t == focal_timestep
        )
      
      
      # Select the damage states using h_max.
      
      damage_states_at_h_max <-
        data_strategy_for_phase_and_time %>%
        filter(
          hormone == h_max
        )
      
      
      # Find the lowest damage at which h_max is used.
      #
      # If h_max is not used, switching_damage will be NA.
      
      switching_damage_for_phase_and_time <-
        find_lowest_damage_using_hormone(
          data_strategy_subset =
            data_strategy_for_phase_and_time,
          
          target_hormone =
            h_max
        )
      
      
      # ----------------------------------------------------------------------
      # Assemble the result for this phase and time
      
      result_for_focal_phase_and_time <- tibble(
        
        breeding_phase =
          focal_phase,
        
        time_since_stressor =
          focal_timestep,
        
        h_max =
          h_max,
        
        switching_damage =
          switching_damage_for_phase_and_time,
        
        # TRUE when h_max does not occur anywhere in this damage profile.
        
        no_h_max_state =
          nrow(damage_states_at_h_max) == 0,
        
        # TRUE when every damage state in the profile uses h_max.
        
        all_states_at_h_max =
          all(
            data_strategy_for_phase_and_time$hormone ==
              h_max
          ),
        
        # Number of damage states represented in this profile.
        
        number_of_damage_states =
          nrow(data_strategy_for_phase_and_time),
        
        # Percentage of damage states using h_max.
        
        percentage_at_h_max =
          100 *
          nrow(damage_states_at_h_max) /
          nrow(data_strategy_for_phase_and_time)
      )
      
      
      # Add the filename and model parameters.
      
      result_for_focal_phase_and_time <- bind_cols(
        tibble(
          filename = basename(filename)
        ),
        params,
        result_for_focal_phase_and_time
      )
      
      
      # Add this result to the accumulating table.
      
      switching_summary <- bind_rows(
        switching_summary,
        result_for_focal_phase_and_time
      )
      
    } # End of time-step loop
    
  } # End of breeding-phase loop
  
  
  return(switching_summary)
  
} # End of summarise_switching_damage_by_time()


# ==============================================================================
# COORDINATING FUNCTION FOR ONE STRATEGY FILE
# ==============================================================================

# This function handles one complete model-output file.
#
# It:
#
#   1. extracts the parameters and strategy data;
#   2. passes those objects to each of the four table functions;
#   3. returns the four resulting tables together in a named list.
#
# A list is useful because the four tables have different numbers of rows and
# therefore cannot be joined directly into one data frame.

summarise_single_file <- function(filename) {
  
  message(
    "Reading: ",
    basename(filename)
  )
  
  
  # --------------------------------------------------------------------------
  # Extract the parameters and strategy data
  
  extracted_file <-
    extract_strategy_and_params_single_file(
      filename = filename
    )
  
  
  # --------------------------------------------------------------------------
  # Create Table 1
  
  # Table 1 returns one row for this model run.
  
  overall_result <-
    summarise_overall_strategy(
      data_strategy =
        extracted_file$data_strategy,
      
      params =
        extracted_file$params,
      
      filename =
        extracted_file$filename
    )
  
  
  # --------------------------------------------------------------------------
  # Create Table 2
  
  # Table 2 returns one row for every breeding phase.
  
  phase_result <-
    summarise_strategy_by_breeding_phase(
      data_strategy =
        extracted_file$data_strategy,
      
      params =
        extracted_file$params,
      
      filename =
        extracted_file$filename
    )
  
  
  # --------------------------------------------------------------------------
  # Create Table 3
  
  # Table 3 returns one row for every observed hormone level within:
  #
  #   - each individual breeding phase; and
  #   - the complete strategy.
  
  frequency_result <-
    summarise_hormone_frequencies(
      data_strategy =
        extracted_file$data_strategy,
      
      params =
        extracted_file$params,
      
      filename =
        extracted_file$filename
    )
  
  
  # --------------------------------------------------------------------------
  # Create Table 4
  
  # Table 4 returns one row for every breeding phase and time combination.
  
  switching_by_time_result <-
    summarise_switching_damage_by_time(
      data_strategy =
        extracted_file$data_strategy,
      
      params =
        extracted_file$params,
      
      filename =
        extracted_file$filename
    )
  
  
  # --------------------------------------------------------------------------
  # Store the four results together
  
  result <- list(
    
    overall =
      overall_result,
    
    by_phase =
      phase_result,
    
    frequencies =
      frequency_result,
    
    switching_by_time =
      switching_by_time_result
  )
  
  
  return(result)
  
} # End of summarise_single_file()


# ==============================================================================
# APPLY THE COORDINATING FUNCTION TO EVERY STRATEGY FILE
# ==============================================================================

# strategy_files is a data frame containing a column called file.
#
# strategy_files$file extracts the vector of filenames.
#
# lapply() passes those filenames through summarise_single_file() one at a time.
#
# The result is a list containing one set of four tables for every file.

results_by_file <-
  lapply(
    strategy_files$file,
    summarise_single_file
  )


# ==============================================================================
# COMBINE TABLE 1 RESULTS FROM ALL MODEL RUNS
# ==============================================================================

# Each element of results_by_file contains an object called $overall.
#
# bind_rows() places those one-row tables underneath one another.

overall_results <-
  bind_rows(
    lapply(
      results_by_file,
      function(file_result) {
        file_result$overall
      }
    )
  )


# ==============================================================================
# COMBINE TABLE 2 RESULTS FROM ALL MODEL RUNS
# ==============================================================================

# Each by-phase result contains one row per breeding phase.
#
# bind_rows() combines those phase rows across all model runs.

phase_results <-
  bind_rows(
    lapply(
      results_by_file,
      function(file_result) {
        file_result$by_phase
      }
    )
  )


# ==============================================================================
# COMBINE TABLE 3 RESULTS FROM ALL MODEL RUNS
# ==============================================================================

# Each frequency result contains phase-specific and overall hormone-frequency
# rows for one model run.

frequency_results <-
  bind_rows(
    lapply(
      results_by_file,
      function(file_result) {
        file_result$frequencies
      }
    )
  )


# ==============================================================================
# COMBINE TABLE 4 RESULTS FROM ALL MODEL RUNS
# ==============================================================================

# Each switching-by-time result contains one row for every:
#
#   breeding phase * time since stressor

switching_by_time_results <-
  bind_rows(
    lapply(
      results_by_file,
      function(file_result) {
        file_result$switching_by_time
      }
    )
  )


# ==============================================================================
# WRITE THE FOUR COMPLETED TABLES
# ==============================================================================

# write_delim() saves the tables as semicolon-separated CSV files.
#
# All four tables are written into data_dir.

write_delim(
  x = overall_results,
  file = file.path(
    data_dir,
    "strategy_summary_overall.csv"
  ),
  delim = ";"
)


write_delim(
  x = phase_results,
  file = file.path(
    data_dir,
    "strategy_summary_by_phase.csv"
  ),
  delim = ";"
)


write_delim(
  x = frequency_results,
  file = file.path(
    data_dir,
    "strategy_hormone_frequencies.csv"
  ),
  delim = ";"
)


write_delim(
  x = switching_by_time_results,
  file = file.path(
    data_dir,
    "strategy_switching_damage_by_time.csv"
  ),
  delim = ";"
)


# Print a message after all four files have been written.

message(
  "Finished writing all four detailed strategy-summary tables to: ",
  data_dir
)
