#!/usr/bin/env Rscript

# ==============================================================================
# PURPOSE AND TERMS
# ==============================================================================
# This script reads the four detailed strategy-summary tables (produced by 
# summarise_strategies.R) and produces four smaller tables for quicker 
# comparisons among model runs and breeding phases.
#
# Naming abbreviations used in the compact tables:
#   min = minimum
#   max = maximum
#   obs = observed
#   h_max = highest hormone permitted by the model (maxH - 1)

rm(list = ls())
suppressPackageStartupMessages(library("tidyverse"))

# ==============================================================================
# SCRIPT INPUTS
# ==============================================================================
# The script only requires the directory containing the detailed tables:
#
# Rscript condensed_strategy_summaries.R DATA_DIR
#
# The varied parameters are read from the detailed tables.

script_arguments <- commandArgs(trailingOnly = TRUE)

if (length(script_arguments) == 0) {
  stop(
    paste0(
      "No arguments specified. Please run the script using:\n\n",
      "Rscript condensed_strategy_summaries.R ",
      "\"DATA_DIR\"\n\n",
      "Example:\n",
      "Rscript condensed_strategy_summaries.R ",
      "\"~/Documents/StressModel/data/mu_sweep_01\""
    )
  )
}

if (length(script_arguments) != 1) {
  stop(
    paste0(
      "Incorrect number of arguments supplied. The script requires exactly ",
      "one argument:\n\n",
      "Rscript condensed_strategy_summaries.R ",
      "\"DATA_DIR\""
    )
  )
}

data_dir <- path.expand(script_arguments[1])

if (!dir.exists(data_dir)) {
  stop("The data directory does not exist: ", data_dir)
}

message("Data directory: ", data_dir)

# ==============================================================================
# READ THE DETAILED TABLES
# ==============================================================================
required_files <- c(
  "strategy_summary_overall.csv",
  "strategy_summary_by_phase.csv",
  "strategy_hormone_frequencies.csv",
  "strategy_switching_damage_by_time.csv"
)

missing_files <- required_files[
  !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Run summarise_strategies.R first. Missing files: ",
    paste(missing_files, collapse = ", ")
  )
}

overall_detailed <- read_delim(
  file.path(data_dir, "strategy_summary_overall.csv"),
  delim = ";",
  show_col_types = FALSE
)

phase_detailed <- read_delim(
  file.path(data_dir, "strategy_summary_by_phase.csv"),
  delim = ";",
  show_col_types = FALSE
)

hormone_frequencies <- read_delim(
  file.path(data_dir, "strategy_hormone_frequencies.csv"),
  delim = ";",
  show_col_types = FALSE
)

switching_by_time <- read_delim(
  file.path(data_dir, "strategy_switching_damage_by_time.csv"),
  delim = ";",
  show_col_types = FALSE
)

# Read the varied-parameter definition from the detailed tables.
detailed_tables <- list(
  overall = overall_detailed,
  by_phase = phase_detailed,
  hormone_frequencies = hormone_frequencies,
  switching_by_time = switching_by_time
)

tables_missing_varied_parameters <- names(
  detailed_tables
)[
  !vapply(
    detailed_tables,
    function(detailed_table) {
      "varied_parameters" %in% names(detailed_table)
    },
    logical(1)
  )
]

if (length(tables_missing_varied_parameters) > 0) {
  stop(
    "These detailed tables do not contain varied_parameters: ",
    paste(tables_missing_varied_parameters, collapse = ", "),
    ". Run summarise_strategies.R again first."
  )
}

varied_parameter_definitions <- unique(
  c(
    overall_detailed$varied_parameters,
    phase_detailed$varied_parameters,
    hormone_frequencies$varied_parameters,
    switching_by_time$varied_parameters
  )
)

varied_parameter_definitions <- varied_parameter_definitions[
  !is.na(varied_parameter_definitions) &
    varied_parameter_definitions != ""
]

if (length(varied_parameter_definitions) != 1) {
  stop(
    "The detailed tables do not contain one consistent varied-parameter ",
    "definition. Found: ",
    paste(varied_parameter_definitions, collapse = " | ")
  )
}

varied_parameters <- trimws(
  unlist(
    strsplit(
      varied_parameter_definitions[1],
      split = ",",
      fixed = TRUE
    )
  )
)

varied_parameters <- varied_parameters[varied_parameters != ""]

missing_parameter_columns <- unique(
  unlist(
    lapply(
      detailed_tables,
      function(detailed_table) {
        setdiff(varied_parameters, names(detailed_table))
      }
    )
  )
)

if (length(missing_parameter_columns) > 0) {
  stop(
    "These varied-parameter value columns are missing from the detailed ",
    "tables: ",
    paste(missing_parameter_columns, collapse = ", ")
  )
}

message(
  "Varied parameters read from detailed tables: ",
  paste(varied_parameters, collapse = ", ")
)

# Create one name column and one value column for each varied parameter.
for (parameter_number in seq_along(varied_parameters)) {
  parameter_name_column <- paste0(
    "varied_parameter_",
    parameter_number
  )

  parameter_value_column <- paste0(
    "parameter_",
    parameter_number,
    "_value"
  )

  parameter_name <- varied_parameters[parameter_number]

  overall_detailed[[parameter_name_column]] <- parameter_name
  overall_detailed[[parameter_value_column]] <- overall_detailed[[parameter_name]]

  phase_detailed[[parameter_name_column]] <- parameter_name
  phase_detailed[[parameter_value_column]] <- phase_detailed[[parameter_name]]

  hormone_frequencies[[parameter_name_column]] <- parameter_name
  hormone_frequencies[[parameter_value_column]] <- hormone_frequencies[[parameter_name]]

  switching_by_time[[parameter_name_column]] <- parameter_name
  switching_by_time[[parameter_value_column]] <- switching_by_time[[parameter_name]]
}

compact_identifier_columns <- c(
  "filename",
  unlist(
    lapply(
      seq_along(varied_parameters),
      function(parameter_number) {
        c(
          paste0("varied_parameter_", parameter_number),
          paste0("parameter_", parameter_number, "_value")
        )
      }
    )
  )
)

parameter_value_columns <- paste0(
  "parameter_",
  seq_along(varied_parameters),
  "_value"
)

# ==============================================================================
# TABLE 1: strategy_summary_overall_compact.csv
# ==============================================================================
# One row per model run.

overall_compact <- overall_detailed %>%
  rename(
    post_stressor_mean_hormone = overall_post_stressor_mean_hormone,
    far_from_stressor_mean_hormone = overall_far_from_stressor_mean_hormone,
    response_amplitude = overall_response_amplitude,
    min_obs_hormone = minimum_observed_hormone,
    max_obs_hormone = maximum_observed_hormone,
    max_obs_below_h_max = maximum_observed_below_h_max,
    percentage_at_obs_min = percentage_at_observed_minimum,
    distinct_hormone_levels = number_of_distinct_hormone_levels
  ) %>%
  mutate(
    min_equals_h0 = min_obs_hormone == h0,
    max_equals_h_max = max_obs_hormone == h_max,
    percentage_intermediate_h = if_else(
      min_obs_hormone < h_max,
      percentage_below_h_max - percentage_at_obs_min,
      0
    )
  ) %>%
  select(
    all_of(compact_identifier_columns),
    distinct_hormone_levels,
    overall_mean_hormone,
    post_stressor_mean_hormone,
    far_from_stressor_mean_hormone,
    response_amplitude,
    min_equals_h0,
    min_obs_hormone,
    percentage_at_obs_min,
    max_equals_h_max,
    max_obs_hormone,
    max_obs_below_h_max,
    percentage_at_h_max,
    percentage_below_h_max,
    percentage_intermediate_h
  ) %>%
  mutate(
    across(
      where(is.numeric) & !all_of(parameter_value_columns),
      ~ round(.x, 2)
    )
  )

# ==============================================================================
# TABLE 2: strategy_summary_by_phase_compact.csv
# ==============================================================================
# One row per model run and breeding phase.

phase_compact <- phase_detailed %>%
  rename(
    min_obs_hormone = minimum_observed_hormone,
    max_obs_hormone = maximum_observed_hormone,
    max_obs_below_h_max = maximum_observed_below_h_max,
    percentage_at_obs_min = percentage_at_observed_minimum,
    distinct_hormone_levels = number_of_distinct_hormone_levels
  ) %>%
  mutate(
    min_equals_h0 = min_obs_hormone == h0,
    max_equals_h_max = max_obs_hormone == h_max,
    percentage_intermediate_h = if_else(
      min_obs_hormone < h_max,
      percentage_below_h_max - percentage_at_obs_min,
      0
    )
  ) %>%
  select(
    all_of(compact_identifier_columns),
    breeding_phase,
    distinct_hormone_levels,
    overall_mean_hormone_in_phase,
    post_stressor_mean_hormone,
    far_from_stressor_mean_hormone,
    response_amplitude,
    min_equals_h0,
    min_obs_hormone,
    percentage_at_obs_min,
    max_equals_h_max,
    max_obs_hormone,
    max_obs_below_h_max,
    percentage_at_h_max,
    percentage_below_h_max,
    percentage_intermediate_h
  ) %>%
  mutate(
    across(
      where(is.numeric) & !all_of(parameter_value_columns),
      ~ round(.x, 2)
    )
  )

# ==============================================================================
# TABLE 3: switching_damage_by_phase_compact.csv
# ==============================================================================
# One row per model run and breeding phase.

switching_compact <- switching_by_time %>%
  arrange(
    across(all_of(compact_identifier_columns)),
    breeding_phase,
    time_since_stressor
  ) %>%
  group_by(
    across(all_of(compact_identifier_columns)),
    breeding_phase
  ) %>%
  mutate(
    far_switch_damage = last(switching_damage),
    matches_far_switch_damage =
      !is.na(far_switch_damage) &
      !is.na(switching_damage) &
      switching_damage == far_switch_damage,
    remains_at_far_switch_damage =
      rev(cumall(rev(matches_far_switch_damage)))
  ) %>%
  summarise(
    number_switch_damage_values = n_distinct(
      switching_damage,
      na.rm = TRUE
    ),
    switch_damage_at_t0 = first(
      switching_damage[time_since_stressor == 0],
      default = NA_real_
    ),
    switch_damage_t0_is_min = if (
      all(is.na(switching_damage))
    ) NA else switch_damage_at_t0 == min(switching_damage, na.rm = TRUE),
    switch_damage_at_t1 = first(
      switching_damage[time_since_stressor == 1],
      default = NA_real_
    ),
    switch_damage_at_t_far = first(far_switch_damage),
    switch_damage_t_far_is_max = if (
      all(is.na(switching_damage))
    ) NA else switch_damage_at_t_far == max(switching_damage, na.rm = TRUE),
    switch_damage_range = if (
      all(is.na(switching_damage))
    ) NA_real_ else
      max(switching_damage, na.rm = TRUE) -
      min(switching_damage, na.rm = TRUE),
    switch_damage_stabilisation_time = if (
      any(remains_at_far_switch_damage)
    ) min(time_since_stressor[remains_at_far_switch_damage]) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric) & !all_of(parameter_value_columns),
      ~ round(.x, 2)
    )
  )

# ==============================================================================
# TABLE 4: hormone_percentage_by_phase.csv
# ==============================================================================
# One row per model run and observed hormone level.

hormone_percentage_by_phase <- hormone_frequencies %>%
  select(
    all_of(compact_identifier_columns),
    breeding_phase,
    hormone,
    percentage_of_states
  ) %>%
  mutate(
    percentage_column = if_else(
      breeding_phase == "all",
      "percentage_all_phases",
      paste0("percentage_phase_", breeding_phase)
    )
  ) %>%
  select(-breeding_phase) %>%
  pivot_wider(
    names_from = percentage_column,
    values_from = percentage_of_states,
    values_fill = 0
  ) %>%
  select(
    all_of(compact_identifier_columns),
    hormone,
    starts_with("percentage_phase_"),
    percentage_all_phases
  ) %>%
  arrange(
    across(all_of(compact_identifier_columns)),
    hormone
  ) %>%
  mutate(
    across(
      where(is.numeric) & !all_of(parameter_value_columns),
      ~ round(.x, 2)
    )
  )

# ==============================================================================
# WRITE THE FOUR COMPACT TABLES
# ==============================================================================
write_delim(
  overall_compact,
  file.path(data_dir, "strategy_summary_overall_compact.csv"),
  delim = ";"
)

write_delim(
  phase_compact,
  file.path(data_dir, "strategy_summary_by_phase_compact.csv"),
  delim = ";"
)

write_delim(
  switching_compact,
  file.path(data_dir, "switching_damage_by_phase_compact.csv"),
  delim = ";"
)

write_delim(
  hormone_percentage_by_phase,
  file.path(data_dir, "hormone_percentage_by_phase.csv"),
  delim = ";"
)

message(
  "Finished writing all four compact summary tables to: ",
  data_dir
)
