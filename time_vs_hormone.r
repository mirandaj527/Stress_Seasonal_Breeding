#!/usr/bin/env Rscript
# loop over a bunch of strategy files and
# summarize relationship between season and
# hormone level when it comes to time since last
# attack vs damage

# get rid of everything there
rm(list=ls())

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("filesstrings"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("tidyverse"))

# please adjust the number of cores to what you have
# on your computer. If you run this on Rstudio, you
# could go v high, assuming that it will not take mega long
# try it out!
registerDoParallel(cores = 5)

# load file with general functions
# to deal with the files coming out of the DP programme
script.dir <- here()

source(file.path(script.dir,
                "stress_file_functions.r"))

# create a data.frame with all the file names that contain the
# stress strategy data
strategy_summary <- data.frame(file = list.files(path=".",
        pattern="^sim_stress_life"))

# function to summarise a single file
# we'll later call this function on each file in the data.frame above
summarize_single_file <- function(filename) {
    
    # get the parameters of this run
    params <- read.parameters(filename)

    params$autocorrelation = 1.0 - params$lambdaA - params$lambdaL
    params$risk = params$lambdaA / (params$lambdaL + params$lambdaA)

    # get the strategy data
    data.strategy <- as_tibble(read.stress.file(filename)) 

    t_attack = 0
    t_baseline = params$maxT - 1
    
    # now calculate averages per season, we average over all levels of 
    # damage just to make start
    data.strategy.per.s.attack <- data.strategy %>% filter(t == 0) %>% group_by(s) %>%
        summarize(mean_hormone = mean(hormone)) %>%
        pivot_wider(names_from = s, values_from = mean_hormone, names_prefix="season_attack_")

    data.strategy.per.s.baseline <- data.strategy %>% filter(t == t_baseline) %>% group_by(s) %>%
        summarize(mean_hormone = mean(hormone)) %>%
        pivot_wider(names_from = s, values_from = mean_hormone, names_prefix="season_baseline_")

    # combine parameters and season result columns
    result <- bind_cols(params, data.strategy.per.s.attack, data.strategy.per.s.baseline)

    # return result
    return(result)
} # end function summarize_single_file

# now run the summarize single file function and automagically
# it will generate a data.frame with parameters and results
# for each simulation run and put it all in a data.frame or tibble or whatevz

# notice that foreach is parallel, using the %dopar% keyword means that
# multiple cores are processing things simultaneously
outcome <- foreach (strategy_file_idx=seq(1,nrow(strategy_summary),1), 
        .combine=bind_rows) %dopar% {

    the_file <- strategy_summary[strategy_file_idx,"file"]

    if (!is.na(the_file))
    {
        summarize_single_file(filename = the_file)
    }
}

write_delim(x = outcome, file="summary_seasonal.csv", delim=";")
