#!/usr/bin/env Rscript

# This script creates a heat map for one stress DP output file.
# It displays one panel for each breeding phase. Damage is on the
# y-axis, time since the last attack (tau) is on the x-axis, and
# optimal hormone level is represented by colour from dark purple
# (low hormone) to yellow (high hormone).

# First, we delete all existing variables and objects from the current R session so no old objects can interfere with the script)
rm(list=ls())

# Loading packages
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("filesstrings"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("here"))

# load file with general functions
# to deal with the files coming out of the DP programme
# (There's a file in the repository which already contains lots of useful functions for handling the stress data, so we load it so we can use those functions)
script.dir <- here()

source(file.path(script.dir,
                "stress_file_functions.r"))


#### WORK OUT THE FILENAMES ####

if (!exists("strategy.file.name"))
{
    # get the file containing the data which we want
    # to plot from the command line
    args = commandArgs(trailingOnly=TRUE)
    strategy.file.name <- args[[1]]
}

# Remove the input file's final extension to create the PDF filename.
# using filesstrings::before_last_dot()
strategy.file.base <- before_last_dot(strategy.file.name) 

plot.file.name <- paste0("plot_",basename(strategy.file.base),".pdf")



### GET THE ACTUAL DATA ###

# get the parameters of this run
params <- read.parameters(strategy.file.name)

# get the strategy data
data.strategy <- read.stress.file(strategy.file.name)

### PLOTTING ###
list.plots <- list()

# make annotation in the plot

# first provide a map of variables we want to plot
# (Select the model parameters to display in the title above the plot.)
params.to.list <- c("lambdaL"
                    ,"lambdaA"
                    ,"rho"
                    ,"mu"
                    ,"gamma"
                    ,"pAtt"
                    ,"maxT")

param.str <- ""
for (param.idx in 1:length(params.to.list))
{
  param.i <- params.to.list[[param.idx]]
  
  if (param.i %in% names(params))
  {
    param.str <- paste0(param.str
                        ,ifelse(str_length(param.str) > 0, ", ", "")
                        ,param.i
                        ,": "
                        ,params[param.i])
  }
}

# This code creates a duplicate data table of the strategy (there's already the data.strategy data table object), and it outputs this table as a file, always called test.csv (so it will be overwritten if we run this plotting script again on another strategy file, for example).
# I'm going to comment this code out for now, since I don't need to produce this csv file every time I need a heat map made.

#write.table(x = data.strategy,
#        file="test.csv",
#        sep=";",
#        row.names=F)


# Create a heat map of the optimal hormone strategy.
list.plots <- c(list.plots, list(ggplot(data=data.strategy) +
        geom_tile(aes(x = t, y = d, fill = hormone)) +
          # Note on scaling the axes:
            # Tile centres run from 0 to maxT/maxD - 1, so half-unit
            # boundaries show the complete first and last tiles, including d=0.
          scale_x_continuous(
            limits = c(-0.5, as.numeric(params["maxT"]) - 0.5),
            expand = c(0,0)
          ) +
          scale_y_continuous(
            limits = c(-0.5, as.numeric(params["maxD"]) - 0.5),
            expand = c(0,0)
          ) +
        scale_fill_viridis(option="plasma", limits=c(0,as.numeric(params["maxH"]))) +
        xlab("Time since last attack, tau") +
        theme_classic() + 
        facet_grid(~s) +
        ggtitle("Hormone strategy over tau and damage")))

# Convert lambdaA and lambdaL back into the more intuitive
# environmental risk and autocorrelation measures.
autocorr <- 1.0 - (params["lambdaL"] + params["lambdaA"])
risk <- params["lambdaA"] / (params["lambdaL"] + params["lambdaA"])

param.str <- paste0(param.str,
        " risk: ",risk,",",
        " autocorr: ",autocorr,",")

# Compatibility with older output files: if a parameter called
# n_repro_bout is present, add it to the plot annotation.
# Current model outputs use maxS and do not include n_repro_bout,
# so this block is normally skipped.
if ("n_repro_bout" %in% names(params))
{
    param.str <- paste0(param.str,
            " bout: ",params["n_repro_bout"])
}

# Compatibility with older output files: add any fecundity parameters
# whose names begin with "K" to the plot annotation.
# Current model outputs use fec_scale and fec_baseline, so this block
# normally finds no matching parameters and does nothing.
# find fecundity values
Kcolumns <- grepl("^K",names(params))
fec.cols <- names(params)[Kcolumns]

for (col in fec.cols) 
{
    param.str <- paste0(param.str,", ",col,": ",params[col])
}


all.plots <- wrap_plots(list.plots, ncol=1) + 
    plot_annotation(title=param.str,
            caption=strategy.file.name,
            theme=theme(text=element_text(size=8)
                    )) 
#    plot_layout(heights = unit(c(rep(1,times=4),8),rep("null",times=5)))


# Save the completed heat map and its parameter annotation as a PDF.
ggsave(filename=plot.file.name,height=5,width=15)
