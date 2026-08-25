#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("filesstrings"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("here"))

# load file with general functions
# to deal with the files coming out of the DP programme
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

# obtain file name without final extension (.txt)
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

write.table(x = data.strategy,
        file="test.csv",
        sep=";",
        row.names=F)

# now a levelplot showing the stress levels
list.plots <- c(list.plots, list(ggplot(data=data.strategy) +
        geom_tile(aes(x = t, y = d, fill = hormone)) +
        scale_x_continuous(limits = c(-0.5,as.numeric(params["maxT"])), expand=c(0,0)) +
        scale_y_continuous(limits = c(0,as.numeric(params["maxD"])), expand=c(0,0)) +
        scale_fill_viridis(option="plasma", limits=c(0,as.numeric(params["maxH"]))) +
        xlab("Time since last attack, tau") +
        theme_classic() + 
        facet_grid(~s) +
        ggtitle("Hormone strategy over tau and damage")))


autocorr <- 1.0 - (params["lambdaL"] + params["lambdaA"])
risk <- params["lambdaA"] / (params["lambdaL"] + params["lambdaA"])

param.str <- paste0(param.str,
        " risk: ",risk,",",
        " autocorr: ",autocorr,",")

if ("n_repro_bout" %in% names(params))
{
    param.str <- paste0(param.str,
            " bout: ",params["n_repro_bout"])
}



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

ggsave(filename=plot.file.name,height=5,width=15)
