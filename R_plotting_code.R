library(ggplot2)
library(gridExtra)
library(cowplot)

### STEP 1: MAKING DATA FRAMES #########################################################################################################################################################################

# Creating function to split data produced by simulation so that it 'tracks' individuals based on what breeding phase they began the simulation in (and so what phase they were attacked at time step 10)
remap_scenario <- function(fileName, scenarioName, maxS,
                           closePhases, farPhases)
{
  # Read in data, skipping first line to get to headers
  df <- read.table(fileName, header=TRUE, skip=1, stringsAsFactors=FALSE)
  # Columns in data tables: time, s, nInd, meanD, sdD, meanH, sdH
  
  # Get max time (50 timesteps)
  simTime <- max(df$time)
  
  outDF <- data.frame()
  
  for (sVal in 0:(maxS-1)) {
    
    # Find group's attack phase at t=10
    atkPhase <- (sVal + 10) %% maxS
    
    # "close" vs "far" attack phase groups
    if (sVal %in% closePhases) {
      atkType <- "close"
    } else {
      atkType <- "far"
    }
    
    timeVec  <- 0:simTime
    meanHVec <- numeric(length(timeVec))
    
    for (i in seq_along(timeVec)) {
      t       <- timeVec[i]
      sWanted <- (sVal + t) %% maxS
      rowSub  <- subset(df, time==t & s==sWanted)
      meanHVec[i] <- if(nrow(rowSub)==1) rowSub$meanH else NA
    }
    
    tmp <- data.frame(
      time        = timeVec,
      startPhase  = sVal,
      meanH       = meanHVec,
      scenario    = scenarioName,
      attackPhase = atkPhase,
      attackType  = atkType
    )
    outDF <- rbind(outDF, tmp)
  }
  
  # Build combo factor with type of attack and starting phase
  outDF$attackCombo <- paste0(outDF$attackType, "_s", outDF$attackPhase)
  
  # Identify all unique combos in cycle scenario
  combos <- sort(unique(outDF$attackCombo))
  
  # Separate them into "close_" combos vs "far_" combos
  closeCombos <- combos[grepl("^close_", combos)]
  farCombos   <- combos[grepl("^far_",   combos)]
  
  # Colour ramps
  colClose <- colorRampPalette(c("plum","purple"))(length(closeCombos))
  colFar   <- colorRampPalette(c("powderblue","turquoise4"))(length(farCombos))
  
  comboColors <- c(
    setNames(colClose, closeCombos),
    setNames(colFar,   farCombos)
  )
  
  # Specifying colours for each combo
  outDF$comboColor <- comboColors[outDF$attackCombo]
  
  # Return the final, sorted dataframe
  outDF
}

# S=1
close1 <- c()       # none
far1   <- c(0)
dfS1 <- remap_scenario("SimAttack_SimS1.txt", "S=1", maxS=1,
                       closePhases=close1, farPhases=far1)

# S=5
close5 <- c(3,4)
far5   <- c(0,1,2)
dfS5 <- remap_scenario("SimAttack_SimS5.txt", "S=5", maxS=5,
                       closePhases=close5, farPhases=far5)

# S=10
close10 <- c(7,8,9)
far10   <- c(0,1,2,3,4,5,6)
dfS10 <- remap_scenario("SimAttack_SimS10.txt", "S=10", maxS=10,
                        closePhases=close10, farPhases=far10)

# S=20
close20 <- c(9,8,7,6,5)
far20   <- c(19,18,17,16,15,14,13,12,11,10,4,3,2,1,0)
dfS20 <- remap_scenario("SimAttack_SimS20.txt", "S=20", maxS=20,
                        closePhases=close20, farPhases=far20)

# CSV files
# write.csv(dfS1, file = "dfS1.csv")
# write.csv(dfS5, file = "dfS5.csv")
# write.csv(dfS10, file = "dfS10.csv")
# write.csv(dfS20, file = "dfS20.csv")

### STEP 2: PLOTTING GRAPHS WITH ALL HORMONE LINES, COLOUR CODED BY ATTACK TYPE #########################################################################################

# Finding the different attack phase combos
build_combo_colors <- function(df) {
  combos <- sort(unique(df$attackCombo))
  closeCombos <- combos[grepl("^close_", combos)]
  farCombos   <- combos[grepl("^far_",   combos)]
  
  # If there's only 1 close combo, just use darker colour
  if (length(closeCombos) == 1) {
    colClose <- "purple"
  } else {
    # multiple close combos = colour ramp
    colClose <- colorRampPalette(c("plum","purple"))(length(closeCombos))
  }
  
  # If there's only 1 far combo, just use darker colour
  if (length(farCombos) == 1) {
    colFar <- "turquoise4"
  } else {
    # multiple far combos = colour ramp
    colFar <- colorRampPalette(c("powderblue","turquoise4"))(length(farCombos))
  }
  
  c(
    setNames(colClose, closeCombos),
    setNames(colFar,   farCombos)
  )
}

# Plotting function
plot_scenario_nolegend <- function(dfScenario, plotTitle) {
  
  # maximum time
  maxTime <- max(dfScenario$time)
  
  # build color map for combos
  colMap <- build_combo_colors(dfScenario)
  
  # creating a data frame for the red rectangle from t=9.5..10.5 (predator attack)
  attackDF <- data.frame(
    xmin = 9.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf
  )
  
  # Create the plot
  p <- ggplot() +
    
    # red rectangle around t=10
    geom_rect(
      data = attackDF,
      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
      fill="red",
      alpha=0.2,
      inherit.aes=FALSE,  
      color=NA
    ) +
    
    # dotted grey baseline at y=0.2
    geom_segment(
      aes(x=0, xend=maxTime, y=0.2, yend=0.2),
      color="grey60",
      linetype="dotted",
      size=1
    ) +
    
    # scenario lines
    geom_line(
      data=dfScenario,
      aes(
        x=time,
        y=meanH,
        color=attackCombo,
        group=interaction(startPhase, attackCombo)
      ),
      size=1
    ) +
    
    # color scale for combos (but no legend)
    scale_color_manual(values=colMap, guide="none") +
    
    labs(
      title = plotTitle,
      x     = "Timestep",
      y     = "Mean Hormone Level"
    ) +
    theme_minimal() +
    
    # add margin on right for manual
    theme(
      plot.margin     = margin(t=5, r=100, b=5, l=5, unit="pt"),
      legend.position = "none"
    )
  
  p
}

# Build each scenario's plot
p1  <- plot_scenario_nolegend(dfS1,  "a) Breeding cycle length S=1")
p5  <- plot_scenario_nolegend(dfS5,  "b) Breeding cycle length S=5")
p10 <- plot_scenario_nolegend(dfS10, "c) Breeding cycle length S=10")
p20 <- plot_scenario_nolegend(dfS20, "d) Breeding cycle length S=20")

# Arrange them in 2x2 grid
grid.arrange(p1, p5, p10, p20, ncol=2)

### STEP 3: PREDATOR RISK CALCULATIONS ################################################################################################################

# Define parameters
lambdaL <- 0.065   # Leaving probability
lambdaA <- 0.035   # Arriving probability
pAtt    <- 0.5     # Probability predator attacks if present

# Define the time range (max 24)
tau_values <- 0:24

# Initialize a vector to store P_pred for each tau
pPred <- numeric(length(tau_values))  # numeric vector length 24

# Compute P_pred(1)
pPred[1] <- 1
pPred[2] <- 1 - lambdaL

# For tau > 1, use equation
for (t in 3:25) {
  pPred[t] <- ( pPred[t-1]*(1 - pAtt)*(1 - lambdaL) + (1 - pPred[t-1])*lambdaA ) /
    ( 1 - pPred[t-1]*pAtt )
}

# Combine into a data frame
riskDF <- data.frame(
  tau    = tau_values,
  P_pred = pPred
)

## STEP 4: SCENARIO DATAFRAMES WITH RISK ###########################################################################################################################

## dfFinal5 #################################################################

# Subset your data for S=5, starting phase=3
subdat5 <- subset(dfS5, startPhase == 3)  # group data

# Background phase data
maxTime5 <- max(subdat5$time)
bgDF5 <- data.frame(t = 0:maxTime5)
bgDF5$phase <- (3 + bgDF5$t) %% 5
bgDF5$xmin  <- bgDF5$t
bgDF5$xmax  <- bgDF5$t + 1

# Matching timestep to tua value
timeVals5 <- 0:maxTime5
tauVals5  <- numeric(length(timeVals5))

for (i in seq_along(timeVals5)) {
  t <- timeVals5[i]
  if (t < 10) {
    tauVals5[i] <- 24
  } else {
    # how many steps since attack at t=10?
    elapsed5 <- t - 10
    # limit of 24
    tauVals5[i] <- min(elapsed5, 24)
  }
}

timeToTau5 <- data.frame(time = timeVals5, tau = tauVals5)

# Merge into data frame
timeRisk5 <- merge(timeToTau5, riskDF, by="tau", all.x=TRUE)

timeRisk5 <- timeRisk5[order(timeRisk5$time),]

# Now merge scenario hormone subdat (time, meanH) with timeRisk (time, tau, P_pred)
dfFinal5 <- merge(subdat5, timeRisk5, by="time", all.x=TRUE)

## dfFinal10 #################################################################

# Subset for S=10, starting phase=7
subdat10 <- subset(dfS10, startPhase == 7)  # group data

# Background phase data
maxTime10 <- max(subdat10$time)
bgDF10 <- data.frame(t = 0:maxTime10)
bgDF10$phase <- (7 + bgDF10$t) %% 10
bgDF10$xmin  <- bgDF10$t
bgDF10$xmax  <- bgDF10$t + 1

timeVals10 <- 0:maxTime10
tauVals10  <- numeric(length(timeVals10))

for (i in seq_along(timeVals10)) {
  t <- timeVals10[i]
  if (t < 10) {
    tauVals10[i] <- 24
  } else {
    # how many steps since attack at t=10?
    elapsed10 <- t - 10
    tauVals10[i] <- min(elapsed10, 24)
  }
}

timeToTau10 <- data.frame(time = timeVals10, tau = tauVals10)
timeRisk10 <- merge(timeToTau10, riskDF, by="tau", all.x=TRUE)
timeRisk10 <- timeRisk10[order(timeRisk10$time),]
dfFinal10 <- merge(subdat10, timeRisk10, by="time", all.x=TRUE)

## dfFinal20 #################################################################

# Subset for S=20, starting phase=6
subdat20 <- subset(dfS20, startPhase == 6)  # group data

# Background phase data
maxTime20 <- max(subdat20$time)
bgDF20 <- data.frame(t = 0:maxTime20)
bgDF20$phase <- (6 + bgDF20$t) %% 20
bgDF20$xmin  <- bgDF20$t
bgDF20$xmax  <- bgDF20$t + 1

timeVals20 <- 0:maxTime20
tauVals20  <- numeric(length(timeVals20))

for (i in seq_along(timeVals20)) {
  t <- timeVals20[i]
  if (t < 10) {
    tauVals20[i] <- 24
  } else {
    elapsed20 <- t - 10
    tauVals20[i] <- min(elapsed20, 24)
  }
}

timeToTau20 <- data.frame(time = timeVals20, tau = tauVals20)
timeRisk20 <- merge(timeToTau20, riskDF, by="tau", all.x=TRUE)
timeRisk20 <- timeRisk20[order(timeRisk20$time),]
dfFinal20 <- merge(subdat20, timeRisk20, by="time", all.x=TRUE)


### STEP 5: PLOTTING WITH RISK DATA ##################################################################################################

#######################
# PLOT 1 (S=5)
#######################
# colour ramp
phaseLevels5 <- 0:4
phaseColors5 <- colorRampPalette(c("turquoise3","purple"))(length(phaseLevels5))
names(phaseColors5) <- as.character(phaseLevels5)

p1 <- ggplot() +
  # background rectangles (phases)
  geom_rect(
    data = bgDF5,
    aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=factor(phase)),
    alpha=0.2, color=NA
  ) +
  # Predator attack risk line = red
  geom_line(
    data = dfFinal5,
    aes(
      x=time, 
      y=P_pred, 
      color="Predator attack risk", 
      linetype="Predator attack risk"
    ),
    size=0.8
  ) +
  # Mean hormone line = black
  geom_line(
    data = dfFinal5,
    aes(
      x=time, 
      y=meanH,
      color="Mean hormone level",
      linetype="Mean hormone level"
    ),
    size=1
  ) +
  # Baseline hormone line = grey
  geom_segment(
    aes(
      x=0, xend=50,
      y=0.2, yend=0.2,
      color="Baseline hormone level",
      linetype="Baseline hormone level"
    ),
    size=0.8
  ) +
  # fill scale for phases
  scale_fill_manual(
    name   = "Breeding Phase",
    values = phaseColors5
  ) +
  # colour scale for lines
  scale_color_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "black",
      "Baseline hormone level" = "grey70",
      "Predator attack risk"   = "red"
    )
  ) +
  # linetype scale
  scale_linetype_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "solid",
      "Baseline hormone level" = "dotted",
      "Predator attack risk"   = "dashed"
    )
  ) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title = "a) Stress response during phase s=3, breeding cycle length S=5",
    x     = "Time",
    y     = "Hormone / Risk"
  ) +
  theme_minimal() +
  # Hide legend
  theme(legend.position="none")


#######################
# PLOT 2 (S=10)
#######################
phaseLevels10 <- 0:9
phaseColors10 <- colorRampPalette(c("turquoise3","purple"))(length(phaseLevels10))
names(phaseColors10) <- as.character(phaseLevels10)

p2 <- ggplot() +
  geom_rect(
    data = bgDF10,
    aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=factor(phase)),
    alpha=0.2, color=NA
  ) +
  geom_line(
    data = dfFinal10,
    aes(
      x=time, y=P_pred, 
      color="Predator attack risk", 
      linetype="Predator attack risk"
    ),
    size=0.8
  ) +
  geom_line(
    data = dfFinal10,
    aes(
      x=time, y=meanH,
      color="Mean hormone level",
      linetype="Mean hormone level"
    ),
    size=1
  ) +
  geom_segment(
    aes(
      x=0, xend=50,
      y=0.2, yend=0.2,
      color="Baseline hormone level",
      linetype="Baseline hormone level"
    ),
    size=0.8
  ) +
  scale_fill_manual(
    name="Breeding Phase",
    values=phaseColors10
  ) +
  scale_color_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "black",
      "Baseline hormone level" = "grey70",
      "Predator attack risk"   = "red"
    )
  ) +
  scale_linetype_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "solid",
      "Baseline hormone level" = "dotted",
      "Predator attack risk"   = "dashed"
    )
  ) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title="b) Stress response during phase s=7, breeding cycle length S=10",
    x="Time",
    y="Hormone / Risk"
  ) +
  theme_minimal() +
  theme(legend.position="none")


#######################
# PLOT 3 (S=20)
#######################
phaseLevels20 <- 0:19
phaseColors20 <- colorRampPalette(c("turquoise3","purple"))(length(phaseLevels20))
names(phaseColors20) <- as.character(phaseLevels20)

p3 <- ggplot() +
  geom_rect(
    data = bgDF20,
    aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=factor(phase)),
    alpha=0.2, color=NA
  ) +
  geom_line(
    data = dfFinal20,
    aes(
      x=time, y=P_pred,
      color="Predator attack risk", 
      linetype="Predator attack risk"
    ),
    size=0.8
  ) +
  geom_line(
    data = dfFinal20,
    aes(
      x=time, y=meanH,
      color="Mean hormone level",
      linetype="Mean hormone level"
    ),
    size=1
  ) +
  geom_segment(
    aes(
      x=0, xend=50,
      y=0.2, yend=0.2,
      color="Baseline hormone level",
      linetype="Baseline hormone level"
    ),
    size=0.8
  ) +
  scale_fill_manual(
    name="Breeding Phase",
    values=phaseColors20
  ) +
  scale_color_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "black",
      "Baseline hormone level" = "grey70",
      "Predator attack risk"   = "red"
    )
  ) +
  scale_linetype_manual(
    name="Line",
    values=c(
      "Mean hormone level"     = "solid",
      "Baseline hormone level" = "dotted",
      "Predator attack risk"   = "dashed"
    )
  ) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title="c) Stress response during phase s=16, breeding cycle length S=20",
    x="Time",
    y="Hormone / Risk"
  ) +
  theme_minimal() +
  theme(legend.position="none")


##############################################################################
# Tiny legend plot just for legend basically
##############################################################################
dummy_df <- data.frame(
  x=1:4,
  y=1:4,
  lineType=c(
    "Mean hormone level",
    "Baseline hormone level",
    "Predator attack risk",
    "Predator attack risk"
  ),
  phase=c("0","0","final","final")
)

legend_plot <- ggplot(dummy_df, aes(x=x, y=y)) +
  # lines (3 categories):
  geom_line(aes(color=lineType, linetype=lineType), size=1) +
  # squares for phases:
  geom_point(aes(fill=phase), shape=22, size=8, alpha=0.2, color=NA) +
  
  # colour scale
  scale_color_manual(
    name="Line",
    breaks=c("Mean hormone level", "Baseline hormone level", "Predator attack risk"),
    labels=c(
      "Mean hormone level"     = "Mean hormone level",
      "Baseline hormone level" = "Baseline hormone level",
      "Predator attack risk"   = "Predator attack risk"
    ),
    values=c(
      "Mean hormone level"     = "black",
      "Baseline hormone level" = "grey70",
      "Predator attack risk"   = "red"
    )
  ) +
  # linetype scale
  scale_linetype_manual(
    name="Line",
    breaks=c("Mean hormone level", "Baseline hormone level", "Predator attack risk"),
    labels=c(
      "Mean hormone level"     = "Mean hormone level",
      "Baseline hormone level" = "Baseline hormone level",
      "Predator attack risk"   = "Predator attack risk"
    ),
    values=c(
      "Mean hormone level"     = "solid",
      "Baseline hormone level" = "dotted",
      "Predator attack risk"   = "dashed"
    )
  ) +
  # fill scale for phases = just 2 examples: "0", "final"
  scale_fill_manual(
    name="Phase in breeding cycle",
    breaks=c("0","final"),
    labels=c("0"="Reproductive phase, s = 0","final"="Final non-reproductive phase, s = S-1"),
    values=c("0"="turquoise3","final"="purple")
  ) +
  
  theme_void() +
  theme(
    legend.box = "vertical",
    legend.position = "right",
    legend.key.size = unit(1.2,"cm"),
    legend.text     = element_text(size=12),
    legend.title    = element_text(size=14, face="bold")
  )

# Extract only the legend as a grob
legend_only <- get_legend(legend_plot)

##############################################################################
# Arrange plots and legend
##############################################################################
layout_mat <- rbind(
  c(1,2),
  c(3,4)
)

grid.arrange(
  p1, p2, p3, legend_only,
  layout_matrix = layout_mat,
  widths  = c(1,1),
  heights = c(1,1)
)


### STEP 6: PEAK HORMONE DATABASE ###########################################################################################

extract_peaks_df <- function(df, scenarioName) {
  # Identify unique (attackType, attackPhase) pairs
  uniquePairs <- unique(df[, c("attackType", "attackPhase")])
  
  # Create an empty data frame to store results
  peak_df <- data.frame(
    Scenario                    = character(),
    Attack_type                 = character(),
    Attack_phase               = numeric(),
    Peak_Hormone_level_t10     = numeric(),
    Peak_Hormone_level_t_over10 = numeric(), 
    stringsAsFactors = FALSE
  )
  
  # Loop over each unique pair
  for (i in seq_len(nrow(uniquePairs))) {
    atType  <- uniquePairs$attackType[i]
    atPhase <- uniquePairs$attackPhase[i]
    
    # Subset rows for this (attackType, attackPhase)
    subData <- subset(df, attackType==atType & attackPhase==atPhase)
    
    #'Peak_Hormone_level_t10' => meanH at time=10
    sub_t10 <- subset(subData, time==10)
    peakH_t10 <- if (nrow(sub_t10)>0) sub_t10$meanH[1] else NA_real_
    
    #'Peak_Hormone_level_t_over10' => max of meanH from time<10
    sub_after10 <- subset(subData, time>=11 & time<=50)
    peakH_after10 <- if (nrow(sub_after10)>0) max(sub_after10$meanH) else NA_real_
    
    # Add a single row to the result
    peak_df <- rbind(
      peak_df,
      data.frame(
        Scenario                     = scenarioName,
        Attack_type                  = atType,
        Attack_phase                = atPhase,
        Peak_Hormone_level_t10      = peakH_t10,
        Peak_Hormone_level_t_over10= peakH_after10,
        stringsAsFactors = FALSE
      )
    )
  }
  
  peak_df
}


# For dfS1
peak_dfS1 <- extract_peaks_df(dfS1, "S=1")
# For dfS5
peak_dfS5 <- extract_peaks_df(dfS5, "S=5")
# For dfS10
peak_dfS10 <- extract_peaks_df(dfS10, "S=10")
# For dfS20
peak_dfS20 <- extract_peaks_df(dfS20, "S=20")

peak_df_all <- rbind(peak_dfS1, peak_dfS5, peak_dfS10, peak_dfS20)
head(peak_df_all)


### STEP 7: PEAK BAR CHART ###############################################################################################

# Finda far and close attack types
create_three_cat <- function(df) {
  out <- data.frame()  # empty
  
  for(i in seq_len(nrow(df))) {
    rowi <- df[i,]
    sc   <- rowi$Scenario
    at   <- rowi$Attack_type
    
    if(at == "far") {
      # only one bar: far_t10
      out <- rbind(out, data.frame(
        Scenario = sc,
        cat3     = "far_t10",
        # want rowi$Peak_Hormone_level_t10
        peakVal  = rowi$Peak_Hormone_level_t10
      ))
    } else if(at == "close") {
      # 2 bars for "close": one at t10, one tOver10
      out <- rbind(out,
                   # (1) close_t10
                   data.frame(
                     Scenario = sc,
                     cat3     = "close_t10",
                     peakVal  = rowi$Peak_Hormone_level_t10
                   ),
                   # (2) close_tUnder10
                   data.frame(
                     Scenario = sc,
                     cat3     = "close_tOver10",
                     peakVal  = rowi$Peak_Hormone_level_t_over10
                   )
      )
    }
  }
  out
}

# Apply to peak database
peakBarsAll <- create_three_cat(peak_df_all)

# Subset for only S=5, S=10, S=20
subsetBars <- peakBarsAll %>%
  filter(Scenario %in% c("S=5","S=10","S=20"))

# Convert 'Scenario' into a factor 
subsetBars$Scenario <- factor(
  subsetBars$Scenario,
  levels = c("S=5","S=10","S=20")
)

# Ensure cat3 is also a factor with chosen order
subsetBars$cat3 <- factor(
  subsetBars$cat3,
  levels = c("far_t10","close_t10","close_tOver10")
)

# Bar chart
ggplot(subsetBars, aes(x=Scenario, y=peakVal, fill=cat3, group=cat3)) +
  stat_summary(
    fun="mean",
    geom="col",
    width=0.5,
    position=position_dodge(width=0.7)
  ) +

  # Place numeric text on top of each bar
  # stat_summary(
  #  fun = mean,
  #  aes(label = round(after_stat(y), 2)),   # round to 2 decimal places
  #  geom = "text",
  #  position = position_dodge(width=0.7),
  #  vjust = -0.5,     # move text above the bar
  #  color = "black",  # text color
  #  size = 4          # text size
  #) +
  
  scale_fill_manual(
    name = "Peak Category",      # legend title
    # Match the breaks to factor levels
    breaks = c("far_t10", "close_t10", "close_tOver10"),
    # Provide corresponding labels
    labels = c(
      "Attacked far from cycle end",             # 1) far_t10
      "Attacked close to cycle end: initial peak",    # 2) close_t10
      "Attacked close to cycle end: sendondry peak"  # 3)close_tOver10
    ),
    # define colour values
    values = c(
      "far_t10"       = "turquoise3",
      "close_t10"     = "mediumpurple3",
      "close_tOver10" = "plum"
    )
  ) +
  labs(
    title = "",
    x     = "Breeding Cycle Length",
    y     = "Mean Peak Hormone Level"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(
      size = 14,                 # bigger font
      margin = margin(t = 10)    # extra space on top
    ),
    axis.title.y = element_text(
      size = 14,                 # bigger font
      margin = margin(r = 10)    # extra space on the right
    )
  )

### STEP 8: HORMONE & DAMAGE DATAFRAMES ######################################################################################

# Creating function to split data produced by simulation. Extracts both meanH and meanD.
remap_scenario <- function(fileName, scenarioName, maxS,
                           closePhases, farPhases)
{
  #Read the data
  df <- read.table(fileName, header=TRUE, skip=1, stringsAsFactors=FALSE)

  simTime <- max(df$time)
  
  outDF <- data.frame()
  
  # For each possible starting breeding phase sVal
  for (sVal in 0:(maxS-1)) {
    
    # group's attack phase at t=10
    atkPhase <- (sVal + 10) %% maxS
    
    # "close" vs "far" based on whether sVal is in closePhases
    if (sVal %in% closePhases) {
      atkType <- "close"
    } else {
      atkType <- "far"
    }
    
    # vectors for meanH and meanD over time=0 to simTime
    timeVec    <- 0:simTime
    meanHVec   <- numeric(length(timeVec))
    meanDVec   <- numeric(length(timeVec))
    
    # Reconstruct group's time series
    for (i in seq_along(timeVec)) {
      t       <- timeVec[i]
      sWanted <- (sVal + t) %% maxS
      rowSub  <- subset(df, time==t & s==sWanted)
      
      if (nrow(rowSub)==1) {
        meanHVec[i] <- rowSub$meanH
        meanDVec[i] <- rowSub$meanD
      } else {
        meanHVec[i] <- NA
        meanDVec[i] <- NA
      }
    }
    
    # Build a partial data frame for group
    tmp <- data.frame(
      time        = timeVec,
      startPhase  = sVal,
      meanH       = meanHVec,
      meanD       = meanDVec,
      scenario    = scenarioName,
      attackPhase = atkPhase,
      attackType  = atkType,
      stringsAsFactors = FALSE
    )
    
    outDF <- rbind(outDF, tmp)
  }
  
  # Add a combined column "attackCombo"
  outDF$attackCombo <- paste0(outDF$attackType, "_s", outDF$attackPhase)
  
  # Return the final data frame
  outDF
}

# S=1
close1 <- c()   # none
far1   <- c(0)
dfS1_dh <- remap_scenario("SimAttack_SimS1.txt", "S=1", maxS=1,
                          closePhases=close1, farPhases=far1)

# S=5
close5 <- c(3,4)
far5   <- c(0,1,2)
dfS5_dh <- remap_scenario("SimAttack_SimS5.txt", "S=5", maxS=5,
                          closePhases=close5, farPhases=far5)

# S=10
close10 <- c(7,8,9)
far10   <- c(0,1,2,3,4,5,6)
dfS10_dh <- remap_scenario("SimAttack_SimS10.txt", "S=10", maxS=10,
                           closePhases=close10, farPhases=far10)

# S=20
close20 <- c(9,8,7,6,5)
far20   <- c(19,18,17,16,15,14,13,12,11,10,4,3,2,1,0)
dfS20_dh <- remap_scenario("SimAttack_SimS20.txt", "S=20", maxS=20,
                           closePhases=close20, farPhases=far20)

# Mean dataframe for S5
meanS5 <- dfS5_dh %>%
  group_by(time, attackType) %>%
  summarise(
    meanH = mean(meanH, na.rm=TRUE),
    meanD = mean(meanD, na.rm=TRUE),
    .groups="drop"
  )

# Scale damage by 1/100 
meanS5$damageScaled <- meanS5$meanD / 100

### STEP 9: HORMONE DAMAGE PLOTS ##########################################################################################################

# Red attack bar data
attackDF <- data.frame(
  xmin = 9.5,
  xmax = 10.5,
  ymin = -Inf,   
  ymax = Inf,
  event= "Predator Attack"  # used for fill legend
)

ggplot(meanS5, aes(x=time)) +
  
  # red rectangle around t=10
  geom_rect(
    data = attackDF,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=event),
    inherit.aes = FALSE,   
    alpha=0.2,
    color=NA
  ) +
  
  # Close hormone
  geom_line(
    data = subset(meanS5, attackType=="close"),
    aes(
      y=meanH,
      color="Close Hormone",
      linetype="Close Hormone"
    ),
    size=1
  ) +
  
  # Close damage
  geom_line(
    data = subset(meanS5, attackType=="close"),
    aes(
      y=damageScaled,
      color="Close Damage",
      linetype="Close Damage"
    ),
    size=1
  ) +
  
  # Far hormone
  geom_line(
    data = subset(meanS5, attackType=="far"),
    aes(
      y=meanH,
      color="Far Hormone",
      linetype="Far Hormone"
    ),
    size=1
  ) +
  
  # Far damage 
  geom_line(
    data = subset(meanS5, attackType=="far"),
    aes(
      y=damageScaled,
      color="Far Damage",
      linetype="Far Damage"
    ),
    size=1
  ) +
  
  # Now define single legend for all four lines
  scale_color_manual(
    name   = "Line",
    breaks = c("Close Hormone","Close Damage","Far Hormone","Far Damage"),
    labels = c(
      "Close Hormone" = "Attack close to repro: Hormone",
      "Close Damage"  = "Attack close to repro: Damage",
      "Far Hormone"   = "Attack far from repro: Hormone",
      "Far Damage"    = "Attack far from repro: Damage"
    ),
    values = c(
      "Close Hormone"="purple",
      "Close Damage" ="purple4",
      "Far Hormone"  ="turquoise2",
      "Far Damage"   ="turquoise4"
    )
  ) +
  scale_linetype_manual(
    name   = "Line",
    breaks = c("Close Hormone","Close Damage","Far Hormone","Far Damage"),
    labels = c(
      "Close Hormone" = "Attack close to repro: Hormone",
      "Close Damage"  = "Attack close to repro: Damage",
      "Far Hormone"   = "Attack far from repro: Hormone",
      "Far Damage"    = "Attack far from repro: Damage"
    ),
    values = c(
      "Close Hormone"="solid",
      "Close Damage" ="dashed",
      "Far Hormone"  ="solid",
      "Far Damage"   ="dashed"
    )
  ) +
  
  # Fill scale for the rectangle
  scale_fill_manual(
    name   = "Event",                    # legend title
    breaks = "Predator Attack",          # factor level
    values = c("Predator Attack"="red"), # colour for rectangle
    labels = c("Predator Attack"="Predator Attack") 
  ) +
  
  labs(
    title = "",
    x="Time",
    y="Mean Hormone / Damage Level"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(
      size = 14,                 # bigger font
      margin = margin(t = 10)    # extra space on top
    ),
    axis.title.y = element_text(
      size = 14,                 # bigger font
      margin = margin(r = 10)    # extra space on the right
    )
  )



