#!/usr/bin/env Rscript
library("lubridate")
library("datetime")

# risk and autocorrelation determining stressor strength
risk <- c(0.1,0.5,0.9)
autocorr <- seq(0.51,0.99,0.01)

# mortality probability, 1/mu is the lifespan
mu <- c(0.001,0.01,0.05,0.1,0.5)

# effects of deviations from naturally selected optimum
# hormone level on damage build up
#omega <- c(0.1,0.01,0.5)
omega <- c(0.1)

# cost of damage on reproductive output
gamma_g <- c(0.1,0.5, 0.01)

# maximum time steps until (and including) next reproductive episode 
maxS <- c(5,10,20,30)

# fixed rate of repair
rho <- c(0.1, 0.5, 0.01, 1.0)

# current_time, formatted_time for the basename
file_prefix <- "sim_stress_life"
current_time <- now()
formatted_time <- format(current_time, "%Y%m%d_%H%M%S")
basename_dp <- paste0(file_prefix, "_", formatted_time)
basename_sim <- paste0(file_prefix, "_", formatted_time)

ctr <- 0

to_screen <- 0

exe <- "./stress_lh.exe"

for (risk_i in risk)
{
    for (autocorr_i in autocorr)
    {
        lambdaL = (1.0 - autocorr_i) / (1.0 + (risk_i/(1.0 - risk_i)))
        lambdaA = 1.0 - lambdaL - autocorr_i   
        for (mu_i in mu)
        {
            for (omega_i in omega)
            {
                for (gamma_g_i in gamma_g)
                {
                    for (maxS_i in maxS)
                    {
                        for (rho_i in rho)
                        {
                            writeLines(text = paste0(
                                            exe," ",
                                            "lambdaA=",lambdaA," ",
                                            "lambdaL=",lambdaL," ",
                                            "mu=",mu_i," ",
                                            "rho=",rho_i," ",
                                            "omega=",omega_i," ",
                                            "gamma_g=",gamma_g_i," ",
                                            "maxS=",maxS_i," ",
                                            "stdoutput=",to_screen," ",
                                            "dpFile=",basename_dp,"_",ctr,".csv ",
                                            "simBase=",basename_sim,"_",ctr,".csv"))
                            ctr <- ctr + 1
                        }
                    }
                }
            }
        }
    }
}
