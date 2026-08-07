library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(truncnorm)
library(cmdstanr)

optimise_model <- "
data {
  int<lower=0> nTimesteps;
  real<lower=1> N0w;
  real<lower=0> N0c;
  real<lower=0> base_lambdaw;
  real<lower=0> flu_lambdaw;
  real<lower=0> lambdac;
  real<lower=0.001, upper=0.999> phi;
  real<lower=0.001, upper=0.999> psi;
  array[nTimesteps-1] int<lower=0, upper=1> flu;
}

parameters {
  array[nTimesteps-1] real<lower=0, upper=1> gamma;
  array[nTimesteps-1] real<lower=0, upper=1> epsilon;
}

model {
  vector[nTimesteps] Nw;
  vector[nTimesteps] Nc;
  
  Nw[1] = N0w;
  Nc[1] = N0c;
  for (t in 2:nTimesteps){
    Nw[t] = (Nw[t-1]*flu[t-1]*flu_lambdaw*(1-epsilon[t-1]))+
    (Nw[t-1]*(1-flu[t-1])*base_lambdaw*(1-epsilon[t-1]))+
    (Nc[t-1]*(gamma[t-1])*lambdac*psi);
    Nc[t] = (Nc[t-1]*(1-gamma[t-1])*lambdac)+
    (Nw[t-1]*flu[t-1]*flu_lambdaw*epsilon[t-1]*phi)+
    (Nw[t-1]*(1-flu[t-1])*base_lambdaw*epsilon[t-1]*phi);
  }
  
  target += Nw[nTimesteps];
}
"
optimise_model_path <- write_stan_file(optimise_model)
optimise_model_compiled <- cmdstan_model(optimise_model_path)

optimise_params <- function(model_data, stanModel, noptreps, nsimreps, updateProgress){
  
  all_data <- data.frame()
  for (i in 1:noptreps){
    n0_w <- rtruncnorm(1, 0, Inf,
                       mean = model_data$n_0_mu, 
                       sd = model_data$n_0_sd)
    n0_c <- model_data$n_0_c
    base_lambdaw <- rtruncnorm(1, 0, 1,
                               mean = model_data$mu_base_lambdaw,
                               sd = model_data$sd_base_lambdaw)
    flu_lambdaw <- rtruncnorm(1, 0, 1,
                              mean = model_data$mu_flu_lambdaw,
                              sd = model_data$sd_flu_lambdaw)
    lambdac <- rtruncnorm(1, 0, 1,
                          mean = model_data$captive_lambda_mu,
                          sd = model_data$captive_lambda_sd)
    phi <- rtruncnorm(1, 0, 1,
                      mean = model_data$phi_mu,
                      sd = model_data$phi_sd)
    psi <- rtruncnorm(1, 0, 1,
                      mean = model_data$psi_mu,
                      sd = model_data$psi_sd)
    flu<-vector()
    for (t in 1:9){
      flu[t] <- as.numeric(runif(1,0,1) < model_data[["pr_flu"]][t])
    }
    
    optData <- list(
      nTimesteps = 10,
      N0w = n0_w,
      N0c = n0_c,
      base_lambdaw = base_lambdaw,
      flu_lambdaw = flu_lambdaw,
      lambdac = lambdac,
      phi = phi,
      psi = psi,
      flu = flu
    )
    
    try({
      optModel <- optimise_model_compiled$optimize(
        data = optData,
        init = list(list(gamma = rep(0.1, 9),
                         epsilon = rep(0.1, 9),
                         N0w = n0_w,
                         base_lambdaw = base_lambdaw,
                         flu_lambdaw = flu_lambdaw,
                         lambdac = lambdac,
                         phi = phi,
                         psi = psi,
                         flu = flu))
      )
      
      rowSummary <- optModel$summary() |>
        filter_out(variable == "lp__") |> 
        separate_wider_delim(variable, "[", names = c("var", "timestep")) |>
        mutate(timestep = as.factor(str_remove(timestep, "\\]"))) |>
        mutate(id = i) 
      all_data <- bind_rows(rowSummary, all_data)
    })
    updateProgress(detail = paste0(as.character(i), "/1000"))
  }
  
  outPlot <- all_data |>
    ggplot(aes(x = timestep, y = estimate, colour = var)) +
    geom_boxplot()
  
  medianVals <- all_data |>
    summarise(estimate = median(estimate),
              .by = c(var, timestep))
  
  model_data$ctw <- medianVals$estimate[1:9]
  model_data$wtc <- medianVals$estimate[10:18]
  
  all_data <- data.frame()
  for (i in 1:nsimreps){
    n0_w <- rtruncnorm(1,
                       mean = model_data$n_0_mu, 
                       sd = model_data$n_0_sd,
                       a = 0, b = Inf)
    n0_c <- model_data$n_0_c
    base_lambdaw <- rtruncnorm(1,
                               mean = model_data$mu_base_lambdaw,
                               sd = model_data$sd_base_lambdaw,
                               a = 0, b = 1)
    flu_lambdaw <- rtruncnorm(1,
                              mean = model_data$mu_flu_lambdaw,
                              sd = model_data$sd_flu_lambdaw,
                              a = 0, b = 1)
    lambdac <- rtruncnorm(1,
                          mean = model_data$captive_lambda_mu,
                          sd = model_data$captive_lambda_sd,
                          a = 0, b = 1)
    phi <- rtruncnorm(1,
                      mean = model_data$phi_mu,
                      sd = model_data$phi_sd,
                      a = 0, b = 1)
    psi <- rtruncnorm(1,
                      mean = model_data$psi_mu,
                      sd = model_data$psi_sd,
                      a = 0, b = 1)
    flu<-vector()
    for (t in 1:9){
      flu[t] <- as.numeric(runif(1,0,1) < model_data[["pr_flu"]][t])
    }
    row_data <- data.frame(id = i,
                           t = 1:10, 
                           nw = NA,
                           nc = NA,
                           n_base = NA)
    row_data[row_data$t == 1, "nw"] <- n0_w
    row_data[row_data$t == 1, "n_base"] <- n0_w
    row_data[row_data$t == 1, "nc"] <- n0_c
    for (t in 2:10){
      row_data$nw[row_data$t==t] <- 
        row_data$nw[row_data$t==t-1]*base_lambdaw*(1-flu[t-1])*(1-model_data[["wtc"]][t-1]) +
        row_data$nw[row_data$t==t-1]*flu_lambdaw*(flu[t-1])*(1-model_data[["wtc"]][t-1]) +
        row_data$nc[row_data$t==t-1]*lambdac*(model_data[["ctw"]][t-1])*psi 
      
      row_data$nc[row_data$t==t] <- 
        row_data$nc[row_data$t==t-1]*lambdac*(1-model_data[["ctw"]][t-1]) +
        row_data$nw[row_data$t==t-1]*base_lambdaw*(1-flu[t-1])*(model_data[["wtc"]][t-1]*phi) +
        row_data$nw[row_data$t==t-1]*flu_lambdaw*(flu[t-1])*(model_data[["wtc"]][t-1]*phi) 
      
      row_data$n_base[row_data$t==t] <- 
        row_data$n_base[row_data$t==t-1]*base_lambdaw*(1-flu[t-1]) +
        row_data$n_base[row_data$t==t-1]*flu_lambdaw*(flu[t-1])
      
      if (row_data$nw[row_data$t==t] <= 10){
        row_data$nw[row_data$t==t] <- 0
      }
      if (row_data$nc[row_data$t==t] <= 10){
        row_data$nc[row_data$t==t] <- 0
      }
      if (row_data$n_base[row_data$t==t] <= 10){
        row_data$n_base[row_data$t==t] <- 0
      }
    }
    all_data <- rbind(all_data, row_data)
  }
  
  
  
  avg_trendw <- all_data |>
    summarise(n = mean(nw),
              .by = t)
  avg_trendc <- all_data |>
    summarise(n = mean(nc),
              .by = t)
  avg_trendb <- all_data |>
    summarise(n = mean(n_base),
              .by = t)
  
  predictPlot <- all_data |>
    pivot_longer(cols = c(nw, nc),
                 names_to = "groupv",
                 values_to = "n") |>
    ggplot() +
    geom_line(data = all_data, aes(x = t, y = nw, group = id), 
              colour = "sienna1", alpha = 0.05, linewidth = 1) +
    geom_line(data = all_data, aes(x = t, y = nc, group = id), 
              colour = "aquamarine", alpha = 0.05, linewidth = 1) +
    geom_line(data = avg_trendw, aes(x = t, y = n), 
              colour = "sienna4", linewidth = 2, alpha = 1) +
    geom_line(data = avg_trendc, aes(x = t, y = n),
              colour = "aquamarine4", linewidth = 2, alpha = 1) +
    scale_colour_manual(limits = c("nw", "nc"),
                        labels = c("wild", "captive"),
                        values = c("sienna2", "aquamarine4")) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Timestep", y = "Population", 
         title = "Wild (orange) and captive (green) populations") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
  
  base_plot <- ggplot() +
    geom_line(data = all_data,
              aes(x = t, y = n_base, group = id), colour = "gray40", alpha = 0.05, linewidth = 1) +
    geom_line(data = avg_trendb, aes(x = t, y = n), 
              colour = "gray20", linewidth = 2, alpha = 1) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Timestep", y = "Population", 
         title = "Baseline projection") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
  
  survival <- all_data |>
    filter(t == 10) |>
    mutate(surv = nw >= 10) |>
    pull(surv) |>
    sum()/nsimreps
  
  base_survival <- all_data |>
    filter(t == 10) |>
    mutate(surv = n_base >= 10) |>
    pull(surv) |>
    sum()/nsimreps
  
  
  outList <- list(plot = outPlot,
                  predictPlot = predictPlot,
                  base_plot = base_plot,
                  survival = survival,
                  base_survival = base_survival)
  
  return(outList)
}