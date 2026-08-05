# Server logic for H5N1 simulator

library(tidyverse)
library(shiny)
library(ggtext)
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

# Function to make predictions
simulate_baseline <- function(model_data, nreps){
  all_data <- data.frame()
  for (i in 1:nreps){
    n0 <- rnorm(1,
                mean = model_data$n_0_mu, 
                sd = model_data$n_0_sd)
    base_lambdaw <- rnorm(1,
                          mean = model_data$mu_base_lambdaw,
                          sd = model_data$sd_base_lambdaw)
    flu_lambdaw <- rnorm(1,
                         mean = model_data$mu_flu_lambdaw,
                         sd = model_data$sd_flu_lambdaw)
    flu<-vector()
    for (t in 1:9){
      flu[t] <- as.numeric(runif(1,0,1) < model_data[["pr_flu"]][t])
    }
    row_data <- data.frame(i = i,
                           t = 1:10, 
                           n = NA)
    row_data[row_data$t == 1, "n"] <- n0
    for (t in 2:10){
      row_data$n[row_data$t==t] <- row_data$n[row_data$t==t-1]*base_lambdaw*(1-flu[t-1])+
        row_data$n[row_data$t==t-1]*flu_lambdaw*(flu[t-1])
      
      if (row_data$n[row_data$t==t] <= 10){
        row_data$n[row_data$t==t] <- 0
      }
    }
    
    all_data <- rbind(all_data, row_data)
  }
  
  avg_trend <- all_data |>
    summarise(n = mean(n), .by = t)
  
  plot <- all_data |>
    ggplot(aes(x = t, y = n)) +
    geom_line(aes(group = i), alpha = 0.05, colour = "gray20", linewidth = 1) +
    geom_line(data = avg_trend, colour = "sienna2", linewidth = 2) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Timestep", y = "Population") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
  
  survival <- all_data |>
    filter(t == 10) |>
    mutate(surv = n >= 10) |>
    pull(surv) |>
    sum()/nreps
  
  outList <- list(plot = plot,
                  survival = survival)
  
  return(outList)
}

simulate_translocation <- function(model_data, nreps){
  all_data <- data.frame()
  for (i in 1:nreps){
    n0_w <- rnorm(1,
                mean = model_data$n_0_mu, 
                sd = model_data$n_0_sd)
    n0_c <- model_data$n_0_c
    base_lambdaw <- rnorm(1,
                          mean = model_data$mu_base_lambdaw,
                          sd = model_data$sd_base_lambdaw)
    flu_lambdaw <- rnorm(1,
                         mean = model_data$mu_flu_lambdaw,
                         sd = model_data$sd_flu_lambdaw)
    lambdac <- rnorm(1,
                     mean = model_data$captive_lambda_mu,
                     sd = model_data$captive_lambda_sd)
    phi <- rnorm(1,
                 mean = model_data$phi_mu,
                 sd = model_data$phi_sd)
    psi <- rnorm(1,
                 mean = model_data$psi_mu,
                 sd = model_data$psi_sd)
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
  
  plot <- all_data |>
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
         title = "<span style='color: orange;'>Wild</span> and <span style='color: green;'>captive</span> populations") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
          plot.title = element_markdown())
  
  base_plot <- ggplot() +
    geom_line(data = all_data,
              aes(x = t, y = n_base, group = id), colour = "gray40", alpha = 0.05, linewidth = 1) +
    geom_line(data = avg_trendb, aes(x = t, y = n), 
              colour = "gray20", linewidth = 2, alpha = 1) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Timestep", y = "Population", 
         title = "Baseline projection") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
          plot.title = element_markdown())
  
  survival <- all_data |>
    filter(t == 10) |>
    mutate(surv = nw >= 10) |>
    pull(surv) |>
    sum()/nreps
  
  base_survival <- all_data |>
    filter(t == 10) |>
    mutate(surv = n_base >= 10) |>
    pull(surv) |>
    sum()/nreps
  
  outList <- list(plot = plot,
                  base_plot = base_plot,
                  survival = survival,
                  base_survival = base_survival)
  
  return(outList)
}

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
    n0_w <- rnorm(1,
                  mean = model_data$n_0_mu, 
                  sd = model_data$n_0_sd)
    n0_c <- model_data$n_0_c
    base_lambdaw <- rnorm(1,
                          mean = model_data$mu_base_lambdaw,
                          sd = model_data$sd_base_lambdaw)
    flu_lambdaw <- rnorm(1,
                         mean = model_data$mu_flu_lambdaw,
                         sd = model_data$sd_flu_lambdaw)
    lambdac <- rnorm(1,
                     mean = model_data$captive_lambda_mu,
                     sd = model_data$captive_lambda_sd)
    phi <- rnorm(1,
                 mean = model_data$phi_mu,
                 sd = model_data$phi_sd)
    psi <- rnorm(1,
                 mean = model_data$psi_mu,
                 sd = model_data$psi_sd)
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
         title = "<span style='color: orange;'>Wild</span> and <span style='color: green;'>captive</span> populations") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
          plot.title = element_markdown())
  
  base_plot <- ggplot() +
    geom_line(data = all_data,
              aes(x = t, y = n_base, group = id), colour = "gray40", alpha = 0.05, linewidth = 1) +
    geom_line(data = avg_trendb, aes(x = t, y = n), 
              colour = "gray20", linewidth = 2, alpha = 1) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Timestep", y = "Population", 
         title = "Baseline projection") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
          plot.title = element_markdown())
  
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


# Define server logic
function(input, output, session) {

    output$n_w_proj <- renderPlot({

      data.frame(x = (input$n_0_w_mu-(3*input$n_0_w_sd)):(input$n_0_w_mu+(3*input$n_0_w_sd))) |>
        mutate(pr = dnorm(x, mean = input$n_0_w_mu, sd = input$n_0_w_sd)) |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr)) +
        geom_line(colour = "gray40", 
                  alpha = 0.5, linewidth = 3) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())

    })
    
    output$n_proj2 <- renderPlot({
      
      data.frame(x = (input$n_0_w_mu-(3*input$n_0_w_sd)):(input$n_0_w_mu+(3*input$n_0_w_sd))) |>
        mutate(pr = dnorm(x, mean = input$n_0_w_mu, sd = input$n_0_w_sd)) |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr)) +
        geom_line(colour = "gray40", 
                  alpha = 0.5, linewidth = 3) +
        geom_vline(aes(xintercept = input$n_0_c), colour = "aquamarine3", alpha = 0.7, linewidth = 3) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())
      
    })
    
    
    output$lambda_proj <- renderPlot({
      
      data.frame(x = c(seq((input$base_lambda_mu-(3*input$base_lambda_sd)),
                           (input$base_lambda_mu+(3*input$base_lambda_sd)),
                           by = 0.05),
                       seq((input$flu_lambda_mu-(3*input$flu_lambda_sd)),
                           (input$flu_lambda_mu+(3*input$flu_lambda_sd)),
                           by = 0.05)
                 )) |>
        mutate(basepr = dnorm(x, mean = input$base_lambda_mu, sd = input$base_lambda_sd),
               flupr = dnorm(x, mean = input$flu_lambda_mu, sd = input$flu_lambda_sd)) |>
        pivot_longer(cols = c(basepr, flupr),
                     names_to = "group", values_to = "pr") |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr, colour = group)) +
        geom_line(alpha = 0.5, linewidth = 3) +
        scale_colour_manual("",
                            limits = c("basepr", "flupr"),
                            labels = c("Base", "Flu"),
                            values = c("gray40", "sienna2")) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())
      
    })
    
    output$lambda_proj2 <- renderPlot({
      
      data.frame(x = c(seq((input$base_lambda_mu-(3*input$base_lambda_sd)),
                           (input$base_lambda_mu+(3*input$base_lambda_sd)),
                           by = 0.05),
                       seq((input$flu_lambda_mu-(3*input$flu_lambda_sd)),
                           (input$flu_lambda_mu+(3*input$flu_lambda_sd)),
                           by = 0.05),
                       seq((input$captive_lambda_mu-(3*input$captive_lambda_sd)),
                           (input$captive_lambda_mu+(3*input$captive_lambda_sd)),
                           by = 0.05)
      )) |>
        mutate(basepr = dnorm(x, mean = input$base_lambda_mu, sd = input$base_lambda_sd),
               flupr = dnorm(x, mean = input$flu_lambda_mu, sd = input$flu_lambda_sd),
               captivepr = dnorm(x, mean = input$captive_lambda_mu, sd = input$captive_lambda_sd)) |>
        pivot_longer(cols = c(basepr, flupr, captivepr),
                     names_to = "group", values_to = "pr") |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr, colour = group)) +
        geom_line(alpha = 0.5, linewidth = 3) +
        scale_colour_manual("",
                            limits = c("basepr", "flupr", "captivepr"),
                            labels = c("Base", "Flu", "Captive"),
                            values = c("gray40", "sienna2", "aquamarine3")) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())
      
    })
    
    output$translocate_proj <- renderPlot({
      
      data.frame(x = c(seq((input$phi_mu-(3*input$phi_sd)),
                           (input$phi_mu+(3*input$phi_sd)),
                           by = 0.05),
                       seq((input$psi_mu-(3*input$psi_sd)),
                           (input$psi_mu+(3*input$psi_sd)),
                           by = 0.05))) |>
        mutate(phipr = dnorm(x, mean = input$phi_mu, sd = input$phi_sd),
               psipr = dnorm(x, mean = input$psi_mu, sd = input$psi_sd)) |>
        pivot_longer(cols = c(psipr, phipr),
                     names_to = "group", values_to = "pr") |>
        filter(x >= 0 & x <= 1) |>
        ggplot(aes(x = x, y = pr, colour = group)) +
        geom_line(alpha = 0.5, linewidth = 3) +
        scale_colour_manual("",
                            limits = c("phipr", "psipr"),
                            labels = c("wild->captive", "captive->wild"),
                            values = c("sienna2", "aquamarine3")) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())
      
    })
    
    
    baseline_sims <- reactive({
      input$sim_button
      
      model_data <- list(
        nTimesteps = 10,
        n_0_mu = input$n_0_w_mu,
        n_0_sd = input$n_0_w_sd,
        mu_base_lambdaw = input$base_lambda_mu,
        sd_base_lambdaw = input$base_lambda_sd,
        mu_flu_lambdaw = input$flu_lambda_mu,
        sd_flu_lambdaw = input$flu_lambda_sd,
        pr_flu = c(input$t1_flu, input$t2_flu, input$t3_flu, input$t4_flu,
                   input$t5_flu, input$t6_flu, input$t7_flu, input$t8_flu,
                   input$t9_flu)
      )
       
      output <- simulate_baseline(model_data, nreps = 1000) 
      
      return(output)
    }) |>
      bindEvent(input$sim_button)
    
    output$baseline_sims <- renderPlot(baseline_sims()$plot)
    output$baseline_survival <- renderText(baseline_sims()$survival)
    
    translocate_sims <- reactive({
        input$sim_translocations
      
      model_data <- list(
        nTimesteps = 10,
        n_0_mu = input$n_0_w_mu,
        n_0_sd = input$n_0_w_sd,
        n_0_c = input$n_0_c,
        mu_base_lambdaw = input$base_lambda_mu,
        sd_base_lambdaw = input$base_lambda_sd,
        mu_flu_lambdaw = input$flu_lambda_mu,
        sd_flu_lambdaw = input$flu_lambda_sd,
        captive_lambda_mu = input$captive_lambda_mu,
        captive_lambda_sd = input$captive_lambda_sd,
        phi_mu = input$phi_mu,
        phi_sd = input$phi_sd,
        psi_mu = input$psi_mu,
        psi_sd = input$psi_sd,
        pr_flu = c(input$t1_flu, input$t2_flu, input$t3_flu, input$t4_flu,
                   input$t5_flu, input$t6_flu, input$t7_flu, input$t8_flu,
                   input$t9_flu),
        ctw = c(input$ctw_1, input$ctw_2, input$ctw_3, input$ctw_4,
                input$ctw_5, input$ctw_6, input$ctw_7, input$ctw_8,
                input$ctw_9),
        wtc = c(input$wtc_1, input$wtc_2, input$wtc_3, input$wtc_4,
                input$wtc_5, input$wtc_6, input$wtc_7, input$wtc_8,
                input$wtc_9)
      )
      
      output <- simulate_translocation(model_data, nreps = 1000) 
      
      return(output)
    }) |>
      bindEvent(input$sim_translocations)
    
    output$translocate_sims <- renderPlot(translocate_sims()$plot)
    output$translocate_survival <- renderText(translocate_sims()$survival)
    output$translocate_sims_base <- renderPlot(translocate_sims()$base_plot)
    output$translocate_survival_base <- renderText(translocate_sims()$base_survival)
    
    optimised <- reactive({
    
      input$run_optimiser
      
      model_data <- list(
        nTimesteps = 10,
        n_0_mu = input$n_0_w_mu,
        n_0_sd = input$n_0_w_sd,
        n_0_c = input$n_0_c,
        mu_base_lambdaw = input$base_lambda_mu,
        sd_base_lambdaw = input$base_lambda_sd,
        mu_flu_lambdaw = input$flu_lambda_mu,
        sd_flu_lambdaw = input$flu_lambda_sd,
        captive_lambda_mu = input$captive_lambda_mu,
        captive_lambda_sd = input$captive_lambda_sd,
        phi_mu = input$phi_mu,
        phi_sd = input$phi_sd,
        psi_mu = input$psi_mu,
        psi_sd = input$psi_sd,
        pr_flu = c(input$t1_flu, input$t2_flu, input$t3_flu, input$t4_flu,
                   input$t5_flu, input$t6_flu, input$t7_flu, input$t8_flu,
                   input$t9_flu)
      )
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Optimising...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a callback function to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      
      output <- optimise_params(model_data, optimise_model_compiled, noptreps = 1000, nsimreps = 10000, updateProgress) 
      
      return(output)
    }) |>
      bindEvent(input$run_optimiser)
    
      output$optimised_translocations <- renderPlot(optimised()$plot)
      output$optimised_translocations_preds <- renderPlot(optimised()$predictPlot)
      output$optimised_base <- renderPlot(optimised()$base_plot)
      output$optimised_survival <- renderText(optimised()$survival)
      output$optimised_survival_base <- renderText(optimised()$base_survival)
}
