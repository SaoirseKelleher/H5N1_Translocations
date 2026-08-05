# Server logic for H5N1 simulator

library(tidyverse)
library(shiny)
#library(cmdstanr)
library(ggtext)

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


# Define server logic
function(input, output, session) {

    output$n_w_proj <- renderPlot({

      data.frame(x = (input$n_0_mu-(3*input$n_0_sd)):(input$n_0_mu+(3*input$n_0_sd))) |>
        mutate(pr = dnorm(x, mean = input$n_0_mu, sd = input$n_0_sd)) |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr)) +
        geom_line(colour = "gray40", 
                  alpha = 0.5, linewidth = 3) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())

    })
    
    output$n_proj2 <- renderPlot({
      
      data.frame(x = (input$n_0_mu2-(3*input$n_0_sd2)):(input$n_0_mu2+(3*input$n_0_sd2))) |>
        mutate(pr = dnorm(x, mean = input$n_0_mu2, sd = input$n_0_sd2)) |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr)) +
        geom_line(colour = "gray40", 
                  alpha = 0.5, linewidth = 3) +
        geom_vline(aes(xintercept = input$n_0_c2), colour = "aquamarine3", alpha = 0.7, linewidth = 3) +
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
      
      data.frame(x = c(seq((input$base_lambda_mu2-(3*input$base_lambda_sd2)),
                           (input$base_lambda_mu2+(3*input$base_lambda_sd2)),
                           by = 0.05),
                       seq((input$flu_lambda_mu2-(3*input$flu_lambda_sd2)),
                           (input$flu_lambda_mu2+(3*input$flu_lambda_sd2)),
                           by = 0.05),
                       seq((input$captive_lambda_mu2-(3*input$captive_lambda_sd2)),
                           (input$captive_lambda_mu2+(3*input$captive_lambda_sd2)),
                           by = 0.05)
      )) |>
        mutate(basepr = dnorm(x, mean = input$base_lambda_mu2, sd = input$base_lambda_sd2),
               flupr = dnorm(x, mean = input$flu_lambda_mu2, sd = input$flu_lambda_sd2),
               captivepr = dnorm(x, mean = input$captive_lambda_mu2, sd = input$captive_lambda_sd2)) |>
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
        n_0_mu = input$n_0_mu,
        n_0_sd = input$n_0_sd,
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
        n_0_mu = input$n_0_mu2,
        n_0_sd = input$n_0_sd2,
        n_0_c = input$n_0_c2,
        mu_base_lambdaw = input$base_lambda_mu2,
        sd_base_lambdaw = input$base_lambda_sd2,
        mu_flu_lambdaw = input$flu_lambda_mu2,
        sd_flu_lambdaw = input$flu_lambda_sd2,
        captive_lambda_mu = input$captive_lambda_mu2,
        captive_lambda_sd = input$captive_lambda_sd2,
        phi_mu = input$phi_mu,
        phi_sd = input$phi_sd,
        psi_mu = input$psi_mu,
        psi_sd = input$psi_sd,
        pr_flu = c(input$t1_flu2, input$t2_flu2, input$t3_flu2, input$t4_flu2,
                   input$t5_flu2, input$t6_flu2, input$t7_flu2, input$t8_flu2,
                   input$t9_flu2),
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
}
