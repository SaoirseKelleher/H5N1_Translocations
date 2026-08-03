# Server logic for H5N1 simulator

library(tidyverse)
library(shiny)
library(cmdstanr)

# Function to make predictions
plot_baseline <- function(model_data, nreps){
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
      
      if (row_data$n[row_data$t==t] <= 0){
        row_data$n[row_data$t==t] <- 0
      }
    }
    
    all_data <- rbind(all_data, row_data)
  }
  
  all_data |>
    ggplot(aes(x = t, y = n, group = i)) +
    geom_line(alpha = 0.1)
}
get_baseline_surv <- function(model_data, nreps){
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
      
      if (row_data$n[row_data$t==t] <= 0){
        row_data$n[row_data$t==t] <- 0
      }
    }
    
    all_data <- rbind(all_data, row_data)
  }
  
  outSurv <- all_data |>
    filter(t == 10) |>
    mutate(surv = n > 0) |>
    pull(surv)
  
  sum(outSurv)/nreps
    
}

# Define server logic required to draw a histogram
function(input, output, session) {

    output$n_w_proj <- renderPlot({

      data.frame(x = (input$n_0_mu-(3*input$n_0_sd)):(input$n_0_mu+(3*input$n_0_sd))) |>
        mutate(pr = dnorm(x, mean = input$n_0_mu, sd = input$n_0_sd)) |>
        filter(x > 0) |>
        ggplot(aes(x = x, y = pr)) +
        geom_line(colour = "cadetblue4", 
                  alpha = 0.5, linewidth = 3) +
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
                            values = c("cadetblue4", "firebrick")) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.title = element_blank())
      
    })
    
    output$baseline_sims <- renderPlot({
      model_data <- list(
        nTimesteps = 10,
        n_0_mu = input$n_0_mu,
        n_0_sd = input$n_0_sd,
        mu_base_lambdaw = input$base_lambda_mu,
        sd_base_lambdaw = input$base_lambda_sd,
        mu_flu_lambdaw = input$flu_lambda_mu,
        sd_flu_lambdaw = input$flu_lambda_mu,
        pr_flu = c(input$t1_flu, input$t2_flu, input$t3_flu, input$t4_flu,
                   input$t5_flu, input$t6_flu, input$t7_flu, input$t8_flu,
                   input$t9_flu)
      )
      plot_baseline(model_data, nreps = 1000)
    })
    
    output$baseline_survival <- renderText({
      model_data <- list(
        nTimesteps = 10,
        n_0_mu = input$n_0_mu,
        n_0_sd = input$n_0_sd,
        mu_base_lambdaw = input$base_lambda_mu,
        sd_base_lambdaw = input$base_lambda_sd,
        mu_flu_lambdaw = input$flu_lambda_mu,
        sd_flu_lambdaw = input$flu_lambda_mu,
        pr_flu = c(input$t1_flu, input$t2_flu, input$t3_flu, input$t4_flu,
                   input$t5_flu, input$t6_flu, input$t7_flu, input$t8_flu,
                   input$t9_flu)
      )
      get_baseline_surv(model_data, nreps = 1000)
      
    })
}
