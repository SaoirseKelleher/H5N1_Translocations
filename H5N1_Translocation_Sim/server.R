# Server logic for H5N1 simulator
library(dplyr)
library(tidyr)
library(ggplot2)
library(shiny)
library(stringr)
library(readr)
library(truncnorm)
library(cmdstanr)

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
