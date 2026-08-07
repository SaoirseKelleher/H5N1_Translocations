simulate_translocation <- function(model_data, nreps){
  all_data <- data.frame()
  for (i in 1:nreps){
    n0_w <- rtruncnorm(1,
                  mean = model_data$n_0_mu, 
                  sd = model_data$n_0_sd,
                  a = 0, b = 1)
    n0_c <- model_data$n_0_c
    base_lambdaw <- rnorm(1,
                          mean = model_data$mu_base_lambdaw,
                          sd = model_data$sd_base_lambdaw,
                          a = 0, b = 1)
    flu_lambdaw <- rnorm(1,
                         mean = model_data$mu_flu_lambdaw,
                         sd = model_data$sd_flu_lambdaw,
                         a = 0, b = 1)
    lambdac <- rnorm(1,
                     mean = model_data$captive_lambda_mu,
                     sd = model_data$captive_lambda_sd,
                     a = 0, b = 1)
    phi <- rnorm(1,
                 mean = model_data$phi_mu,
                 sd = model_data$phi_sd,
                 a = 0, b = 1)
    psi <- rnorm(1,
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
