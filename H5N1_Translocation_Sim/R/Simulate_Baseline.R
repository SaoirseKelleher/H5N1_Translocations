simulate_baseline <- function(model_data, nreps){
  all_data <- data.frame()
  for (i in 1:nreps){
    n0 <- rnorm(1,
                mean = model_data$n_0_mu, 
                sd = model_data$n_0_sd)
    base_lambdaw <- rtruncnorm(1,
                               mean = model_data$mu_base_lambdaw,
                               sd = model_data$sd_base_lambdaw,
                               a = 0, b = 1)
    flu_lambdaw <- rtruncnorm(1,
                              mean = model_data$mu_flu_lambdaw,
                              sd = model_data$sd_flu_lambdaw,
                              a = 0, b = 1)
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