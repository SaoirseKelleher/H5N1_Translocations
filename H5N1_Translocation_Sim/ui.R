# UI for H5N1 shiny app

library(shiny)
library(bslib)
library(shinythemes)

# Baseline page
baseline_page <- layout_sidebar(
  sidebar = sidebar(title = "Set Parameters",
                    width = 500,
                    card(
                      card_header("Initial population"),
                      sliderInput(
                        "n_0_mu",
                        label = "Mean",
                        min = 0,
                        max = 10000,
                        value = 1000
                      ),
                      sliderInput(
                        "n_0_sd",
                        label = "Std. Deviation",
                        min = 0,
                        max = 10000,
                        value = 100
                      ),
                      plotOutput("n_w_proj", height = "100px")
                    ),
                    card(
                      card_header("Growth rate"),
                      sliderInput(
                        "base_lambda_mu",
                        label = "Mean (base)",
                        min = 0,
                        max = 3,
                        value = 1,
                        step = 0.01
                      ),
                      sliderInput(
                        "base_lambda_sd",
                        label = "Std. Deviation (base)",
                        min = 0,
                        max = 1,
                        value = .1,
                        step = 0.01
                      ),
                      sliderInput(
                        "flu_lambda_mu",
                        label = "Mean (flu)",
                        min = 0,
                        max = 3,
                        value = 0.6, 
                        step = 0.01
                      ),
                      sliderInput(
                        "flu_lambda_sd",
                        label = "Std. Deviation (flu)",
                        min = 0,
                        max = 1,
                        value = .1,
                        step = 0.01
                      ),
                      plotOutput("lambda_proj", height = "100px")
                    ),
                    card(
                      card_header("Probability of flu"),
                      numericInput("t1_flu",
                                   "t=1",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t2_flu",
                                   "t=2",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t3_flu",
                                   "t=3",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t4_flu",
                                   "t=4",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t5_flu",
                                   "t=5",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t6_flu",
                                   "t=6",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t6_flu",
                                   "t=6",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t6_flu",
                                   "t=6",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t7_flu",
                                   "t=7",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t8_flu",
                                   "t=8",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                      numericInput("t9_flu",
                                   "t=9",
                                   value = 0.3,
                                   min = 0,
                                   max = 1
                      ),
                    )
  ),
  card(
    card_header("Simulated populations"),
    actionButton("sim_button", "Run simulation"),
    card(
      plotOutput("baseline_sims")
    ),
    value_box(
      title = "Proportion of simulations where population persists:",
      value = textOutput("baseline_survival"),
      theme = "primary"
    ),
  )
)

translocation_page <- layout_sidebar(
  sidebar = sidebar(title = "Set Parameters",
                    width = 500,
                    card(
                      card_header("Initial population"),
                      sliderInput(
                        "n_0_w_mu",
                        label = " Mean (Wild)",
                        min = 0,
                        max = 10000,
                        value = 1000
                      ),
                      sliderInput(
                        "n_0_w_sd",
                        label = "Std. Deviation (Wild)",
                        min = 0,
                        max = 10000,
                        value = 100
                      ),
                      sliderInput(
                        "n_0_c",
                        label = "Captive",
                        min = 0,
                        max = 10000,
                        value = 0
                      )
                    )),
  card(
    card_header("Simulate translocations"),
    actionButton("sim_translocations", "Simulate"),
    card(
      "Placeholder"
    )
  )
  )


# Define UI for application that draws a histogram
ui <- page_navbar(theme = bs_theme(preset = "superhero"),
  title = "H5N1 Simulator",
  nav_panel(
    "Baseline simulations",
    baseline_page
  ), 
  nav_panel(
    "Translocation simulations",
    translocation_page
  )
)
