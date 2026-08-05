# UI for H5N1 shiny app

library(shiny)
library(bslib)
library(shinythemes)

# Baseline page
baseline_page <- card(
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

translocation_page <- card(
  card_header("Simulated translocations"),
  actionButton("sim_translocations", "Run simulation"),
  card(
    plotOutput("translocate_sims")
  ),
  card(
    plotOutput("translocate_sims_base")
  ),
  value_box(
    title = "Proportion of simulations where population persists (with translocations):",
    value = textOutput("translocate_survival"),
    theme = "primary"
  ),
  value_box(
    title = "Proportion of simulations where population persists (no translocations):",
    value = textOutput("translocate_survival_base"),
    theme = "secondary"
  ),
)

optimisation_page <- card(
  card_header("Optimise"),
  actionButton("run_optimiser", "Run optimisation"),
  card(
    plotOutput("optimised_translocations"),
    plotOutput("optimised_translocations_preds"),
    plotOutput("optimised_base"),
    value_box(
      title = "Proportion of simulations where population persists (with translocations):",
      value = textOutput("optimised_survival"),
      theme = "primary"
    ),
    value_box(
      title = "Proportion of simulations where population persists (no translocations):",
      value = textOutput("optimised_survival_base"),
      theme = "secondary"
    )
  )
)


# Define UI for application that draws a histogram
ui <- page_navbar(theme = bs_theme(preset = "superhero"),
                  title = "H5N1 Simulator",
                  id = "nav",
                  sidebar = sidebar(
                    title = "Set Parameters",
                    width = 500,
                    card(
                      card_header("Initial wild population"),
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
                      )
                    ),
                    card(
                      card_header("Wild growth rates"),
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
                      )),
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
                    ),
                    conditionalPanel(
                      "input.nav === 'Baseline simulations'",
                      card(
                        card_header("Estimated initial population"),
                        plotOutput("n_w_proj", height = "100px")
                      ),
                      card(
                        card_header("Estimated wild growth rates"),
                        plotOutput("lambda_proj", height = "100px")
                      ),
                    ),
                  conditionalPanel(
                    "['Translocation simulations', 'Translocation optimisation'].includes(input.nav)",
                    card(
                      card_header("Initial captive population"),
                      sliderInput(
                        "n_0_c",
                        label = "Captive",
                        min = 0,
                        max = 10000,
                        value = 0
                      ),
                      plotOutput("n_proj2", height = "100px")),
                    card(
                      card_header("Captive growth rate"),
                      sliderInput(
                        "captive_lambda_mu",
                        label = "Mean (captive)",
                        min = 0,
                        max = 3,
                        value = 0.9, 
                        step = 0.01
                      ),
                      sliderInput(
                        "captive_lambda_sd",
                        label = "Std. Deviation (captive)",
                        min = 0,
                        max = 1,
                        value = 0.1,
                        step = 0.01
                      ),
                      plotOutput("lambda_proj2", height = "100px")
                    ),
                    card(
                      card_header("Translocation success rates"),
                      sliderInput(
                        "phi_mu",
                        label = "Mean (wild->captive)",
                        min = 0,
                        max = 1,
                        value = 0.9,
                        step = 0.01
                      ),
                      sliderInput(
                        "phi_sd",
                        label = "Std. Deviation (wild->captive)",
                        min = 0,
                        max = 1,
                        value = .1,
                        step = 0.01
                      ),
                      sliderInput(
                        "psi_mu",
                        label = "Mean (captive->wild)",
                        min = 0,
                        max = 1,
                        value = 0.8,
                        step = 0.01
                      ),
                      sliderInput(
                        "psi_sd",
                        label = "Std. Deviation (captive->wild)",
                        min = 0,
                        max = 1,
                        value = .1,
                        step = 0.01
                      ),
                      plotOutput("translocate_proj", height = "100px")
                    ),
                    card(
                      card_header("Translocations"),
                      sliderInput("wtc_1",
                                  "Prop. of wild -> captive at t=1",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_2",
                                  "Prop. of wild -> captive at t=2",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_3",
                                  "Prop. of wild -> captive at t=3",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_4",
                                  "Prop. of wild -> captive at t=4",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_5",
                                  "Prop. of wild -> captive at t=5",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_6",
                                  "Prop. of wild -> captive at t=6",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_7",
                                  "Prop. of wild -> captive at t=7",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_8",
                                  "Prop. of wild -> captive at t=8",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("wtc_9",
                                  "Prop. of wild -> captive at t=9",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_1",
                                  "Prop. of captive -> wild at t=1",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_2",
                                  "Prop. of captive -> wild at t=2",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_3",
                                  "Prop. of captive -> wild at t=3",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_4",
                                  "Prop. of captive -> wild at t=4",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_5",
                                  "Prop. of captive -> wild at t=5",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_6",
                                  "Prop. of captive -> wild at t=6",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_7",
                                  "Prop. of captive -> wild at t=7",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_8",
                                  "Prop. of captive -> wild at t=8",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01),
                      sliderInput("ctw_9",
                                  "Prop. of captive -> wild at t=9",
                                  min = 0,
                                  max = 1,
                                  value = 0,
                                  step = 0.01)
                    )
                  )
),
nav_panel(
  "Baseline simulations",
  baseline_page
), 
nav_panel(
  "Translocation simulations",
  translocation_page
),
nav_panel(
  "Translocation optimisation",
  optimisation_page
)
)
