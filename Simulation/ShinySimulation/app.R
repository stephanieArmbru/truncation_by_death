
# Shiny App to explore data simulation

# Library -----------------------------------------------------------------
library(shiny)
library(tidyverse)

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Simulation: Longitudinal Marker and Mortality"),

  sidebarLayout(
    sidebarPanel(
      h1("Parameter Input"),
      tags$p("Any vectors should be comma-separated."),

      h3("Global Parameters"),
      numericInput("n", "Sample size (n):", value = 250000, min = 1),
      textInput("times", "Discrete Time Partition:",
                  value = paste(seq(0, 7, by = 0.5),
                                collapse = ",")),

      h3("Mortality Parameters"),
      textInput("lambda_t", "Logit baseline mortality probability (Lambda.t):",
                value = paste(rep(logit(0.025), 14), collapse = ",")),
      numericInput("theta_t", "Global marker effect on mortality probability (Theta.t):",
                   value = 0.05),
      textInput("varphi_t", "Local marker effect on mortality probability (Varphi.t):",
                value = paste(c(0.15, 0.15, 0.1, 0.1, 0.15, 0.2, 0.4, 1.1, 0.4, 0.2, 0.5, 0.2, 0.1, 0.1),
                              collapse = ",")),
      textInput("xi_t", "Time-independent covariate effect on mortality probability (Xi.t):",
                value = paste(c(0, 0, 0), collapse = ",")),

      h3("Longitudinal Marker Parameters"),
      textInput("zeta_long", "Baseline marker trend (Zeta.long):",
                value = paste(sort(c(runif(n = 15, min = 2, max = 2.5)), decreasing = F), collapse = ",")),
      textInput("eta_long", "Local autoregressive component (Eta.long):",
                value = paste(sort(runif(n = 14, min = 0.7, max = 1.8), decreasing = T),
                              collapse = ",")),
      textInput("beta_long", "Time-independent covariate effect on marker (Beta.long):",
                value = paste(c(0, 0, 0), collapse = ",")),
      textInput("sd_long_trajectory", "Standard deviation for longitudinal marker (sigma.t):",
                value = "0.25"),

      actionButton("simulate", "Let's go!")
    ),

    mainPanel(tabsetPanel(
      # Tab for Parameters Table
      tabPanel(
        title = "Parameters",
        h3("Check the parameters!"),
        tableOutput("parameters_table")
      ),

      # Tab for Patients Who Died
      tabPanel(
        title = "Patients Who Died",
        h3("Number of patients who died:"),
        tableOutput("died_table")
      ),

      # Tab for Baseline Mortality Probability
      tabPanel(
        title = "Baseline Mortality",
        h3("Baseline mortality probability"),
        plotOutput("baseline_mortality")
      ),

      # Tab for Baseline Marker Trend
      tabPanel(
        title = "Baseline Marker Trend",
        h3("Baseline marker trend"),
        plotOutput("baseline_long_trend")
      ),

      # Tab for Patients at Risk
      tabPanel(
        title = "Patients at Risk",
        h3("Number of patients who remain at risk over time:"),
        plotOutput("at_risk")
      ),

      # Tab for Heterogeneity in marker trajectory
      tabPanel(
        title = "Patient-specific marker trajectories",
        h3("Heterogeneity in marker trajectory between example patients:"),
        plotOutput("heterogeneity")
      ),

      # Tab for patient-specific mortality probability
      tabPanel(
        title = "Patient-specific mortality probability",
        h3("Mortality probability for example patients; conditional on surviving to previous time partition:"),
        plotOutput("mortality_prob"),

        h3("Cumulative mortality probability for example patients:"),
        plotOutput("cum_mortality_prob")
      ),

      # Tab for Mortal and Immortal Trajectory
      tabPanel(
        title = "Mortal vs Immortal Trajectory",
        h3("Longitudinal trajectory over time for mortal and immortal cohort"),
        plotOutput("immortal_mortal_marker"),

        h3("Difference in Mortal and Immortal Trajectory"),
        plotOutput("diff_marker")
      )


    )
    )
  )
)



# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive expression triggered by the action button
  simulation_results <- reactive({
    # Parse the inputs into R objects
    times <- as.numeric(unlist(strsplit(input$times, ",")))
    lambda_t <- as.numeric(unlist(strsplit(input$lambda_t, ",")))
    theta_t <- as.numeric(input$theta_t)
    varphi_t <- as.numeric(unlist(strsplit(input$varphi_t, ",")))
    xi_t <- as.numeric(unlist(strsplit(input$xi_t, ",")))
    zeta_long <- as.numeric(unlist(strsplit(input$zeta_long, ",")))
    eta_long <- as.numeric(unlist(strsplit(input$eta_long, ",")))
    beta_long <- as.numeric(unlist(strsplit(input$beta_long, ",")))
    sd_long_trajectory <- as.numeric(input$sd_long_trajectory)


    # Simulate the data
    sim_dat <- SimulateJointLongData(
      n.sample = input$n,
      times = times,
      p = length(beta_long),
      lambda.t = lambda_t,
      theta.t = theta_t,
      varphi.t = varphi_t,
      xi.t = xi_t,
      zeta.long = zeta_long,
      eta.long = eta_long,
      beta.long = beta_long,
      sd.long.trajectory = sd_long_trajectory
    )

    list(
      parameters = data.frame(
        Parameter = c("Sample Size (n)", "Times", "Lambda.t", "Theta.t",
                      "Varphi.t", "Xi.t", "Zeta.long", "Eta.long", "Beta.long"),
        Value = c(
          input$n,
          paste(times, collapse = ", "),
          paste(lambda_t, collapse = ", "),
          theta_t,
          paste(varphi_t, collapse = ", "),
          paste(xi_t, collapse = ", "),
          paste(zeta_long, collapse = ", "),
          paste(eta_long, collapse = ", "),
          paste(beta_long, collapse = ", ")
        )
      ),
      sim_data = sim_dat$df.return %>% ungroup()
    )
  }) %>%
    bindEvent(input$simulate)

  # Render parameter table
  output$parameters_table <- renderTable({
    req(simulation_results())  # Ensure the simulation is triggered
    simulation_results()$parameters
  })

  # Render plots and tables triggered by the action button
  output$baseline_mortality <- renderPlot({
    req(simulation_results())  # Ensure the simulation is triggered
    plot_baseline_mortality_risk(
      times = as.numeric(unlist(strsplit(input$times, ","))),
      lambda.t = as.numeric(unlist(strsplit(input$lambda_t, ",")))
    )
  })

  output$baseline_long_trend <- renderPlot({
    req(simulation_results())
    plot_longitudinal_trend(
      times = as.numeric(unlist(strsplit(input$times, ","))),
      zeta.long = as.numeric(unlist(strsplit(input$zeta_long, ",")))
    )
  })

  output$died_table <- renderTable({
    req(simulation_results())
    simulation_results()$sim_data %>%
      group_by(ID) %>%
      reframe(dies = any(Death == 1)) %>%
      pull(dies) %>%
      table()
  })

  output$at_risk <- renderPlot({
    req(simulation_results())
    plot_patients_at_risk(simulation_results()$sim_data)
  })

  output$immortal_mortal_marker <- renderPlot({
    req(simulation_results())
    plot_immortal_mortal_longitudinal_trajectory(simulation_results()$sim_data)
  })

  output$heterogeneity <- renderPlot({
    req(simulation_results())
    plot_heterogeneity(simulation_results()$sim_data)
  })

  output$mortality_prob <- renderPlot({
    req(simulation_results())
    plot_mortality_prob(simulation_results()$sim_data)
  })

  output$cum_mortality_prob <- renderPlot({
    req(simulation_results())
    plot_cum_mortality_prob(simulation_results()$sim_data)
  })

  output$diff_marker <- renderPlot({
    req(simulation_results())
    plot_diff_moral_immortal(simulation_results()$sim_data)
  })
}


# Run ---------------------------------------------------------------------

shinyApp(ui = ui, server = server)








