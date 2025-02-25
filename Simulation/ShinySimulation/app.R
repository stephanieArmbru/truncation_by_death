
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
                  value = paste(seq(0, 10, by = 0.25),
                                collapse = ",")),

      h3("Mortality Parameters"),
      numericInput("slope_t_threshold", "Threshold for marker change (Slope.t.threshold):",
                   value = -3),
      numericInput("long_t_threshold", "Threshold for marker value (Long.t.threshold):",
                   value = 15),
      textInput("lambda_t", "Logit baseline mortality probability (Lambda.t):",
                value = paste(rep(logit(0.02), 40), collapse = ",")),
      numericInput("theta_t", "Global marker effect on mortality probability (Theta.t):",
                   value = 0.02),
      textInput("varphi_t", "Local marker effect on mortality probability (Varphi.t):",
                value = paste(rep(0.05, 40),
                              collapse = ",")),
      textInput("varphi_ReLU_t", "Local threshold marker effect on mortality probability (Varphi.ReLU.t):",
                value = paste(rep(1.5, 40),
                              collapse = ",")),
      textInput("varphi_slope_t", "Local threshold marker slope effect on mortality probability (Varphi.slope.t):",
                value = paste(rep(2, 40),
                              collapse = ",")),
      textInput("xi_t", "Time-independent covariate effect on mortality probability (Xi.t):",
                value = paste(c(1, -2, 3), collapse = ",")),

      h3("Longitudinal Marker Parameters"),
      numericInput("long_threshold", "Threshold for marker change (Long.threshold):",
                   value = 24),
      textInput("zeta_long", "Baseline marker trend (Zeta.long):",

                # pweibull(seq(0, 2.5, length.out = 41),
                #          shape = 2.5, scale = 1) * 2 + 23

                # value = paste(sort(25 - 2 / (1 + exp(-2.5 * (seq(0, 2.5, length.out = 41) - 2))),
                #                    decreasing = T),
                #               collapse = ",")),

                value = paste(25 - 0.015 * time^2,
                              collapse = ",")),



      textInput("zeta_ReLU_long", "Threshold trend effect (Zeta.ReLU.long):",
                value = paste(rep(-1.5, 41),
                              collapse = ",")),

      textInput("eta_long", "Local autoregressive component (Eta.long):",
                value = paste(rep(0.9, 40),
                  # runif(n = 40, min = 0.85, max = 0.95),
                  collapse = ",")),
      numericInput("slope_threshold", "Threshold for slope effect (Slope.threshold):",
                   value = -1),
      textInput("eta_slope_long", "Local autoregressive slope effect (Eta.slope.long):",
                value = paste(rep(-1.5, 40),
                              collapse = ",")),

      textInput("beta_long", "Time-independent covariate effect on marker (Beta.long):",
                value = paste(c(0.5, 0.2, -0.1), collapse = ",")),
      textInput("sd_long_trajectory", "Standard deviation for longitudinal marker (Sd.long.trajectory):",
                value = 0.25),

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
        actionButton("resample", "Resample!"),

        h3("Heterogeneity in marker trajectory between example patients:"),
        plotOutput("heterogeneity")
      ),

      # Tab for patient-specific mortality probability
      tabPanel(
        title = "Patient-specific mortality probability",
        actionButton("resample", "Resample!"),
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
      ),

      # Tab for calibration metrics (sanity check of my simulation)
      tabPanel(
        title = "Calibration",
        h3("Calibration plot"),
        plotOutput("calibration"),

        h3("D-calibration"),
        tableOutput("calibration_II")
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

    slope_t_threshold <- as.numeric(input$slope_t_threshold)
    long_t_threshold <- as.numeric(input$long_t_threshold)
    lambda_t <- as.numeric(unlist(strsplit(input$lambda_t, ",")))
    theta_t <- as.numeric(input$theta_t)
    varphi_t <- as.numeric(unlist(strsplit(input$varphi_t, ",")))
    varphi_ReLU_t <- as.numeric(unlist(strsplit(input$varphi_ReLU_t, ",")))
    varphi_slope_t <- as.numeric(unlist(strsplit(input$varphi_slope_t, ",")))
    xi_t <- as.numeric(unlist(strsplit(input$xi_t, ",")))


    long_threshold <- as.numeric(input$long_threshold)
    zeta_long <- as.numeric(unlist(strsplit(input$zeta_long, ",")))
    zeta_ReLU_long <- as.numeric(unlist(strsplit(input$zeta_ReLU_long, ",")))
    eta_long <- as.numeric(unlist(strsplit(input$eta_long, ",")))
    slope_threshold <- as.numeric(input$slope_threshold)
    eta_slope_long <- as.numeric(unlist(strsplit(input$eta_slope_long, ",")))
    beta_long <- as.numeric(unlist(strsplit(input$beta_long, ",")))
    sd_long_trajectory <- as.numeric(input$sd_long_trajectory)


    # Simulate the data
    sim_dat <- SimulateJointLongData(
      n.sample = input$n,
      times = times,
      p = length(beta_long),

      slope.t.threshold = slope_t_threshold,
      long.t.threshold = long_t_threshold,
      lambda.t = lambda_t,
      theta.t = theta_t,
      varphi.t = varphi_t,
      varphi.ReLU.t = varphi_ReLU_t,
      varphi.slope.t = varphi_slope_t,
      xi.t = xi_t,

      long.threshold = long_threshold,
      zeta.long = zeta_long,
      zeta.ReLU.long = zeta_ReLU_long,
      eta.long = eta_long,
      slope.threshold = slope_threshold,
      eta.slope.long = eta_slope_long,
      beta.long = beta_long,
      sd.long.trajectory = sd_long_trajectory
    )

    list(
      parameters = data.frame(
        Parameter = c("Sample Size (n)", "Times",

                      "Slope Threshold.t", "Long. Threshold.t",
                      "Lambda.t", "Theta.t",
                      "Varphi.t", "Varphi.ReLU.t", "Varphi.slope.t",
                      "Xi.t",

                      "Long. Threshold",
                      "Zeta.long", "Zeta.ReLU.long",
                      "Eta.long", "Slope Threshold", "Eta.slope.long",
                      "Beta.long"),
        Value = c(
          input$n,
          paste(times, collapse = ", "),

          slope_t_threshold,
          long_t_threshold,
          paste(lambda_t, collapse = ", "),
          theta_t,
          paste(varphi_t, collapse = ", "),
          paste(varphi_ReLU_t, collapse = ", "),
          paste(varphi_slope_t, collapse = ", "),
          paste(xi_t, collapse = ", "),

          long_threshold,
          paste(zeta_long, collapse = ", "),
          paste(zeta_ReLU_long, collapse = ", "),
          paste(eta_long, collapse = ", "),
          slope_threshold,
          paste(eta_slope_long, collapse = ", "),
          paste(beta_long, collapse = ", ")
        )
      ),
      sim_data = sim_dat$df.return.IM %>% ungroup()
    )
  }) %>%
    bindEvent(input$simulate)

  # reactive for resampling
  resample <- reactive({
    IDs <- sample(seq(1, input$n),
           10)

    list(IDs = IDs)
  })  %>%
    bindEvent(input$resample)

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


  output$calibration <- renderPlot({
    req(simulation_results())
    plot_calibration(simulation_results()$sim_data)
  })

  output$calibration_II <- renderTable({
    req(simulation_results())
    D_res <- estimate_d_calibration(simulation_results()$sim_data)

    D_res$contingency_table
  })

  output$immortal_mortal_marker <- renderPlot({
    req(simulation_results())
    plot_immortal_mortal_longitudinal_trajectory(simulation_results()$sim_data)
  })

  output$heterogeneity <- renderPlot({
    req(simulation_results())
    req(resample())
    plot_heterogeneity(simulation_results()$sim_data,
                       sample_IDs = resample()$IDs)
  })

  output$mortality_prob <- renderPlot({
    req(simulation_results())
    req(resample())
    plot_mortality_prob(simulation_results()$sim_data,
                        sample_IDs = resample()$IDs)
  })

  output$cum_mortality_prob <- renderPlot({
    req(simulation_results())
    req(resample())
    plot_cum_mortality_prob(simulation_results()$sim_data,
                            sample_IDs = resample()$IDs)
  })

  output$diff_marker <- renderPlot({
    req(simulation_results())
    plot_diff_moral_immortal(simulation_results()$sim_data)
  })
}


# Run ---------------------------------------------------------------------

shinyApp(ui = ui, server = server)








