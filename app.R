# Thomas M. Massie, 03.10.2017, Zurich
rm(list = ls())



# Load libraries.
library(shiny)
library(shinythemes)
library(tidyverse)
library(ggplot2)
library(ggpubr)   # devtools::install_github("kassambara/ggpubr")
library(tweenr)




# ------------------------------
# Colors
palette.1 <- rev(c("#FF3A21", "#FF2C8C", "#990CE8", "#0004FF", "#0C9BE8", "#72C258"))
palette.2 <- c("#3D010A", "#4B000C", "#61000F", "#930017", "#C7001F", "#C70000")
palette.3 <- c("#032EB0", "#159C00", "#899400", "#F2B800", "#F26700", "#C70000")
# palette.4 <- c("#", "#", "#", "#", "#", "#")
col.pal <- palette.3




# shinyApp(
ui <- fluidPage(theme = shinytheme("simplex"),
                titlePanel("Time-descrete logistic population growth"),
                sidebarLayout(
                  sidebarPanel(
                    # textInput("txt", "Text input:", "text here"),
                    sliderInput("lambda", "Growth rate:", 
                                min = 0, 
                                max = 3, 
                                value = 1,
                                step = 0.01),
                    sliderInput("K", "Carrying capacity:", 
                                min = 0, 
                                max = 100, 
                                value = 50,
                                step = 1),
                    sliderInput("N_0", "Initial population density:", 
                                min = 0, 
                                max = 100, 
                                value = 1,
                                step = 0.1),
                    sliderInput("t_stop", "Maximum simulation time:", 
                                min = 1, 
                                max = 500, 
                                value = 50),
                    sliderInput("z.val", "Stochasticity level:",
                                min = 0,
                                max = 3,
                                value = 0.005,
                                step = 0.005)
                  ),
                  mainPanel(
                    plotOutput(
                      outputId = "figure",
                      width = "100%",
                      height = "90vh"
                      )
                  )
                )
)





server <- function(input, output) {
  
  dataInput <- reactive({
    
    # ------------------------------
    # Parameters
    
    # Variable model parameter values.
    lambda <- input$lambda
    K      <- input$K
    N_0    <- input$N_0
    t_stop <- input$t_stop
    z.val  <- input$z.val
    
    # Fixed model parameter values.
    t_start <- 0     # starting time for simulation
    t_range <- seq(t_start, t_stop, 1)
    
    # ...for demographic stochasticity.
    set.seed(1)
    zeta <- rnorm(as.numeric(t_stop), 0, as.numeric(z.val))
    
    # Parameter sets
    parameters <- c(lambda, K, N_0, t_start, t_stop)
    parameters.noise <- c(lambda, K, N_0, t_start, t_stop, zeta)
    
    
    # ------------------------------
    # Growth equations/functions.
    
    # Without noise.
    dlogistic <- function(parameters) {
      N <- c(N_0, numeric(t_stop))
      for (i in 1:t_stop) N[i + 1] <- {
        N[i] + lambda * N[i] * (1 - N[i] / K)
      }
      return(N)
    }
    
    # With demographic noise.
    dlogistic.demo.noise <- function(parameters.noise) {
      N <- c(N_0, numeric(t_stop))
      for (i in 1:t_stop) N[i + 1] <- {
        N[i] + lambda * N[i] * (1 - N[i] / K) + zeta[i]
      }
      return(N)
    }
    
    
    # ------------------------------
    # Run the simulation with function dlogistic.
    temp.Dens <- dlogistic()
    Density   <- data.frame(Time = t_range, PopDens = temp.Dens)
    
    
    # Run the simulation with function dlogistic.demo.noise.
    temp.Dens.noise <- dlogistic.demo.noise()
    Density.noise   <- data.frame(Time = t_range, PopDens = temp.Dens.noise)
    
    
    # -----------------
    # Cobweb data (-1 in length...)
    j <- seq(1, t_stop+1,1)
    
    NN <- seq(1,2*(t_stop+1),1)
    NN[j+j-1] <- Density$PopDens
    NN[j+j] <- Density$PopDens
    Cobweb <- data.frame(x = NN[1:length(NN)-1],
                         y = NN[2:length(NN)])
    
    NN.noise <- seq(1,2*(t_stop+1),1)
    NN.noise[j+j-1] <- Density.noise$PopDens
    NN.noise[j+j] <- Density.noise$PopDens
    Cobweb.noise <- data.frame(x = NN.noise[1:length(NN.noise)-1],
                               y = NN.noise[2:length(NN.noise)])
    
    
    # ---------------
    # Making plots for latter...
    ts.simple <- ggplot(data = Density,
                        aes(x = Time, y = PopDens)) +
      geom_line(alpha = 0.8, 
                color = col.pal[1],
                size = 0.9) +
      geom_point(size = 0.5, 
                 alpha = 0.8, 
                 color = col.pal[1]) +
      geom_hline(yintercept = input$K,
                 color = "#A9A9A9",
                 size = 0.2) +
      theme_classic() +
      labs(x = "Time", y = "Population density",
           title = "Simple logistic growth")
    
    
    ts.noise <- ggplot(data = Density.noise,
                       aes(x = Time, y = PopDens)) +
      geom_line(alpha = 0.8, 
                color = col.pal[6]) +
      geom_point(size = 0.3, 
                 alpha = 0.8,
                 color = col.pal[6]) +
      geom_hline(yintercept = input$K,
                 color = "#A9A9A9",
                 size = 0.2) +
      theme_classic() +
      labs(x = "Time", y = "Population density",
           title = "Logistic growth with demographic stochasticity")
    
    
    cw.simple <- ggplot(data = Cobweb,
                        aes(x = x, y = y)) +
      geom_line(alpha = 0.8, 
                color = col.pal[1]) +
      geom_point(size = 0.3, 
                 alpha = 0.8,
                 color = col.pal[1]) +
      theme_classic() +
      labs(x = "N_{t}", y = "N_{t+1}",
           title = "Cobweb of population dynamics")
    
    
    cw.noise <- ggplot(data = Cobweb.noise,
                       aes(x = x, y = y)) +
      geom_line(alpha = 0.8, 
                color = col.pal[1]) +
      geom_point(size = 0.3, 
                 alpha = 0.8,
                 color = col.pal[1]) +
      theme_classic() +
      labs(x = "N_{t}", y = "N_{t+1}",
           title = "Cobweb of population dynamics")
    
    
    # list(Density = Density)
    list(Density = Density, 
         Density.noise = Density.noise,
         parameters.noise = parameters.noise,
         Cobweb = Cobweb,
         Cobweb.noise = Cobweb.noise,
         ts.simple,
         ts.noise,
         cw.simple,
         cw.noise
    )
    
  })
  
  
  output$figure <- renderPlot({
    
    # Combine these 4 plots:
    
    # ggarrange(dataInput()[["ts.simple"]], 
    #           dataInput()[["ts.noise"]],  
    #           ggarrange(dataInput()[["cw.noise"]], 
    #                     dataInput()[["cw.simple"]],
    #                     labels = c("c", "d")),
    #           labels = c("a", "b"),
    #           nrow = 3) 
    
    ggarrange(ggplot(data = dataInput()[["Density"]],
                     aes(x = Time, y = PopDens)) +
                geom_line(alpha = 0.8, 
                          color = col.pal[1],
                          size = 0.9) +
                geom_point(size = 0.5, 
                           alpha = 0.8, 
                           color = col.pal[1]) +
                geom_hline(yintercept = input$K,
                           color = "#A9A9A9",
                           size = 0.2) +
                theme_classic() +
                labs(x = "Time", y = "Population density",
                     title = "Simple logistic growth"), 
              ggplot(data = dataInput()[["Density.noise"]],
                     aes(x = Time, y = PopDens)) +
                geom_line(alpha = 0.8, 
                          color = col.pal[6]) +
                geom_point(size = 0.3, 
                           alpha = 0.8,
                           color = col.pal[6]) +
                geom_hline(yintercept = input$K,
                           color = "#A9A9A9",
                           size = 0.2) +
                theme_classic() +
                labs(x = "Time", y = "Population density",
                     title = "Logistic growth with demographic stochasticity"),  
              ggarrange(ggplot(data = dataInput()[["Cobweb"]],
                               aes(x = x, y = y)) +
                          geom_hline(yintercept = input$K,
                                     color = "#A9A9A9",
                                     size = 0.2) +
                          geom_vline(xintercept = input$K,
                                     color = "#A9A9A9",
                                     size = 0.2) +
                          geom_abline(slope = 1,
                                      color = "#A9A9A9",
                                      size = 0.2) +
                          geom_line(alpha = 0.8, 
                                    color = col.pal[1]) +
                          geom_path(size = 0.3, 
                                     alpha = 0.8,
                                     color = col.pal[1]) +
                          theme_classic() +
                          labs(x = "N_{t}", y = "N_{t+1}",
                               title = "Cobweb, simple dynamics"), 
                        ggplot(data = dataInput()[["Cobweb.noise"]],
                               aes(x = x, y = y)) +
                          geom_hline(yintercept = input$K,
                                     color = "#A9A9A9",
                                     size = 0.2) +
                          geom_vline(xintercept = input$K,
                                     color = "#A9A9A9",
                                     size = 0.2) +
                          geom_abline(slope = 1,
                                      color = "#A9A9A9",
                                      size = 0.2) +
                          geom_path(alpha = 0.8, 
                                    color = col.pal[6]) +
                          geom_point(size = 0.3, 
                                     alpha = 0.8,
                                     color = col.pal[6]) +
                          theme_classic() +
                          labs(x = "N_{t}", y = "N_{t+1}",
                               title = "Cobweb, stochastic dynamics"),
                        labels = c("c", "d")),
              labels = c("a", "b"),
              nrow = 3) 
    
  })
  
  
  
  # output$figure1 <- renderPlot({
  #   # shiny::validate(
  #   
  #   # ---------------
  #   # Plotting growth curves
  #   ggplot(data = dataInput()[["Density"]],
  #          aes(x = Time, y = PopDens)) +
  #     geom_line(alpha = 0.8, 
  #               color = col.pal[1],
  #               size = 0.9) +
  #     geom_point(size = 0.5, 
  #                alpha = 0.8, 
  #                color = col.pal[1]) +
  #     geom_hline(yintercept = input$K,
  #                color = "#A9A9A9",
  #                size = 0.2) +
  #     # scale_y_continuous(labels = scales::percent) +
  #     # theme_bw() +
  #     theme_classic() +
  #     # facet_grid(. ~ as.factor(PopRate)) +
  #     # labs(x = "Time", y = "Population density") +
  #     labs(x = "Time", y = "Population density",
  #          title = "Simple logistic growth")
  #   # )
  # })
  
  # output$figure2 <- renderPlot({
  #   # shiny::validate(
  #   # ---------------
  #   # Plotting growth curves
  #   dd <- dataInput()[["Density.noise"]]
  #   ggplot(data = dd,
  #          aes(x = Time, y = PopDens)) +
  #     geom_line(alpha = 0.8, 
  #               color = col.pal[6]) +
  #     geom_point(size = 0.3, 
  #                alpha = 0.8,
  #                color = col.pal[6]) +
  #     geom_hline(yintercept = input$K,
  #                color = "#A9A9A9",
  #                size = 0.2) +
  #     # theme_bw() +
  #     theme_classic() +
  #     labs(x = "Time", y = "Population density",
  #          title = "Logistic growth with demographic stochasticity")
  #   # )
  # })
  
  
  # output$figure3 <- renderPlot({
  #   # shiny::validate(
  #   # ---------------
  #   # Plotting growth curves
  #   dd <- dataInput()[["Cobweb"]]
  #   ggplot(data = dd,
  #          aes(x = N_t, y = N_t1)) +
  #     geom_line(alpha = 0.8, 
  #               color = col.pal[2]) +
  #     geom_point(size = 0.3, 
  #                alpha = 0.8,
  #                color = col.pal[2]) +
  #     # theme_bw() +
  #     theme_classic() +
  #     labs(x = "N_{t}", y = "N_{t+1}",
  #          title = "Cob web of population dynamics")
  #   # )
  # })
  
}



# Run the application 
shinyApp(ui = ui, server = server)











