#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(ggplot2)

# Define server logic required to draw a histogram
shinyServer(
    function(session, input, output) {
      data <- reactive({req(input$fichier1)
        read.csv(input$fichier1$datapath, header = TRUE, sep = ",", quote = "\"'")
        })
 
      ## Premier fichier
      output$contents1 <- renderTable({
        data()
      })

  }
)
  

