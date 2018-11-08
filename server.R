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
      
      ## Choix des gènes
      observe({
        genes <- colnames(data())
        gene_choice <- updateSelectInput(session, "Gene", choices = genes)
      })
      ## Premier fichier
      
      output$contents1 <- renderTable({
        data()
      })
      
      ## Deuxième fichier
      output$contents2 <- renderTable({
        Fichier2 <- input$fichier2
        read.csv(Fichier2$datapath, header = TRUE, sep = ",")
      })
      
      ## Boxplot non fonctionnel
      output$boxplot <- renderPlot({
       data2<- read.csv(input$fichier2$datapath, header = TRUE, sep = ",")
       data_gene <- data()
       boxplot(data2,input$Gene)
      })
  
  }
)
  

