#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  fluidRow(
    column(12,
           
           titlePanel("C1C2 classification explorer"),
           verbatimTextOutput("lib.msg"),
           verbatimTextOutput("data.msg"),
           fileInput(inputId= "expression.data", label= NULL),
           verbatimTextOutput("exp.msg"),
           verbatimTextOutput("gene.presence.msg"),
           actionButton("classify", "Run"),
           uiOutput("download.table"),
           tableOutput("class.props"),
           plotOutput("plot", width="800px"),
           textOutput("test1"),
           textOutput("test2")
    )
  ),
  
  mainPanel(
    
  )
))
