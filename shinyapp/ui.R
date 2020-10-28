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
  
  titlePanel(
    h1("C1C2 classification explorer", align = "center")
  ),
  
  hr(),

  
  fluidRow(
    
    column(6,
           
           h3("Upload an .xlsx file with subjects in columns and genes in rows"),
           fileInput(inputId= "expression.data", label= NULL),

           p("The first element of the first row must be the word 'SYMBOL' or 'ENTREZID', 
             and the first column must contain either SYMBOL or ENTREZID gene identifiers accordingly.
             The remaining elements of the first row are the sample identifiers and must be unique.
             The expression data must be properly normalized within samples and expressed in a logarithimc scale.
             An example_expression.xlsx file can be found in the data directory."),
           
           hr(),
           
           h3("Run the app"),
           p("After having uploaded a valid expression file, the app can be executed"),
           uiOutput("run.button"),
           
           hr(),
           
           h3("Download the classification and signature results"),
           p("Once the processing is done, the results can be downloaded as an .xlsx file"),
           uiOutput("download.table")
    ),
    
    column(6,
           
           hr(),
           
           verbatimTextOutput("lib.msg"),
           verbatimTextOutput("data.msg"),
           verbatimTextOutput("exp.msg"),
           verbatimTextOutput("gene.presence.msg")
    )
  ),
    
  mainPanel(width=12,
    hr(),
    h2("Results", align = "center"),
    
    
    column(8, align="center",
           tableOutput("class.props"),
           plotOutput("plot", width="800px")  
    )
  )
))
