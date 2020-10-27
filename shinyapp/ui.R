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
           p("Upload an .xlsx file with subjects in columns and genes in rows</br>
             The first element of the first row must be either SYMBOL or ENTREZID, 
             and the first column must contain either SYMBOL or ENTREZID gene identifiers accordingly</br>
             The remaining elements of the first row are the sample identifiers and must be unique</br>
             The expression data must be properly normalized within samples and expressed in a logarithimc scale"),
           uiOutput("run.button"),
           verbatimTextOutput("lib.msg"),
           verbatimTextOutput("data.msg"),
           fileInput(inputId= "expression.data", label= NULL),
           verbatimTextOutput("exp.msg"),
           verbatimTextOutput("gene.presence.msg"),

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
