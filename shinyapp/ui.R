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
           p("Given a gene expression matrix, survival data, an user-provided transcriptomic signature and an user-provided classification,
             this program will apply the C1/C2 classification, analize it according to six predefined transcriptomic signatures and the user-provided molecular signature,
             perform a survival analysis for the C1/C2 classification and the user-provided classification, and study the prognostic potential of combining the C1/C2 classification with the user-provided classification."),
           
           hr(),
           
           h1("Load data"),
           hr(),
           
           h2("Expression data"),
           p("An .xlsx file with samples in columns and genes in rows.
             First row must contain unique sample identifiers.
             First column must contain gene labels as either SYMBOL or ENTREZID.
             ENTREZIDs are preffered, SYMBOL will be converted to ENTREZID with org.Hs.eg.db package and genes without a valid ENTREZID will be discarded.
             The expression data should be properly normalized between samples and expressed in a logarithm scale.
             Genes with duplicated ENTREZID will be averaged.You can load the example expression data with SYMBOL as gene identifier from the 'example_expression.xlsx' file found in the 'data' directory."),
           fileInput(inputId= "expression.data", label= NULL),
           
           p("Indicate whether the gene annotation provided in the first column of the expression file is SYMBOL or ENTREZID"),
           radioButtons(inputId = "expression.gene.annotation", label = NULL, choices = c("SYMBOL", "ENTREZID")),
           
           h2("Survival data (optional)"),
           p("An .xlsx file with three mandatory columns named 'identifier', 'survival_time' and 'survival_status'.
             These three columns will be used in the survival analysis, if this file is ommited, no survival analysis will be performed.
             A fourth column named 'custom_marker' can be included. This must be a variable with two categories and it will be used to test the interaction with the C1/C2 classification. You can load the survival example data with a custom classification from the 'example_survival.xlsx' file in the 'data' directory."),
           fileInput(inputId= "survival.data", label= NULL),
           
           h2("Custom signature (optional)"),
           p("A custom signature as a .txt file with one gene per line can be included. 
             The average of the z-scores of the genes in this signatures will be analyzed along with the other signatures.
             Genes can be provided either as SYMBOL or ENTREZID, but SYMBOL will be converted to ENTREZID using the org.Hs.eg.db package, and genes without a valid ENTREZID will be removed. You can load a custom signature with ENTREZID as gene identifiers from the 'example_custom_signature.txt' file in the 'data' directory."),
           fileInput(inputId= "custom.signature.data", label= NULL),
           
           p("Indicate whether the gene annotation provided the in custom signature file is SYMBOL or ENTREZID. If no custom signature was provided, choose any option."),
           radioButtons(inputId = "signature.gene.annotation", label = NULL, choices = c("SYMBOL", "ENTREZID")),
           
           hr(),
           
           h1("Load packages"),
           p("Check and load the needed libraries. If the needed libraries are not properly installed and loaded, the app will fail. You may need to manually install the libraries listed in the readme."),
           actionButton("check.libraries", "Check and load packages"),
           textOutput("check.libraries"),
           
           hr(),
           
           h1("Run program"),
           
           p("Processing time may take a long time depending on the number of samples, twenty minutes for 300 samples is expectable."),
           actionButton("run", "Run")
           
           )
           ),
  
  hr(),
  
  mainPanel(
    
    h1("Results"),
    
    h2("Classification results"),
    
    textOutput("missing.genes"),
    textOutput("performance"),
    textOutput("classres"),
    h3("Download the classification and signature results:"),
    p("Signature scores are calculated only for assigned samples"),
    uiOutput("download.table"),
    
    h2("Signature plots"),

    plotOutput("signature.plots", width = "800px", height="1000px"),
    
    h2("Survival plots"),
    
    plotOutput("survplots", width = "800", height="800"),
    
    h2("Model results"),
    tableOutput("deviance")
  )
  ))
