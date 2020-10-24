#SHINY APP----
library(shiny)

options(shiny.maxRequestSize=100*1024^2)

shinyServer(function(input, output) {
  source("code.R")$value
  
  #PARAMETERS----
  #TRUE if gene annotation is symbol, FALSE if gene annotation is entrezid
  annot.symbol<- reactive({
    ifelse(input$expression.gene.annotation == "SYMBOL", TRUE, FALSE)
  })
  
  #TRUE if survival data was provided
  survdata<- reactive({
    !is.null(input$survival.data)
  })
  
  #TRUE if custom marker was provided
  custom.marker<- reactive({
    "custom_marker" %in% colnames()
  })
  
  #TRUE if a custom signature was provided
  custom.sig<- reactive({
    !is.null(input$custom.signature.data)
  })
  
  #TRUE if custom signature is symbol
  custom.is.symbol<- reactive({
    if(custom.sig()){
      ifelse(input$signature.gene.annotation == "SYMBOL", TRUE, FALSE)
    } else {
      NULL
    }
  })
  
  #DATA----
  #expression
  exp<- reactive({
    read.xlsx(input$expression.data$datapath)
  })
  
  #survival
  surv<- reactive({
    if(survdata()){
      surv<- read.xlsx(input$survival.data$datapath)
    } else {
      surv<- NULL
    }
    return(surv)
  })
  
  #custom signature
  custom.sig.genes<- reactive({
    if(custom.sig()){
      as.character(read.table(input$custom.signature.data$datapath, stringsAsFactors = FALSE)$V1)
    } else {
      NULL
    }
  })
  
  #CHECK LIBRARIES----
  check.libraries <- eventReactive(input$check.libraries, {

    showModal(modalDialog("Checking, installing and loading libraries", footer = NULL))

    check.libraries <- f_check_libraries()

    removeModal()

    return(check.libraries)
  })

  output$check.libraries<- renderText(check.libraries())

  #RUN CODE----
  results <- eventReactive(input$run, {
    
    showModal(modalDialog("Running, may take up to 30 minutes depending on the number of samples", footer=NULL))
    
    results <- f_code(a.exp = exp(),
                      a.surv = surv(),
                      a.custom.sig = custom.sig.genes(),
                      a.annot.is.symbol= annot.symbol(),
                      a.custom.sig.symbol= custom.is.symbol())
    
    removeModal()
    return(results)
    
  })
  
  #C)OUTPUTS----

  output$missing.genes<- renderText(results()$missing.genes)
  output$performance<-  renderText(results()$performance)
  output$classres<- renderText(results()$missingclass)
  
  output$classtable<- renderTable(results()$classification.table)
  
  output$signature.plots<- renderPlot(results()$signature.plots, width=800, heigh=1000)
  output$survplots<- renderPlot(results()$survplots, width=800, height = 800)
  output$deviance<- renderTable(results()$deviance)
  
  #download button
  output$classification.xlsx <- downloadHandler(
    filename = function() {
      "classification.xlsx"
    },
    content = function(file) {
      write.xlsx(results()$classification.table, file, row.names = FALSE)
    }
  )
  
  output$download.table <- renderUI({
    req(results())
    downloadButton("classification.xlsx")
  })
  
})
