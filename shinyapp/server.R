#SHINY APP----
library(shiny)

options(shiny.maxRequestSize=100*1024^2)

shinyServer(function(input, output) {
  source("code.R")$value
  
  #A) PARAMETERS----
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
  
  #B)LIBRARIES----
  #define libraries to be loaded
  b.libraries<- c("openxlsx",
                  "ggplot2",
                  "ggpubr",
                  "plyr",
                  "patchwork", 
                  "AnnotationDbi",
                  "limma",
                  "org.Hs.eg.db",
                  "BiocParallel",
                  "pbcmc")
  
  #try to load libraries
  b.lib.res<- f_load_libs(b.libraries)
  
  #output libraries----
  b.lib.msg<- do.call(paste, c(b.lib.res, sep = "\n"))
  
  output$lib.msg<- renderText(b.lib.msg)
  
  #C)DATA----

  #pam50----
  c.pam50<- tryCatch({
    read.xlsx("data/pam50.xlsx")
  }, error = function(e){
    "error"
  })
  
  c.pam50.msg<- "PAM50 signature correctly loaded"
  
  if(c.pam50[[1]][1] == "error"){
    c.pam50.msg <- renderText("Could not load pam50.xlsx file, make sure it exists in the data directory")
  }
  
  #signatures----
  c.signatures<- tryCatch({
    list("Proliferation"= read.table("data/proliferation.txt")$V1,
         "CA20"= read.table("data/ca20.txt")$V1,
         "RB"= read.table("data/rb.txt")$V1,
         "TP53"= read.table("data/tp53.txt")$V1,
         "Differentiation"= read.table("data/differentiation.txt")$V1,
         "Core_95"= read.table("data/95_core.txt")$V1)
  }, error = function(e){
    "error"
  })
  
  c.sig.msg<- "Signatures correctly loaded"
  
  if(c.signatures[[1]][1] == "error"){
    c.pam50.msg <- renderText("Could not load signature files, make sure they exist in the data directory")
  }
  
  #output signatures----
  c.data.msg<- paste(c.pam50.msg, c.sig.msg, sep="\n")
  output$data.msg<- renderText(c.data.msg)
  
  #expression----
  c.exp<- reactive({
    f_check_exp(read.xlsx(input$expression.data$datapath))
  })
  
  #check signatures in expression
  c.sig.condition<- reactive(c.exp.input.done() & c.pam50[[1]][1] != "error" & c.signatures[[1]][1] != "error")
  
  c.sig.exp<- reactive({
    
    if(c.sig.condition()){
      
      f_check_signatures(exp.data = c.exp()$exp, 
                                   signatures.list = c(list("PAM50"= c.pam50$EntrezGene.ID), c.signatures))
    }
    
  })
  
  #output expression----
  #pam50 and sig load
  c.data.msg<- paste(c.pam50.msg, c.sig.msg, sep="\n")
  output$data.msg<- renderText(c.data.msg)
  
  #expression load
  c.exp.input.done<- reactive({
    !is.null(input$expression.data)
  })
  
  output$exp.msg<- reactive({
    
    if(c.exp.input.done()){
      c.exp()$msg
    }
  })
  
  #gene presence check
  output$gene.presence.msg<-  renderText(c.sig.exp()) 
  
  #D)CLASSIFY----
  d.class.table<- eventReactive(input$classify, {
    
    showModal(modalDialog("Running, may take a long time depending on the number of samples. For the example data it should take less than 3 minutes.", footer=NULL))
    
    d.class.table <- f_classify(exp.data = c.exp()$exp)
    
    removeModal()
    return(d.class.table)
    
  })
  
  # inform classification results
  d.class.sum<- reactive({
    f_class_sum(d.class.table()$class)
  })
  
  #output class----
  output$class.props<- renderTable(d.class.sum())
  
  #E)SIGNATURE SCORES----
  e.new.targets<- reactive({
    req(d.class.table())
    f_sigscore_wrap(expression.data = c.exp()$exp,
                    targets.data = d.class.table(),
                    signature.list = c.signatures)
    
  })
  
  #output download table----
  #download button
  output$classification.xlsx <- downloadHandler(
    filename = function() {
      "classification.xlsx"
    },
    content = function(file) {
      write.xlsx(e.new.targets(), file, row.names = FALSE)
    }
  )
  
  output$download.table <- renderUI({
    req(e.new.targets())
    downloadButton("classification.xlsx")
  })
  
  #F)PLOTS----
  f.plot<- reactive({
    
    req(e.new.targets())
    f_plot_wrap(sig.list = c.signatures, targets = e.new.targets())

  })
  
  #output plot----
  output$plot<- renderPlot(f.plot())


  #test output----
  output$test1<- reactive({
    req(e.new.targets(), f.plot())
    class(f.plot())
    
    })
  
})
