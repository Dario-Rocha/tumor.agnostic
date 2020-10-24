#INSTALL LIBRARIES----

#function to install libraries
f_libraries<- function(){
  #cran packages
  cran.packages<- c("openxlsx",
                    "plyr",
                    "ggplot2",
                    "ggpubr",
                    "patchwork",
                    "survival",
                    "survminer")
  
  for(aux.pack in cran.packages){
    if(!require(aux.pack, character.only = TRUE)) install.packages(aux.pack)
    library(aux.pack,character.only = TRUE)
  }
  rm(aux.pack)

  #bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  bioc.packages<- c("AnnotationDbi",
                    "org.Hs.eg.db",
                    "limma",
                    "BiocParallel")
  for(aux.pack in bioc.packages){
    if(!require(aux.pack, character.only = TRUE)){
      BiocManager::install(bioc.packages)
    } 
    library(aux.pack,character.only = TRUE)
  }
  rm(aux.pack)

  #pbcmc package
  if(!require("pbcmc", character.only = TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pbcmc")
  }
  
  message<- "All libraries seem to be present, properly loaded and installed"
  
  return(message)
}

#function to check if installation was successful
f_check_libraries<- function(){
  message<- tryCatch({
    f_libraries()
  }, error = function(e) {
    "There was an error when installing or loading libraries, please manually check if they are all properly installed"
  })
  
  return(message)
}


f_code<- function(a.exp, a.surv, a.custom.sig, a.annot.is.symbol, a.custom.sig.symbol){

  a.results<- list()
  
  #A)DATA----
  f_annot_act<- function(symbol.vec){
    entrezids1 <- mapIds(org.Hs.eg.db, keys=symbol.vec, column="ENTREZID", keytype="SYMBOL", multiVals="first")
    entrezids2<- mapIds(org.Hs.eg.db, keys=symbol.vec, column="ENTREZID", keytype="ALIAS", multiVals="first")
    entrez<- ifelse(is.na(entrezids1), entrezids2, entrezids1)
    entrez[sapply(entrez, is.null)]<- NA
    entrez<- unlist(entrez)
    return(entrez)
  }
  
  #a.load----
  a.pam50<- read.xlsx("data/pam50.xlsx")
  
  #check if surv data was provided
  a.surv.provided<- ifelse(!is.null(a.surv), TRUE, FALSE)
  
  #if surv was not loaded, create fake object
  if(is.null(a.surv)){
    a.surv<- data.frame("identifier"= colnames(a.exp)[-1],
                        stringsAsFactors = FALSE)
  }
  
  #if custom marker was provided, create indicator
  a.custom.marker<- "custom_marker" %in% colnames(a.surv)
  
  #a.annot----
  a.annot<- as.character(a.exp[,1])
  
  #if symbols were provided, convert to entrez
  if(a.annot.is.symbol){
    a.annot<- f_annot_act(a.annot)
  }
  
  #a.elist----
  a.expression<- as.matrix(a.exp[,-1])
  storage.mode(a.expression) <- "numeric"
  
  a.elist<- new("EList", list(E= a.expression, 
                              targets= a.surv, 
                              genes= data.frame("entrezid"= a.annot,
                                                "dummy"= NA,
                                                stringsAsFactors=FALSE)))
  
  #remove rows with missing annotation
  a.elist<- a.elist[!is.na(a.elist$genes$entrezid),]
  
  #a.avereps----
  a.elist<- avereps(a.elist, ID = a.elist$genes$entrezid)
  
  #pbcmc format
  a.elist$genes$EntrezGene.ID<- a.elist$genes$entrezid
  a.elist$genes<- join(a.elist$genes, a.pam50, by="EntrezGene.ID")
  
  row.names(a.elist$E)<- a.elist$genes$EntrezGene.ID
  row.names(a.elist$genes)<- a.elist$genes$EntrezGene.ID
  
  #a.signatures----
  a.signatures<- list("proliferation"= read.table("data/proliferation.txt")$V1,
                      "ca20"= read.table("data/ca20.txt")$V1,
                      "rb"= read.table("data/rb.txt")$V1,
                      "tp53"= read.table("data/tp53.txt")$V1,
                      "differentiation"= read.table("data/differentiation.txt")$V1,
                      "core_95"= read.table("data/95_core.txt")$V1)
  
  a.signatures<- lapply(a.signatures, as.character)
  
  #a.custom signature----
  if(!is.null(a.custom.sig)){
    
    if(a.custom.sig.symbol){
      a.custom_signature<- mapIds(x = org.Hs.eg.db,
                                  keys = a.custom.sig,
                                  column = "ENTREZID",
                                  keytype = "SYMBOL",
                                  multiVals = "first")
    } else {
      a.custom_signature<- a.custom.sig
    }
    
    a.custom_signature<- a.custom_signature[!is.na(a.custom_signature)]
    
    a.signatures$custom_signature<- a.custom_signature
  }
  
  #B)CLASSIFY----
  #b.check genes----
  b.missing<- which(!(a.pam50$EntrezGene.ID %in% a.elist$genes$entrezid))
  b.missing<- a.pam50$NCBI.gene.symbol[b.missing]
  
  if(isEmpty(b.missing)){
    b.missing<- "no"
  }
  
  b.missing<- paste0(b.missing, " PAM50 genes are missing in the expression matrix")
  
  a.results$missing.genes<- b.missing
  
  #b.pbcmc----
  b.pam50object<- PAM50(exprs = a.elist$E, annotation = a.elist$genes)
  b.filtratedobject<- filtrate(b.pam50object)
  b.classifiedobject<- classify(b.filtratedobject, std = "median")
  b.permutatedobject<- permutate(b.classifiedobject, nPerm=100, pCutoff=0.01, where="fdr",
                                 corCutoff=0.1, keep=FALSE, verbose=TRUE,
                                 BPPARAM=bpparam())
  
  #add pam50 subtypes
  a.elist$targets$pam50<- as.character(b.permutatedobject@permutation$subtype$PAM50)
  a.elist$targets$pam50[a.elist$targets$pam50=="Her2"]<- "Her2e"
  
  #add class
  a.elist$targets$class<- as.character(b.permutatedobject@permutation$subtype$Class)
  a.elist$targets$class[a.elist$targets$class == "LumA"]<- "C1"
  a.elist$targets$class[which(a.elist$targets$class %in% c("Basal", "Her2", "LumB"))]<- "C2"
  a.elist$targets$class[a.elist$targets$class == "Normal"]<- "unassigned"
  a.elist$targets$class[is.na(a.elist$targets$class)]<- "unassigned"
  
  #return summary
  b.performance<- c("C1"= sum(a.elist$targets$class=="C1"),
                    "C2"= sum(a.elist$targets$class=="C2"),
                    "unassigned"= sum(a.elist$targets$class=="unassigned"))
  
  b.somesamples<- b.performance["C1"]== 0 | b.performance["C2"]== 0
  b.nosamples<- b.performance["C1"]== 0 & b.performance["C2"]== 0
  
  b.performance<- round( 100*b.performance/ncol(a.elist), 2)
  
  b.performance<- paste0(b.performance, "%")
  b.performance<- paste0(b.performance[1], " C1, ", 
                         b.performance[2], " C2 and ",
                         b.performance[3], " unassigned samples")
  
  a.results$performance<- b.performance
  
  #message if one class missing or confirming that it's all fine
  a.results$missingclass<- NULL
  
  if(b.somesamples){
    a.results$missingclass<- "one class has no samples: signature scores will still be calculated for assiged samples, but no further analysis will be performed"
  }
  
  if(b.nosamples){
    a.results$missingclass<- "no samples could be confidently assigned: no further analysis will be performed"
  }
  
  #C)SIGNATURE SCORES----
  #c.functions----
  f_sigscores_ca20<- function(aux.sig, aux.cancer){
    
    #centrar y escalar la base
    #para que sea comparable entre bases
    aux.expr<- aux.cancer$E-median(aux.cancer$E, na.rm=TRUE)
    aux.expr<- aux.expr/sd(aux.expr, na.rm=TRUE)
    
    #recortar genes
    row.names(aux.expr)<- aux.cancer$genes$entrezid
    aux.genes<- intersect(row.names(aux.expr), aux.sig)
    aux.expr<- aux.expr[aux.genes,]
    
    #calcular los scores por suma
    aux.scores<- colSums(aux.expr, na.rm=TRUE)
    summary(aux.scores)
    
    #objeto salida
    aux.salida<- data.frame("identifier"= colnames(aux.cancer$E),
                            "sigscore"= aux.scores, stringsAsFactors = FALSE)
    
    return(aux.salida)
  }
  
  f_sigscores<- function(aux.sig, aux.cancer){
    
    #extraer expresión
    aux.expr<- aux.cancer$E
    
    #recortar genes
    row.names(aux.expr)<- aux.cancer$genes$entrezid
    aux.genes<- intersect(row.names(aux.expr), aux.sig)
    aux.expr<- aux.expr[aux.genes,]
    
    #estandarizar genes
    aux.expr<- t(scale(t(aux.expr)))
    
    #calcular score por promedio
    aux.scores<- colMeans(aux.expr, na.rm=TRUE)
    
    #objeto salida
    aux.salida<- data.frame("identifier"= colnames(aux.cancer$E),
                            "sigscore"= aux.scores, stringsAsFactors = FALSE)
    
    return(aux.salida)
  }
  
  #c.sigscores----
  if(!b.nosamples){
    c.ca20<- f_sigscores_ca20(a.signatures$ca20, aux.cancer = a.elist[,a.elist$targets$class != "unassigned"])
    
    c.other<- lapply(a.signatures[ !names(a.signatures)=="ca20" ], function(one.signature){
      f_sigscores(one.signature, a.elist[, a.elist$targets$class != "unassigned"])
    })
    
    c.sigscores<- c(list("ca20"=c.ca20), c.other)
    
    c.sigscores<- sapply(c.sigscores, "[", i=, j="sigscore")
    c.sigscores<- cbind.data.frame("identifier"= c.ca20$identifier, c.sigscores)
  }
  
  #c.save data----
  if(!b.nosamples){
    a.elist$targets<- join(a.elist$targets, c.sigscores, by="identifier")
  }
  
  a.results$classification.table<- a.elist$targets
  
  #D)PLOT SIGNATURES----
  d.targets<- a.elist$targets[a.elist$targets$class != "unassigned",]
  
  d.sigplots<- list()
  
  if(!b.somesamples){
    for(aux.sig.name in names(a.signatures)){
      d.sigplots[[aux.sig.name]]<- ggboxplot(d.targets, x = "class", y = aux.sig.name,
                                             color = "class", 
                                             palette = c("C1"="#005397", "C2"="#E58031"),
                                             add = "jitter") +
        stat_compare_means(method = "wilcox.test")
    }
    rm(aux.sig.name)
    
    d.sigplots<- wrap_plots(d.sigplots, ncol=3)
  }
  
  #print sig plot
  a.results$signature.plots<- d.sigplots
  
  #E)SURVIVAL----
  e.survplots<- list()
  
  #error text if ggsurv fails
  e.ggsurv.error<- ggplot() + 
    annotate(geom = 'text', 
             x=1, 
             y=1, 
             label="Something went wrong when trying to generate the survival plot.\n Maybe there are not enough assigned samples and/or events", 
             size = 4) + 
    theme_void()
  
  #e.c1/c2----
  #if survival data was provided, make km plot of c1/c2
  f_auxiliar_1<- function(survdata){
    
    aux.fit<- do.call(survfit, list(formula = Surv(survival_time, survival_status) ~ class, data = survdata))
    
    aux.ggsurv<- ggsurvplot(aux.fit,
                            pval = TRUE, risk.table = TRUE, palette = c("#005397", "#E58031"),
                            title = "C1/C2 classification")
    
    aux.plot<- wrap_plots(aux.ggsurv$plot, aux.ggsurv$table, ncol=1, heights= c(5, 1))
    
    return(aux.plot)
  }
  
  if(a.surv.provided & !b.somesamples){
    e.survplots[["C1C2"]]<- tryCatch({
      f_auxiliar_1(d.targets)
    }, error= function(e){
      e.ggsurv.error
    })
  } else {
    e.survplots[["C1C2"]]<- ggplot() + 
      annotate(geom = 'text', 
               x=1, 
               y=1, 
               label="Insufficient survival data for the C1/C2 model", 
               size = 4) + 
      theme_void()
  }
  
  #e.custom----
  #if custom marker was provided, make km plot of custom marker
  f_auxiliar_2<- function(survdata){
    
    aux.fit<- do.call(survfit, list(formula = Surv(survival_time, survival_status) ~ custom_marker, data = survdata))
    
    aux.ggsurv<- ggsurvplot(aux.fit,
                            pval = TRUE, risk.table = TRUE, palette = c("#005397", "#E58031"),
                            title = "Custom classification")
    
    aux.plot<- wrap_plots(aux.ggsurv$plot, aux.ggsurv$table, ncol=1, heights= c(5, 1))
    
    return(aux.plot)
  }
  
  if(a.surv.provided & !b.somesamples & a.custom.marker){
    e.survplots[["custom"]]<- tryCatch({
      f_auxiliar_2(d.targets)
    }, error= function(e){
      e.ggsurv.error
    })
  } else {
    e.survplots[["custom"]]<- ggplot() + 
      annotate(geom = 'text', 
               x=1, 
               y=1, 
               label="Insufficient survival data for the custom marker model", 
               size = 4) + 
      theme_void()
  }
  
  #e.compare----
  #data for comparison must have complete cases
  e.complete.index<- !is.na(d.targets$class) & 
    !is.na(d.targets$custom_marker) &
    !is.na(d.targets$survival_time) &
    !is.na(d.targets$survival_status)
  
  e.comparedata<- d.targets[e.complete.index,]
  
  #if custom marker was provided, compare models
  f_compare_models<- function(){
    e.modcustom<- coxph(Surv(survival_time, survival_status) ~ custom_marker, data=e.comparedata)
    e.modfull<- coxph(Surv(survival_time, survival_status) ~ custom_marker * class, data=e.comparedata)
    e.modsum<- summary(e.modfull)
    
    e.comparison<- anova(e.modcustom, e.modfull)
    
    e.compar.res<- c("model.p" = e.modsum$sctest[3],
                     "interaction.p"= e.modsum$coefficients[3,5],
                     "custom.loglik"= e.comparison$loglik[1],
                     "full.loglik"= e.comparison$loglik[2],
                     "deviance.p"= e.comparison$`P(>|Chi|)`[2])
    
    return(e.compar.res)
  }
  
  if(a.surv.provided & a.custom.marker & !b.somesamples){
    e.compar.res<- tryCatch({
      f_compare_models()
    }, error = function(e) {
      e.compar.res<- "The models could not be properly compared"
    })
  } else {
    e.compar.res<- "Insufficient survival data for the custom marker model"
  }
  
  #if a custom marker was provided, plot combined model
  f_auxiliar_3<- function(survdata){
    
    aux.fit<- do.call(survfit, list(formula = Surv(survival_time, survival_status) ~ class + custom_marker, data = survdata))
    
    aux.ggsurv<- ggsurvplot(aux.fit,
                            pval = TRUE, risk.table = TRUE, palette = c("#005397", "#96BEDF", "#E58031", "#FFC08F"),
                            title = "Combined model")
    
    aux.plot<- wrap_plots(aux.ggsurv$plot, aux.ggsurv$table, ncol=1, heights= c(4, 1))
    
    return(aux.plot)
  }
  
  if(a.surv.provided & !b.somesamples & a.custom.marker){
    e.survplots[["combined"]]<- tryCatch({
      f_auxiliar_3(d.targets)
    }, error= function(e){
      e.ggsurv.error
    })
  } else {
    e.survplots[["combined"]]<- ggplot() + 
      annotate(geom = 'text', 
               x=1, 
               y=1, 
               label="Insufficient survival data for the combined model", 
               size = 4) + 
      theme_void()
  }
  
  #e.save----
  a.results$survplots<- wrap_plots(e.survplots, ncol = 2)
  
  #inform comparison results
  if(a.surv.provided & !b.somesamples & a.custom.marker){
    #create table output
    e.compar.inform<- data.frame(c("Custom marker univariate log-likelihood",
                                   "Combined markers log-likelihood",
                                   "Interaction coefficient p-value",
                                   "Analysis of Deviance p-value"
    ),
    c(e.compar.res[c(3, 4, 1, 5)]))
    names(e.compar.inform)<- NULL
    
    a.results$deviance<- e.compar.inform
  } else {
    a.results$deviance<- "Insufficient survival data"
  }
  
  return(a.results)
}
