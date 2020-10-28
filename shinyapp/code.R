#A)LOAD LIBRARIES----
f_load_libs<- function(libraries){
  load.msg<- list()
  for(aux.lib in libraries){
    
    #load and print message
    aux.msg<- tryCatch({
      library(aux.lib, character.only = TRUE)
      paste0(aux.lib, " was correctly loaded")
    }, error = function(e){
      paste0("there was an error when loading ", aux.lib)
    })
    
    load.msg[[aux.lib]]<- aux.msg
  }
  
  return(load.msg)
}

#B)CHECK EXPRESSION----
f_annot_act<- function(symbol.vec){
  entrezids1 <- mapIds(org.Hs.eg.db, keys=symbol.vec, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  entrezids2<- mapIds(org.Hs.eg.db, keys=symbol.vec, column="ENTREZID", keytype="ALIAS", multiVals="first")
  entrez<- ifelse(is.na(entrezids1), entrezids2, entrezids1)
  entrez[sapply(entrez, is.null)]<- NA
  entrez<- unlist(entrez)
  return(entrez)
}

#check, convert and avereps
#returns formatted matrix and action log
f_check_exp<- function(exp.xlsx){
  #save object
  aux.results<- list()
  
  #format expression
  aux.results$exp<- as.matrix(exp.xlsx[,-1])
  storage.mode(aux.results$exp) <- "numeric"
  
  #create gene annotation

  #if symbol, convert to entrezid
  if(colnames(exp.xlsx)[1] == "SYMBOL"){
    aux.annot<- tryCatch({
      f_annot_act(exp.xlsx$SYMBOL)
    }, error = function(e){
      "update.error"
    })
      
  } else if(colnames(exp.xlsx)[1] == "ENTREZID"){
    aux.annot<- exp.xlsx$ENTREZID
  }
  
  #remove missing genes
  present.index<- !is.na(aux.annot)
  
  aux.results$exp<- aux.results$exp[present.index,]
  aux.annot<- aux.annot[present.index]
  
  #avereps
  row.names(aux.results$exp)<- aux.annot
  aux.results$exp<- tryCatch({
    avereps(aux.results$exp, ID = aux.annot)
  }, error = function(e){
    "error.avereps"
  })
    
  #write action log
  #check if the column was provided as asked
  aux.results$annot.found<- "Could not find gene annotation, make sure the input file is properly formated"
  
  if(colnames(exp.xlsx)[1] == "SYMBOL"){
    aux.results$annot.found <- "SYMBOL gene annotation detected"
  }
  
  if(colnames(exp.xlsx)[1] == "ENTREZID"){
    aux.results$annot.found <- "ENTREZID gene annotation detected"
  }
  
  #inform annotation conversion
  if(colnames(exp.xlsx)[1] == "SYMBOL" & aux.annot[1] == "update.error"){
    aux.results$annot.convert<- "There was an error when trying to convert SYMBOL to ENTREZID"
  }

  if(colnames(exp.xlsx)[1] == "SYMBOL" & aux.annot[1] != "update.error"){
    aux.results$annot.convert<- "SYMBOL converted to ENTREZID using org.Hs.eg.db"
  }
  
  #inform postprocessing
  aux.results$gene.delete<- "Genes with missing ENTREZID were removed"
  
  #inform avereps
  aux.results$avereps<- "Genes with duplicated identifiers were averaged"
  
  if(aux.results$exp[1] == "error.avereps"){
    aux.results$avereps<- "There was an error when trying to average repeated genes, maybe there weren not enough valid genes"
  }
  
  #combine text text
  aux.results$msg<- do.call(paste, c(aux.results[-1], sep="\n"))

  return(aux.results)
}

#function to check the number of signature genes found in the expression
#returns a message
f_check_signatures<- function(exp.data, signatures.list){

  aux.n<- list()
  for(aux.sig in names(signatures.list)){
    present<- sum(signatures.list[[aux.sig]] %in% row.names(exp.data))
    total<- length(signatures.list[[aux.sig]])
    
    aux.n[[aux.sig]]<- paste(present,
                             "of",
                             total,
                             aux.sig,
                             "genes present in the expression data")

  }
  
  aux.msg<- do.call(paste, c(aux.n, sep="\n"))
  
  return(aux.msg)
}

#C)CLASSIFY AND ANALYZE----
#function to generate the classification
#returns targets
f_classify<- function(exp.data){
  
  #create pbcmc formatted annotation
  pbcmc.annot<- data.frame("probe"= rownames(exp.data),
                           "EntrezGene.ID"= rownames(exp.data),
                           "NCBI.gene.symbol"= NA,
                           stringsAsFactors = FALSE)
  
  row.names(pbcmc.annot)<- pbcmc.annot$probe
  
  #apply pbcmc
  pam50object<- PAM50(exprs = exp.data, annotation = pbcmc.annot)
  filtratedobject<- filtrate(pam50object)
  classifiedobject<- classify(filtratedobject, std = "median")
  permutatedobject<- permutate(classifiedobject, nPerm=100, pCutoff=0.01, where="fdr",
                                 corCutoff=0.1, keep=FALSE, verbose=TRUE,
                                 BPPARAM=bpparam())
  
  #create targets
  targets<- data.frame("identifier"= colnames(exp.data),
                       "pam50"= as.character(permutatedobject@permutation$subtype$PAM50),
                       "class"= as.character(permutatedobject@permutation$subtype$Class),
                       stringsAsFactors = FALSE)
  
  #format targets
  targets$class[targets$class == "LumA"]<- "C1"
  targets$class[which(targets$class %in% c("Basal", "Her2", "LumB"))]<- "C2"
  targets$class[targets$class == "Normal"]<- "unassigned"
  targets$class[is.na(targets$class)]<- "unassigned"
  
  targets$pam50[targets$pam50=="Her2"]<- "Her2e"

  return(targets)
}

#classification summary
f_class_sum<- function(class.vector){
  
  #counts and props
  counts<-  c("C1"= sum(class.vector == "C1"),
              "C2"= sum(class.vector == "C2"),
              "Unassigned"= sum(class.vector == "unassigned"))
  
  props<- round(100* counts/length(class.vector), 2)
  
  #make table
  text.content<- paste0(counts, " (", props, "%)")
  
  table.res<- data.frame("Class"= c("C1", "C2", "Unassigned"),
                         "Frequency"= text.content)

  
  return(table.res)
}

#D)SIGNATURES----
f_sigscores_ca20<- function(gene.vector, expression.matrix){
  
  #centrar y escalar la base
  aux.expr<- expression.matrix-median(expression.matrix, na.rm=TRUE)
  aux.expr<- aux.expr/sd(aux.expr, na.rm=TRUE)
  
  #recortar genes
  aux.genes<- as.character(intersect(row.names(aux.expr), gene.vector))
  aux.expr<- aux.expr[aux.genes,]
  
  #calcular los scores por suma
  aux.scores<- colSums(aux.expr, na.rm=TRUE)
  summary(aux.scores)
  
  #objeto salida
  aux.salida<- data.frame("identifier"= colnames(expression.matrix),
                          "sigscore"= aux.scores, stringsAsFactors = FALSE)
  
  return(aux.salida)
}

f_sigscores<- function(gene.vector, expression.matrix){
  
  #recortar genes
  aux.genes<- as.character(intersect(row.names(expression.matrix), gene.vector))
  aux.expr<- expression.matrix[aux.genes,]
  
  #estandarizar genes
  aux.expr<- t(scale(t(aux.expr)))
  
  #calcular score por promedio
  aux.scores<- colMeans(aux.expr, na.rm=TRUE)
  
  #objeto salida
  aux.salida<- data.frame("identifier"= colnames(expression.matrix),
                          "sigscore"= aux.scores, stringsAsFactors = FALSE)
  
  return(aux.salida)
}

#expression has rownames as entrezid
#expression columns match target identifiers
#targets has colnames identifier and class
#signatures is a list of character vectors with entrezid
#CA20 is one of the signatures
f_sigscore_wrap<- function(expression.data, targets.data, signature.list){
  
  #select classified only
  aux.index<- which(targets.data$class %in% c("C1", "C2"))
  aux.exp<- expression.data[, aux.index]
  
  #calculate generic sigscores
  scores.list<- lapply(signature.list, function(one.signature){
    f_sigscores(gene.vector = one.signature, expression.matrix = aux.exp)
  })
  
  #calculate CA20 separately
  scores.list$CA20<- f_sigscores_ca20(gene.vector = signature.list$CA20, expression.matrix = aux.exp)
  
  #rename columns
  for(aux.name in names(scores.list)){
    colnames(scores.list[[aux.name]])[2]<- aux.name
  }

  #combine with targets
  new.targets<- plyr::join_all(scores.list, by="identifier", type="full")
  new.targets<- join(targets.data, new.targets, by="identifier", type="full")
  
  return(new.targets)
}

#E)PLOTS----
#can only handle two classes
f_plot<- function(sig.name, targets.table){
  error.plot<- ggplot() + 
    annotate(geom = 'text', 
             x=1, 
             y=1, 
             label= paste0("Could not analyze ", sig.name, "\nmaybe there are insufficient classified samples"), 
             size = 4) + 
    theme_void()
  
  el.plot<- tryCatch({
    
    ggboxplot(targets.table, x = "class", y = sig.name,
              color = "class", 
              palette = c("C1"="#005397", "C2"="#E58031"),
              add = "jitter") + stat_compare_means(method = "wilcox.test")
    
  }, error = function(e){
    error.plot
  })
  
  return(el.plot)
}

f_plot_wrap<- function(sig.list, targets){
  
  this.targets<- targets[which(targets$class %in% c("C1", "C2")),]
  
  signames<- intersect(names(sig.list), colnames(this.targets))
  
  plotlist<- lapply(signames, function(one.sig){
    f_plot(sig.name = one.sig, targets.table = this.targets)
  })
  
  plotfin<- wrap_plots(plotlist, ncol=3, guides = "collect")
  
  return(plotfin)
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
  b.permutatedobject<- permutate(b.classifiedobject, nPerm=1000, pCutoff=0.01, where="fdr",
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
    
    aux.plot<- wrap_plots(aux.ggsurv$plot, aux.ggsurv$table, ncol=1, heights= c(3.5, 1))
    
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
    e.compar.inform<- data.frame(c("Custom classification univariate log-likelihood",
                                   "Combined classifications log-likelihood",
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
