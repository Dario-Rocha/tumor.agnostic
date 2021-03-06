# C1C2 classification explorer Shiny app

## Installation

### Download the app
![](https://github.com/Dario-Rocha/tumor.agnostic/blob/main/readme_images/download.jpg?raw=true)

### Unzip the downloaded file and verify that the directory structure is like this
![](https://github.com/Dario-Rocha/tumor.agnostic/blob/main/readme_images/structure.jpg?raw=true)

### Open RStudio and verify that the following packages are properly installed
* openxlsx 
* plyr 
* ggplot2 
* ggpubr
* patchwork
* AnnotationDbi
* org.Hs.eg.db
* limma
* BiocParallel 
* pbcmc (*you may need to manually download and install this package from [PBCMC Bioconductor page](https://bioconductor.riken.jp/packages/3.3/bioc/html/pbcmc.html)*)

#### If some of the packages are not installed, try to paste and run the following code in R or RStudio
````
#install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#install required packages
aux.list<- c("openxlsx",
             "ggplot2",
             "ggpubr",
             "plyr",
             "patchwork", 
             "AnnotationDbi",
             "limma",
             "org.Hs.eg.db",
             "BiocParallel")

for(aux.pack in aux.list){
  if(!require(aux.pack, character.only = TRUE)){
    
    BiocManager::install(aux.pack, dependencies = TRUE)
    
  }
}

rm(aux.pack, aux.list)
````
### Open the ui.R file in RStudio and click on "Run app" to start the application
A new window should pop up with further instructions
![](https://github.com/Dario-Rocha/tumor.agnostic/blob/main/readme_images/runapp.jpg)

## Authors

* **Elmer Andrés Fernández** - *Original Idea* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CIDIE]- [CONICET](http://www.conicet.gov.ar) - [UCC](http://www.ucc.edu.ar)
* **Dario Rocha** - *Developer and maintener" - [UNC-FCeFyN](https://fcefyn.unc.edu.ar/)
