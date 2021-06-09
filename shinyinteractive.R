rm(list = ls(all.names = TRUE))
libDir = "Y:/Home/roeder/R-4.0.5.2/library"
appDir = "Y:/Home/roeder/daten/Forschungsplanung/EpiVisR/"
baseDir = "Z:/EWAS_MercapturicAcid/"
dataDir = "Z:/EWAS_MercapturicAcid/female/"
baseDir = "F:/roeder/EWAS_DrugMetabolites/"
dataDir = "F:/roeder/EWAS_DrugMetabolites/all/"
baseDir = "Z:/EWAS_VOC/"
dataDir = "Z:/EWAS_VOC/all/"
dataDir = "Z:/EWAS_VOC/female/"

setwd(dataDir)
source(paste0(appDir,"shiny_functions.R"))

#gc(verbose=T)

#PHENOFileName<-"y:/home/roeder/daten/KOLLEGEN/Wissenbach/UrineScreening/Urinescreening_short_KindID.csv"
#baseURL_EWASDataHub = "https://bigd.big.ac.cn/ewas/datahub/probe/" #https://bigd.big.ac.cn/ewas/datahub/probe/cg16867657
#baseURL_MRCEWASCatalog = "http://ewascatalog.org/?query=" #http://ewascatalog.org/?query=cg00029284
#installLibraries()
loadLibraries()
setwd(appDir)
configApp <- config::get()
setwd(dataDir)
rm(configData)
configData <- config::get()
config <- config::merge(configApp, configData)
loadObjects()

setwd(appDir)
setwd(dataDir)
conflicts(detail=TRUE)

# tell shiny to log all reactivity
#reactlog_enable()
deepDebugMode = FALSE
deepDebugMode = TRUE
debugMode = FALSE
debugMode = TRUE
debugMode = config$debugMode

#setwd(appDir)
source(paste0(appDir,"util.R"))
source(paste0(appDir,"ui.R"))
source(paste0(appDir,"server.R"))
shiny::shinyApp(ui,server)

########################Testing###########################


# once app has closed, display reactlog from shiny
shiny::reactlogShow()

options(scipen=0, digits=7)
#reactiveConsole(TRUE)
debugMode=TRUE #FALSE #TRUE
#gctorture(on = TRUE)
#dev.off()
#shiny::shinyApp(ui_short,server_short)
shiny::shinyApp(ui,server)
#shiny::shinyApp(ui,server_structured)

########################Testing###########################

devtools::install_github('yonicd/shinyHeatmaply')

#clustering
X <- matrix(runif(197*780), ncol=780)
d <- stats::dist(X)
197*780
h <- hclust(d)
plot(h)


#Clustering using k-means
# Data
x <- rbind(matrix(rnorm(70000, sd = 0.3), ncol = 2),
           matrix(rnorm(70000, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")

# CAH without kmeans : doesn't work necessarily
library(FactoMineR)
cah.test <- HCPC(x, graph=FALSE, nb.clust=-1)

# CAH with kmeans : works more quickly
cl <- kmeans(x, 1000, iter.max=20)
cah <- HCPC(cl$centers, graph=FALSE, nb.clust=-1)
plot.HCPC(cah, choice="tree")