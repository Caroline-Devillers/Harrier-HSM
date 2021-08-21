# title: "CODE_MEMOIRE"
# author: "Caroline Devillers"
# date: "18/08/2021"


# Script master thesis
#Required packages 
library(rgdal)
library(raster)
library(sp)
library(biomod2)
library(gplots)
library(Factoshiny)
library(FactoMineR)
library(factoextra)
library(readr)


# Working directory

setwd("~/R/HSM")


## Variables environnementales
# Loading data
#CLC = "Corin Land Cover"
#DP  = "Population dDnsity"
#ELEV= "Elevation"
#N2K = "Natura2000"
#WAW = "Water And Wetness"

CLC <- raster("environnement/CLC.tif")
DP <- raster("environnement/DP.tif")
DP[is.na(DP[])] <- 0
ELEV <- raster("environnement/ELEV.tif")
N2K <- raster("environnement/N2K.tif")
N2K[is.na(N2K[])] <- 0 
WAW <- raster("environnement/WAW.tif")


# Check variables

CLC
projection(CLC)
bbox(CLC) ; ncol(CLC) ; nrow(CLC) ; res(CLC)

DP
projection(DP)
bbox(DP) ; ncol(DP) ; nrow(DP) ; res(DP)

ELEV
projection(ELEV)
bbox(ELEV) ; ncol(ELEV) ; nrow(ELEV) ; res(ELEV)

N2K
projection(N2K)
bbox(N2K) ; ncol(N2K) ; nrow(N2K) ; res(N2K)

WAW
projection(WAW)
bbox(WAW) ; ncol(WAW) ; nrow(WAW) ; res(WAW)


## Modifier les variables

# CLC
#Resolution (100 -> 1000) + extent (pour correspondre à DP)

CLC1 <- projectRaster(CLC,DP,method="ngb")
#check
projection(CLC1)
bbox(CLC1) ; ncol(CLC1) ; nrow(CLC1) ; res(CLC1)


# DP 
#projection = LAEA
#résolution = pixel de 1km^2
#étendue de référence 


# ELEV
#Changement de projection

ELEV <- projectRaster(ELEV,crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
#check
projection(ELEV)
bbox(ELEV) ; ncol(ELEV) ; nrow(ELEV) ; res(ELEV)

#Changement de résolution et d'étendue

ELEV2<-resample(ELEV,DP, method="bilinear")
#check
projection(ELEV2)
bbox(ELEV2) ; ncol(ELEV2) ; nrow(ELEV2) ; res(ELEV2)


# N2K
#Changement d'étendue

N2K1<-projectRaster(N2K,DP,method="ngb")
#check
projection(N2K1)
bbox(N2K1) ; ncol(N2K1) ; nrow(N2K1) ; res(N2K1)


# WAW
#Changement de résolution et d'étendue

WAW1<- projectRaster(WAW,DP,method="ngb")
#check
projection(WAW1)
bbox(WAW1) ; ncol(WAW1) ; nrow(WAW1) ; res(WAW1)


# Visualisation

plot(CLC1, main="Corin Land Cover")
plot(DP, main="Population density")
plot(ELEV2, main="Elevation")
plot(N2K1, main="Natura2000")
plot(WAW1, main="Water and wetness")


# Empiler

variables.stack <- stack(CLC1,DP,ELEV2,WAW1,N2K1)


# Données d'espèces
# Importer le dataset

circus <- read_csv("circus/test.csv", col_types = cols(CA = col_double(), CC = col_double(), CP = col_double(), DATE = col_number()))
summary(circus)


# Definir un data set par espèce
CA = *Circus aeruginosus*
  CC = *Circus cyaneus*
  CP = *Circus pygargus*
  
CA1 <- circus[c(1:533),1:2]
CC1 <- circus[c(534:1663),1:2]
CP1 <- circus[c(1664:4252),1:2]
CA <- CA1
CC <- CC1
CP <- CP1
#check
summary(CA)
summary(CC)
summary(CP)

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------##
##Pour le reste du scripte le C. aeruginosus servira d'exemple (pour obtenir les infos correspondantes des autres espèces, il suffit de remplacer "CA" par "CC" ou "CP").##
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------##

# Définir l'objects comme objet géographique

colnames(CA) <- c("CAX","CAY")
coordinates(CA)<-c("CAX","CAY")


# Definir un SCR 

projection(CA) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


# Exctraction des donnée environnementales

pts.env.CA <- extract(variables.stack, CA)
CA.env <- na.omit(data.frame(cbind(coordinates(CA),pts.env.CA)))


# enregistrer les résultats

write.table(CA.env,file="CAenv.csv",col.names = TRUE ,sep = ";")


################################################################
## Analyses descriptives de variables environnemantales mixtes ##

# Préparation

v1 <- subset(variables.stack,1)
e <- extent(v1)
p <- as(e, 'SpatialPolygons')


*Choose randomly 4000 through study zone

backgroundPoint <- coordinates(spsample(p,4000,type = "random"))


*Extract environmental values from random point : point values (PV)

PV <- extract(variables.stack,backgroundPoint)
PV <- na.omit(data.frame(cbind(coordinates(backgroundPoint),PV)))


## Analyse univariée 

PV1<- PV
PV1$CLC <- as.character(PV1$CLC)
PV1$WAW <- as.character(PV1$WAW)
PV1$N2K <- as.character(PV1$N2K)
summary(PV1)

#DP
hist(PV1$DP,main="Histogramme de la densité de population humaine", xlab="Densité de population (individus/km²)")
hist(log10(PV1$DP),main="Histogramme de la densité 
de population humaine", xlab="Densité de population (individus/km²)")

#ELEV
hist(PV1$ELEV,main="Histogramme de l'altitude", xlab="Altitude (m)")

#CLC
sort(table(PV1$CLC), decreasing = TRUE)
barplot(sort(table(PV1$CLC), decreasing = TRUE),main="Histogramme de l'occupation du sol", xlab="Catégorie d'occupation du sol")

#WAW
sort(table(PV1$WAW), decreasing = TRUE)
barplot(sort(table(PV1$WAW), decreasing = TRUE),main="Histogramme des niveaux d'humidité", xlab="Niveaux d'humidité")

#N2K
table(PV1$N2K)


## Analyse bivariée

library("dplyr")
library(tidyverse)
require(tidyverse)
require(rcompanion)


# Trié les variables

data.cor <- PV
data.cor$CLC <- as.factor(data.cor$CLC)
data.cor$WAW <- as.factor(data.cor$WAW)
data.cor$N2K <- as.factor(data.cor$N2K)
data.cor$DP  <- as.numeric(data.cor$DP)
data.cor$ELEV  <- as.numeric(data.cor$ELEV)


# Calculate a pairwise association between all variables in a data-frame
#In particular nominal vs nominal with Chi-square, 
#numeric vs numeric with Pearson correlation, 
#and nominal vs numeric with ANOVA.

mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  
                        stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      cv = cramerV(as.character(x), as.character(y), 
                   bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), 
                      complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}



# test de corrélation
correlations_test <- mixed_assoc(data.cor[,-c(1:2)])
write.table(correlations_test, file="correlation-test.csv",col.names = TRUE , 
            row.names = TRUE, sep = ";",dec=".")

###################################################################################
## Analyses descriptives des variables environnementales en fonction des espèces ##

## Aanalyses univariées

CA.env$CLC <- as.character(CA.env$CLC)
CA.env$WAW <- as.character(CA.env$WAW)
CA.env$N2K <- as.character(CA.env$N2K)
summary(CA.env)

#DP
hist(CA.env$DP,main="Histogramme de la densité de population humaine", xlab="Densité de population (individus/km²)")
hist(log(CA.env$DP))

#ELEV
hist(CA.env$ELEV,main="Histogramme de l'altitude", xlab="Altitude (m)")

#CLC
sort(table(CA.env$CLC), decreasing = TRUE)
barplot(sort(table(CA.env$CLC), decreasing = TRUE),main="Histogramme de l'occupation du sol", xlab="Catégorie d'occupation du sol")

#WAW
sort(table(CA.env$WAW), decreasing = TRUE)
barplot(sort(table(CA.env$WAW), decreasing = TRUE),main="Histogramme des niveaux d'humidité", xlab="Niveaux d'humidité")

#N2K
table(CA.env$N2K)


## Analyses bivariées

shapiro.test(PV$DP)
shapiro.test(CA.env$DP)
shapiro.test(CC.env$DP)
shapiro.test(CP.env$DP)

st_DP_PV_CA<-t.test(CA.env$DP,PV$DP,alternative = "less",
                    conf.level = 0.95,var.equal = TRUE)
st_DP_CA_CC<-t.test(CA.env$DP,CC.env$DP,alternative = "greater",
                    conf.level = 0.95,var.equal = TRUE)
st_DP_CP_CA<-t.test(CP.env$DP,CA.env$DP,alternative = "less",
                    conf.level = 0.95,var.equal = TRUE)

pvaleursDP<- c(st_DP_PV_CA$p.value, st_DP_CA_CC$p.value,st_DP_CP_CA$p.value) 
names(pvaleursDP) <- c(" PV-CA ", "CA-CC","CP-CA") 
pvaleurs_corrigéesDP <- p.adjust(pvaleursDP, method="bonferroni") 
cbind(pvaleursDP, pvaleurs_corrigéesDP) 

shapiro.test(PV$ELEV)
shapiro.test(CA.env$ELEV)
shapiro.test(CC.env$ELEV)
shapiro.test(CP.env$ELEV)

st_ELV_PV_CA<-t.test(CA.env$ELEV,PV$ELEV,alternative = "less",
                     conf.level = 0.95,var.equal = TRUE)
st_ELV_CA_CC<-t.test(CA.env$ELEV,CC.env$ELEV,alternative = "less",
                     conf.level = 0.95,var.equal = TRUE)
st_ELV_CP_CA<-t.test(CP.env$ELEV,CA.env$ELEV,alternative = "greater",
                     conf.level = 0.95,var.equal = TRUE)

pvaleursELV<- c(st_ELV_PV_CA$p.value,st_ELV_CA_CC$p.value,st_ELV_CP_CA$p.value) 
names(pvaleursELV) <- c(" PV-CA ","CA-CC","CP-CA") 
pvaleurs_corrigéesELV <- p.adjust(pvaleursELV, method="bonferroni") 
cbind(pvaleursELV, pvaleurs_corrigéesELV)

############################
## Analyses exploratoires ##

# Preparation

CA.env$CLC <- as.character(CA.env$CLC)
CA.env$WAW <- as.character(CA.env$WAW)
CA.env$N2K <- as.character(CA.env$N2K)
summary(CA.env)

# AFM
newDF.CA <- CA.env[,c("CAX","CAY","DP","ELEV","CLC","WAW","N2K")]
res.MFA.CA<-MFA(newDF.CA,group=c(2,1,1,1,1,1),
                type=c("c","s","c","n","n","n"),
                name.group=c("COOR","DP","ELEV","CLC","WAW","N2K"),
                num.group.sup=c(1),graph=FALSE)
summary(res.MFA.CA)


# Plot
plot.MFA(res.MFA.CA,choix="ind",lab.par=FALSE,invisible= c('quali','quali.sup'),
         select='cos2 0.999999',title="Graphe des individus",
         cex=1.5,cex.main=1.3,cex.axis=1.5)
plot.MFA(res.MFA.CA, choix="var",habillage='group',
         title="Cercle des corrélations",invisible= c('quanti.sup'))
plot.MFA(res.MFA.CA, choix="group",invisible = c('quanti.sup'),
         title="Graphe des groupes")
plot.MFA(res.MFA.CA, choix="ind",lab.par=FALSE,invisible= c('ind','quali.sup'),
         title="Graphe des variables qualitatives",
         cex=0.5,cex.main=1,cex.axis=1,habillage=c('WAW'))
plot.MFA(res.MFA.CA, choix="ind",lab.par=FALSE,invisible= c('ind','quali.sup'),
         title="Graphe des variables qualitatives",
         cex=0.5,cex.main=1,cex.axis=1,habillage=c('N2K'))


##########################
## Analyses predictives ##
Préparation des variables environnementales 
* Charger les raster

CLC.m <- raster("env/CLC.tif")
DP.m <- raster("env/DP.tif")
ELEV.m <- raster("env/ELEV.tif")
N2K.m <- raster("env/N2K.tif")
WAW.m <- raster("env/WAW.tif")


# Modifier les catégories selon l'annexe ??? (Quelques exemples parmis les 47 modifications apportées)

CLC.m[values(CLC.m)==2]<-1
CLC.m[values(CLC.m)==10]<-1
CLC.m[values(CLC.m)==6]<-2
CLC.m[values(CLC.m)==7]<-3
WAW.m[values(WAW.m)==4]<-1


# Assigner les même emprises et les mêmes SCR

ELEV.m1<- projectRaster(ELEV.m,crs = crs(DP.m))
ELEV.m2<- resample(ELEV.m1,DP.m, method="bilinear")
N2K.m1<-resample(N2K.m,DP.m,method="ngb")
WAW.m1<- resample(WAW.m,DP.m,method="ngb")
CLC.m1 <- resample(CLC.m,DP.m,method="ngb")


# Signaler les variables catégorielles et empiler les cartes

values(N2K.m1) = as.factor(round(values(N2K.m1))) 
values(WAW.m1) = as.factor(round(values(WAW.m1)))
values(CLC.m2) = as.factor(round(values(CLC.m2)))
variables.stack.m <- stack(CLC.m1,DP.m,ELEV.m2,WAW.m1,N2K.m1)


# Préparation des données d'espèce et de biomod2

circus <- read.csv("circus/test.csv", h=T)
DataSpecies_CA <- as.data.frame(circus[c(1:533),c(1:2,5)])
head(DataSpecies_CA) # nom de l especee 
myRespName_CA <- 'CA' # DonnÃ©es présences/ absences de l espéce (0/1) 
myResp_CA <- as.numeric(DataSpecies_CA[,myRespName_CA]) #Les coordonnÃ©es des occurences 
myRespXY_CA <- DataSpecies_CA[,c("X","Y")]
myExpl = variables.stack.m


# Formater les données

myBiomodData_CA <- BIOMOD_FormatingData(resp.var = myResp_CA, 
                                        expl.var = myExpl, # le stack des couches de variables 
                                        PA.nb.rep = 1, # pseudo-absences? (O=non, 1=oui) 
                                        PA.nb.absences = 2000,# nombre de pseudo-absences 
                                        PA.strategy = 'random', 
                                        resp.xy = myRespXY_CA, # les donnéees d occurences 
                                        resp.name = myRespName_CA)


# Biomod option : par défaut

myBiomodOption <- BIOMOD_ModelingOptions() 


# Modelisation

myBiomodModelOut_CA <- BIOMOD_Modeling( myBiomodData_CA, 
                                        models = c('GLM', 'RF','MAXENT.Phillips'), 
                                        models.options = myBiomodOption, 
                                        NbRunEval=3,# le nombre de cross-validations 
                                        DataSplit=70, # pourcentage de splitting 
                                        VarImport=3, 
                                        models.eval.meth = c('TSS','ROC'), 
                                        SaveObj = TRUE, 
                                        rescal.all.models = TRUE, 
                                        do.full.models = FALSE, 
                                        modeling.id = paste(myRespName_CA,"firstmodeling",sep=""))


# Evaluation outputs

TSS_CA<-get_evaluations(myBiomodModelOut_CA)["TSS",,,,]
TSS_CA
ROC_CA<-get_evaluations(myBiomodModelOut_CA)["ROC",,,,]
ROC_CA
myModelsVarImport_CA<-get_variables_importance(myBiomodModelOut_CA)
VarImport_CA<-apply(myModelsVarImport, c(1,2), mean)

myBiomodEM_CA <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut_CA,
                                         em.by = 'all',
                                         eval.metric = 'TSS', #avec CP on remplace TSS par ROC
                                         eval.metric.quality.threshold = 0.7,
                                         prob.mean = F,
                                         prob.mean.weight = T,
                                         prob.mean.weight.decay = 'proportional' )

evalAss_CA <- get_evaluations(myBiomodEM_CA)[[1]]
evalAss_CA #permet de connaître le cutoff pour binariser les résultat finaux


# Store the evaluation outputs

spName <- "Circus_aeruginosus"
write.table(TSS_CA,paste0(spName,"_TSS.txt"),sep="\t",row.names=T)
write.table(ROC_CA,paste0(spName,"_ROC.txt"),sep="\t",row.names=T)
write.table(evalAss_CA,paste0(spName,"_eval_assemble.txt"),sep="\t",row.names=T)
write.table(VarImport_CA,paste0(spName,"_varImport.txt"),sep="\t",row.names=T)


# Courbes de réponses des variables environnementales par modèle

CA_glm <- BIOMOD_LoadModels(myBiomodModelOut_CA,models='GLM')
CA_rf <- BIOMOD_LoadModels(myBiomodModelOut_CA,models='RF')
CA_max <- BIOMOD_LoadModels(myBiomodModelOut_CA,models='MAXENT.Phillips')

glm_eval_strip_CA <- biomod2::response.plot2(
  models = CA_glm,
  Data=get_formal_data(myBiomodModelOut_CA,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut_CA,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric='median',
  legend=F,
  display_title=F,
  data_species=get_formal_data(myBiomodModelOut_CA,'resp.var'))

rf_eval_strip_CA <- biomod2::response.plot2(
  models = CA_rf,
  Data=get_formal_data(myBiomodModelOut_CA,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut_CA,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric='median',
  legend=F,
  display_title=F,
  data_species=get_formal_data(myBiomodModelOut_CA,'resp.var'))

max_eval_strip_CA <- biomod2::response.plot2(
  models = CA_max,
  Data=get_formal_data(myBiomodModelOut_CA,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut_CA,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric='median',
  legend=F,
  display_title=F,
  data_species=get_formal_data(myBiomodModelOut_CA,'resp.var'))


# Projection

myBiomodProj_CA <- BIOMOD_Projection(modeling.output = myBiomodModelOut_CA, 
                                     new.env = variables.stack.m, 
                                     proj.name = 'current',
                                     binary.meth = 'TSS',
                                     build.clamping.mask = F) 


# Projection of the consensus 

PredEF_CA <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_CA,
                                        projection.output = myBiomodProj_CA,
                                        binary.meth = "TSS")

PredEF_CA<-raster(paste0("CA/proj_current/proj_current_CA","_ensemble.grd"))


# Store the rasters

writeRaster(PredEF_CA, 
            paste0("CA/proj_current/proj_current_CA","_ensemble.tif"), 
            datatype='INT2S', 
            options="COMPRESS=LZW", 
            overwrite=TRUE)

threshold  <-  398.0
RMatrix<- c(-Inf,threshold,0,threshold,Inf,1)
map_CA<- reclassify(PredEF_CA, RMatrix,right=T)

writeRaster(map_CA, 
            paste0("CA/proj_current/proj_current_CA","_ensemble_bin.tif"), 
            datatype='INT2S', 
            options="COMPRESS=LZW", 
            overwrite=TRUE)


# Autres plot

plot(myBiomodData_CA)
models_scores_graph(myBiomodModelOut_CA,by="models",metrics = c("ROC",'TSS'), 
                    xlim=c(.5,1), ylim=c(.5,1))
models_scores_graph(myBiomodModelOut_CA,by="cv_run",metrics = c("ROC",'TSS'), 
                    xlim=c(.5,1), ylim=c(.5,1))

