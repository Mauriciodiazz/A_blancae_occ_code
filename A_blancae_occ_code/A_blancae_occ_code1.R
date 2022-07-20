#Modelo de ocupacion para Atlapetes blancae
#Mauricio Díaz

setwd("C:/Users/Mauricio Diaz/Documents/U/Tesis/Modelo de ocupacion/Datos y script")
#upload detection histories and variables
hd<- read.table('Historias de deteccion.txt', header=T, sep = "\t")
#occupancy
ocu<-read.table('ocupacion.txt', header=T, sep = "\t")
#detection
det_ACI<-read.table('det_ACI.txt', header=T, sep = "\t")
det_hora<-read.table('det_hora.txt', header=T, sep = "\t")
det_ruido<-read.table('det_ruido.txt', header=T, sep = "\t")

##detection histories
abla_hd<-hd[,2:8]

##occupancy variables (site)
ocu_num<- scale(ocu[,c(-1,-9)]) #standardize
id<- ocu[,1] #recorders ID
site.cov<- data.frame(ocu_num,Cob_ROD=as.factor(ocu$Cob_Rod), A_vsbpe_50C=scale((ocu$A_vsbpe_50)^2))#characters to factors
head(site.cov)


##detection variables (obs)
obs.cov<-list(hora=matrix(scale(as.numeric(as.matrix(det_hora))),nrow=80, ncol=7), ruido=matrix(scale(as.numeric(as.matrix(det_ruido))),nrow=80, ncol=7), ACI=matrix(scale(as.numeric(as.matrix(det_ACI))),nrow=80, ncol=7))


#correlation
library(corrplot)
vble.ocu<-site.cov[,c(2,4,5,8:28)]
cor_ocu2<- cor(vble.ocu, method = "spearman")
cor_ocu<- cor(site.cov[,1:19], method = "spearman")
corrplot.mixed(cor_ocu2, lower = "number", 
               lower.col ="black", number.cex=.7,
               upper="square", order="FPC", title=" ",
               tl.cex=.7, tl.pos="lt", diag = "u")

#Occupancy models
#cargar paqueteria
library(unmarked) #Occupancy models package
library(dplyr) #Data management
library(kableExtra) #tables management
library(ggplot2) #figure management
library(gridExtra) #figure management
library(AICcmodavg) #to AICc

abla_umf<- unmarkedFrameOccu(y=abla_hd, siteCovs = site.cov, obsCovs= obs.cov)

#detection submodels
ab_det1<- occu(~1~1, start=c(1,1), data=abla_umf) #null model
ab_det2<- occu(~hora~1, data=abla_umf)
ab_det3<- occu(~ACI~1, data=abla_umf) 
ab_det4<- occu(~ruido~1, data=abla_umf)
ab_det5<- occu(~hora+ACI+ruido~1, data=abla_umf)

#list of models
fms_detn <- fitList ("p(.)psi(.)"               =ab_det1,   
                     "p(hora)psi(.)"            =ab_det2)
(ms_detn <- modSel(fms_detn))

#AICc values
detsim1<-AICc(ab_det1, return.K = FALSE)
detsim2<-AICc(ab_det2, return.K = FALSE)

detsimDF<- data.frame(model=ms_detn@Full$model, 
                      nPars=ms_detn@Full$nPars ,
                      AIC= ms_detn@Full$AIC, 
                      delta=ms_detn@Full$delta, 
                      AICwt=ms_detn@Full$AICwt, 
                      cumltvWt=ms_detn@Full$cumltvWt, 
                      "AICc"= c(detsim2,detsim1),
                      delta_AICc=NA)
detsimDF[1,8]<-detsimDF[1,7]-detsimDF[1,7]
detsimDF[2,8]<-detsimDF[2,7]-detsimDF[1,7]

#detection sub models table
kbl(detsimDF[,c(1,7,8)], caption = "Group Rows")%>% 
  kable_paper("striped", full_width = T) %>%
  pack_rows("Submodelos de detección",1, 2, bold=T) 

#complete models
ab_ocsim1 <- occu(~hora~Rios_b50 + Rios_b500 + Sdt_b50 + Sdt_b500 + A_vsbpe_50+ slop_50 + slop_500 +CH_GEDI_50 , data=abla_umf) 
ab_ocsim2 <- occu(~hora~A_vsbpe_50+Rios_b50, data=abla_umf)
ab_ocsim3 <- occu(~hora~A_vsbpe_50+Rios_b500, data=abla_umf)
ab_ocsim4 <- occu(~hora~A_vsbpe_50+Sdt_b50, data=abla_umf)
ab_ocsim5 <- occu(~hora~A_vsbpe_50+Sdt_b500, data=abla_umf)
ab_ocsim6 <- occu(~hora~A_vsbpe_50+slop_50, data=abla_umf)
ab_ocsim7 <- occu(~hora~A_vsbpe_50+slop_500, data=abla_umf)
ab_ocsim8 <- occu(~hora~A_vsbpe_50+CH_GEDI_50, data=abla_umf)
ab_ocsim9 <- occu(~hora~A_vsbpe_50 + A_vsbpe_500 + Rios_b50 + Rios_b500 + CH_GEDI_50, data=abla_umf)
ab_ocsim10 <- occu(~hora~A_vsbpe_50 + A_vsbpe_500 + Sdt_b50 + Sdt_b500 + slop_50 + slop_500 + CH_GEDI_50, data=abla_umf)
ab_ocsim11 <- occu(~hora~A_vsbpe_50 + A_vsbpe_500 + CH_GEDI_50, data=abla_umf)
ab_ocsim12 <- occu(~hora~A_vsbpe_50 + Sdt_b50 + slop_50 + CH_GEDI_50, data=abla_umf)
ab_ocsim13 <- occu(~hora~A_vsbpe_500 + Sdt_b500 + slop_500 + CH_GEDI_50, data=abla_umf)
ab_ocsim14 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2) + A_vsbpe_500  + Rios_b50 + Rios_b500 + CH_GEDI_50, data=abla_umf)
ab_ocsim15 <-occu(~hora~A_vsbpe_50 + Rios_b50 + Rios_b500 + CH_GEDI_50 + I(CH_GEDI_50), data=abla_umf)
ab_ocsim16 <-occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2) + CH_GEDI_50 , data=abla_umf)
ab_ocsim17 <-occu(~hora~A_vsbpe_50 + A_vsbpe_500 + CH_GEDI_50 + I(CH_GEDI_50), data=abla_umf)
ab_ocsim18 <-occu(~hora~A_vsbpe_50 + CH_GEDI_50 + I(CH_GEDI_50) + Sdt_b50, data=abla_umf)
ab_ocsim19 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+Rios_b50, data=abla_umf)
ab_ocsim20 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+Rios_b500, data=abla_umf)
ab_ocsim21 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+Sdt_b50, data=abla_umf)
ab_ocsim22 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+Sdt_b500, data=abla_umf)
ab_ocsim23 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+slop_50, data=abla_umf)
ab_ocsim24 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2)+slop_500, data=abla_umf)
ab_ocsim27 <- occu(~hora~A_vsbpe_50 + A_vsbpe_500+I(A_vsbpe_500^2) + Rios_b50 + Rios_b500 + CH_GEDI_50, data=abla_umf)
ab_ocsim28 <- occu(~hora~A_vsbpe_50+I(A_vsbpe_50^2) + A_vsbpe_500 + CH_GEDI_50, data=abla_umf)
ab_ocsim29 <- occu(~hora~A_vsbpe_50+A_vsbpe_500 + I(A_vsbpe_500^2) + CH_GEDI_50, data=abla_umf)

#occupancy models table
fms_ocsim <- fitList ("p(.)psi(.)"               = ab_det1,
                      "p(hora)psi(Rios_b50 + Rios_b500 + Sdt_b50 + Sdt_b500 + A_vsbpe_50+ slop_50 + slop_500 +CH_GEDI_50)"= ab_ocsim1,
                      "p(hora)psi(A_vsbpe_50+Rios_b50)"         = ab_ocsim2,
                      "p(hora)psi(A_vsbpe_50+Rios_b500)"         = ab_ocsim3,
                      "p(hora)psi(A_vsbpe_50+Sdt_b50)"         = ab_ocsim4,
                      "p(hora)psi(A_vsbpe_50+Sdt_b500)"          = ab_ocsim5,
                      "p(hora)psi(A_vsbpe_50+slop_50)"       = ab_ocsim6,
                      "p(hora)psi(A_vsbpe_50+slop_500)"       = ab_ocsim7,
                      "p(hora)psi(A_vsbpe_50+CH_GEDI_50)"       = ab_ocsim8,
                      "p(hora)psi(A_vsbpe_50 + Rios_b50 + Rios_b500 + CH_GEDI_50)"        = ab_ocsim9,
                      "p(hora)psi(A_vsbpe_50 + Sdt_b50 + Sdt_b500 + slop_50 + slop_500 + CH_GEDI_50)"  = ab_ocsim10,
                      "p(hora)psi(A_vsbpe_50 + A_vsbpe_500 + CH_GEDI_50)"          = ab_ocsim11,
                      "p(hora)psi(A_vsbpe_50 + Sdt_b50 + slop_50 + CH_GEDI_50)"          =ab_ocsim12,
                      "p(hora)psi(A_vsbpe_50^2  + Rios_b50 + Rios_b500 + CH_GEDI_50)"        =ab_ocsim14,
                      "p(hora)psi(A_vsbpe_50 + Rios_b50 + Rios_b500 +CH_GEDI_50 ^2)"        =ab_ocsim15,
                      "p(hora)psi(A_vsbpe_50^2 + CH_GEDI_50)"        =ab_ocsim16,
                      "p(hora)psi(A_vsbpe_50 +CH_GEDI_50 ^2)"        =ab_ocsim17,
                      "p(hora)psi(A_vsbpe_50 +CH_GEDI_50 ^2 + Sdt_b50)"        =ab_ocsim18,
                      "p(hora)psi(A_vsbpe_50^2+Rios_b50)"         = ab_ocsim19,
                      "p(hora)psi(A_vsbpe_50^2+Rios_b500)"         = ab_ocsim20,
                      "p(hora)psi(A_vsbpe_50^2+Sdt_b50)"         = ab_ocsim21,
                      "p(hora)psi(A_vsbpe_50^2+Sdt_b500)"          = ab_ocsim22,
                      "p(hora)psi(A_vsbpe_50^2+slop_50)"       = ab_ocsim23,
                      "p(hora)psi(A_vsbpe_50^2+slop_500)"       = ab_ocsim24,
                      "p(hora)psi(A_vsbpe_50 + A_vsbpe_500^2 + Rios_b50 + Rios_b500 + CH_GEDI_50)"  = ab_ocsim27,
                      "p(hora)psi(A_vsbpe_50^2 + A_vsbpe_500 + CH_GEDI_50)"       = ab_ocsim28,
                      "p(hora)psi(A_vsbpe_50 + A_vsbpe_500^2 + CH_GEDI_50)"       = ab_ocsim29
)
(ms_ocsim <- modSel(fms_ocsim))

#AICc values for each model - to build the kabletable
nulle<-AICc(ab_det1, return.K = FALSE)
ocsim1<-AICc(ab_ocsim1, return.K = FALSE)
ocsim2<-AICc(ab_ocsim2, return.K = FALSE)
ocsim3<-AICc(ab_ocsim3, return.K = FALSE)
ocsim4<-AICc(ab_ocsim4, return.K = FALSE)
ocsim5<-AICc(ab_ocsim5, return.K = FALSE)
ocsim6<-AICc(ab_ocsim6, return.K = FALSE)
ocsim7<-AICc(ab_ocsim7, return.K = FALSE)
ocsim8<-AICc(ab_ocsim8, return.K = FALSE)
ocsim9<-AICc(ab_ocsim9, return.K = FALSE)
ocsim10<-AICc(ab_ocsim10, return.K = FALSE)
ocsim11<-AICc(ab_ocsim11, return.K = FALSE)
ocsim12<-AICc(ab_ocsim12, return.K = FALSE)
ocsim13<-AICc(ab_ocsim13, return.K = FALSE)
ocsim14<-AICc(ab_ocsim14, return.K = FALSE)
ocsim15<-AICc(ab_ocsim15, return.K = FALSE)
ocsim16<-AICc(ab_ocsim16, return.K = FALSE)
ocsim17<-AICc(ab_ocsim17, return.K = FALSE)
ocsim18<-AICc(ab_ocsim18, return.K = FALSE)
ocsim19<-AICc(ab_ocsim19, return.K = FALSE)
ocsim20<-AICc(ab_ocsim20, return.K = FALSE)
ocsim21<-AICc(ab_ocsim21, return.K = FALSE)
ocsim22<-AICc(ab_ocsim22, return.K = FALSE)
ocsim23<-AICc(ab_ocsim23, return.K = FALSE)
ocsim24<-AICc(ab_ocsim24, return.K = FALSE)
ocsim27<-AICc(ab_ocsim27, return.K = FALSE)
ocsim28<-AICc(ab_ocsim28, return.K = FALSE)
ocsim29<-AICc(ab_ocsim29, return.K = FALSE)

Snulle<-summary(ab_det1)
Scompl<-summary(ab_ocsim1)
Socsim2<-summary(ab_ocsim2)
Socsim3<-summary(ab_ocsim3)
Socsim4<-summary(ab_ocsim4)
Socsim5<-summary(ab_ocsim5)
Socsim6<-summary(ab_ocsim6)
Socsim7<-summary(ab_ocsim7)
Socsim8<-summary(ab_ocsim8)
Socsim9<-summary(ab_ocsim9)
Socsim10<-summary(ab_ocsim10)
Socsim11<-summary(ab_ocsim11)
Socsim12<-summary(ab_ocsim12)
Socsim13<-summary(ab_ocsim13)
Socsim14<-summary(ab_ocsim14)
Socsim15<-summary(ab_ocsim15)
Socsim16<-summary(ab_ocsim16)
Socsim17<-summary(ab_ocsim17)
Socsim18<-summary(ab_ocsim18)
Socsim19<-summary(ab_ocsim19)
Socsim20<-summary(ab_ocsim20)
Socsim21<-summary(ab_ocsim21)
Socsim22<-summary(ab_ocsim22)
Socsim23<-summary(ab_ocsim23)
Socsim24<-summary(ab_ocsim24)
Socsim27<-summary(ab_ocsim27)
Socsim28<-summary(ab_ocsim28)
Socsim29<-summary(ab_ocsim29)

ocsim1DF<- data.frame(model=ms_ocsim@Full$model, 
                    nPars=ms_ocsim@Full$nPars ,
                    AIC= ms_ocsim@Full$AIC, 
                    delta=ms_ocsim@Full$delta, 
                    AICwt=ms_ocsim@Full$AICwt, 
                    cumltvWt=ms_ocsim@Full$cumltvWt, 
                    "AICc"= c(ocsim22,                              
                              ocsim16,                              ocsim10,
                              ocsim21,                              ocsim5,
                              ocsim28,                              ocsim14,
                              ocsim12,                              ocsim1,
                              ocsim18,                              ocsim8,
                              ocsim24,                              ocsim4,
                              ocsim23,                              ocsim20,
                              ocsim19,                              ocsim11,
                              ocsim29,                              ocsim17,
                              ocsim9,                              ocsim7,
                              ocsim15,                              ocsim6,
                              ocsim2,                              ocsim27,
                              ocsim3,                              nulle),
                    delta_AICc=NA) 

#AICc delta
for (x in 1:length(ocsim1DF$model)) {
  ocsim1DF[x,8]<-ocsim1DF[x,7]-ocsim1DF[1,7]
}

ocsim1DF <- ocsim1DF[order(ocsim1DF$AICc),]

#table for best occupancy model
estimates_state<-data.frame(ab_ocsim22@estimates@estimates$state@estimates)
names_state<-row.names(estimates_state)
summ<-summary(ab_ocsim22)

estimates_det<-data.frame(ab_ocsim22@estimates@estimates$det@estimates)
names_det<-row.names(estimates_det)

tab_state<- data.frame(names=names_state, Estimate=estimates_state[,1], SE=summ$state$SE, z=summ$state$z, "P(>|z|)"=summ$state$`P(>|z|)`, check.names = F)
tab_det<- data.frame(names=names_det, Estimate= estimates_det[,1], SE=summ$det$SE, z=summ$det$z, "P(>|z|)"=summ$det$`P(>|z|)`, check.names = F)

tab_mod<-rbind(tab_state, tab_det)
kbl(tab_mod, caption = "Group Rows")%>% 
  kable_paper("striped", full_width = T) %>%
  pack_rows("Occupancy",1, 4, bold=T) %>%
  pack_rows("Detection",5, 6, bold=T)


####Variables prediction####

#for occupancy
ocu.mean<-predict(ab_ocsim22, type="state")
det.mean<-predict(ab_ocsim22, type="det")

Cov_Sdt<- seq(range(abla_umf@siteCovs$Sdt_b500)[1]-0.25,
             range(abla_umf@siteCovs$Sdt_b500)[2]+0.25,
             length.out=100)

Cov_vsbpe <- seq(range(abla_umf@siteCovs$A_vsbpe_50)[1]-0.25,
                 range(abla_umf@siteCovs$A_vsbpe_50)[2]+0.25,
                 length.out=100)

#Terrain second derivative (sdt)
newData_sdt<- data.frame(Sdt_b500=Cov_Sdt, A_vsbpe_50=0)
preddata_sdt <- predict(ab_ocsim22, newdata=newData_sdt, type="state")
preddata_sdt$Sdt_b500<- Cov_Sdt
preddata_sdt$Sdt_b500nst<- (preddata_sdt$Sdt_b500*sd(ocu$Sdt_b500))+mean(ocu$Sdt_b500) #"des-standardize"
  
#to build the plot
predplot_sdt<-ggplot(preddata_sdt,aes(x=Sdt_b500nst,y=Predicted)) +
  labs(x="Terrain second derivative (rad/m)", y=expression(paste('Pred. prob. ', psi))) +
  geom_ribbon(data = preddata_sdt,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf") +
  geom_line(data = preddata_sdt, colour="#000000", size=0.5) +
  annotate(geom="text", x=0.0002, y=1, label="A",
           color="black", size=5) +
  theme_classic()
predplot_sdt

#Shrubs and herbaceous vegetation
newData_vsbpe<- data.frame(Sdt_b500=0, A_vsbpe_50=Cov_vsbpe)

preddata_vsbpe <- predict(ab_ocsim22, newdata=newData_vsbpe, type="state")
preddata_vsbpe$A_vsbpe_50<- Cov_vsbpe
preddata_vsbpe$A_vsbpe_50nst<- (preddata_vsbpe$A_vsbpe_50*sd(ocu$A_vsbpe_50))+mean(ocu$A_vsbpe_50) #"des-standardize"
preddata_vsbpe$A_vsbpe_50porc<- (preddata_vsbpe$A_vsbpe_50nst*100/8247.656853) #percentage of SHV

#to build the plot
predplot_vsbpe<-ggplot(preddata_vsbpe,aes(x=A_vsbpe_50porc,y=Predicted))+
  labs(x=expression(paste('% SHV')), y=expression(paste('Pred. prob. ', psi)))+
  geom_ribbon(data = preddata_vsbpe,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf")+
  geom_line(data = preddata_vsbpe, colour="#000000", size=0.5)+
  annotate(geom="text", x=97, y=0.85, label="B",
           color="black", size=5) +
  theme_classic()
predplot_vsbpe

#detection variable (hour)
newData_hora <- data.frame(hora=seq(1, 51, length=51))
preddata_statehora <-predict(ab_ocsim22, type="det", newdata=newData_hora, appendData=TRUE)
hora_real<- read.table("C:/Users/Mauricio Diaz/Documents/U/Tesis/Modelo de ocupacion/Datos y script/horas.txt", sep="\t")
preddata_statehora$hora_real<-hora_real[,1]

#to build the plot
predplot_hora <- ggplot(preddata_statehora,aes(x=hora,y=Predicted)) +
  labs(x="Hour",y="Pred. prob. p") +
  geom_ribbon(data = preddata_statehora,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf") +
  geom_line(data = preddata_statehora, colour="#000000", size=0.5) +
  scale_x_discrete(limits = preddata_statehora$hora_real, guide = guide_axis(angle = 90)) +
  annotate(geom="text", x=48, y=0.17, label="C",
           color="black", size=5) +
  theme_classic() 
predplot_hora

grid.arrange(predplot_sdt,  predplot_vsbpe, predplot_hora, nrow = 3, ncol = 1) 

#model chi-squared

oc_m22_ba<-mb.gof.test(ab_ocsim22, nsim = 500, plot.hist = TRUE)
oc_m22_ba
hist(oc_m22_ba$t.star, xlab=expression(paste("Xi"^"2")), ylab="Frecuencia", col="lightgrey",font.lab=2, cex.lab=0.9, main="Prueba bondad de ajuste modelo", sub="P-value = 0.342 / c-hat = 0.99")
abline(v=oc_m22_ba$chi.square, lty=2, lwd=3, col="red")

###Spatial projection####
library(sf) #shapes
library(rgdal) #raster
library(fasterize) #raster
library(raster) #raster
library(terra) #raster

cobs<-rast("C:/Users/Mauricio Diaz/Documents/U/Tesis/Becas/San pedro/Buffer/project_vsbpe/cobs_vsb_peWGS84.tif")
sdt <- rast("C:/Users/Mauricio Diaz/Documents/U/Tesis/Becas/San pedro/Buffer/project_sdt/dxy_53m.tif")

cobs.ras <- aggregate(cobs, fact=53, fun="sum") 
rast53_vsbpe<-cobs.ras

rast53_sdt<- sdt
ext(rast53_vsbpe)<- ext(rast53_sdt)
rast53_vsbpe <- resample(rast53_vsbpe, rast53_sdt) 

x.sdt.st<- (rast53_sdt - mean(ocu$Sdt_b500))/sd(ocu$Sdt_b500)
x.cob.st<- (rast53_vsbpe - mean(ocu$A_vsbpe_50))/sd(ocu$A_vsbpe_50)

A_vsbpe_50<-x.cob.st
Sdt_b500<-x.sdt.st

x.stack<-c(A_vsbpe_50, Sdt_b500) 
x.stack<- stack(x.stack)
names(x.stack)<- c("A_vsbpe_50", "Sdt_b500") 
x.stack

E.Psi<- predict(ab_ocsim22, type="state", newdata= x.stack)
raster::crs(E.Psi)<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 
