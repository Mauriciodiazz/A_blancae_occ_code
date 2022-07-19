#Modelo de ocupacion para Atlapetes blancaeeee
#Mauricio Díaz

setwd("C:/Users/Mauricio Diaz/Documents/U/Tesis/Modelo de ocupacion/Datos y script")
#cargar historias de deteccion y vbles
hd<- read.table('Historias de deteccion.txt', header=T, sep = "\t")
head(hd)
#ocupacion
ocu<-read.table('ocupacion.txt', header=T, sep = "\t")
head(ocu)
#deteccion
det_ACI<-read.table('det_ACI.txt', header=T, sep = "\t")
det_hora<-read.table('det_hora.txt', header=T, sep = "\t")
det_ruido<-read.table('det_ruido.txt', header=T, sep = "\t")

#organizar historias deteccion
abla_hd<-hd[,2:8]
str(abla_hd)
names(abla_hd)

#organizar variables de ocupacion (site)
head(ocu)
names(ocu)
ocu_num<- scale(ocu[,c(-1,-9)])
id<- ocu[,1] #ID de grabadoras
site.cov<- data.frame(ocu_num,Cob_ROD=as.factor(ocu$Cob_Rod), A_vsbpe_50C=scale((ocu$A_vsbpe_50)^2))#debo convertir los characters in factors
head(site.cov)


##organizar variables de deteccion (obs)
#las variables de deteccion estan asociadas a cada sitio y cada repeticion, por ejemplo sitio 0 dia 1, sitio 0 dia 2 y as? suscesivamente 

#debo convertir las variuables de deteccion en una lista
obs.cov<-list(hora=matrix(scale(as.numeric(as.matrix(det_hora))),nrow=80, ncol=7), ruido=matrix(scale(as.numeric(as.matrix(det_ruido))),nrow=80, ncol=7), ACI=matrix(scale(as.numeric(as.matrix(det_ACI))),nrow=80, ncol=7)) #las variables ya estan estandarizadas, se debe estandarizar como si fueran una matriz, de modo contrario R asume que cada columna es una variable, de esta manera asume que toda la matriz es una sola variable


#correlacion
library(corrplot)
vble.ocu<-site.cov[,c(2,4,5,8:28)]
cor_ocu2<- cor(vble.ocu, method = "spearman")
cor_ocu<- cor(site.cov[,1:19], method = "spearman")
corrplot.mixed(cor_ocu2, lower = "number", 
               lower.col ="black", number.cex=.7,
               upper="square", order="FPC", title=" ",
               tl.cex=.7, tl.pos="lt", diag = "u")

#ahora el modelo...
#cargar paqueteria
library(unmarked)
library(dplyr) #manejo de datos
library(kableExtra) #tablas
library(ggplot2)
library(gridExtra)
library(AICcmodavg)

abla_umf<- unmarkedFrameOccu(y=abla_hd, siteCovs = site.cov, obsCovs= obs.cov)

plot(abla_umf, ylab="Sitios", main="Historias de detecci?n Atlapetes blancae",xlab="D?as") #plot de las historias de derecci?n

#vamos pues con los modelos

#veamos el modelo de ocupaci?n simple

ab_det1 <- occu(~1~1, start=c(1,1), data=abla_umf) #modelo nulo, deteccion, ocupacion, poner  start=c(1,1) son valores de inicio en los que se acota o se da una franja donde el modelo hace las iteraciones,
#cambiar esos valores no cambia realmente la probabilidad de ocupaci?n. el primer numero es de la deteccion y segundo ocupacion

ab_det2<- occu(~hora~1, data=abla_umf)
ab_det3<- occu(~ACI~1, data=abla_umf) 
ab_det4<- occu(~ruido~1, data=abla_umf)
ab_det5<- occu(~hora+ACI+ruido~1, data=abla_umf)

fms_detn <- fitList ("p(.)psi(.)"               =ab_det1,   
                     "p(hora)psi(.)"            =ab_det2)
(ms_detn <- modSel(fms_detn))

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

#tabla submodelos de deteccion
kbl(detsimDF[,c(1,7,8)], caption = "Group Rows")%>% 
  kable_paper("striped", full_width = T) %>%
  pack_rows("Submodelos de detección",1, 2, bold=T) 

#modelo simple
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

AICc(ab_ocsim1, return.K = FALSE)
AICc(ab_ocsim2, return.K = FALSE)
AICc(ab_ocsim3, return.K = FALSE)
AICc(ab_ocsim4, return.K = FALSE)
AICc(ab_ocsim5, return.K = FALSE)
AICc(ab_ocsim6, return.K = FALSE)
AICc(ab_ocsim7, return.K = FALSE)
AICc(ab_ocsim8, return.K = FALSE)
AICc(ab_ocsim9, return.K = FALSE)
AICc(ab_ocsim10, return.K = FALSE)
AICc(ab_ocsim11, return.K = FALSE)
AICc(ab_ocsim12, return.K = FALSE)
AICc(ab_ocsim13, return.K = FALSE)
AICc(ab_ocsim14, return.K = FALSE)
AICc(ab_ocsim15, return.K = FALSE)
AICc(ab_ocsim16, return.K = FALSE)
AICc(ab_ocsim17, return.K = FALSE)
AICc(ab_ocsim18, return.K = FALSE)
AICc(ab_ocsim19, return.K = FALSE)
AICc(ab_ocsim20, return.K = FALSE)
AICc(ab_ocsim21, return.K = FALSE)
AICc(ab_ocsim22, return.K = FALSE) #mejor modelo
AICc(ab_ocsim23, return.K = FALSE)
AICc(ab_ocsim24, return.K = FALSE)
AICc(ab_ocsim27, return.K = FALSE)
AICc(ab_ocsim28, return.K = FALSE)
AICc(ab_ocsim29, return.K = FALSE)

#esto es para la tabla bonita - AICc
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
ocsim22<-AICc(ab_ocsim22, return.K = FALSE) #mejor modelo
ocsim23<-AICc(ab_ocsim23, return.K = FALSE)
ocsim24<-AICc(ab_ocsim24, return.K = FALSE)
ocsim27<-AICc(ab_ocsim27, return.K = FALSE)
ocsim28<-AICc(ab_ocsim28, return.K = FALSE)
ocsim29<-AICc(ab_ocsim29, return.K = FALSE)

#esto es para la tabla bonita - estos vectores contienen el resumen de cada modelo
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
Socsim22<-summary(ab_ocsim22) #mejor modelo
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
                    delta_AICc=NA) #primero debo hacer un dataframe porque msdetn es un objeto de unmarked

#calculo del delta de AICc
for (x in 1:length(ocsim1DF$model)) {
  ocsim1DF[x,8]<-ocsim1DF[x,7]-ocsim1DF[1,7]
}

ocsim1DF <- ocsim1DF[order(ocsim1DF$AICc),]
write.table(ocsim1DF, "modelos11072022.txt", sep="\t")

#ahora voy a hacer una tabla con el modelo
ab_ocsim22
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


# Predicci?n de las variables
########################################d###
#primero los promedios
#para ocupacion
ocu.mean<-predict(ab_ocsim22, type="state")
mean(ocu.mean$Predicted) #ocupaci?n
mean(ocu.mean$SE) #error estandar asociado
#para deteccion
det.mean<-predict(ab_ocsim22, type="det")
mean(det.mean$Predicted) #ocupaci?n
mean(det.mean$SE) #error estandar asociado

ab_ocsim22
Cov_Sdt<- seq(range(abla_umf@siteCovs$Sdt_b500)[1]-0.25,
             range(abla_umf@siteCovs$Sdt_b500)[2]+0.25,
             length.out=100)

Cov_vsbpe <- seq(range(abla_umf@siteCovs$A_vsbpe_50)[1]-0.25,
                 range(abla_umf@siteCovs$A_vsbpe_50)[2]+0.25,
                 length.out=100)

#para terrain second derivative (sdt)
newData_sdt<- data.frame(Sdt_b500=Cov_Sdt, A_vsbpe_50=0)
preddata_sdt <- predict(ab_ocsim22, newdata=newData_sdt, type="state")
preddata_sdt$Sdt_b500<- Cov_Sdt
preddata_sdt$Sdt_b500nst<- (preddata_sdt$Sdt_b500*sd(ocu$Sdt_b500))+mean(ocu$Sdt_b500) #"des estandarizar"
  
head(preddata_sdt)
predplot_sdt<-ggplot(preddata_sdt,aes(x=Sdt_b500nst,y=Predicted)) +
  labs(x="Terrain second derivative (rad/m)", y=expression(paste('Pred. prob. ', psi))) +
  geom_ribbon(data = preddata_sdt,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf") +
  geom_line(data = preddata_sdt, colour="#000000", size=0.5) +
  annotate(geom="text", x=0.0002, y=1, label="A",
           color="black", size=5) +
  theme_classic()
predplot_sdt

#para área de vsb y pe
newData_vsbpe<- data.frame(Sdt_b500=0, A_vsbpe_50=Cov_vsbpe)

preddata_vsbpe <- predict(ab_ocsim22, newdata=newData_vsbpe, type="state")
preddata_vsbpe$A_vsbpe_50<- Cov_vsbpe
preddata_vsbpe$A_vsbpe_50nst<- (preddata_vsbpe$A_vsbpe_50*sd(ocu$A_vsbpe_50))+mean(ocu$A_vsbpe_50) ##"des estandarizar"
preddata_vsbpe$A_vsbpe_50porc<- (preddata_vsbpe$A_vsbpe_50nst*100/8247.656853)
head(preddata_vsbpe)

predplot_vsbpe<-ggplot(preddata_vsbpe,aes(x=A_vsbpe_50porc,y=Predicted))+
  labs(x=expression(paste('% SHV')), y=expression(paste('Pred. prob. ', psi)))+
  geom_ribbon(data = preddata_vsbpe,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf")+
  geom_line(data = preddata_vsbpe, colour="#000000", size=0.5)+
  annotate(geom="text", x=97, y=0.85, label="B",
           color="black", size=5) +
  theme_classic()
predplot_vsbpe

#deteccion para hora
newData_hora <- data.frame(hora=seq(1, 51, length=51))
preddata_statehora <-predict(ab_ocsim22, type="det", newdata=newData_hora, appendData=TRUE)
hora_real<- read.table("C:/Users/Mauricio Diaz/Documents/U/Tesis/Modelo de ocupacion/Datos y script/horas.txt", sep="\t")
preddata_statehora$hora_real<-hora_real[,1]

predplot_hora <- ggplot(preddata_statehora,aes(x=hora,y=Predicted)) +
  labs(x="Hour",y="Pred. prob. p") +
  geom_ribbon(data = preddata_statehora,aes(ymin=lower,ymax=upper),alpha=0.8, fill = "#dfdfdf") +
  geom_line(data = preddata_statehora, colour="#000000", size=0.5) +
  scale_x_discrete(limits = preddata_statehora$hora_real, guide = guide_axis(angle = 90)) +
  annotate(geom="text", x=48, y=0.17, label="C",
           color="black", size=5) +
  theme_classic() 

predplot_hora

grid.arrange(predplot_sdt,  predplot_vsbpe, predplot_hora, nrow = 3, ncol = 1) #todos en una sola tabla

########estos vectores contienen los valores no estandarizados de las variables para ambas predicciones##
vsbpe_nsta<-(preddata_vsbpe$A_vsbpe_50*sd(ocu$A_vsbpe_50))+mean(ocu$A_vsbpe_50)
sdt_nsta<-(preddata_sdt$Sdt_b200*sd(ocu$Sdt_b200))+mean(ocu$Sdt_b200)

#¿cuales son los promedios de vegetación para los sitios con presencia?
#12, 33, 38, 51, 52, 69, 70
vsb.pre<- ocu[c(13, 34, 39, 52, 53, 70, 71),c(1,29)]
mean(vsb.pre[,2])
sd(vsb.pre[,2])

#chi cuadrado para el modelos

library(AICcmodavg)
#Ocu
oc_m22_ba<-mb.gof.test(ab_ocsim22, nsim = 500, plot.hist = TRUE)
oc_m22_ba
hist(oc_m22_ba$t.star, xlab=expression(paste("Xi"^"2")), ylab="Frecuencia", col="lightgrey",font.lab=2, cex.lab=0.9, main="Prueba bondad de ajuste modelo", sub="P-value = 0.342 / c-hat = 0.99")
abline(v=oc_m22_ba$chi.square, lty=2, lwd=3, col="red")

########################################d###
#¿cuál es el resultado de la ocupación y la detección?
backTransform(linearComb(ab_ocsim27, coefficients=c(1,0,0,0), type="state"))
#estimado es el lambda, es decir un promedio de las abundancias locales de los sitios de muestreo. Se peude interpretar como una abundancia relativa. indice de abundancia por unidad de muestreo, 1.95 individuos por sitio. con esta funci?n voy a poner los valores en "coefficients" de los valores que quiero para cada variable, si son 4 variables, deben de ir 5 numeros; el primero para el intercepto y los otros para el valor de cada variable
#no da un numero de individuos, pero si da un ?ndice de la intensidad de la abundancia en el ?rea de interes. Intensidad de puntos con los que se tiene una aproximacion de que te tan abundante es la especie en el sitio. Se puede obtener una abundancia local para cada grabadora, pero el que arroja el modelo de primera es el promedio de todas las abundancias locales de los sitios definidos

backTransform(linearComb(ablamod7, coefficients=c(2,1), type="det"))
#la probabilidad de deteccion se hace mas baja porque los modelos RN se usan para evaluar heterogeniedad entre sitios y como hace eso, entonces ajusta los valores, por eso a veces se usa para corregir. Si es muy baja, se debe mirar  
#intervalos de confianza


#vamos con la proyección espacial pues, respira

#vsb y pe: para esta capa hay que transformar el shapefile y unir las dos categorías; vegetación secundaria baja y pastos enmalezados
library(sf) #tratamiento de shapes
library(rgdal) #tratamiento de raster
library(fasterize)
library(terra)

#voy a sacar las variables para la proyecci?n espacial del modelo

#cargar archivo raster con variable de interes
cobs<-rast("C:/Users/Asus/Documents/U/Tesis/Becas/San pedro/Buffer/project_vsbpe/cobs_vsb_peWGS84.tif")
sdt <- rast("C:/Users/Asus/Documents/U/Tesis/Becas/San pedro/Buffer/project_sdt/dxy_53m.tif")
#sdt<- project(sdt, y = "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") #OJO reproyectar raster puede cambiar los valores de las celdas


#si la variable se quiere remuestrear a un tamaño similar al de la celda,#entonces se puede hacer lo siguiente tratando de que el area de la celda#sea lo mas similar al area del buffer que se quiere representar, en este #caso particular, el buffer es de 30m 
cobs.ras <- aggregate(cobs, fact=53, fun="sum") #con esto sale la de 
rast53_vsbpe<-cobs.ras
#rast53 <- aggregate(cobs, fact=53, fun="sum") #este es el raster "base"
rast53_sdt<- sdt
ext(rast53_vsbpe)<- ext(rast53_sdt)
rast53_vsbpe <- resample(rast53_vsbpe, rast53_sdt) 

plot(rast53_sdt)
plot(rast53_vsbpe, add=T)



#Si el interes es representar un buffer extenso pero para una celda de menor 
#tamano, por ejemplo, si queremos representar lo presente en un buffer de 500m
#de radio pero en una celda de 53m de lado, hay dos vias:

#1 VIA: esta seria la via estrictamente correcta PERO al menos como lo estoy 
#haciendo es inviable por tiempo
#convierto raster a puntos
#rast53_p <- as.points(rast53)

#realizo el buffer y extraigo el promedio para el buffer exacto para cada punto
#start_time <- Sys.time()
#out <- extract(x=rast10_sdt, y=rast53_p)
#values(rast53) <- out

#end_time <- Sys.time()

#end_time-start_time
#en mi computador, se demora 5.21 minutos en hacer 100 buffer. Son 84970, 
#osea que, en mi compu, seria un poco mas de 70 horas!

#2 VIA: en esta calculamos el buffer como si fuera un cuadrado con los mismos 
#lineamientos de antes (por ejemplo, el area de un cuadrado de 886m de lado 
#es similar al de un buffer de 500m de radio)
#rast_agg_1000 <- aggregate(sdt, fact=10, fun=mean) #este es el raster de la variable grande

#extraemos los valores de TODAS las celdas incluyendo las que no tienen valor asignado del raster
#rast53_p_all <- values(rast53) 
#averiguo cuales son las que tienen valores asignados
#cuales <- which(is.na(rast53_p_all[,1]) != TRUE, arr.ind=TRUE)
#extraigo para cada punto el valor del raster pero haciendo interpolacion (esto ultimo es importante porque a cada punto no simplemente repite el valor de la celda) 
#rast53_p_valores <- extract(rast_agg_1000, rast53_p, method="bilinear")

#values(rast53)[cuales] <- rast53_p_valores
#rast53_p_all[cuales,1] <- rast53_p_valores$sdt_recortado
#Luego le asigno esos valores al raster original
#values(rast53) <- rast53_p_all
#rast53_sdt<-rast53
#rast53_vsbpe

terra::writeRaster(x=rast53_sdt, "rast53_sdt.tif", overwrite=TRUE)
terra::writeRaster(x=rast53_vsbpe, "rast53_vsbpe.tif", overwrite=TRUE)

cob.shp<-st_read("./fwdcoberturasshape/Cob_1600_2_WGS84.shp") 
#2. cargar el shape que va a servir de máscara
mpios<- st_read("./municipios//Municipios_100_Febrero_2012_3116.shp")

#voy a hacer una seleccion de los municipios al rededor de san pedro para el mapa
which(mpios$NOMBRE_ENT == "SAN PEDRO")
which(mpios$NOMBRE_ENT == "BELLO")
which(mpios$NOMBRE_ENT == "SAN JERÓNIMO")
which(mpios$NOMBRE_ENT == "SOPETRÁN")
which(mpios$NOMBRE_ENT == "BELMIRA")
which(mpios$NOMBRE_ENT == "ENTRERRIOS")
which(mpios$NOMBRE_ENT == "SANTA ROSA DE OSOS")
which(mpios$NOMBRE_ENT == "DON MATÍAS")
which(mpios$NOMBRE_ENT == "GIRARDOTA")
which(mpios$NOMBRE_ENT == "COPACABANA")
which(mpios$NOMBRE_ENT == "BARBOSA")
which(mpios$NOMBRE_ENT == "MARINILLA")
which(mpios$NOMBRE_ENT == "RIONEGRO")
which(mpios$NOMBRE_ENT == "MEDELLÍN")
which(mpios$NOMBRE_ENT == "EBÉJICO")


sanp<-mpios[631,]
#Aqui hago la selección
area<-mpios[c(631,612,1053,651,676,658,688,642,613,1050,624,584,581,1124,610),]

#sanp<- st_transform(sanp, "+proj=tmerc +lat_0=4.59620041666667 +lon_0=-74.0775079166667 +k=1 +x_0=1000000+y_0=1000000 +ellps=GRS80 +units=m +no_defs") con esta función puedo reproyectar un shape!!!

#debo estandarizar los raster?????????????????????? Si, pero la media y la sd debe ser la de los datos que le meti al modelo para ajustar los datos del raster a los datos predichos inicialmente

x.sdt.st<- (rast53_sdt - mean(ocu$Sdt_b500))/sd(ocu$Sdt_b500)
x.cob.st<- (rast53_vsbpe - mean(ocu$A_vsbpe_50))/sd(ocu$A_vsbpe_50)

hist(x.cob.st)
hist(x.sdt.st)

#para correr predict() deben de tener los mismos nombres de las variables que use en el modelo

A_vsbpe_50<-x.cob.st
Sdt_b500<-x.sdt.st

#se crea un stack con ambos raters
x.stack<-c(A_vsbpe_50, Sdt_b500) #en el paquete terra la funci?n de stack es c
x.stack<- stack(x.stack)
plot(x.stack)
names(x.stack)<- c("A_vsbpe_50", "Sdt_b500") #recuerda: mismos nombres que el modelo!
x.stack

#se corre el predict
E.Psi<- predict(ab_ocsim22, type="state", newdata= x.stack)
#hay que asignarle una proyección
raster::crs(E.Psi)<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 

#guardar el raster
raster::writeRaster(x= E.Psi, filename = "Stack_modoc_abla_5001.tif", overwrite=TRUE)
stackSave(x= E.Psi, filename = "Stack_modoc_abla_400")


a<-data.frame(E.Psi$Predicted)
a<-na.omit(a)
hist(a$Predicted)
summary(a)