#Session<set Working Directory<To source file location, per introduïr les dades
setwd("C:/Users/Usuario/Desktop/Xavi/UDG/4rt/TFG/R")
dades_inicial<-read.csv(file = "dadess_balancejat.csv", header = T, sep = ";",dec=",")
#O bé copiant les dades de l'excel amb:
dades_inicial<-read.table("clipboard",header=T,sep="\t")
names(dades_inicial)
summary(dades_inicial)


#Transformem a factor les variables que els hi calgui.
dades_inicial$Mesocosmos<-factor(dades_inicial$Mesocosmos)
dades_inicial$tractament<-factor(dades_inicial$tractament)
dades_inicial$data<-factor(dades_inicial$data)
dades_inicial$sampleID<-factor(dades_inicial$sampleID)
dades_inicial$grup<-factor(dades_inicial$grup)
dades_inicial$codi<-factor(dades_inicial$codi)
dades_inicial$taxon<-factor(dades_inicial$taxon)
summary(dades_inicial)


#calcul de la biomassa individual (en microg per individu)
dades_inicial$biom_ind<-ifelse(dades_inicial$grup=="CHI",   5.1*((dades_inicial$a/1000)^2.32),
                ifelse(dades_inicial$grup=="AMP",  5.7*((dades_inicial$a/1000)^2.721),
                ifelse(dades_inicial$grup=="GAS",  27.3*((dades_inicial$a/1000)^2.4992),
                ifelse(dades_inicial$grup=="OST",  28.42*((dades_inicial$a/1000)^2.8),
                ifelse(dades_inicial$grup=="CAL",  2.9946*((dades_inicial$a/1000)^2.1951),
                ifelse(dades_inicial$grup=="CYC",  1.8504*((dades_inicial$a/1000)^2.034),
                ifelse(dades_inicial$grup=="HAR",  1.8504*((dades_inicial$a/1000)^2.034),
                ifelse(dades_inicial$grup=="ROT",  36.4*(dades_inicial$a/1000)*((dades_inicial$b/1000)^2),
                ifelse(dades_inicial$grup=="CIL",  98.8*(dades_inicial$a/1000)*((dades_inicial$b/1000)^2),
               NA)))))))))


#Fem transformacións necesàries a les variables per poder treballar amb elles:
dades_inicial$countsxdilucioperL<-((dades_inicial$counts/dades_inicial$mesures)*dades_inicial$dilucio)/dades_inicial$volumf
dades_inicial$biom_total<-dades_inicial$biom_ind*dades_inicial$countsxdilucioperL
summary(dades_inicial)


#Calculem la biomassa total per cada taxó i tractament.
library(dplyr)
.dades_grup_taxo_replica_tract_data_sample<-group_by(dades_inicial,Mesocosmos,taxon,tractament,replica,data,sampleID)
.dades_grup_taxo_replica_tract_data_sample
summarise(.dades_grup_taxo_replica_tract_data_sample)
dades_grup_biom<-summarise(.dades_grup_taxo_replica_tract_data_sample,suma_biomasa=sum(biom_total,na.rm=T))
dades_grup_biom

#Calculem els countsxdilucioperL totals per cada taxo i tractament (tenint en compte la dilució, mesures i volum)
dades_grup_counts<-summarise(.dades_grup_taxo_replica_tract_data_sample,counts_total=sum(countsxdilucioperL))
dades_grup_counts


#reorganitzem la taula de biomassa per poder fer desprès l'anàlisis d'ordenació.
library(reshape2)
dades_grup_biom_dcast<-dcast(dades_grup_biom,sampleID~taxon,value.var="suma_biomasa")
#Convertim la SampleID en rowname
row.names(dades_grup_biom_dcast) <- dades_grup_biom_dcast$sampleID
#Eliminem la primera columna (si no no podrem fer PCA, perquè totes les variables haurien de ser numèriques)
dades_grup_biom_dcast[1] <- NULL
#Substituïm els NA per zeros
dades_grup_biom_dcast[is.na(dades_grup_biom_dcast)] <- 0

#Fem el mateix per counts
dades_grup_counts_dcast<-dcast(dades_grup_counts,sampleID~taxon,value.var="counts_total")
#Convertim la SampleID en rowname
row.names(dades_grup_counts_dcast) <- dades_grup_counts_dcast$sampleID
#Eliminem la primera columna (si no no podrem fer PCA, perquè totes les variables haurien de ser numèriques)
dades_grup_counts_dcast[1] <- NULL
#Substituïm els NA per zeros
dades_grup_counts_dcast[is.na(dades_grup_counts_dcast)] <- 0


#Abans hem de fer una transformació Numerical ecology cap.3 (https://we.tl/t-hsBeZmOHrp). lg+1 o hellinger. 
#fem anàlisis d'ordenació PCA, per dades_grup_biomasa_dcast
library(ade4)
library(vegan)  
library(gclus) 
library(cluster)
library(RColorBrewer)  
library(labdsv)
library(ape)
#S'ha de transformar amb hellinger
dades_grup_biom_dcast_hellinger<-decostand(dades_grup_biom_dcast,"hellinger")
str(dades_grup_biom_dcast_hellinger)
dades_grup_biom_dcast_hellinger_PCA <- rda(dades_grup_biom_dcast_hellinger, scale = TRUE)
summary(dades_grup_biom_dcast_hellinger_PCA)
par(mfrow=c(1,2))
barplot(dades_grup_biom_dcast_hellinger_PCA$CA$eig, ylab='Eigenvalues', las=1, main='Kaiser-Gutman')
abline(h=mean(dades_grup_biom_dcast_hellinger_PCA$CA$eig), lty=2, col='dark red')
screeplot(dades_grup_biom_dcast_hellinger_PCA, bstick = TRUE, npcs = length(dades_grup_biom_dcast_hellinger_PCA$CA$eig), ylab='Eigenvalues', main='Broken-Stick')
#segons kaiser-Guntam hauriem de triar les 6 primeres dimensions, segons Broken-Stick hauriem de triar les domensions 1, 2, 3, 6, i 7.
#Grafiquem els eixos PC1 i PC2:
par(mfrow=c(1,2))
biplot(dades_grup_biom_dcast_hellinger_PCA, scaling = 1, main = "Biomass PCA - scaling 1", xlab = 'PC1 (23.9%)',  ylab = 'PC2 (15.7%)',cex.main=2.2,cex.lab=1.52,cex.axis=1.25, cex.text=1)
biplot(dades_grup_biom_dcast_hellinger_PCA, main = "PCA biomassa - scaling 2", xlab = 'PC1 (23.9%)',  ylab = 'PC2 (15.7%)',cex.main=2,cex.lab=1.5)
#Scalling 1 per veure la diferència entre objectes, i el 2 per veure com es relacionen les diferents variables.
par(mfrow=c(1,1))


#Fem l'anàlisis de PCA per dades_grup_counts_dcast
dades_grup_counts_dcast_hellinger<-decostand(dades_grup_counts_dcast,"hellinger")
str(dades_grup_counts_dcast_hellinger)
dades_grup_counts_dcast_hellinger_PCA <- rda(dades_grup_counts_dcast_hellinger, scale = TRUE)
summary(dades_grup_counts_dcast_hellinger_PCA)
par(mfrow=c(1,2))
barplot(dades_grup_counts_dcast_hellinger_PCA$CA$eig, ylab='Eigenvalues', las=1, main='Kaiser-Gutman')
abline(h=mean(dades_grup_counts_dcast_hellinger_PCA$CA$eig), lty=2, col='dark red')
screeplot(dades_grup_counts_dcast_hellinger_PCA, bstick = TRUE, npcs = length(dades_grup_counts_dcast_hellinger_PCA$CA$eig), ylab='Eigenvalues', main='Broken-Stick')
#segons kaiser-Guntam hauriem de triar les 6 primeres dimensions, segons Broken-Stick hauriem de triar les domensions 1, 2, 3, i 6.
#Grafiquem els eixos PC1 i PC2:
par(mfrow=c(1,2))
biplot(dades_grup_counts_dcast_hellinger_PCA, scaling = 1, main = "PCAcounts - scaling 1", xlab = 'PC1 (25.9%)',  ylab = 'PC2 (16.6%)')
biplot(dades_grup_counts_dcast_hellinger_PCA, main = "PCAcounts - scaling 2", xlab = 'PC1 (25.9%)',  ylab = 'PC2 (16.6%)')
#Scalling 1 per veure la diferència entre objectes, i el 2 per veure com es relacionen les diferents variables.
par(mfrow=c(1,1))


#Comparem els PCA 1 (diferències entre objectes) de biomassa i counts:
par(mfrow=c(1,2))
biplot(dades_grup_biom_dcast_hellinger_PCA, scaling = 1, main = "PCAbiom - scaling 1", xlab = 'PC1 (23.9%)',  ylab = 'PC2 (15.7%)')
biplot(dades_grup_counts_dcast_hellinger_PCA, scaling = 1, main = "PCAcounts - scaling 1", xlab = 'PC1 (25.9%)',  ylab = 'PC2 (16.6%)')


#Comparem els PCA2 (relacions entre variables) de biomassa i counts:
biplot(dades_grup_biom_dcast_hellinger_PCA, main = "PCAbiom - scaling 2", xlab = 'PC1 (23.9%)',  ylab = 'PC2 (15.7%)')
biplot(dades_grup_counts_dcast_hellinger_PCA, main = "PCAcounts - scaling 2", xlab = 'PC1 (25.9%)',  ylab = 'PC2 (16.6%)')
par(mfrow=c(1,1))



#Provem de fer un PRC (Principal Response Curves), primer per biomassa
#Principal response cool (segon link(http://edild.github.io/prc1/). Recodificar data 1,2, cada corva és un tractament, que es compara amb els controls. Haig de dir el control és tractament =o (relevel o algu aixi¡)).
#colocar una llista d'espècies de manera que els que estigion cap a baix siguin les més abundants en el control, i les que estàn cap a dalt les més característiques dels altres tractaments.
dades_grup_biom_dcast_prepararPRC<-dcast(dades_grup_biom,sampleID+Mesocosmos+data+tractament~taxon,value.var="suma_biomasa")
#Convertim la SampleID en rowname
row.names(dades_grup_biom_dcast_prepararPRC) <- dades_grup_biom_dcast_prepararPRC$sampleID
#Eliminem la primera columna (si no no podrem fer PCA, perquè totes les variables haurien de ser numèriques)
dades_grup_biom_dcast_prepararPRC[1] <- NULL
#Substituïm els NA per zeros
dades_grup_biom_dcast_prepararPRC[is.na(dades_grup_biom_dcast_prepararPRC)] <- 0
#Eliminem package dplyr perquè fa interferència amb el pack car (si no no podem fer recode).
detach("package:dplyr", unload = TRUE)
library(car)
dades_grup_biom_dcast_prepararPRC$data<-recode(dades_grup_biom_dcast_prepararPRC$data,"'29/03/2020'=2;'02/03/2020'=1")
dades_grup_biom_dcast_prepararPRC$data
dades_grup_biom_dcast_prepararPRC$tractament<-relevel(dades_grup_biom_dcast_prepararPRC$tractament,ref = "O")
dades_grup_biom_dcast_prepararPRC$tractament

#Fem el PRC amb les dades transformades per hellinger, aprofitant la transformació feta en el PCA
dades_grup_biom_dcast_PRC <- prc(response = dades_grup_biom_dcast_hellinger, treatment = dades_grup_biom_dcast_prepararPRC$tractament, time = dades_grup_biom_dcast_prepararPRC$data)
dades_grup_biom_dcast_PRC
plot(dades_grup_biom_dcast_PRC, scaling = 1)
#sembla ser que el tractament P és el que més s'allunya del control, per tant és el que més pertorbacions crea a la xarxa de comunitats aqüàtiques, el segueix el Continu (C) i per últim el PC.
.dades_grup_biom_dcast_PRC_sum <- summary(dades_grup_biom_dcast_PRC)
.dades_grup_biom_dcast_PRC_sum
.dades_grup_biom_dcast_PRC_sum$sp
head(.dades_grup_biom_dcast_PRC_sum$sp)

#es fiquen al gràfic només les espècies amb scores majors a 0,5.
###############################################################################################################################
##########Al transformar les dades per hellinger, la següent comanda em dona un missatge d'advertència que no se com arreglar (sense transformar les dades no em sortia l'error)
plot(dades_grup_biom_dcast_PRC, select = abs(.dades_grup_biom_dcast_PRC_sum$sp) > 0.5, scaling = 1)
# PRC plot; at the right of it, only species with large total
# (log-transformed) abundances are reported
dades_grup_biom_dcast_log<-colSums(dades_grup_biom_dcast_hellinger)
dades_grup_biom_dcast_log
plot(dades_grup_biom_dcast_PRC,select = dades_grup_biom_dcast_log>200)
# Statistical test
# Ditches are randomized, we have a time series, and are only
# interested in the first axis
.PRC_biom_permutations <-
  how(plots = Plots(strata = dades_grup_biom_dcast_prepararPRC$Mesocosmos, type = "free"), 
      within = Within(type = "series"), nperm = 999)
.PRC_biom_permutations
anova(dades_grup_biom_dcast_PRC, permutations = .PRC_biom_permutations, first = TRUE)
#Hi ha més diferències entre tractaments que dins de cada tractament.



#Provem de fer un PRC (Principal Response Curves), ara per counts
dades_grup_counts_dcast_prepararPRC<-dcast(dades_grup_counts,sampleID+Mesocosmos+data+tractament~taxon,value.var="counts_total")
#Convertim la SampleID en rowname
row.names(dades_grup_counts_dcast_prepararPRC) <- dades_grup_counts_dcast_prepararPRC$sampleID
#Eliminem la primera columna (si no no podrem fer PCA, perquè totes les variables haurien de ser numèriques)
dades_grup_counts_dcast_prepararPRC[1] <- NULL
#Substituïm els NA per zeros
dades_grup_counts_dcast_prepararPRC[is.na(dades_grup_counts_dcast_prepararPRC)] <- 0
#Eliminem package dplyr perquè fa interferència amb el pack car (si no no podem fer recode).
detach("package:dplyr", unload = TRUE)
library(car)
dades_grup_counts_dcast_prepararPRC$data<-recode(dades_grup_counts_dcast_prepararPRC$data,"'29/03/2020'=2;'02/03/2020'=1")
dades_grup_counts_dcast_prepararPRC$data
dades_grup_counts_dcast_prepararPRC$tractament<-relevel(dades_grup_counts_dcast_prepararPRC$tractament,ref = "O")
dades_grup_counts_dcast_prepararPRC$tractament

#Fem el PRC amb les dades transformades per hellinger, aprofitant la transformació feta en el PCA
dades_grup_counts_dcast_PRC <- prc(response = dades_grup_counts_dcast_hellinger, treatment = dades_grup_counts_dcast_prepararPRC$tractament, time = dades_grup_counts_dcast_prepararPRC$data)
dades_grup_counts_dcast_PRC
plot(dades_grup_counts_dcast_PRC, scaling = 1)
#sembla ser que el tractament P és el que més s'allunya del control, per tant és el que més pertorbacións crea a la xarxa de comunitats aqüàtiques, el segueix el (PC) i per últim el C.
.dades_grup_counts_dcast_PRC_sum <- summary(dades_grup_counts_dcast_PRC)
.dades_grup_counts_dcast_PRC_sum
.dades_grup_counts_dcast_PRC_sum$sp
head(.dades_grup_counts_dcast_PRC_sum$sp)

#es fiquen al gràfic només les espècies amb scores majors a 0,5.
plot(dades_grup_counts_dcast_PRC, select = abs(.dades_grup_counts_dcast_PRC_sum$sp) > 0.5, scaling = 1)
# PRC plot; at the right of it, only species with large total
# (log-transformed) abundances are reported
.dades_grup_counts_dcast_log<-colSums(dades_grup_counts_dcast_hellinger)
.dades_grup_counts_dcast_log
plot(dades_grup_counts_dcast_PRC,select = .dades_grup_counts_dcast_log>200)
# Statistical test
# Ditches are randomized, we have a time series, and are only
# interested in the first axis
.PRC_counts_permutations <-
  how(plots = Plots(strata = dades_grup_counts_dcast_prepararPRC$Mesocosmos, type = "free"), 
      within = Within(type = "series"), nperm = 999)
.PRC_counts_permutations
anova(dades_grup_counts_dcast_PRC, permutations = .PRC_counts_permutations, first = TRUE)
#Hi ha més diferències entre tractaments que dins de cada tractament.



#comparem els PRC biomasa i counts
par(mfrow=c(1,2))
###############################################################################################################################
##########No se perquè, però no es veuen els 2 gràfics alhora
plot(dades_grup_biom_dcast_PRC, scaling = 1)
plot(dades_grup_counts_dcast_PRC, scaling = 1)

plot(dades_grup_biom_dcast_PRC, select = abs(.dades_grup_biom_dcast_PRC_sum$sp) > 0.5, scaling = 1)
plot(dades_grup_counts_dcast_PRC, select = abs(.dades_grup_counts_dcast_PRC_sum$sp) > 0.5, scaling = 1)

plot(dades_grup_biom_dcast_PRC,select = dades_grup_biom_dcast_log>200)
plot(dades_grup_counts_dcast_PRC,select = .dades_grup_counts_dcast_log>200)
par(mfrow=c(1,1))





#Per avaluar l'estructura de la comunitat:
#preparem la taula de dades
library(dplyr)
.dades_grup_taxo_tract_data<-group_by(dades_inicial,taxon,tractament,data,sampleID,replica)
.dades_grup_taxo_tract_data
summarise(.dades_grup_taxo_tract_data)
dades_grup_counts_estructura<-summarise(.dades_grup_taxo_tract_data,counts_tot=sum(countsxdilucioperL,na.rm=T))
dades_grup_counts_estructura
dades_grup_counts_estructura_dcast<-dcast(dades_grup_counts_estructura,tractament+data+replica~taxon,value.var="counts_tot")


Cop.Calanoida<-dades_grup_counts_estructura_dcast$`Copepoda Calanoida`+dades_grup_counts_estructura_dcast$`Copepoda Calanoida + ous`+dades_grup_counts_estructura_dcast$`Nauplius Calanoida`
dades_grup_counts_estructura_dcast[21]<-Cop.Calanoida
dades_grup_counts_estructura_dcast[6]<-NULL
dades_grup_counts_estructura_dcast[6]<-NULL
dades_grup_counts_estructura_dcast[13]<-NULL
dades_grup_counts_estructura_dcast<-rename(dades_grup_counts_estructura_dcast,"Calanoida"=V21)

Cop.Cyclopoida<-dades_grup_counts_estructura_dcast$`Copepoda Cyclopoida`+dades_grup_counts_estructura_dcast$`Copepoda Cyclopoida + ous`+dades_grup_counts_estructura_dcast$`Nauplius Cyclopoida`
dades_grup_counts_estructura_dcast[19]<-Cop.Cyclopoida
dades_grup_counts_estructura_dcast[6]<-NULL
dades_grup_counts_estructura_dcast[6]<-NULL
dades_grup_counts_estructura_dcast[11]<-NULL
dades_grup_counts_estructura_dcast<-rename(dades_grup_counts_estructura_dcast,"Cyclopoida"=V19)

Ostrac<-dades_grup_counts_estructura_dcast$`Larva Ostracoda`+dades_grup_counts_estructura_dcast$Ostracoda
dades_grup_counts_estructura_dcast[17]<-Ostrac
dades_grup_counts_estructura_dcast[10]<-NULL
dades_grup_counts_estructura_dcast[11]<-NULL
dades_grup_counts_estructura_dcast<-rename(dades_grup_counts_estructura_dcast,"Ostracoda"=V17)

Harpact<-dades_grup_counts_estructura_dcast$`Copepoda Harpacticoida`+dades_grup_counts_estructura_dcast$`Nauplius Harpacticoida`
dades_grup_counts_estructura_dcast[16]<-Harpact
dades_grup_counts_estructura_dcast[6]<-NULL
dades_grup_counts_estructura_dcast[9]<-NULL
dades_grup_counts_estructura_dcast<-rename(dades_grup_counts_estructura_dcast,"Harpacticoida"=V16)


dades_grup_counts_estructura_dcast[is.na(dades_grup_counts_estructura_dcast)] <- 0
#Afegim columna d'abundància d'espècies
.dades_grup_counts_estructura_dcast_abundanciaa<-rowSums(dades_grup_counts_estructura_dcast[4:14])
.dades_grup_counts_estructura_dcast_abundanciaa
dades_grup_counts_estructura_dcast_final<-dades_grup_counts_estructura_dcast
dades_grup_counts_estructura_dcast_final[4:14]<-NULL
dades_grup_counts_estructura_dcast_final[4]<-.dades_grup_counts_estructura_dcast_abundanciaa
dades_grup_counts_estructura_dcast_final<-rename(dades_grup_counts_estructura_dcast_final,abundancia=V4)
#Afegim columna de riquesa
.dades_grup_counts_estructura_dcast_riquesaa<-dades_grup_counts_estructura_dcast
.dades_grup_counts_estructura_dcast_riquesaa2<-.dades_grup_counts_estructura_dcast_riquesaa[4:14]>0
detach("package:dplyr", unload = TRUE)
library(car)
.dades_grup_counts_estructura_dcast_riquesaa2<-recode(.dades_grup_counts_estructura_dcast_riquesaa2,"'FALSE'=0;'TRUE'=1")
library(dplyr)
.dades_grup_counts_estructura_dcast_riquesaa2
.dades_grup_counts_estructura_dcast_riquesa<-rowSums(.dades_grup_counts_estructura_dcast_riquesaa2)
.dades_grup_counts_estructura_dcast_riquesa
dades_grup_counts_estructura_dcast_final[5]<-.dades_grup_counts_estructura_dcast_riquesa
dades_grup_counts_estructura_dcast_final<-rename(dades_grup_counts_estructura_dcast_final,Riquesa=V5)
#Afegim columna de shannon
.dades_grup_counts_estructura_dcast_shannon<-diversity(dades_grup_counts_estructura_dcast[4:14],index = "shannon")
.dades_grup_counts_estructura_dcast_shannon
dades_grup_counts_estructura_dcast_final[6]<-.dades_grup_counts_estructura_dcast_shannon
dades_grup_counts_estructura_dcast_final<-rename(dades_grup_counts_estructura_dcast_final,Shannon=V6)
summary(dades_grup_counts_estructura_dcast_final)
dades_grup_counts_estructura_dcast_final$tractament_dia<-dades_grup_counts_estructura_dcast_final$tractament
dades_grup_counts_estructura_dcast_final$tractament_dia<-c(rep("C-i",3),rep("C-f",3),rep("O-i",5),rep("O-f",5),rep("P-i",4),rep("P-f",4),rep("PC-i",3),rep("PC-f",3))
dades_grup_counts_estructura_dcast_final$data_tract<-dades_grup_counts_estructura_dcast_final$data
dades_grup_counts_estructura_dcast_final$data_tract<-c(rep("02/03/2020-nutrients",3),rep("29/03/2020-nutrients",3),rep("02/03/2020-control",5),rep("29/03/2020-control",5),rep("02/03/2020-nutrients",4),rep("29/03/2020-control",4),rep("02/03/2020-nutrients",3),rep("29/03/2020-control",3))
dades_grup_counts_estructura_dcast_final$data_tract<-factor(dades_grup_counts_estructura_dcast_final$data_tract,levels=c("02/03/2020-control","29/03/2020-control","02/03/2020-nutrients","29/03/2020-nutrients"))
dades_grup_counts_estructura_dcast_final$tractament<-factor(dades_grup_counts_estructura_dcast_final$tractament,levels=c("O","C","P","PC"))
dades_grup_counts_estructura_dcast_final$tractament_dia<-factor(dades_grup_counts_estructura_dcast_final$tractament_dia,levels=c("O-i","O-f","PC-i","PC-f","C-i","C-f","P-i","P-f"))

#ANOVA  de la abundancia
ANOVA_abundancia<-lm(dades_grup_counts_estructura_dcast_final$abundancia~dades_grup_counts_estructura_dcast_final$tractament*dades_grup_counts_estructura_dcast_final$data)
library(car)
leveneTest(ANOVA_abundancia)#ok homogeneïtat
shapiro.test(ANOVA_abundancia$res)#no normalitat
logANOVA_abundancia<-lm(log(dades_grup_counts_estructura_dcast_final$abundancia)~dades_grup_counts_estructura_dcast_final$tractament*dades_grup_counts_estructura_dcast_final$data)
leveneTest(logANOVA_abundancia)
shapiro.test(logANOVA_abundancia$res)
anova(logANOVA_abundancia)
library(car)
#GRÀFIC COMPARANT abundància PER tractament
summary(dades_grup_counts_estructura_dcast_final)
plot(dades_grup_counts_estructura_dcast_final$tractament,dades_grup_counts_estructura_dcast_final$abundancia, ylab="Abundància",xlab = "Tractament",cex.lab=1.5, col=c("skyblue1","seagreen1", "lightcoral","plum1"))
#Tenim 4 tractaments diferents, quin és diferent de qun?
#grafic comparant abundància per data
plot(dades_grup_counts_estructura_dcast_final$data,dades_grup_counts_estructura_dcast_final$abundancia, ylab="Abundància",xlab="Data",main="Diagrama de caixes de l'abundància respecte la data de mostreig",cex.main=2.2,cex.lab=1.55,cex.axis=1.25, col=c("paleturquoise1","tan1"))
boxplot(dades_grup_counts_estructura_dcast_final$abundancia~dades_grup_counts_estructura_dcast_final$tractament_dia,col=c("cornflowerblue","green2","cornflowerblue","orange","cornflowerblue","orangered","cornflowerblue","red3"), at=c(1,2,4,5,7,8,10,11),main="Boxplot of individuals abundance (initial and final) with respect to each treatment used",xlab = "Treatment and sampling day",ylab="Individuals abundance",cex.main=2.2,cex.lab=1.55,cex.axis=1.25)
boxplot(dades_grup_counts_estructura_dcast_final$abundancia~dades_grup_counts_estructura_dcast_final$data_tract, ylab="Abundància d'espècies",xlab="Control-nutrients segons la data de mostreig",main="Diagrama de caixes de l'abundància respecte la data de mostreig i presència / absència de nutrients",cex.main=2.2,cex.lab=1.55,cex.axis=1.25,col=c("paleturquoise1","green2","cornflowerblue","orangered1"),at=c(1,2,4,5))
par(mfrow=c(1,1))



#ANOVA Shannon
ANOVA_shannon<-lm(dades_grup_counts_estructura_dcast_final$Shannon~dades_grup_counts_estructura_dcast_final$tractament*dades_grup_counts_estructura_dcast_final$data)
leveneTest(ANOVA_shannon)#ok homogeneitat
shapiro.test(ANOVA_shannon$res)#ok normalitat
anova(ANOVA_shannon)
plot(dades_grup_counts_estructura_dcast_final$tractament,dades_grup_counts_estructura_dcast_final$Shannon, ylab="Shannon", xlab="Tractament",cex.lab=1.5, col=c("skyblue1","seagreen1", "lightcoral","plum1"))
#Tenim 4 tractaments diferents, quin és diferent de qun?
#Grafic comparant shannon per data
plot(dades_grup_counts_estructura_dcast_final$data,dades_grup_counts_estructura_dcast_final$Shannon, ylab="Shannon",xlab="Data",main="Diagrama de caixes de la diversitat específica respecte la data de mostreig",cex.main=2.2,cex.lab=1.55,cex.axis=1.25, col=c("paleturquoise1","tan1"))
#Eliminem la interacció tractament*data
boxplot(dades_grup_counts_estructura_dcast_final$Shannon~dades_grup_counts_estructura_dcast_final$tractament_dia,col=c("cornflowerblue","green2","cornflowerblue","orange","cornflowerblue","orangered","cornflowerblue","red3"), at=c(1,2,4,5,7,8,10,11),main="Boxplot of taxonomic diversity (initial and final) with respect to each treatment used",xlab = "Treatment and sampling day",ylab="Taxonomic diversity (according to the Shannon-Weaver index)",cex.main=2.2,cex.lab=1.55,cex.axis=1.25)
boxplot(dades_grup_counts_estructura_dcast_final$Shannon~dades_grup_counts_estructura_dcast_final$data_tract, ylab="Diversitat taxonòmica (segons l'índex Shannon-Weaver)",xlab="Control-nutrients segons la data de mostreig",main="Diagrama de caixes de la diversitat respecte la data de mostreig i presència / absència de nutrients",cex.main=2.2,cex.lab=1.55,cex.axis=1.25,col=c("paleturquoise1","green2","cornflowerblue","orangered1"),at=c(1,2,4,5))

#ANOVA  de la riquesa
ANOVA_riquesa<-lm(dades_grup_counts_estructura_dcast_final$Riquesa~dades_grup_counts_estructura_dcast_final$tractament*dades_grup_counts_estructura_dcast_final$data)
library(car)
leveneTest(ANOVA_riquesa)#ok homogeneïtat
shapiro.test(ANOVA_riquesa$res)#ok normalitat
logANOVA_riquesa<-lm(log(dades_grup_counts_estructura_dcast_final$Riquesa)~dades_grup_counts_estructura_dcast_final$tractament*dades_grup_counts_estructura_dcast_final$data)
leveneTest(logANOVA_riquesa)
shapiro.test(logANOVA_riquesa$res)
anova(logANOVA_riquesa)
#Totes les interaccións són significatives (tractament, data i tractament:data), tractament és molt significatiu
library(effects)
###############################################################################################################################
##########No em surt el gràfic, haig de mirar el perquè
plot(allEffects(ANOVA_riquesa), ask=FALSE)#FEM EL GRÀFIC DE LA INTERACCIÓ
#GRÀFIC COMPARANT abundància PER tractament
summary(dades_grup_counts_estructura_dcast_final)
plot(dades_grup_counts_estructura_dcast_final$tractament,dades_grup_counts_estructura_dcast_final$Riquesa, ylab="Riquesa",xlab = "Tractament",cex.lab=1.5, col=c("skyblue1","seagreen1", "lightcoral","plum1"))
#Tenim 4 tractaments diferents, quin és diferent de qun?
#grafic comparant abundància per data
plot(dades_grup_counts_estructura_dcast_final$data,dades_grup_counts_estructura_dcast_final$Riquesa, ylab="Riquesa",xlab="Data",main="Diagrama de caixes de l'abundància respecte la data de mostreig",cex.main=2.2,cex.lab=1.55,cex.axis=1.25, col=c("paleturquoise1","tan1"))
boxplot(dades_grup_counts_estructura_dcast_final$Riquesa~dades_grup_counts_estructura_dcast_final$tractament_dia,col=c("cornflowerblue","green2","cornflowerblue","orange","cornflowerblue","orangered","cornflowerblue","red3"), at=c(1,2,4,5,7,8,10,11),main="Boxplot of taxonomic richness (initial and final) with respect to each treatment used",xlab = "Treatment and sampling day",ylab="Taxonomic richness",cex.main=2.2,cex.lab=1.55,cex.axis=1.25)
boxplot(dades_grup_counts_estructura_dcast_final$Riquesa~dades_grup_counts_estructura_dcast_final$data_tract, ylab="Riquesa d'espècies",xlab="Control-nutrients segons la data de mostreig",main="Diagrama de caixes de la riquesa respecte la data de mostreig i presència / absència de nutrients",cex.main=2.2,cex.lab=1.55,cex.axis=1.25,col=c("paleturquoise1","green2","cornflowerblue","orangered1"),at=c(1,2,4,5))



dades_grup_counts_estructura_dcast_final_i<-dades_grup_counts_estructura_dcast_final
dades_grup_counts_estructura_dcast_final_i<-dades_grup_counts_estructura_dcast_final_i[dades_grup_counts_estructura_dcast_final_i$data=="02/03/2020",]
dades_grup_counts_estructura_dcast_final_f<-dades_grup_counts_estructura_dcast_final
dades_grup_counts_estructura_dcast_final_f<-dades_grup_counts_estructura_dcast_final_f[dades_grup_counts_estructura_dcast_final_f$data=="29/03/2020",]

summary(dades_grup_counts_estructura_dcast_final_i)
summary(dades_grup_counts_estructura_dcast_final_f)

dades_grup_counts_estructura_dcast_final_i_O<-dades_grup_counts_estructura_dcast_final_i[dades_grup_counts_estructura_dcast_final_i$tractament=="O",]
dades_grup_counts_estructura_dcast_final_f_O<-dades_grup_counts_estructura_dcast_final_f[dades_grup_counts_estructura_dcast_final_f$tractament=="O",]

summary(dades_grup_counts_estructura_dcast_final_i_O)
summary(dades_grup_counts_estructura_dcast_final_f_O)

dades_grup_counts_estructura_dcast_final_i_C<-dades_grup_counts_estructura_dcast_final_i[dades_grup_counts_estructura_dcast_final_i$tractament=="C",]
dades_grup_counts_estructura_dcast_final_f_C<-dades_grup_counts_estructura_dcast_final_f[dades_grup_counts_estructura_dcast_final_f$tractament=="C",]

summary(dades_grup_counts_estructura_dcast_final_i_C)
summary(dades_grup_counts_estructura_dcast_final_f_C)

dades_grup_counts_estructura_dcast_final_i_P<-dades_grup_counts_estructura_dcast_final_i[dades_grup_counts_estructura_dcast_final_i$tractament=="P",]
dades_grup_counts_estructura_dcast_final_f_P<-dades_grup_counts_estructura_dcast_final_f[dades_grup_counts_estructura_dcast_final_f$tractament=="P",]

summary(dades_grup_counts_estructura_dcast_final_i_P)
summary(dades_grup_counts_estructura_dcast_final_f_P)

dades_grup_counts_estructura_dcast_final_i_PC<-dades_grup_counts_estructura_dcast_final_i[dades_grup_counts_estructura_dcast_final_i$tractament=="PC",]
dades_grup_counts_estructura_dcast_final_f_PC<-dades_grup_counts_estructura_dcast_final_f[dades_grup_counts_estructura_dcast_final_f$tractament=="PC",]

summary(dades_grup_counts_estructura_dcast_final_i_PC)
summary(dades_grup_counts_estructura_dcast_final_f_PC)
###############################################################################################################################
##########Material i métodes
##########Resultats #estudi site, disseny experimental
#Introducció: parlar de pertorbacions pulse i continous, cites bibliogràfiques, etc.
