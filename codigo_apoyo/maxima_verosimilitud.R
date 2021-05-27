library(readr)
library(stringr)
library(tidyverse)
library(survival)
library(KMsurv)
library(actuar)
library(BGPhazard)
library(rmutil)
library(survminer)
library(fitdistrplus)
library(logspline)

df<-read_csv('DIG.csv')
colnames(df)<-tolower(colnames(df))

df_death<-df[c("id","trtmt","death","deathday","reason", "diabetes","hyperten", "angina", "prevmi", "sex", "bmi")]
df_death$reason[is.na(df$reason)]<-0

df_hosp<-df[c("id","trtmt","dwhf","dwhfdays")]


plot.new()
descdist(df_death$deathday,boot=10000, discrete = FALSE)

fitw <- fitdist(df_death$deathday, "weibull")
fitu <- fitdist(df_death$deathday, "norm")
fitln <- fitdist(df_death$deathday, "lnorm")
fitg <- fitdist(df_death$deathday, "gamma")

plot(fitw)
plot(fitu)
plot(fitln)
plot(fitg)

# hacer el plot de todos juntos
par(mfrow=c(2,2))
plot.legend <- c("weibull", "norm", "lnorm", "gamma")
denscomp(list(fitw, fitu, fitln, fitg), addlegend = FALSE)
qqcomp(list(fitw, fitu, fitln, fitg), legendtext = plot.legend)
cdfcomp(list(fitw, fitu, fitln, fitg), addlegend = FALSE)
ppcomp(list(fitw, fitu, fitln, fitg), addlegend = FALSE)

# con base en las gráficas anteriores, vemos que la distribución que más se 
# acerca a nuestros datos es la weibull


# ---------------------------- M U E R T E S -----------------------------------

par(mfrow=c(1,2), oma=c(0,0,1.5,0)) 


#en general
regw <- survreg( Surv(deathday, death) ~ 1, data = df_death, dist="weibull")
fitw <- survfit(Surv(deathday, death) ~ 1, data = df_death)
plot(fitw, xlab = "Días desde inicio del estudio", ylab = "Prob. de supervivencia", conf.int = TRUE)

pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(regw,  type='quantile',
                 p=pct, se=TRUE)
lines(x=ptime$fit[1,], y=1-pct, col="purple")


# shape es el recíproco de scale
shape <- 1/regw$scale #a
# scale es el exponencial del intercepto
scale <- exp( coef(r) ) #b


# separado por treatment
regw <- survreg( Surv(deathday, death) ~ trtmt, data = df_death, dist="weibull")
regwp <- survreg( Surv(deathday, death) ~ trtmt, data = subset(df_death, trtmt == 0), dist="weibull")
regwt <- survreg( Surv(deathday, death) ~ trtmt, data = subset(df_death, trtmt == 1), dist="weibull")

fitw <- survfit(Surv(deathday, death) ~ trtmt, data = df_death)

plot(fitw, xlab = "Días desde inicio del estudio", ylab = "Prob. de supervivencia", conf.int = TRUE)
lines(predict(regw, newdata=list(trtmt=0),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red")
lines(predict(regw, newdata=list(trtmt=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green")

title("Efecto de la digoxina en la mortalidad", outer = TRUE)




# ---------------------- H O S P I T A L I Z A C I Ó N -------------------------

par(mfrow=c(1,2), oma=c(0,0,1.5,0)) 
#en general
regw <- survreg( Surv(dwhfdays, dwhf) ~ 1, data = df_hosp, dist="weibull")
regwp <- survreg( Surv(dwhfdays, dwhf) ~ 1, data = subset(df_hosp, trtmt == 0), dist="weibull")
regwt <- survreg( Surv(dwhfdays, dwhf) ~ 1, data = subset(df_hosp, trtmt == 1), dist="weibull")

fitw <- survfit(Surv(dwhfdays, dwhf) ~ 1, data = df_hosp)
plot(fitw, xlab = "Días desde inicio del estudio", ylab = "Prob. de supervivencia", conf.int = TRUE)

pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(regw,  type='quantile',
                 p=pct, se=TRUE)
lines(x=ptime$fit[1,], y=1-pct, col="purple")

# separado por treatment
regw <- survreg( Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp, dist="weibull")
fitw <- survfit(Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp)

plot(fitw, xlab = "Días desde inicio del estudio", ylab = "Prob. de supervivencia", conf.int = TRUE)
lines(predict(regw, newdata=list(trtmt=0),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red")
lines(predict(regw, newdata=list(trtmt=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green")

title("Efecto de la digoxina en la hospitalización por falla cardiaca", outer = TRUE)


# ----------------------- L O G - R A N K     T E S T --------------------------
#                              M U E R T E S 

# muertes normales
survdiff(Surv(deathday, death) ~ trtmt, data = df_death)

# muertes por falla cardiaca
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1))

# muertes por falla cardiaca con cormobilidades

survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & diabetes == 1))
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & hyperten == 1))
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & angina == 1))
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & prevmi == 1))
df_death$bmio <- ifelse(df_death$bmi < 25, 0, 1)
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & bmio == 1))

# muertes por falla cardiaca estratificadas por demográficos
survdiff(Surv(deathday, death) ~ trtmt +  strata(sex), data = subset(df_death, reason == 1))
# bmi > 25 es sobrepeso
df_death$bmio <- ifelse(df_death$bmi < 25, 0, 1)
survdiff(Surv(deathday, death) ~ trtmt +  strata(bmio), data = subset(df_death, reason == 1))
survdiff(Surv(deathday, death) ~ trtmt, data = subset(df_death, reason == 1 & bmio == 1))


# ----------------------- L O G - R A N K     T E S T --------------------------
#                      H O S P I T A L I Z A C I O N E S

# hospitalizaciones normales
# no hay

# hospitalizaciones por falla cardiaca
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp)

# hospitalizaciones por falla cardiaca que además tenían otros padecimientos
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = subset(df_hosp, diabetes == 1))
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = subset(df_hosp, hyperten == 1))
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = subset(df_hosp, angina == 1))
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = subset(df_hosp, prevmi == 1))
df_death$bmio <- ifelse(df_death$bmi < 25, 0, 1)
survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = subset(df_hosp, bmio == 1))
