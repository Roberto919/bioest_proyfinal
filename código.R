# Cargamos las librerías
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

# leer los datos
df<-read_csv('DIG.csv')
colnames(df)<-tolower(colnames(df))

control<-c( 
  "id","trtmt","age","race","sex","ejf_per","ejfmeth","chestx","bmi","klevel","creat","digdoser","chfdur","rales",  
  "elevjvp","pedema","restdys","exertdys","actlimit","s3","pulcong","nsym","heartrte","diabp","sysbp","functcls","chfetiol","prevmi", 
  "angina","diabetes","hyperten","diguse","diuretk","diuret","ksupp","aceinhib","nitrates","hydral","vasod","digdose")

df_control<-df[control]

df_death<-df[c("id","trtmt","death","deathday","reason")]
df_death$reason[is.na(df$reason)]<-0

df_hosp<-df[c("id","trtmt","dwhf","dwhfdays")]


# ----------------- A N Á L I S I S     E X P L O R A T O R I O ----------------
#                              D E     D A T O S

# Distribuciones de la edad
vals <- df %>% 
  select(age) %>% 
  summarise(
    mean = mean(age),
    median = median(age)
  )

vals

df %>% 
  ggplot(
    aes(x = age)
  ) + 
  geom_histogram(fill = "#C3D7A4") +
  theme_bw() +
  geom_vline(xintercept = vals$mean, color = "red") + geom_vline(xintercept = vals$median, color = "blue") +
  geom_text(
    aes(x = round(vals$mean, 2), label = paste0("Media = ", round(vals$mean, 2)), y = 650),
    colour="red", angle=90, vjust = -1
  ) +
  geom_text(
    aes(x = round(vals$median, 2), label = paste0("\nMediana = ", round(vals$median, 2)), y = 650),
    colour="blue", angle=90
  ) +
  xlab("Edad") + ylab("Conteo") + ggtitle("Distribución de la edad")


### Conteos de sexo segmentados por raza
df %>% 
  select(c("race", "sex")) %>% 
  mutate(
    race = ifelse(race == 1, "Blanco/a", "Otra"),
    sex = ifelse(sex == 1, "Hombre", "Mujer"),
  ) %>% 
  group_by(sex, race) %>%
  summarise(Conteo = n()) %>% 
  ggplot(
    aes(fill = race, y = Conteo, x = sex)
  ) + 
  theme_bw() +
  geom_bar(position = "stack", stat = "identity") +
  xlab(" ") + ggtitle("Participantes segmentados por sexo y raza") + labs(fill = "Raza")



### Enfermedades crónicas
# Infarto y angina de pecho
plt_ma <- df %>%
  select(c("prevmi", "angina")) %>%
  filter(
    !is.na(prevmi),
    !is.na(angina),
  ) %>%
  mutate(
    prevmi = ifelse(prevmi == 0, "Sin infarto", "Con infarto"),
    angina = ifelse(angina == 0, "Sin angina", "Con angina")
  ) %>%
  group_by(prevmi, angina) %>%
  summarise(n = n()) %>%
  ggplot(
    aes(x = prevmi, y = angina, fill = n)
  ) +  
  theme_minimal() +
  theme(legend.position = "none") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = round(n, 2))) +
  xlab("") + ylab("") + labs(fill = "Conteo")


## Diabetes e hipertensión
plt_dh <- df %>% 
  select(c("diabetes", "hyperten")) %>% 
  filter(
    !is.na(diabetes),
    !is.na(hyperten),
  ) %>% 
  mutate(
    dia = ifelse(diabetes == 0, "Sin diabetes", "Con diabetes"),
    hyp = ifelse(hyperten == 0, "Sin hipertensión", "Con hipertensión")
  ) %>% 
  group_by(dia, hyp) %>% 
  summarise(Conteo = n()) %>% 
  ggplot(
    aes(x = dia, y = hyp, fill = Conteo)
  ) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = round(Conteo, 2))) + 
  xlab("") + ylab("") + labs(fill = "Conteo")


## Juntando las dos gráficas
ggarrange(
  plt_dh, plt_ma,
  nrow = 1
) %>% 
  annotate_figure(
    top = text_grob("Enfermedades crónicas de los participantes previo al estudio", size = 16)
  )


### MBI
df %>% 
  select("bmi") %>% 
  ggplot(
    aes(x = bmi)
  ) +
  geom_boxplot(fill = "#56B4E9") +
  theme_bw() +
  theme(
    axis.title.y=element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
  ) +
  xlab(" ") + ggtitle("Distribución del índice de masa corporal (BMI)")

# ¿Qué porcentaje de la gente que tomo el medicamento se murió?
df_death %>% 
  group_by(trtmt) %>% 
  summarise(
    total = n(),
    deaths = sum(death)
  ) %>% 
  mutate(
    prop = round(deaths/total, 3),
    trtmt = ifelse(trtmt == 0, "Sin Tratamiento", "Con Tratamiento")
  ) %>% 
  ggplot(
    aes(x = trtmt, y = prop, fill = trtmt)
  ) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) + 
  ylab("Proporción de muertes") + xlab("") + labs(fill = "") +
  ggtitle("Proporción de muertes diferenciando por tratamiento")


# ---------------------------- E S T I M A D O R -------------------------------
#                         K A P L A N   -   M E I E R

### Muertes
ggsurvplot(
  fit = survfit(Surv(deathday, death) ~ trtmt, data = df_death,conf.type = "log-log"), 
  xlab = "Días", 
  ylab = "Probabilidad de Supervivencia",
  title = "Efecto de la Digoxina en Mortalidad",
  legend.labs = c("Control",
                  "Tratamiento"),
  conf.int = TRUE,
  pval = TRUE,
  size =.0,
  linetype = "strata",
  times=seq(0,1781+30,120),
  xlim=c(0,1800),
  censor.shape = "|",
  censor.size=2.5)
mort<-survdiff(Surv(deathday, death) ~ trtmt,rho=0, data = df_death)

mortfit<-survfit(Surv(deathday, death) ~ trtmt, data = df_death,conf.type = "log-log") 
quantile(mortfit,c(.15,.25,.35,.45,.50,.55))

summary(mortfit)
ci<-cmprsk::cuminc(
  ftime = df_death$deathday, 
  fstatus = df_death$reason, 
  group = df_death$trtmt,
  #strata = df_death$trtmt[df_death$death==1],
  rho = 1,
  cencode = 0,
)

tp<-timepoints(ci,seq(1,1781,1))$est

meantp<-as.vector(rowSums(1-tp,na.rm=TRUE))
sdtp<-c()
for(i in c(1:nrow(tp))){
  sds<-c()
  for (j in c(1:1000)){
    s<-sample(1-tp[i,],replace=TRUE)
    sds<-c(sds,sum(s,na.rm = TRUE))
    
  }
  
  sdtp<-c(sdtp,sd(sds,na.rm = TRUE))
}

ciinf<-meantp-1.96*sdtp
cisup<-meantp+1.96*sdtp




ggcompetingrisks(ci,
                 multiple_panels = TRUE,
                 xlab = "Días",
                 ylab = "Incidencia Acumulada",
                 title = "Efecto de la Digoxina por Tipo de Muerte",
                 
                 legend.title = "Evento",
                 conf.int = TRUE,
                 xlim=c(0,1800)
                 
)

ci$Tests


print(xtable(tibble(x$time,x$n.risk,x$n.event,x$surv,x$std.err,x$lower,x$upper)),include.rownames = FALSE)

# --------------------- E S T I M A C I Ó N     P O R --------------------------
#                M Á X I M A     V E R O S I M I L I T U D

df_death<-df[c("id","trtmt","death","deathday","reason", "diabetes","hyperten", "angina", "prevmi", "sex", "bmi")]
df_death$reason[is.na(df$reason)]<-0

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

# coeficiente de información de Akaike
extractAIC(regw)

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



# --------------------- M O D E L O S     D E    R E G R E S I Ó N -------------

#Modelo 1

df$reason[is.na(df$reason)]<-0
df$race[df$race==1]<-0
df$race[df$race==2]<-1

xfit1<-coxph(Surv(dwhfdays,dwhf)~trtmt+age+race+sex+bmi+hyperten+diabetes+angina+prevmi,data=df)


xfit2<-coxph(Surv(deathday,death)~trtmt+age+race+sex+bmi+hyperten+diabetes+angina+prevmi,data=df)


xfit3<-coxph(Surv(death,as.factor(reason),type="mstate")~trtmt+age+race+sex+bmi+hyperten+diabetes+angina+prevmi,data=df,id=id)


xfit3<-coxph(Surv(death,as.factor(reason),type="mstate")~trtmt+age+race+sex+bmi+hyperten+diabetes+angina+prevmi,data=df,id=id)
xfit4<-coxph(Surv(death,as.factor(reason),type="mstate")~trtmt++age+race+sex+bmi+hyperten+diabetes+angina+prevmi+diabetes*prevmi*angina*hyperten*bmi,data=df,id=id)      




ggsurvplot(
  fit = survfit(Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp,conf.type = "log-log"), 
  xlab = "Días", 
  ylab = "Probabilidad de Supervivencia",
  title = "Efecto de la Digoxina en Hospitalización por Falla Cardiaca",
  subtitles="por Falla Cardiaca",
  legend.labs = c("Control",
                  "Tratamiento"),
  conf.int = TRUE,
  pval = TRUE,
  size =.5,
  linetype = "strata",
  times=seq(0,1781+30,120),
  xlim=c(0,1800))

hosp<-survdiff(Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp) 

hospfit<-survfit(Surv(dwhfdays, dwhf) ~ trtmt, data = df_hosp,conf.type = "log-log") 
quantile(hospfit,c(.15,.25,.35,.45,.50,.55))


x<-summary(hospfit,times = seq(0,1781+30,120))
print(xtable(tibble(x$time,x$n.risk,x$n.event,x$surv,x$std.err,x$lower,x$upper)),include.rownames = FALSE)
