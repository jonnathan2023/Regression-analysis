#####################################################################################
# Tercer Modelo v2
#####################################################################################

datos <- Conjunto_2[,-1]
str(datos)
library(dplyr)
me<- rename(datos, "Y" ="Compressive Strength (28-day)(Mpa)")

# Se usara los datos me
m <- lm(Y~.,data=me[,-c(8,9)])
summary(modelo)
shapiro.test(m$residuals)
library(lmtest)
bptest(m)
###################################################################
# Transformación                                               ####
###################################################################
#--------------------------------------------------------------
# Renombrar variables con espacios o puntos 
names(me) <- gsub(" ", "_", names(me))        # Reemplaza espacios por guiones bajos
names(me) <- gsub("\\.", "", names(me))       # Elimina puntos (opcional)

# O manualmente:
names(me)[names(me) == "Fly ash"] <- "Fly_ash"
names(me)[names(me) == "Fine Aggr."] <- "Fine_Aggr"
# Ejecutar
residualPlots(m)
# Regresar a los nombres originales a las variables
#---------------------------------------------------------------------------- 
m3_1 <- lm(Y ~ Cement + Slag +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`+ `Coarse Aggr.`, data = me[,-c(8,9)])
vif(m3_1)
m3 <- lm(Y ~ Cement + Slag +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`, data = me[,-c(8,9)])

# Como no se cumple la normalidad en el modelo m3, se realizara una transformación
library(car)
box<-boxCox(m3, lambda = seq(-2, 2, by = 0.1)) # Prueba de valores
lambda<-box$x[which.max(box$y)] # Valor de λ optimo

ybox<-(me$Y^lambda-1)/lambda 
m3b <- lm(ybox ~ Cement + Slag +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`, data = me[,-c(8,9,10)])
########################################################
# Analisamos las observaciones influyentes y outlier####
########################################################
## MODIFICAR m2b POR m3b
# Outliers 
# Residuos Estandarizados
ri <-abs(rstandard(m3b)) > 2 
ri[ri==TRUE] # Son valores atipicos y pueden ser puntos influyentes:  33, 36
# R. Estunderizados
ti<- abs(rstudent(m3b)) >2
ti[ti == TRUE]  # Son valores atipicos y pueden ser puntos influyentes : 33, 66
# Bonferroni
library(car)
outlierTest(m3b, cutoff=Inf, n.max=5, order=TRUE)

# Influyentes
# Matriz H ####
hii<-hatvalues(m3b)
sum(hii)
2*7/103 > 1
r<- abs(hii) >2*7/103
r[r==TRUE]
# Distancia De Cook ####
cooks.distance(m3b)
di<- cooks.distance(m3b) > 1
di[di==TRUE]
Di <- cooks.distance(m3b) > qf(0.05, 5, 95)
Di[Di==TRUE]
# Puntos influyentes: 31 y 51 para el B globales

#DFBETAS ####
DFB<-dfbeta(m3b)
DFBE<-abs(DFB[,-1]) >2/sqrt(100)
which(apply(DFBE,1,sum)==2) # No hay puntos influyentes para los Bj por individual

# DFFITS ####
DFFI <-abs(dffits(m3b)) > 2*sqrt(5/100)
DFFI[DFFI==TRUE] # Las observaciones son puntos influyentes en la prediccion del modelo, 

# Covratio ####
c1<- covratio(m3b) > 1 + 3*5/100
c2<- covratio(m3b) < 1 - 3*5/100
c1[c1== TRUE] # puntos influyentes en la presicion del modelo
c2[c2== TRUE] # puntos influyentes en la presicion del modelo

########################################
# Seleccionamos Variables ##############
########################################

  ################################
  #Método por pasos (Stepwise)####
  ################################
  
  #Backward
  step(m3b,trace=T,direction="backward") # Cement + Slag + `Fly ash` + + I(Water^2) +SP + `Fine Aggr.`
  
  #Forward
  horizonte <- formula(ybox ~ Cement + Slag +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`)
  modelo0 <- lm(ybox ~ 1,data=me[,-c(6,8,9,10)])
  step(modelo0, trace = T, direction = "forward", scope = horizonte) # Cement + Slag + `Fly ash` + + I(Water^2) +SP + `Fine Aggr.`
  
  #Both
  step(modelo0, trace=T, direction="both", scope=horizonte) # Cement + Slag + `Fly ash` + + I(Water^2) +SP + `Fine Aggr.`
  
  # Como resultado tenemos tres posibles
  m3_1 <- lm(ybox ~ Cement + Slag + `Fly ash` + + I(Water^2) +SP + `Fine Aggr.`,data=me[,-c(6,8,9,10)])

#########################################
# Medir la multicolinealidad #########
######################################
library(car)
# Modelo m2_1
vif(m3_1) # No hay multicolinealidad alta
m3f <- lm(ybox~ Cement + Slag + `Fly ash` + + I(Water^2) +SP + `Fine Aggr.`,data=me[,-c(6,8,9,10)])
################
# Supuestos ####
################
shapiro.test(m3f$residuals) # 0.165
bptest(m3f) # 0.2347
summary(m3f)
_##########################################################################
# VALIDACIÓN CRUZADA REALISTA (1000 repeticiones) PARA EL MODELO m3f
##########################################################################

library(car)
library(lmtest)
library(olsrr)
library(dplyr)

#-------------------------------
# Preparar los datos
#-------------------------------
datos <- Conjunto_2[,-1]
datos <- rename(datos, Y = "Compressive Strength (28-day)(Mpa)")

# Usar mismo subconjunto que tu modelo
datos_modelo <- datos[,-c(8,9)]   
# Columnas disponibles:
# Cement, Slag, Fly ash, Water, SP, Coarse Aggr., Fine Aggr., Y

#-------------------------------
# Calcular lambda del Box–Cox
#-------------------------------
m_box <- lm(Y ~ Cement + Slag + `Fly ash` + Water + SP + `Fine Aggr.`,
            data = datos_modelo)

box <- boxCox(m_box, lambda = seq(-2, 2, by = 0.1))
lambda <- box$x[which.max(box$y)]

#-------------------------------
# Transformación inversa
#-------------------------------
inv_boxcox <- function(y, lambda){
  if(abs(lambda) > 1e-6){
    (lambda*y + 1)^(1/lambda)
  } else {
    exp(y)
  }
}

#-------------------------------
# VALIDACIÓN CRUZADA REALISTA
#-------------------------------
set.seed(40)
R <- 1000
n <- nrow(datos_modelo)
ne <- round(0.8*n)

MSE_vec <- numeric(R)
MAE_vec <- numeric(R)

for(i in 1:R){
  
  # Partición 80/20
  idx <- sample(n, ne)
  me_i <- datos_modelo[idx, ]
  mp_i <- datos_modelo[-idx, ]
  
  # Transformación Box–Cox
  me_i$ybox <- (me_i$Y^lambda - 1)/lambda
  
  # Ajustar modelo final m3f
  mod_i <- lm(
    ybox ~ Cement + Slag + `Fly ash` + I(Water^2) + SP + `Fine Aggr.`,
    data = me_i
  )
  
  # Predicción en escala transformada
  pred_trans <- predict(mod_i, newdata = mp_i)
  
  # Predicciones en escala original
  pred_original <- inv_boxcox(pred_trans, lambda)
  real_original <- mp_i$Y
  
  # Errores
  err <- pred_original - real_original
  
  MSE_vec[i] <- mean(err^2)
  MAE_vec[i] <- mean(abs(err))
}

RMSE_final <- sqrt(mean(MSE_vec))
MAE_final  <- mean(MAE_vec)

cat("==============================================\n")
cat(" VALIDACIÓN CRUZADA REALISTA (1000 repeticiones)\n")
cat("==============================================\n")
cat("RMSE :", round(RMSE_final, 3), "\n")
cat("MAE  :", round(MAE_final, 3), "\n")

#-------------------------------
# Indicadores adicionales
#-------------------------------

# Ajustar modelo final con TODOS los datos
datos_modelo$ybox <- (datos_modelo$Y^lambda - 1)/lambda

m3f <- lm(
  ybox ~ Cement + Slag + `Fly ash` + I(Water^2) + SP + `Fine Aggr.`,
  data = datos_modelo
)

AIC_valor <- AIC(m3f)
BIC_valor <- BIC(m3f)
R2_adj <- summary(m3f)$adj.r.squared
R2_pred <- ols_pred_rsq(m3f)  # PRESS

cat("\n========== INDICADORES ADICIONALES ==========\n")
cat("AIC           :", round(AIC_valor, 3), "\n")
cat("BIC           :", round(BIC_valor, 3), "\n")
cat("R² ajustado   :", round(R2_adj, 3), "\n")
cat("R² predicción :", round(R2_pred, 3), "\n")


