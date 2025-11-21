#####################################################################################
# Primer Modelo (FLOW)
#####################################################################################

datos <- Conjunto_2[,-1]
str(datos)
library(dplyr)
me<- rename(datos, "Y" = "FLOW(cm)")

# Se usara los datos me
m1 <- lm(Y~Cement + Slag +`Fly ash` +Water +SP +`Coarse Aggr.`+`Fine Aggr.`,data=me[,-c(8,10)])
# -------------------------------
summary(m1)
shapiro.test(m1$residuals)
library(lmtest)
bptest(m1)
# ------------------------------

#########################################
# Medir la multicolinealidad #########
######################################
library(car)
vif(m1)
m1x6 <- lm(Y~Cement + Slag +`Fly ash` +Water +SP +`Fine Aggr.`,data=me[,-c(6,8,10)]) #Se usara ese modelo de ahora en adelante
vif(m1x6) # No hay multicolinealidad con: Cement Slag `Fly ash` Water SP `Fine Aggr.`

###################################################################
# Transformación                                               ####
###################################################################

# Como no se cumple la normalidad en el modelo m1x6f, se realizara una transformación
library(car)
box<-boxCox(m1x6, lambda = seq(-2, 2, by = 0.1)) # Prueba de valores
lambda<-box$x[which.max(box$y)] # Valor de λ optimo

ybox<-(me$Y^lambda-1)/lambda 

m1x6 <- lm(ybox~ Cement + Slag +`Fly ash` +Water +SP +`Fine Aggr.`,data=me[,-c(6,8,9,10)])
shapiro.test(m1x6$residuals) # 0.09152
library(lmtest)  
bptest(m1x6) # 0.1987
dwtest(m1x6, alternative="t") # 0.8225
########################################
# Seleccionamos Variables ##############
########################################
  ###################################
  #Método de mejores subconjuntos####
  ###################################
  
  library(olsrr)
  res<-ols_step_all_possible(m1x6)
  #Con el RMSE
  which.min(res$result[,6])
  res$result[63,] # Cement + Slag + `Fly ash` + Water + SP + `Fine Aggr.`
  
  #Con el AIC 
  which.min(res$result[,9])
  res$result[7,] # Slag + Water
  
  #Con el R2 pred 
  which.max(res$result[,7])
  res$result[7,]# Slag + Water
  
  
  ################################
  #Método por pasos (Stepwise)####
  ################################
  
  #Backward
  step(m1x6,trace=T,direction="backward") # Slag + Water
  
  #Forward
  horizonte <- formula(ybox ~ Cement + Slag + `Fly ash` + Water + SP + `Fine Aggr.`)
  modelo0 <- lm(ybox ~ 1, data = me[,-c(6,8,9,10)])
  step(modelo0, trace = T, direction = "forward", scope = horizonte) # Slag + Water
  
  #Both
  step(modelo0, trace=T, direction="both", scope=horizonte) # Slag + Water
  
  # Noa quedamos con el modelo Slag + Water
  m1x6f <- lm(ybox~ Slag + Water,data=me[,-c(6,8,9,10)])

########################################################
# Analisamos las observaciones influyentes y outlier####

# Outliers 
# Residuos Estandarizados
ri <-abs(rstandard(m1x6f)) > 2
ri[ri==TRUE] # Son valores atipicos y pueden ser puntos influyentes:  31, 66 
# R. Estunderizados
ti<- abs(rstudent(m1x6f)) >2
ti[ti == TRUE]  # Son valores atipicos y pueden ser puntos influyentes : 31, 66
# Bonferroni
library(car)
outlierTest(m1x6f, cutoff=Inf, n.max=5, order=TRUE)
# con un alfa de 0.05 no hay observaciones que se puedan considerar outlier estadisticamente

# Influyentes
# Matriz H ####
hii<-hatvalues(m1x6f)
sum(hii)
2*5/100 > 1
r<- abs(hii) >2*7/103
r[r==TRUE]
# Distancia De Cook ####
cooks.distance(m1x6f)
di<- cooks.distance(m1x6f) > 1
di[di==TRUE]
Di <- cooks.distance(m1x6f) > qf(0.05, 5, 95)
Di[Di==TRUE]
# Puntos influyentes: 31 y 51 para el B globales

#DFBETAS ####
DFB<-dfbeta(m1x6f)
DFBE<-abs(DFB[,-1]) >2/sqrt(100)
which(apply(DFBE,1,sum)==2) # No hay puntos influyentes para los Bj por individual

# DFFITS ####
DFFI <-abs(dffits(m1x6f)) > 2*sqrt(5/100)
DFFI[DFFI==TRUE] # Las observaciones son puntos influyentes en la prediccion del modelo, 

# Covratio ####
c1<- covratio(m1x6f) > 1 + 3*5/100
c2<- covratio(m1x6f) < 1 - 3*5/100
c1[c1== TRUE] # puntos influyentes en la presicion del modelo
c2[c2== TRUE] # puntos influyentes en la presicion del modelo


###################################
# Supuestos                    ####
###################################
m1f <- lm(ybox~ Slag + Water,data=me[,-c(6,8,9,10)])
shapiro.test(m1f$residuals) # 0.09768
bptest(m1f) # 0.4088
dwtest(m1f, alternative = "t") # 0.9722
summary(m1f)
###################################################################
# Modelo  final                                      ####
###################################################################

m1f <- lm(ybox~  Slag + Water,data=me[,-c(6,8,9,10)])
###################################################################
# VALIDACIÓN DEL MODELO FLOW (Modelo final: Slag + Water) ####
###################################################################

library(car)
library(lmtest)
library(olsrr)
library(dplyr)

# Renombrar variable respuesta
datos <- Conjunto_2[,-1]
datos <- rename(datos, Y = "FLOW(cm)")

# Eliminar SLUMP y Strength (columnas 8 y 10)
datos_modelo <- datos[,-c(8,10)] 
# Columnas disponibles:
# Cement, Slag, Fly ash, Water, SP, Coarse Aggr., Fine Aggr., Y

###################################################################
# 1. Ajustar modelo inicial para obtener lambda
###################################################################

set.seed(40)
n_total <- nrow(datos_modelo)
ne <- round(0.8 * n_total)
indices_inicial <- sample(n_total, ne)

me <- datos_modelo[indices_inicial, ]

# Modelo sin transformar (para Box-Cox)
m_box <- lm(Y ~ Slag + Water, data = me)

box <- boxCox(m_box, lambda = seq(-2, 2, by = 0.1))
lambda <- box$x[which.max(box$y)]

# Transformación Box–Cox
ybox <- (me$Y^lambda - 1)/lambda
m1f <- lm(ybox ~ Slag + Water, data = me)

###################################################################
# 2. Transformación inversa
###################################################################
transformacion_inversa <- function(y_trans, lambda){
  if (abs(lambda) > 1e-6){
    (lambda*y_trans + 1)^(1/lambda)
  } else {
    exp(y_trans)
  }
}

###################################################################
# 3. Validación Hold-Out Repetida (1000 repeticiones)
###################################################################

set.seed(40)
R <- 1000
MSE_vec <- numeric(R)
MAE_vec <- numeric(R)

for(i in 1:R){
  
  # Nueva partición
  idx <- sample(n_total, ne)
  me_i <- datos_modelo[idx, ]
  mp_i <- datos_modelo[-idx, ]
  
  # Transformación
  y_trans_i <- (me_i$Y^lambda - 1)/lambda
  me_i$ybox <- y_trans_i
  
  # Ajustar modelo final
  mod_i <- lm(ybox ~ Slag + Water, data = me_i)
  
  # Predicción en escala transformada
  pred_trans <- predict(mod_i, newdata = mp_i)
  
  # Volver a escala original
  pred_original <- transformacion_inversa(pred_trans, lambda)
  real_original <- mp_i$Y
  
  errores <- pred_original - real_original
  MSE_vec[i] <- mean(errores^2)
  MAE_vec[i] <- mean(abs(errores))
}

RMSE_final <- sqrt(mean(MSE_vec))
MAE_final  <- mean(MAE_vec)

cat("=== VALIDACIÓN DEL MODELO FLOW (1000 repeticiones) ===\n")
cat("RMSE :", round(RMSE_final, 3), "cm\n")
cat("MAE  :", round(MAE_final, 3), "cm\n")

###################################################################
# 4. Indicadores adicionales del modelo ajustado
###################################################################

AIC_valor <- AIC(m1f)
BIC_valor <- BIC(m1f)
R2_ajustado <- summary(m1f)$adj.r.squared

# R² de predicción (PRESS)
R2_pred <- ols_pred_rsq(m1f)

# Cp de Mallows
ols_all <- ols_step_all_possible(m1f)
Cp_mejor <- min(ols_all$result[,"Mallow's Cp"])

cat("\n=== INDICADORES ADICIONALES ===\n")
cat("AIC            :", round(AIC_valor, 2), "\n")
cat("BIC            :", round(BIC_valor, 2), "\n")
cat("R² ajustado    :", round(R2_ajustado, 3), "\n")
cat("R² predicción  :", round(R2_pred, 3), "\n")
cat("Cp de Mallows  :", round(Cp_mejor, 3), "\n")
