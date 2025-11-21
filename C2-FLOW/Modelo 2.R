#####################################################################################
# Segundo Modelo (Flow)
#####################################################################################
datos <- Conjunto_2[,-1]
str(datos)
library(dplyr)
me<- rename(datos, "Y" = "FLOW(cm)") # Fila 9

# Se usara los datos me
modelo <- lm(Y~Cement + Slag +`Fly ash`+ Water + SP +`Coarse Aggr.` + `Fine Aggr.`,data=me[,-c(8,10)])
# -------------------------------
summary(modelo)
shapiro.test(modelo$residuals)
library(lmtest)
bptest(modelo)

#########################################
# Medir la multicolinealidad #########
######################################
library(car)
vif(modelo)
m <- lm(Y~ Cement + Slag +`Fly ash`+ Water + SP + `Fine Aggr.`,data=me[,-c(6,8,10)])
# -------------------------------
summary(m)
shapiro.test(m$residuals)
library(lmtest)
bptest(m)
# ------------------------------
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
m2 <- lm(Y ~ Cement + I(Slag^2) +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`, data = me[,-c(8,10)])
vif(m2)

# Como no se cumple la normalidad en el modelo m2, se realizara una transformación
library(car)
box<-boxCox(m2, lambda = seq(-2, 2, by = 0.1)) # Prueba de valores
lambda<-box$x[which.max(box$y)] # Valor de λ optimo

ybox<-(me$Y^lambda-1)/lambda 
m2b <- lm(ybox ~ Cement + I(Slag^2) +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`, data = me[,-c(6,8,9,10)])

########################################################
# Analisamos las observaciones influyentes y outlier####
########################################################
## MODIFICAR m2b POR m2
# Outliers 
# Residuos Estandarizados
ri <-abs(rstandard(m2b)) > 2 
ri[ri==TRUE] # Son valores atipicos y pueden ser puntos influyentes:  33, 36
# R. Estunderizados
ti<- abs(rstudent(m2b)) >2
ti[ti == TRUE]  # Son valores atipicos y pueden ser puntos influyentes : 33, 66
# Bonferroni
library(car)
outlierTest(m2b, cutoff=Inf, n.max=5, order=TRUE)

# Influyentes
# Matriz H ####
hii<-hatvalues(m2b)
sum(hii)
2*7/103 > 1
r<- abs(hii) >2*7/103
r[r==TRUE]
# Distancia De Cook ####
cooks.distance(m2b)
di<- cooks.distance(m2b) > 1
di[di==TRUE]
Di <- cooks.distance(m2b) > qf(0.05, 5, 95)
Di[Di==TRUE]
# Puntos influyentes: 31 y 51 para el B globales

#DFBETAS ####
DFB<-dfbeta(m2b)
DFBE<-abs(DFB[,-1]) >2/sqrt(100)
which(apply(DFBE,1,sum)==2) # No hay puntos influyentes para los Bj por individual

# DFFITS ####
DFFI <-abs(dffits(m2b)) > 2*sqrt(5/100)
DFFI[DFFI==TRUE] # Las observaciones son puntos influyentes en la prediccion del modelo, 

# Covratio ####
c1<- covratio(m2b) > 1 + 3*5/100
c2<- covratio(m2b) < 1 - 3*5/100
c1[c1== TRUE] # puntos influyentes en la presicion del modelo
c2[c2== TRUE] # puntos influyentes en la presicion del modelo

########################################
# Seleccionamos Variables ##############
########################################
  
  ################################
  #Método por pasos (Stepwise)####
  ################################
  
  #Backward
  step(m2b,trace=T,direction="backward") # I(Slag^2) + I(Water^2)

  #Forward
  horizonte <- formula(ybox ~ Cement + I(Slag^2) +`Fly ash`+ I(Water^2) + SP+ `Fine Aggr.`)
  modelo0 <- lm(ybox ~ 1,data=me[,-c(6,8,9,10)])
  step(modelo0, trace = T, direction = "forward", scope = horizonte) # I(Slag^2) + I(Water^2)
  
  #Both
  step(modelo0, trace=T, direction="both", scope=horizonte) # I(Slag^2) + I(Water^2)
  
  # Como resultado tenemos tres posibles
  m2_1 <- lm(ybox ~ I(Slag^2) + I(Water^2),data=me[,-c(6,8,9,10)])

#########################################
# Medir la multicolinealidad #########
######################################
library(car)
  # Modelo m2_1
  vif(m2_1) # No hay multicolinealidad alta
  m2f <- lm(ybox~ I(Slag^2) + I(Water^2),data=me[,-c(6,8,9,10)])
################
# Supuestos ####
################
shapiro.test(m2f$residuals) # 0.09608
bptest(m2f) # 0.4947
dwtest(m2f, alternative= "t") # 0.8982
summary(m2f)

###################################################################
# VALIDACIÓN DEL MODELO FLOW (Modelo final: I(Slag^2) + I(Water^2)) ####
###################################################################

library(car)
library(lmtest)
library(olsrr)
library(dplyr)

# Renombrar variable respuesta
datos <- Conjunto_2[,-1]
datos <- rename(datos, Y = "FLOW(cm)")

# Eliminar SLUMP y STRENGTH (columnas 8 y 10)
datos_modelo <- datos[,-c(8,10)]

# Renombrar variables para evitar errores con espacios
names(datos_modelo) <- gsub(" ", "_", names(datos_modelo))
names(datos_modelo) <- gsub("\\.", "", names(datos_modelo))

###################################################################
# 1. Ajustar modelo inicial para obtener lambda de Box-Cox
###################################################################

set.seed(40)
n_total <- nrow(datos_modelo)
ne <- round(0.8 * n_total)
indices_inicial <- sample(n_total, ne)

me <- datos_modelo[indices_inicial, ]

# Modelo preliminar SOLO PARA ESTIMAR LAMBDA
m_box <- lm(Y ~ Slag + Water, data = me)

box <- boxCox(m_box, lambda = seq(-2, 2, by = 0.1))
lambda <- box$x[which.max(box$y)]

# Transformación Box–Cox
ybox <- (me$Y^lambda - 1)/lambda
me$ybox <- ybox

# Modelo final ya transformado
m2f <- lm(ybox ~ I(Slag^2) + I(Water^2), data = me)

###################################################################
# 2. Transformación inversa de Box-Cox
###################################################################
transformacion_inversa <- function(y_trans, lambda){
  if(abs(lambda) > 1e-6){
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
  
  # Transformación Y → ybox
  me_i$ybox <- (me_i$Y^lambda - 1)/lambda
  
  # Ajustar el modelo EXACTO que estás usando
  mod_i <- lm(ybox ~ I(Slag^2) + I(Water^2), data = me_i)
  
  # Predicción en escala transformada
  pred_trans <- predict(mod_i, newdata = mp_i)
  
  # Transformación inversa
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

AIC_valor <- AIC(m2f)
BIC_valor <- BIC(m2f)
R2_ajustado <- summary(m2f)$adj.r.squared

# R² de predicción (PRESS)
R2_pred <- ols_pred_rsq(m2f)

# Cp de Mallows
ols_all <- ols_step_all_possible(m2f)
Cp_mejor <- min(ols_all$result[,"Mallow's Cp"])

cat("\n=== INDICADORES ADICIONALES ===\n")
cat("AIC            :", round(AIC_valor, 2), "\n")
cat("BIC            :", round(BIC_valor, 2), "\n")
cat("R² ajustado    :", round(R2_ajustado, 3), "\n")
cat("R² predicción  :", round(R2_pred, 3), "\n")
cat("Cp de Mallows  :", round(Cp_mejor, 3), "\n")

