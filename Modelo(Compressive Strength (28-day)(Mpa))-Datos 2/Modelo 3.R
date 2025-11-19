#####################################################################################
# Tercer Modelo 
#####################################################################################

datos <- Conjunto_2[,-1]
str(datos)
library(dplyr)
datos<- rename(datos, "Y" ="Compressive Strength (28-day)(Mpa)")

###################################
# Validacion cruzada#####
###################################
#Optimista
set.seed(40)
ne<- round(0.8 *nrow(datos))
indices <- sample(103,ne)

me<- datos[indices,] # Entrenar Modelo
mp<- datos[-indices,] # Probar Modelo

# Se usara los datos me
m1 <- lm(Y~.,data=me[,-c(8,9)])
anova(m1)
shapiro.test(m1$residuals)
library(lmtest)
bptest(m1)
#########################################
# Reespecificación de las varibales #########
######################################
library(GGally)
ggpairs(datos[,-c(8,9)])

# Relación agua/cemento (W/C)
# Porque Water correlaciona con Cement, Slag y Fly as
me$WC <- me$Water / me$Cement

# Binder total (ligante total): Cement + Slag + Fly ash
# Esto reduce la multicolinealidad entre Cement–Slag–Fly ash.
me$Binder <- me$Cement + me$Slag + me$`Fly ash`

# Relación Fine / Coarse Aggregates
# Correlación de –0.602 indica multicolinealidad fuerte entre ambos → reemplázalos.
me$FA_CA <- me$`Fine Aggr.` / me$`Coarse Aggr.`
# Dosis relativa de superplastificante
# Dado que SP correlaciona con Cement (+0.307), lo reemplazas por SP/Cement.
me$SPC <- me$SP / me$Cement

modelo <- lm(Y~ Cement + Slag +`Fly ash`+Water + SP + FA_CA,data=me[,-c(6:9)])
shapiro.test(modelo$residuals)
#########################################
# Medir la multicolinealidad #########
######################################
library(car)
vif(modelo)
m1x6 <- lm(Y~Cement + Slag +`Fly ash`+Water + SP + FA_CA,data=me[,-c(6:9)]) #Se usara ese modelo de ahora en adelante
vif(m1x6) # No hay multicolinealidad con: Cement Slag `Fly ash` Water SP `Fine Aggr.

###################################################################
# Transformación                                               ####
###################################################################

# Como no se cumple la normalidad en el modelo m1x6f, se realizara una transformación
library(car)
box<-boxCox(m1x6, lambda = seq(-2, 2, by = 0.1)) # Prueba de valores
lambda<-box$x[which.max(box$y)] # Valor de λ optimo

ybox<-(me$Y^lambda-1)/lambda # Transformación de Y
m1x6f <- lm(ybox ~ Cement + Slag +`Fly ash`+Water + SP + FA_CA,data=me[,-c(6:10)])


########################################################
# Analisamos las observaciones influyentes y outlier####

# Outliers 
# Residuos Estandarizados
ri <-abs(rstandard(m1x6f)) > 2
ri[ri==TRUE] # Son valores atipicos y pueden ser puntos influyentes:  31, 34, 51, 78 
# R. Estunderizados
ti<- abs(rstudent(m1x6f)) >2
ti[ti == TRUE]  # Son valores atipicos y pueden ser puntos influyentes : 31, 34, 51, 57, 78 
# Bonferroni
library(car)
outlierTest(m1x6f, cutoff=Inf, n.max=5, order=TRUE) # 78
# con un alfa de 0.05 no hay observaciones que se puedan considerar outlier estadisticamente

# Influyentes
# Matriz H ####
hii<-hatvalues(m1x6f)
sum(hii)
2*5/100 > 1
r<- abs(hii) >2*(7)/nrow(me)
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

########################################
# Seleccionamos Variables ##############
########################################


###################################
#Método de mejores subconjuntos####
###################################

library(olsrr)
res<-ols_step_all_possible(m1x6f)
#Con el RMSE
which.min(res$result[,6])
res$result[63,] # Cement + Slag + `Fly ash` + Water + SP + FA_CA

#Con el AIC 
which.min(res$result[,9])
res$result[63,] # Cement + Slag + `Fly ash` + Water + SP + FA_CA

#Con el R2 pred 
which.max(res$result[,7])
res$result[63,]# Cement + Slag + `Fly ash` + Water + SP + FA_CA
m1x6_63<- lm(ybox~ Cement + Slag +`Fly ash` +Water +SP +`Fine Aggr.`,data=me[,-c(6,8,9)])

################################
#Método por pasos (Stepwise)####
################################

#Backward
step(m1x6f,trace=T,direction="backward") # Cement + Slag + `Fly ash` + Water + SP + FA_CA
#Forward
horizonte <- formula(ybox ~  Cement + Slag + `Fly ash` + Water + SP + FA_CA)
modelo0 <- lm(ybox ~ 1, data = me[,-c(6:10)])
step(modelo0, trace = T, direction = "forward", scope = horizonte) #Cement + Slag + `Fly ash` + Water + SP + FA_CA

#Both
step(modelo0, trace=T, direction="both", scope=horizonte) #Cement + Slag + `Fly ash` + Water + SP + FA_CA
# En los 3 metodos por pasos coinciden con el mismo modelo y ademas coincide con el metodo de mejores subconjuntos
m1f <- lm(ybox~ Cement + Slag + `Fly ash` + Water + SP + FA_CA,data=me[,-c(6:10)])
shapiro.test(m1f$residuals)
bptest(m1f)

###################################################################
# Validación del Modelo m1f                                      ####
###################################################################

library(car)
library(lmtest)
library(olsrr)
library(caret)
library(dplyr)

# Renombrar la variable respuesta
datos <- rename(Conjunto_2[,-1], Y = "Compressive Strength (28-day)(Mpa)")

# Crear FA_CA
datos$FA_CA <- datos$`Fine Aggr.` / datos$`Coarse Aggr.`

# Eliminar SLUMP (6), FLOW (8), COARSE (9) y FINE (7, porque lo reemplaza FA_CA)
datos_modelo <- datos[, c("Cement","Slag","Fly ash","Water","SP","FA_CA","Y")]

# ---------------------------------------------------------------------
# AJUSTE INICIAL PARA OBTENER LAMBDA
# ---------------------------------------------------------------------
set.seed(40)
n_total <- nrow(datos_modelo)
ne <- round(0.8 * n_total)
indices_inicial <- sample(n_total, ne)

me <- datos_modelo[indices_inicial, ]

# Ajustar modelo sin transformar para obtener lambda
m1_raw <- lm(Y ~ Cement + Slag + `Fly ash` + Water + SP + FA_CA, data = me)

box <- boxCox(m1_raw, lambda = seq(-2, 2, by = 0.1))
lambda <- box$x[which.max(box$y)]
lambda

# Transformar la variable respuesta Y → ybox
ybox <- (me$Y^lambda - 1) / lambda
me$ybox <- ybox

# Ajustar modelo transformado m1f
m1f <- lm(ybox ~ Cement + Slag + `Fly ash` + Water + SP + FA_CA, data = me)

# ---------------------------------------------------------------------
# FUNCIÓN DE TRANSFORMACIÓN INVERSA
# ---------------------------------------------------------------------
transformacion_inversa <- function(y_trans, lambda) {
  if (abs(lambda) > 1e-6) {
    (lambda * y_trans + 1)^(1 / lambda)
  } else {
    exp(y_trans)
  }
}

# ---------------------------------------------------------------------
# VALIDACIÓN HOLD-OUT REPETIDA (1000 ITERACIONES)
# ---------------------------------------------------------------------
set.seed(40)
R <- 1000
MSE_vec <- numeric(R)
MAE_vec <- numeric(R)

for (i in 1:R) {
  
  # Separación aleatoria
  idx <- sample(n_total, ne)
  me_i <- datos_modelo[idx, ]
  mp_i <- datos_modelo[-idx, ]
  
  # Transformar Y en entrenamiento
  y_trans_i <- if (abs(lambda) > 1e-6) {
    (me_i$Y^lambda - 1) / lambda
  } else {
    log(me_i$Y)
  }
  
  me_i$ybox <- y_trans_i
  
  # Ajustar modelo m1f
  mod_i <- lm(ybox ~ Cement + Slag + `Fly ash` + Water + SP + FA_CA, data = me_i)
  
  # Predicción en escala transformada
  pred_trans <- predict(mod_i, newdata = mp_i)
  
  # Transformar a escala original
  pred_original <- transformacion_inversa(pred_trans, lambda)
  real_original <- mp_i$Y
  
  # Errores
  errores <- pred_original - real_original
  MSE_vec[i] <- mean(errores^2)
  MAE_vec[i] <- mean(abs(errores))
}

# Resultados
RMSE_final <- sqrt(mean(MSE_vec))
MAE_final  <- mean(MAE_vec)

cat("=== RESULTADOS DE VALIDACIÓN (1000 repeticiones) ===\n")
cat("RMSE :", round(RMSE_final, 3), "MPa\n")
cat("MAE  :", round(MAE_final, 3), "MPa\n")

# ---------------------------------------------------------------------
# INDICADORES ADICIONALES DEL MODELO AJUSTADO
# ---------------------------------------------------------------------

AIC_valor <- AIC(m1f)
BIC_valor <- BIC(m1f)
R2_ajustado <- summary(m1f)$adj.r.squared
R2_pred <- ols_pred_rsq(m1f)

ols_all <- ols_step_all_possible(m1f)
Cp_mejor <- min(ols_all$result[,"Mallow's Cp"])

cat("\n=== INDICADORES ADICIONALES ===\n")
cat("AIC            :", round(AIC_valor, 2), "\n")
cat("BIC            :", round(BIC_valor, 2), "\n")
cat("R² ajustado    :", round(R2_ajustado, 3), "\n")
cat("R² predicción  :", round(R2_pred, 3), "\n")
cat("Cp de Mallows  :", round(Cp_mejor, 3), "\n")
