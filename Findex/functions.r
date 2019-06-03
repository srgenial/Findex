# adf test

adf <- function(x, type = c("none", "const", "trend"),
                alternative = c("stationary", "explosive"),
                k = trunc((length(x) - 1)^(1/3))){

    if (NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    if (k < 0)
        stop("k negative")

    type <- match.arg(type)
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))

    k <- k + 1
    x <- as.vector(x, mode = "double")
    y <- diff(x)
    n <- length(y)
    z <- embed(y, k)
    yt <- z[, 1]
    xt1 <- x[k:n]
    tt <- k:n

    if(type == "none"){
      if(k > 1){
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + yt1 - 1)
      }else{
        res <- lm(yt ~ xt1 - 1)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[1, 1]/res.sum$coefficients[1,
                2]
      Table <- adfTable(trend = "nc", statistic = "t")
    }

    if(type == "const"){
      if(k > 1){
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + 1 + yt1)
      }else{
        res <- lm(yt ~ xt1 + 1)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                2]
      Table <- adfTable(trend = "c", statistic = "t")
    }

    if(type == "trend"){
      if (k > 1) {
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + 1 + tt + yt1)
      }else{
        res <- lm(yt ~ xt1 + 1 + tt)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                2]
      Table <- adfTable(trend = "ct", statistic = "t")
    }
    bic <- BIC(res)
    table <- Table$z
    tablen <- dim(table)[2]
    tableT <- Table$x
    tablep <- Table$y
    tableipl <- numeric(tablen)
    for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[,
        i], n, rule = 2)$y
    interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
    if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y))
        if (interpol == min(tablep))
            warning("p-value smaller than printed p-value")
        else warning("p-value greater than printed p-value")
    if (alternative == "stationary")
        PVAL <- interpol
    else if (alternative == "explosive")
        PVAL <- 1 - interpol
    else stop("irregular alternative")
    PARAMETER <- k - 1
    METHOD <- "Augmented Dickey-Fuller Test"
    names(STAT) <- "Dickey-Fuller"
    names(PARAMETER) <- "Lag order"
    structure(list(bic = bic, statistic = STAT, parameter = PARAMETER,
        alternative = alternative, p.value = PVAL, method = METHOD,
        data.name = DNAME), class = "htest")
}

adfTable <-
  function(trend = c("nc", "c", "ct"), statistic = c("t", "n"),
           includeInf = TRUE)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Tables critical values for augmented Dickey-Fuller test.

    # Note:
    #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z

    # Examples:
    #   adfTable()

    # FUNCTION:

    # Match Arguments:
    type = trend = match.arg(trend)
    statistic = match.arg(statistic)

    # Tables:
    if (statistic == "t") {
      # Hamilton Table B.6 - OLS t-Statistic
      if (type == "nc") {
        table = cbind(
          c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
          c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
          c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
      } else if (type == "c") {
        table = cbind(
          c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
          c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
          c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
          c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
          c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
          c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
      } else if (type == "ct") {
        table = cbind(
          c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
          c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
          c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
          c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
          c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
          c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
      } else {
        stop("Invalid type specified")
      }
    } else if (statistic == "z" || statistic == "n") {
      # Hamilton Table B.5 - Based on OLS Autoregressive Coefficient
      if (type == "nc") {
        table = cbind(
          c(-11.9,  -9.3,  -7.3,  -5.3, +1.01, +1.40, +1.79, +2.28),
          c(-12.9,  -9.9,  -7.7,  -5.5, +0.97, +1.35, +1.70, +2.16),
          c(-13.3, -10.2,  -7.9,  -5.6, +0.95, +1.31, +1.65, +2.09),
          c(-13.6, -10.3,  -8.0,  -5.7, +0.93, +1.28, +1.62, +2.04),
          c(-13.7, -10.4,  -8.0,  -5.7, +0.93, +1.28, +1.61, +2.04),
          c(-13.8, -10.5,  -8.1,  -5.7, +0.93, +1.28, +1.60, +2.03))
      } else if (type == "c") {
        table = cbind(
          c(-17.2, -14.6, -12.5, -10.2, -0.76, +0.01, +0.65, +1.40),
          c(-18.9, -15.7, -13.3, -10.7, -0.81, -0.07, +0.53, +1.22),
          c(-19.8, -16.3, -13.7, -11.0, -0.83, -0.10, +0.47, +1.14),
          c(-20.3, -16.6, -14.0, -11.2, -0.84, -0.12, +0.43, +1.09),
          c(-20.5, -16.8, -14.0, -11.2, -0.84, -0.13, +0.42, +1.06),
          c(-20.7, -16.9, -14.1, -11.3, -0.85, -0.13, +0.41, +1.04))
      } else if (type == "ct") {
        table = cbind(
          c(-22.5, -19.9, -17.9, -15.6, -3.66, -2.51, -1.53, -0.43),
          c(-25.7, -22.4, -19.8, -16.8, -3.71, -2.60, -1.66, -0.65),
          c(-27.4, -23.6, -20.7, -17.5, -3.74, -2.62, -1.73, -0.75),
          c(-28.4, -24.4, -21.3, -18.0, -3.75, -2.64, -1.78, -0.82),
          c(-28.9, -24.8, -21.5, -18.1, -3.76, -2.65, -1.78, -0.84),
          c(-29.5, -25.1, -21.8, -18.3, -3.77, -2.66, -1.79, -0.87))
      } else {
        stop("Invalid type specified")
      }
    } else {
      stop("Invalid statistic specified")
    }

    # Transpose:
    Table = t(table)
    colnames(Table) = c("0.010", "0.025", "0.050", "0.100", "0.900",
                        "0.950", "0.975", "0.990")
    rownames(Table) = c(" 25", " 50", "100", "250", "500", "Inf")
    ans = list(
      x = as.numeric(rownames(Table)),
      y = as.numeric(colnames(Table)),
      z = Table)
    class(ans) = "gridData"

    # Exclude Inf:
    if (!includeInf) {
      nX = length(ans$x)
      ans$x = ans$x[-nX]
      ans$z = ans$z[-nX, ]
    }

    # Add Control:
    attr(ans, "control") <-
      c(table = "adf", trend = trend, statistic = statistic)

    # Return Value:
    ans
  }

# function lags
lag.f <- function(x, lags){
  xl <- c(rep(NA, lags), x[1:(length(x)-lags)])
  return(xl)
}

# function MCE bivariate
MCE <- function(regre_lp, lags = 1){

  # dat
  dat <- regre_lp$model

  # diff terms
  dy <- c(NA, diff(dat[,1]))
  dx <- c(NA, diff(dat[,2]))

  # correction error
  ce <- resid(regre_lp)

  if(lags > 0){

    dyl <- dxl <- matrix(0, nrow(dat), lags)
    for(i in 1 : lags){
      dyl[,i] <- lag.f(dy, i)
      dxl[,i] <- lag.f(dx, i)
    }
    # mce
    mce.regre <- lm(dy ~ lag.f(ce, 1) + dx + dxl + dyl)
  }else{
    mce.regre <- lm(dy ~ lag.f(ce, 1) + dx)
  }

  return(mce.regre)
}

# Raíces unitarias (clase 2)

# Fijamos semilla
set.seed(12345)

# PaqueterIa para realizar pruebas de raIces unitarias, simulaciOn de modelos
# y otros estadIsticos
library(tseries)
library(forecast)
library(portes)
library(e1071)  
library(lmtest)

# RaIz unitaria
opp <- par(mfrow = c(2,1))
ts.plot(cumsum(rnorm(100)), ylab = "", xlab = "")
title("Como acumulaci?n de shocks")
ts.plot(varima.sim(model=list(ar=1), n = 100, constant = 100), ylab = "",
        xlab = "")
title("Usando la librer?a portes")
par(opp)


# REplicas
R <- 1000
# Tamanio de la series de tiempo
Tt <- 1000
# Nivel de signficancia
alpha <- 0.05
# Matriz donde guardamos resultados de interEs
stat <- matrix(0, R, 3)
colnames(stat) <- c("t-value", "p-t", "p-df")

# SimulaciOn de caminatas aleatorias, regresiOn y prueba ADF
for(i in 1 : R){
  yt <- cumsum(rnorm(Tt))
  sumregre <- coef(summary(lm(diff(yt) ~ yt[-Tt])))
  stat[i, "t-value"] <- sumregre[2,"t value"]
  stat[i, "p-t"] <- sumregre[2,"Pr(>|t|)"]
  stat[i, "p-df"] <- adf.test(yt)$p.value
}

# Quantual de acuerdo al nivel de significancia
quant <- quantile(stat[, "t-value"], alpha)

# Algunas pruebas
skewness(stat[, "t-value"])
kurtosis(stat[, "t-value"])
shapiro.test(stat[, "t-value"])

# Histograma del estadIstico t
hist(stat[, "t-value"], main = "", xlab = "", col = "blue")
abline(v = quant, col = "red", lwd = 2)

# Plot ilustrativo
plot(density(stat[, "p-t"]), col = "red", main = "", xlab = "")
lines(density(stat[, "p-df"]), col = "blue")
lines(x = rep(alpha, R), y = seq(0,
                                 density(stat[, "p-t"])$y[which(density(stat[, "p-t"])$x > alpha)[1]],
                                 length.out = R), col = "green")

# A cuatro decimales
stat <- round(stat, 4)

# NUmero de veces en que se equivocaron usando t y DF
sum(stat[, "p-t"] < alpha)/R
sum(stat[, "p-df"] < alpha)/R
