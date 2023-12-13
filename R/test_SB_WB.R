
#'  Nonparametric Two-Sample Test for Error Distributions
#'
#' @description This function performs a two-sample test for comparing error distributions
#' in nonparametric regression. It supports two methods: Weighted Bootstrap (WB)
#' and Smooth Bootstrap (SB).
#'
#' @param x1 Numeric vector representing the predictor values for the first sample.
#' @param y1 Numeric vector representing the response values for the first sample.
#' @param x2 Numeric vector representing the predictor values for the second sample.
#' @param y2 Numeric vector representing the response values for the second sample.
#' @param methods Character vector specifying the method(s) to use, either 'SB', 'WB',
#'   or both. Default is c('SB', 'WB').
#'
#' @return A list containing information about the test including sample sizes (n1, n2)
#'   and the p-value.
#'
#' @references
#' Kernel Density Estimation and Its Applications, G. S. Watson, 1983.
#'
#' @export
#' @examples
#' n1 = 39
#' n2 = 49
#' x1 = runif(n1)
#' y1 = runif(n1)
#' x2 = runif(n2)
#' y2 = runif(n2)
#'error.test(x1, y1, x2, y2, methods = 'SB')
#'error.test(x1, y1, x2, y2, methods = 'WB')

error.test = function(x1, y1, x2, y2, methods = c("SB", "WB")) {
    c = 1.5
    a = 0.45
    s = 1234567
    B = 1000
    r = 0.15
    M = 1  #valores cambiables
    nw.mean = function(x, xdata, ydata, h) {
        n = length(xdata)
        h = c * n^(-a)  #parámetro de ventana
        aux <- (xdata - x)/h
        aux <- kernel(aux)
        numerador = sum(aux * ydata)
        denominador = sum(aux)
        salida = numerador/denominador
        return(salida)
    }
    # ========================================================================================
    nw.sd = function(x, xdata, ydata, h) {
        n = length(xdata)
        h = c * n^(-a)
        aux <- (xdata - x)/h
        aux <- kernel(aux)
        denominador = sum(aux)
        zeta = sum(aux * ydata)/denominador
        zeta = (ydata - zeta)^2
        numerador = sum(zeta * aux)
        salida = sqrt(numerador/denominador)
        return(salida)
    }
    # ======================================================================================== Cálculo de residuos residuos en el vector x
    res = function(xdata, ydata, h) {
        n = length(xdata)
        salida = rep(0, n)
        for (i in 1:n) {
            m = nw.mean(xdata[i], xdata, ydata, h)
            d = nw.sd(xdata[i], xdata, ydata, h)
            salida[i] = (ydata[i] - m)/d
        }
        return(salida)
    }
    # ========================================================================================
    n1 = length(x1)
    ee1 = rep(0, n1)
    for (i in 1:n1) {
        h1 = c * n1^(-a)
        m = nw.mean(x1[i], x1, y1, h1)
        d = nw.sd(x1[i], x1, y1, h1)
        ee1[i] = (y1[i] - m)/d
    }
    n2 = length(x2)
    ee2 = rep(0, n2)
    for (i in 1:n2) {
        h2 = c * n2^(-a)
        m = nw.mean(x2[i], x2, y2, h2)
        d = nw.sd(x2[i], x2, y2, h2)
        ee2[i] = (y2[i] - m)/d
    }
    # ========================================================================================
    a11 = outer(ee1, ee1, "-")
    a22 = outer(ee2, ee2, "-")
    a12 = outer(ee1, ee2, "-")
    a21 = outer(ee2, ee1, "-")
    # ========================================================================================

    phi = function(a, r) {
        salida = sqrt(pi/r) * exp(-a^2/(4 * r))
        return(salida)
    }
    phi2 = function(a, r) {
        salida = (1/(2 * r)) * sqrt(pi/r) * exp(-a^2/(4 * r)) * (a^2/(2 * r) - 1)
        return(salida)
    }
    phi1 = function(a, r) {
        salida = -(a/(2 * r)) * sqrt(pi/r) * exp(-a^2/(4 * r))
        return(salida)
    }
    # ========================================================================================
    if (methods == "WB") {
        M1 = (1/(n1^2)) * sum(phi(a11, r))
        M2 = (1/(n2^2)) * sum(phi(a22, r))
        M3 = -(2/(n1 * n2)) * sum(phi(a12, r))
        Tobs = M1 + M2 + M3
        # ======================================================================================== Cálculo del p-valor WB
        N = n1 + n2
        A = (1/N^2) * (sum(phi2(a11, r)) + sum(phi2(a22, r)) + 2 * sum(phi2(a12, r)))
        aux1 = matrix(ee1, n1, n1)
        aux2 = matrix(ee2, n2, n2)
        aux3 = matrix(ee2, n2, n1)
        aux4 = matrix(ee1, n1, n2)
        # ======================================================================================== Cálculo del p-valor WB
        A = (1/N^2) * (sum(phi2(a11, r)) + sum(phi2(a22, r)) + 2 * sum(phi2(a12, r)))
        aux1 = matrix(ee1, n1, n1)
        aux2 = matrix(ee2, n2, n2)
        aux3 = matrix(ee2, n2, n1)
        aux4 = matrix(ee1, n1, n2)
        C = (1/N^2) * (sum(aux1 * phi2(a11, r)) + sum(aux3 * phi2(a21, r)) + sum(aux4 * phi2(a12, r)) + sum(aux2 * phi2(a22, r)))
        D = (1/N^2) * (sum(aux1 * phi1(a11, r)) + sum(t(aux4) * phi1(a21, r)) + sum(t(aux3) * phi1(a12, r)) + sum(aux2 * phi1(a22, r)))
        aux11 = outer(ee1, ee1)
        aux22 = outer(ee2, ee2)
        aux12 = outer(ee1, ee2)
        E = (1/N^2) * (sum(aux11 * phi2(a11, r)) + sum(aux22 * phi2(a22, r)) + 2 * sum(aux12 * phi2(a12, r)))
        F = (1/N^2) * (sum(phi(a11, r)) + sum(phi(a22, r)) + 2 * sum(phi(a12, r)))
        # ======================================================================================== inicio del c?lculo de la matriz M11(jl)
        T1 = phi(a11, r)
        v21 = colSums(phi1(a11, r)) + colSums(phi1(a21, r))
        v22 = -(1/N) * outer(ee1, v21)
        T2 = v22 + t(v22)
        v31 = rowSums(phi1(a11, r) * t(aux1)) + rowSums(phi1(a12, r) * t(aux3))
        aux = (ee1^2 - 1)/2
        v32 = (1/N) * outer(aux, v31)
        T3 = v32 + t(v32)
        v41 = rowSums(phi(a11, r)) + rowSums(phi(a12, r))
        v42 = matrix(v41, n1, n1)
        T4 = -(1/N) * (v42 + t(v42))
        T5 = -outer(ee1, ee1) * A
        v61 = outer(ee1, (ee1^2 - 1)/2)
        T6 = -(v61 + t(v61)) * C
        v71 = outer((ee1^2 - 1)/2, (ee1^2 - 1)/2, "+")
        T7 = -v71 * D
        v81 = outer((ee1^2 - 1)/2, (ee1^2 - 1)/2)
        T8 = -v81 * E
        T9 = matrix(F, n1, n1)
        M11 = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9
        # ======================================================================================== inicio del cálculo de la matriz M22(jl)
        T1 = phi(a22, r)
        v21 = colSums(phi1(a12, r)) + colSums(phi1(a22, r))
        v22 = -(1/N) * outer(ee2, v21)
        T2 = v22 + t(v22)
        v31 = rowSums(phi1(a21, r) * t(aux4)) + rowSums(phi1(a22, r) * t(aux2))
        aux = (ee2^2 - 1)/2
        v32 = (1/N) * outer(aux, v31)
        T3 = v32 + t(v32)
        v41 = rowSums(phi(a21, r)) + rowSums(phi(a22, r))
        v42 = matrix(v41, n2, n2)
        T4 = -(1/N) * (v42 + t(v42))
        T5 = -outer(ee2, ee2) * A
        v61 = outer(ee2, (ee2^2 - 1)/2)
        T6 = -(v61 + t(v61)) * C
        v71 = outer((ee2^2 - 1)/2, (ee2^2 - 1)/2, "+")
        T7 = -v71 * D
        v81 = outer((ee2^2 - 1)/2, (ee2^2 - 1)/2)
        T8 = -v81 * E
        T9 = matrix(F, n2, n2)
        M22 = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9
        # ======================================================================================== inicio del cálculo de la matriz M12(jl)
        T1 = phi(a12, r)
        v21 = colSums(phi1(a12, r)) + colSums(phi1(a22, r))
        v22 = -(1/N) * outer(ee1, v21)
        v211 = colSums(phi1(a11, r)) + colSums(phi1(a21, r))
        v221 = -(1/N) * outer(v211, ee2)
        T2 = v22 + v221
        v31 = rowSums(phi1(a21, r) * t(aux4)) + rowSums(phi1(a22, r) * t(aux2))
        aux = (ee1^2 - 1)/2
        v32 = (1/N) * outer(aux, v31)
        v311 = rowSums(phi1(a11, r) * t(aux1)) + rowSums(phi1(a12, r) * t(aux3))
        aux = (ee2^2 - 1)/2
        v321 = (1/N) * outer(v311, aux)
        T3 = v32 + v321
        v41 = rowSums(phi(a11, r)) + rowSums(phi(a12, r))
        v42 = matrix(v41, n1, n2)
        v411 = rowSums(phi(a21, r)) + rowSums(phi(a22, r))
        v421 = t(matrix(v411, n2, n1))
        T4 = -(1/N) * (v42 + v421)
        T5 = -outer(ee1, ee2) * A
        v61 = outer(ee1, (ee2^2 - 1)/2)
        v62 = outer((ee1^2 - 1)/2, ee2)
        T6 = -(v61 + v62) * C
        v71 = outer((ee1^2 - 1)/2, (ee2^2 - 1)/2, "+")
        T7 = -v71 * D
        v81 = outer((ee1^2 - 1)/2, (ee2^2 - 1)/2)
        T8 = -v81 * E
        T9 = matrix(F, n1, n2)
        M12 = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9
        # ========================================================================================
        set.seed(s)
        auxWB2 = rep(0, B)
        for (b in 1:B) {
            mult1 = stats::rnorm(n1)
            multcent1 = mult1 - mean(mult1)
            mult2 = stats::rnorm(n2)
            multcent2 = mult2 - mean(mult2)
            auxWB2[b] = as.numeric(multcent1 %*% M11 %*% multcent1)/(n1^2) + as.numeric(multcent2 %*% M22 %*% multcent2)/(n2^2) - 2 * as.numeric(multcent1 %*%
                M12 %*% multcent2)/(n1 * n2)
        }
        vale2 = rep(1, B)[auxWB2 > Tobs]
        pvalorWB2 = (1/B) * sum(vale2)
        z1 = pvalorWB2
        # ======================================================================================== cat('\nA two-sample test for the error
        # distribution in nonparametric regression based on WB approximation.')
        # cat('\n--------------------------------------------------------------------------------') cat('\nHypothesis')
        cat("\nHo: The error distributions funcions are equal in both population")
        cat("\nP-value based on Weighted Bootstrap aproximation")
        cat("\n--------------------------------------------------------------------\n")
        M <- list(n1 = n1, n2 = n2, `p-value` = z1)
        suppressWarnings(suppressMessages(M))
        return(M)

    }
    # ========================================================================================
    if (methods == "SB") {
        set.seed(s)
        pvalorB = rep(0, M)
        auxBoot = rep(0, B)
        # cálculo del estadístico observado
        for (m in 1:M) {
            M1 = (1/(n1^2)) * sum(phi(a11, r))
            M2 = (1/(n2^2)) * sum(phi(a22, r))
            M3 = -(2/(n1 * n2)) * sum(phi(a12, r))

            Tobs = M1 + M2 + M3  #valor observado

            # Cálculo del p-valor mediante Bootstrap suavizado (BS)

            # Estandarización de los errores de la muestra conjunta
            ee = c(ee1, ee2)
            ee = (ee - mean(ee))/stats::sd(ee)

            # Cálculo de los datos bootstrap
            yb <- function(xdata, ydata, e, h) {
                n = length(xdata)
                salida = rep(0, n)
                for (i in 1:n) {
                  salida[i] = nw.mean(xdata[i], xdata, ydata, h) + nw.sd(xdata[i], xdata, ydata, h) * e[i]
                }
                return(salida)
            }


            # inicio del bootstrap
            for (b in 1:B) {
                # inicio B PASO 1
                eboot1 = (1 - 4 * n1^(-3/5))^(1/2) * (sample(ee, size = n1, replace = TRUE)) + 2 * n1^(-3/10) * stats::rnorm(n1)
                eboot2 = (1 - 4 * n2^(-3/5))^(1/2) * (sample(ee, size = n2, replace = TRUE)) + 2 * n2^(-3/10) * stats::rnorm(n2)

                # PASO 2
                yboot1 = yb(x1, y1, eboot1, h1)
                yboot2 = yb(x2, y2, eboot2, h2)

                # PASO 3
                eboot1 = res(x1, yboot1, h1)
                eboot2 = res(x2, yboot2, h2)

                a11 = outer(eboot1, eboot1, "-")
                a22 = outer(eboot2, eboot2, "-")
                a12 = outer(eboot1, eboot2, "-")

                M1 = (1/(n1^2)) * sum(phi(a11, r))
                M2 = (1/(n2^2)) * sum(phi(a22, r))
                M3 = -(2/(n1 * n2)) * sum(phi(a12, r))

                auxBoot[b] = M1 + M2 + M3  #valores calculados del estadístico

            }  #fin B

            # Estimación del p-valor

            vale1 = rep(1, B)[auxBoot > Tobs]
            pvalorB = (1/B) * sum(vale1)
        }  #fin M

        z1 = pvalorB
        # ======================================================================================== cat('\nA weighted bootstrap approximation
        # for comparing the error distributions based on SB approximation.')
        # cat('\n----------------------------------------------------------------------------------------\n') cat('\nHypothesis')
        cat("\nHo: The error distributions funcions are equal in both population")
        cat("\nP-value based on Smooth Bootstrap aproximation")
        cat("\n--------------------------------------------------------------------\n")
        y <- list(n1 = n1, n2 = n2, `p-value` = z1)

        return(y)
    }

}



