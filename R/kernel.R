#' @title kernel
#' @param u es un vector con los xs
#' @export kernel
kernel = function(u) {
    0.75 * (1 - u^2) * (u <= 1) * (u >= -1)
}  #Epanechnikov
