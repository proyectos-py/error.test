context("test")
n1 = 39
n2 = 49
x1 = runif(n1)
y1 = runif(n1)
x2 = runif(n2)
y2 = runif(n2)
error.test(x1, y1, x2, y2, methods = "SB")
error.test(x1, y1, x2, y2, methods = "WB")
