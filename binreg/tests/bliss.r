library("SMPracticals")
data(bliss)

source("loadlibs.r")
loadlibs()

y <- rep(0,sum(bliss$m))
x <- rep(0,sum(bliss$m))
k <- 0
for (j in 1:length(bliss$m)) {
  for (i in 1:bliss$m[j]) {
    k <- k+1
    if (i <= bliss$r[j])
      y[k] <- 1
    x[k] <- log(bliss$dose[j])
  }
}

if (TRUE) {
fit.l  <- glm(y ~ x, binomial)
fit.p  <- glm(y ~ x, binomial("probit"))
fit.cl <- glm(y ~ x, binomial("cloglog"))
fit.a  <- binreg(y ~ x, "A")
fit.w  <- binreg(y ~ x, "W")
fit.cw <- binreg(y ~ x, "C")
fit.s  <- binreg(y ~ x, "S")
fit.pr <- binreg(y ~ x, "P")

summary(fit.l)
summary(fit.p)
summary(fit.cl)
summary(fit.a)
summary(fit.w)
summary(fit.cw)
summary(fit.s)
summary(fit.pr)

data     <- dataLD(y ~ x, "cloglog")
b.fit.cl <- LaplacesDemon(ModelCloglog(),  data, GIV(ModelCloglog(), data))
data     <- dataLD(y ~ x, "logit")
b.fit.l  <- LaplacesDemon(ModelLogit(),    data, GIV(ModelLogit(), data))
data     <- dataLD(y ~ x, "probit")
b.fit.p  <- LaplacesDemon(ModelProbit(),   data, GIV(ModelProbit(), data))
data     <- dataLD(y ~ x, "Aranda-Ordaz")
b.fit.a  <- LaplacesDemon(ModelAranda(),   data, GIV(ModelAranda(), data))
data     <- dataLD(y ~ x, "Weibull")
b.fit.w  <- LaplacesDemon(ModelWeibull(),  data, GIV(ModelWeibull(), data))
data     <- dataLD(y ~ x, "CWeibull")
b.fit.cw <- LaplacesDemon(ModelCWeibull(), data, GIV(ModelCWeibull(), data))
data     <- dataLD(y ~ x, "Prentice")
b.fit.pr <- LaplacesDemon(ModelPrentice(), data, GIV(ModelPrentice(), data))
data     <- dataLD(y ~ x, "Stukel")
b.fit.s  <- LaplacesDemon(ModelStukel(),   data, GIV(ModelStukel(), data))

}




