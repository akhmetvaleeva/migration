d <- read.csv("d.csv", encoding = "UTF-8")
di <- read.csv("di.csv", encoding = "UTF-8")
mis <- is.na(d$vrppc[d$year>=2006]) | is.na(d$perc[d$year>=2006])
dimis <- di
dimis[mis,3:ncol(dimis)] <- NA
write.csv(dimis, file="di_noimpute.csv")

y <- "lcoefmigr"
xs <- c("lpop", "ipc", "davto", "dlong", "destprir", "infmort", "dlincome", "goods", "gini", "belowmin", "incgroup",
        "dlbadwater", "dlvrppc", "percmin", "dtotsp", "famsubs", "docload", "unemp")
f1 <- as.formula(paste0(y, "~", paste(xs,collapse="+")))
f2 <- as.formula(paste0(y, "~", paste(xs,collapse="+"), "+as.factor(id)+as.factor(year)"))

# Pairwise plots
for (i in xs) {
  pdf(paste0("img/pw_", i, ".pdf"),7,7)
  plot(di[,i], di$lcoefmigr, xlab=i, ylab="Coef. migr.", cex.lab=1.5)
  dev.off()
}

mod1 <- lm(f1, data = di)
summary(mod1)

library(plm)
mod2 <- plm(f1, data=di, index=c("region", "year"), model="random")
mod3 <- plm(f1, data=di, index=c("region", "year"), model="within", effect="individual")
mod4 <- plm(f1, data=di, index=c("region", "year"), model="within", effect="time")
mod5 <- plm(f1, data=di, index=c("region", "year"), model="within", effect="twoways")

modt1 <- lm(f2, data=di)
modt2 <- lm(f2, data=dimis)
summary(modt1)
summary(modt2)

phtest(mod3, mod2)
phtest(mod4, mod2)
phtest(mod5, mod2)

library(multiwayvcov)
ex <- function(x) sqrt(diag(x))
s1 <- ex(cluster.vcov(mod1, di$region))
s2 <- ex(plm::vcovHC(mod2, cluster="group"))
s3 <- ex(plm::vcovHC(mod3, cluster="group"))
s4 <- ex(plm::vcovHC(mod4, cluster="group"))
s5 <- ex(plm::vcovHC(mod5, cluster="group"))


summary(mod1)
round(1-var(mod1$residuals)/var(di$lcoefmigr),3)
round(1-var(mod2$residuals)/var(di$lcoefmigr),3)
round(1-var(mod3$residuals)/var(di$lcoefmigr),3)
round(1-var(mod4$residuals)/var(di$lcoefmigr),3)
round(1-var(mod5$residuals)/var(di$lcoefmigr),3)


library(stargazer)
stargazer(mod1, mod2, mod3, mod4, mod5, out="naive.tex", report="vc*", se = list(s1, s2, s3, s4, s5))

#############
a <- read.csv("dij2.csv")
a2 <- 1/a
a2[a2==Inf] <- 0
write.csv(a2, "dij3.csv", row.names = FALSE)

aa <- a/1000
summary(aa[aa>0])
a2 <- exp(-a/1000)
diag(a2) <- 0
write.csv(a2, "dij4.csv", row.names = FALSE)

a <- read.table("neighb.csv", sep=",", row.names = 1, header = TRUE, encoding = "UTF-8")
for (i in 1:82) {
  if (sum(a[i,])==0) print(i)
}

a2 <- read.csv("neighb2.csv")
which(diag(as.matrix(a2))==1)

attach(di)

pdf("img/scatterpop.pdf", 7,4.5)
plot(pop, lcoefmigr, main="Pairwise scatterplot without logs for pop")
dev.off()
pdf("img/scatterbadwater.pdf", 7,4.5)
plot(badwater, lcoefmigr, main="Pairwise scatterplot without logs for badwater")
dev.off()
pdf("img/scatterincome.pdf", 7,4.5)
plot(income, lcoefmigr, main="Pairwise scatterplot without logs for income")
dev.off()
pdf("img/scattervrppc.pdf", 7,4.5)
plot(vrppc, lcoefmigr, main="Pairwise scatterplot without logs for vrppc")
dev.off()

pdf("img/scatterlpop.pdf", 7,4.5)
plot(lpop, lcoefmigr, main="Pairwise scatterplot without logs for pop")
dev.off()
pdf("img/scatterlbadwater.pdf", 7,4.5)
plot(lbadwater, lcoefmigr, main="Pairwise scatterplot without logs for badwater")
dev.off()
pdf("img/scatterlincome.pdf", 7,4.5)
plot(lincome, lcoefmigr, main="Pairwise scatterplot without logs for income")
dev.off()
pdf("img/scatterlvrppc.pdf", 7,4.5)
plot(lvrppc, lcoefmigr, main="Pairwise scatterplot without logs for vrppc")
dev.off()

### Moran Index Plot
moran <- read.csv("morans.csv")
m1 <- moran$m[seq(1,nrow(moran), 3)]
m2 <- moran$m[seq(2,nrow(moran), 3)]
m3 <- moran$m[seq(3,nrow(moran), 3)]
s1 <- moran$sd[seq(1,nrow(moran), 3)]
s2 <- moran$sd[seq(2,nrow(moran), 3)]
s3 <- moran$sd[seq(3,nrow(moran), 3)]

pdf("img/moran1.pdf", 7, 5)
plot(2006:2014, m1, main="Moran's I index for 1/d", ylim = range(m1+s1, m1-s1), xlab="Year", ylab="I")
for (i in 1:9) {
  lines(rep(2005+i,2), c(m1[i]-s1[i], m1[i]+s1[i]))
}
abline(h=c(0, moran[1,2]), lty=c(1,2))
dev.off()

pdf("img/moran2.pdf", 7, 5)
plot(2006:2014, m2, main="Moran's I index for exp(-ad)", ylim = range(m2+s2, m2-s2), xlab="Year", ylab="I")
for (i in 1:9) {
  lines(rep(2005+i,2), c(m2[i]-s2[i], m2[i]+s2[i]))
}
abline(h=c(0, moran[2,2]), lty=c(1,2))
dev.off()

pdf("img/moran3.pdf", 7, 5)
plot(2006:2014, m3, main="Moran's I index for d=[0,1]", ylim = range(m3+s3, m3-s3), xlab="Year", ylab="I")
for (i in 1:9) {
  lines(rep(2005+i,2), c(m3[i]-s3[i], m3[i]+s3[i]))
}
abline(h=c(0, moran[3,2]), lty=c(1,2))
dev.off()

m <- t(cbind(m1, m2, m3))
m

pdf("img/d_lcoefmigr.pdf", 7, 4)
plot(density(di$lcoefmigr), main="Log of coefmigr --- dependent variable")
rug(di$lcoefmigr)
dev.off()
pdf("img/d_pop.pdf", 7, 4)
plot(density(di$lpop), main="Log of pop")
dev.off()
pdf("img/d_avto.pdf", 7, 4)
plot(density(di$davto), main="Difference of avto")
dev.off()
pdf("img/d_long.pdf", 7, 4)
plot(density(di$dlong), main="Difference of long")
dev.off()
pdf("img/d_estprir.pdf", 7, 4)
plot(density(di$estprir), main="Difference of estprir")
dev.off()
pdf("img/d_income.pdf", 7, 4)
plot(density(di$dlincome), main="Difference of Log of income")
dev.off()
pdf("img/d_badwater.pdf", 7, 4)
plot(density(di$dlbadwater), main="Difference of Log of badwater")
dev.off()
pdf("img/d_vrppc.pdf", 7, 4)
plot(density(di$dlvrppc), main="Difference of Log of vrppc")
dev.off()
pdf("img/d_totsp.pdf", 7, 4)
plot(density(di$dtotsp), main="Difference of totsp")
dev.off()

# cat(paste0(" - b", vars, "*", "bc(", vars, ", theta)"), sep="\n")
# Selecting the Box-Cox transformation coefficient
y1 <- "lcoefmigr"
f3 <- as.formula(paste0(y1, "~ pop + ipc + avto + long + estprir + infmort + income + goods + gini + belowmin +
                        incgroup + badwater + vrppc + percmin + totsp + famsubs + docload + unemp"))
modstart <- lm(f3, data=di)
summary(modstart)
sc <- modstart$coefficients
names(sc) <- paste0("b", names(sc))
names(sc)[1] <- "b0"
sa <- c(0, sd(modstart$residuals), 0.5)
names(sa) <- c("mu", "sigma", "lambda")
sc <- c(sc, sa)

i <- 0 # For traceback
vals <- list()
LL <- function(b0, bpop, bipc, bavto, blong, bestprir, binfmort, bincome, bgoods, bgini, bbelowmin, bincgroup,
               bbadwater, bvrppc, bpercmin, btotsp, bfamsubs, bdocload, bunemp, mu, sigma, lambda) {
  i <<- i+1
  v <- c(b0, bpop, bipc, bavto, blong, bestprir, binfmort, bincome, bgoods, bgini, bbelowmin, bincgroup,
         bbadwater, bvrppc, bpercmin, btotsp, bfamsubs, bdocload, bunemp, mu, sigma, lambda)
  names(v) <- names(sc)
  vals[[i]] <<- v
  R <- lcoefmigr - b0 - bpop*bc(pop, lambda) - bipc*ipc - bavto*avto - blong*long - bestprir*estprir
  - binfmort*infmort - bincome*income - bgoods*goods - bgini*gini - bbelowmin*belowmin
  - bincgroup*incgroup - bbadwater*badwater - bvrppc*vrppc - bpercmin*percmin
  - btotsp*totsp - bfamsubs*famsubs - bdocload*docload - bunemp*unemp
  RR <- suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
  -sum(RR)
}

library(stats4)
fit <- mle(LL, start=as.list(sc), method="BFGS", control=list(trace=2, maxit=500))
summary(fit)

y1 <- "lcoefmigr"
l <- 1
f3 <- as.formula(paste0(y1, "~ I(pop^l) + ipc + avto + long + estprir + infmort + income + goods + gini + belowmin +
                        incgroup + badwater + vrppc + percmin + totsp + famsubs + docload + unemp"))
Rs <- rep(NA, 100)
for (l in (100:1)/100) {
  f3 <- as.formula(paste0(y1, "~ I(pop^l) + ipc + avto + long + estprir + infmort + income + goods + gini + belowmin +
                        incgroup + badwater + vrppc + percmin + totsp + famsubs + docload + unemp"))
  mod <- lm(f3)
  Rs[i] <- sum(mod$residuals^2)
}


# Methodology change
yy <- di$lcoefmigr
yi <- rep(NA, length(yy))
yt <- rep(NA, length(yy))
for (i in unique(di$id)) {
  yi[di$id==i] <- mean(yy[di$id==i])
}
for (t in unique(di$year)) {
  yt[di$year==t] <- mean(yy[di$year==t])
}
yw <- yy - yi - yt
y10 <- yw[di$year==2010]
y11 <- yw[di$year==2011]
y12 <- yw[di$year==2012]
de1 <- density(y10)
de2 <- density(y11)
de3 <- density(y12)
pdf("img/method.pdf", 9, 2.5)
par(mar=c(2,2,2,1))
plot(de1, xlim=range(de1$x, de2$x, de3$x), lwd=2, lty=1, main="Density of Within(i, t) lcoefmigr")
lines(de2, lwd=2, lty=2)
lines(de3, lwd=2, lty=3)
legend("topright", c("2010", "2011", "2012"), lwd=2, lty=1:3)
dev.off()

pdf("img/method2.pdf", 9, 2.5)
par(mar=c(2,2,2,1))
plot(ecdf(y10), main="ECDF of log(coefmigr)")
lines(ecdf(y11), col="red")
lines(ecdf(y12), col="blue")
legend("topleft", c("2010", "2011", "2012"), lwd=2, col=c("black", "red", "blue"))
dev.off()
ks.test(y10, y11)
ks.test(y11, y12)
ks.test(y10, y12)

res <- read.csv("residuals.csv")$er
r10 <- res[di$year==2010]
r11 <- res[di$year==2011]
r12 <- res[di$year==2012]
de1 <- density(r10)
de2 <- density(r11)
de3 <- density(r12)
pdf("img/method3.pdf", 9, 4.5)
plot(de1, xlim=range(de1$x, de2$x, de3$x), ylim=range(de1$y, de2$y, de3$y), lwd=2, lty=1, main="Density of residuals")
lines(de2, lwd=2, lty=2)
lines(de3, lwd=2, lty=3)
legend("topright", c("2010", "2011", "2012"), lwd=2, lty=1:3)
dev.off()
ks.test(r10, r11)
ks.test(r11, r12)
ks.test(r10, r12)
