allometry<-na.omit(read.csv("orchard_2010_2019 .csv"))
plot(Height~Width1, data=allometry,pch=19,col=allometry$Time)
abline(coef(lm(Height~Width1,data=allometry)))

allometry$Width<-apply(cbind(allometry$Width1,allometry$Width2),1,mean)

allometry<-allometry[which(allometry$Height>0),]

plot(Height~Width, data=allometry,pch=19,col=allometry$Time)
abline(coef(lm(Height~Width,data=allometry)))

library(bbmle)


linmod<-mle2(Width~dnorm(Intercept+Slope*Height,sd=exp(sigma)),data=allometry,
                         start=list(Intercept=30,Slope=1,sigma=1))


powmod<-mle2(Width~dnorm(Intercept*Height^Slope,sd=exp(sigma)),data=allometry,
             start=list(Intercept=1.3,Slope=0.5,sigma=1))


quadmod<-lm(Width~Height+I(Height^2),data=allometry)

AICtab(linmod,powmod,quadmod)

plot(Height~Width, data=allometry,pch=19,col=allometry$Time)
curve(-9.41+1.32*x+-0.0017541*x^2,add=T,col="green",lwd=2)
curve(1.3397*x^0.945062,add=T,col="red3",lwd=2)
abline(coef(lm(Height~Width1,data=allometry)),lwd=2)

mgam<-glm(Height~Width,family=Gamma(link="log"),data=allometry)
mgam2<-glm(Height~Width+I(Width^2),family=Gamma(link="log"),data=allometry)

curve(exp(3.62+1.337e-02*x+-3.604e-05*x^2),add=T,col="brown",lwd=4)
