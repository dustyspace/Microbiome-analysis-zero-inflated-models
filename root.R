root_tidy <- readRDS("/Users/zerui/Desktop/otu_filt_tidy.RDS")
plot(sort(root_tidy$Counts[1:100000]),ylab="Counts")

## Separate otu table based on drought/water control
map <- readRDS("/Users/zerui/Desktop/map.RDS")
otu.filt <- readRDS("/Users/zerui/Desktop/otu_filt_tidy.RDS")
a<-readRDS("/Users/zerui/Desktop/otu_filt.RDS")
ds.sampleid <- map$SampleID[map[,3]=="DS"]
wc.sampleid <- map$SampleID[map[,3]=="WC"]
root.ds <- subset(a,select = ds.sampleid)
root.wc <- subset(a,select = wc.sampleid)

rnamecheck<-function(t1,t2){
  n <- nrow(t1)
  temp <- 0
  for (i in 1:n){
    if( rownames(t1[i,]) != rownames(t2[i,])){
      temp = temp + 1
    }
  }
  return(temp)
}
rnamecheck(root.ds, root.wc)  
## = 0, ds/wc have the same number and order of otu

set.seed(1)
otu.sample <-otu.filt[ sample(nrow(otu.filt),500), ]
otu.sample$index<-as.numeric(gsub("\\D", "", otu.sample[,2]))  ## Extract num from drght.x
otu.sample<-cbind(otu.sample,map[otu.sample$index,3:6])  ## Map the rows to specific index
otu.sample.nona<-na.omit(otu.sample)
saveRDS(otu.sample.nona, file = "otu.sample.nona.rds")

## standard count data models
library(pscl)
otu.sample.nona$OTU_ID <- as.factor(otu.sample.nona$OTU_ID)
otu.sample.nona$Treatment <- as.factor(otu.sample.nona$Treatment)
otu.sample.nona$Compartment <- as.factor(otu.sample.nona$Compartment)
otu.sample.nona$Soil <- as.factor(otu.sample.nona$Soil)
otu.sample.nona$Cultivar <- as.factor(otu.sample.nona$Cultivar)
fm_pois  <- glm(Counts ~ OTU_ID + Treatment + Compartment + Soil + Cultivar,
                data = otu.sample.nona, family = poisson)
fm_qpois <- glm(Counts ~ OTU_ID + Treatment + Compartment + Soil + Cultivar,
                data = otu.sample.nona, family = quasipoisson)
library(MASS)
fm_nb <- glm.nb(Counts ~ OTU_ID + Treatment + Compartment + Soil + Cultivar,
                data = otu.sample.nona)

## with simple inflation 
otu <- a[,1:192]  ##Remove the drght.n with NA cultivar
otu.fit.temp <- data.frame("Counts"=as.numeric(otu[1,]))  ## Create a temp dataframe with Counts from OTU_i
otu.fit.temp$index<-1:192
otu.fit.temp <- cbind(otu.fit.temp, map[otu.fit.temp$index,3:6])
fm_zip  <- zeroinfl(Counts ~ as.factor(Treatment) + as.factor(Compartment) 
                    + as.factor(Soil) + as.factor(Cultivar)| 
                      1, dist="poisson", data = otu.fit.temp)
fm_zinb  <- zeroinfl(Counts ~ as.factor(Treatment) + as.factor(Compartment) 
                    + as.factor(Soil) + as.factor(Cultivar)| 
                      1, dist="negbin", data = otu.fit.temp)
fm_pois  <- glm(Counts ~ as.factor(Treatment) + as.factor(Compartment) 
                + as.factor(Soil) + as.factor(Cultivar),
                data = otu.fit.temp, family = poisson)
fm_nb <- glm.nb(Counts ~ as.factor(Treatment) + as.factor(Compartment) 
                + as.factor(Soil) + as.factor(Cultivar), data = otu.fit.temp)
fm_zip1  <- zeroinfl(Counts ~ as.factor(Soil)+ as.factor(Compartment)+ as.factor(Cultivar)| as.factor(Treatment), dist="poisson", data = otu.fit.temp)
> summary(fm_zip1)

Call:
  zeroinfl(formula = Counts ~ as.factor(Soil) + as.factor(Compartment) + 
             as.factor(Cultivar) | as.factor(Treatment), data = otu.fit.temp, 
           dist = "poisson")

Pearson residuals:
  Min       1Q   Median       3Q      Max 
-1.04202 -0.41267 -0.20713 -0.06055  9.83899 

Count model coefficients (poisson with log link):
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)                0.1449     0.3805   0.381 0.703376    
as.factor(Soil)B          -2.8170     1.0560  -2.667 0.007642 ** 
  as.factor(Soil)D           1.3756     0.3596   3.825 0.000131 ***
  as.factor(Compartment)ES  -2.1463     0.4393  -4.886 1.03e-06 ***
  as.factor(Cultivar)G2     -0.1040     0.3087  -0.337 0.736225    
as.factor(Cultivar)S1     -0.3525     0.2893  -1.219 0.223031    
as.factor(Cultivar)S2     -0.6881     0.3214  -2.141 0.032280 *  
  
  Zero-inflation model coefficients (binomial with logit link):
  Estimate Std. Error z value Pr(>|z|)
(Intercept)             -0.6052     0.5681  -1.065    0.287
as.factor(Treatment)WC   1.0273     0.6309   1.628    0.103
---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Number of iterations in BFGS optimization: 15 
Log-likelihood: -116.7 on 9 Df
> lrtest(fm_zip,fm_zip1)
Likelihood ratio test

Model 1: Counts ~ as.factor(Treatment) + as.factor(Compartment) + as.factor(Soil) + 
  as.factor(Cultivar) | 1
Model 2: Counts ~ as.factor(Soil) + as.factor(Compartment) + as.factor(Cultivar) | 
  as.factor(Treatment)
#Df  LogLik Df  Chisq Pr(>Chisq)    
1   9 -117.27                         
2   9 -116.74  0 1.0566  < 2.2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




#####################################
##       Multiple testing          ##
##   for DS/WC with other fixed    ##
#####################################
map <- readRDS("/Users/zerui/Desktop/map.RDS")
otu.filt<-readRDS("/Users/zerui/Desktop/otu_filt.RDS")
ds.wc <- map$SampleID[map$Compartment=="RS" & map$Soil=="D" & map$Cultivar=="G1"]

library(dplyr)
ds.wc.index <- gsub("\\D","",ds.wc)%>%as.numeric(.)
otu.ds.wc <- otu.filt[,ds.wc]
trt <- colnames(otu.ds.wc)%>%
                          gsub("\\D","",.)%>%
                          as.numeric(.)%>%
                          map[.,]%>%
                          .$Treatment
otu.ds.wc <- rbind(trt,otu.ds.wc)
otu.ds.wc.t<-as.data.frame(t(otu.ds.wc))
colnames(otu.ds.wc.t)[1]<-"Treatment"

fm.pois.1  <- glm(otu.ds.wc.t[,2] ~  as.factor(Treatment),
                data = otu.ds.wc.t, family = poisson)
pvalue.pois<-function(df){
  plist<-c()
  for(i in 2:ncol(df)){
    fm.pois<- glm(df[,i] ~  as.factor(Treatment),
                  data = df, family = poisson)
    result <- anova(fm.pois,test = "Chisq")
    plist<- append(plist, result$`Pr(>Chi)`[2])
  }
  return(plist)
}
p.pois <- pvalue.pois(otu.ds.wc.t)

detect.all.0<-function(df){
  non.0.list<-c()
  for(i in 2:ncol(df)){
    temp <- sum(df[,i])
    if(temp!=0){
      non.0.list <- append(non.0.list, i)
    }
  }
  return(non.0.list)
}
temp <- detect.all.0(otu.ds.wc.t)
temp <- append(1,temp)
otu.ds.wc.t.no0<-otu.ds.wc.t[,temp]
p.pois.no0<- pvalue.pois(otu.ds.wc.t.no0)
hist(p.pois.no0,breaks = 50,col = "black")

library(MASS)
fm.nb.1 <- glm.nb(otu.ds.wc.t[,2] ~  as.factor(Treatment),
                data = otu.ds.wc.t)
pvalue.nb<-function(df){
  plist<-c()
  for(i in 2:ncol(df)){
    fm.pois<- glm.nb(df[,i] ~  as.factor(Treatment),
                  data = df)
    result <- anova(fm.pois)
    plist<- append(plist, result$`Pr(>Chi)`[2])
  }
  return(plist)
}
library(pscl)
fm.zip  <- zeroinfl(otu.ds.wc.t[,2] ~ as.factor(Treatment)| 
                      1, dist="poisson", data = otu.ds.wc.t)
fm.zip.r <- zeroinfl(otu.ds.wc.t[,2] ~ 1| 
                      1, dist="poisson", data = otu.ds.wc.t)
pchisq(-2*(fm.zip.r$loglik-fm.zip$loglik),1)
fm.zip.1  <- zeroinfl(otu.ds.wc.t[,2] ~ as.factor(Treatment)| 
                      1, dist="poisson", data = otu.ds.wc.t)

check.zi<-function(df){
  zi.list<-c()
  for(i in 2:ncol(df)){
    temp <- min(df[,i])
    if(temp==0){
      zi.list <- append(zi.list, i)
    }
  }
  return(zi.list) 
}
temp1<-check.zi(otu.ds.wc.t.no0)
temp1<-append(1,temp1)
otu.ds.wc.t.no0.zi <- otu.ds.wc.t.no0[,temp1]


pvalue.zip<-function(df){
  plist<-c()
  for (i in 2:ncol(df)){
    fm.zip  <- zeroinfl(df[,i] ~ as.factor(Treatment)| 
                          1, dist="poisson", data = df)
    fm.zip.r <- zeroinfl(df[,i] ~ 1| 
                           1, dist="poisson", data = df)
    p.temp <- 1-pchisq(-2*(fm.zip.r$loglik-fm.zip$loglik),1)
    plist<-append(plist,p.temp)
  }
  return(plist)
}

fm.zinb.1  <- zeroinfl(otu.ds.wc.t.no0.zi[,2] ~ as.factor(Treatment)| 
                       1, dist="negbin", data = otu.ds.wc.t)
fm.zinb  <- zeroinfl(otu.ds.wc.t.no0.zi[,3] ~ as.factor(Treatment)| 
                       1, dist="negbin", data = otu.ds.wc.t)
fm.zinb.r <- zeroinfl(otu.ds.wc.t.no0.zi[,3] ~ 1| 
                       1, dist="negbin", data = otu.ds.wc.t)
1-pchisq(-2*(fm.zinb.r$loglik-fm.zinb$loglik),1)
pvalue.zinb<-function(df){
  plist<-c()
  for (i in 2:ncol(df)){
    fm.zinb  <- zeroinfl(df[,i] ~ as.factor(Treatment)| 
                          1, dist="negbin", data = df)
    fm.zinb.r <- zeroinfl(df[,i] ~ 1| 
                           1, dist="negbin", data = df)
    p.temp <- 1-pchisq(-2*(fm.zinb.r$loglik-fm.zinb$loglik),1)
    plist<-append(plist,p.temp)
  }
  return(plist)
}

a.p.pois.no0.zi<-p.adjust(p.pois.no0.zi,method = "BH")
a.p.nb.no0.zi<-p.adjust(p.nb.no0.zi,method = "BH")
a.p.zip.no0.zi<-p.adjust(p.zip.no0.zi,method = "BH")
a.p.zinb.no0.zi<-p.adjust(p.zinb.no0.zi,method = "BH")
p.value <- rbind(a.p.pois.no0.zi,a.p.nb.no0.zi,a.p.zip.no0.zi,a.p.zinb.no0.zi)
row.names(p.value)<-c("pois","nb","zip","zinb")

otu.sig.for.all<-which(apply(p.value,2,max)<=0.05)
otu.sig.for.pois <- which(p.value[1,]<=0.05)
otu.sig.for.nb <- which(p.value[2,]<=0.05)
otu.sig.for.zip <- which(p.value[3,]<=0.05)
otu.sig.for.zinb <- which(p.value[4,]<=0.05)

layout(matrix(c(1, 2,
                3, 4), nrow=2, byrow=TRUE))
hist(p.pois.no0.zi,breaks = 50,col = "black")
#abline(v=0.05,col="red",lwd=3,lty=2)
hist(p.nb.no0.zi,breaks = 50,col = "black")
hist(p.zip.no0.zi,breaks = 50,col = "black")
hist(p.zinb.no0.zi,breaks = 50,col = "black")

estobsk<-function(k){
  cbind(sum(dpois(k,fitted(fm.pois.1))),
        sum(dnbinom(k,mu=fitted(fm.nb.1),size=fm.nb.1$theta)),
        sum(predict(fm.zip.1,type="prob")[,1+k]),
        sum(predict(fm.zinb.1,type="prob")[,1+k]))
}
estobs0to10<-sapply(0:8,estobsk)
estobs11<-length(otu.ds.wc.t.no0.zi[,2])-apply(estobs0to10,1,sum)
estobs<-cbind(estobs0to10,estobs11)
hist(otu.ds.wc.t.no0.zi[,2],xlab="Number of OTU1",ylab="Absolute Frequency",breaks=seq(-0.5,9.5,by=1),main="",cex=1.2)
points(0:9-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
lines(0:9-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
points(0:9-0.05,estobs[2,],col=2,lwd=2,pch=3,lty=3)
lines(0:9-0.05,estobs[2,],col=2,lwd=2,pch=3,lty=3)
lines(0:9+0.05,estobs[3,],col=3,lwd=2,pch=1,lty=2)
points(0:9+0.05,estobs[3,],col=3,lwd=2,pch=1,lty=2)
points(0:9+0.1,estobs[4,],col=4,lwd=2,pch=4,lty=3)
lines(0:9+0.1,estobs[4,],col=4,lwd=2,pch=4,lty=3)

##########################################
##         GLM-Poisson Simulations      ##
##########################################
otu.pois.sig <- otu.ds.wc.t.no0.zi[,-1][otu.sig.for.pois]
otu.pois.sig <- cbind(trt,otu.pois.sig)
colnames(otu.pois.sig)[1]<-"Treatment"
(ncol(otu.pois.sig)-1)/7211
[1] 0.2084316
#fm.pois<- glm(otu.pois.sig[,2] ~  as.factor(Treatment),
#              data = otu.pois.sig, family = poisson)
all <- c(1:(ncol(otu.ds.wc.t.no0.zi[,-1])))
otu.pois.nonsig.index <- all[-otu.sig.for.pois]
otu.pois.nonsig <- otu.ds.wc.t.no0.zi[,-1][otu.pois.nonsig.index]
otu.pois.nonsig <- cbind(trt,otu.pois.nonsig)
colnames(otu.pois.nonsig)[1]<-"Treatment"

estbeta<-function(df){
  beta0list <- c()
  beta1list <- c()
  index <- c()
  for(i in 2:ncol(df)){
    fm.pois<- glm(df[,i] ~  as.factor(Treatment),
                  data = df, family = poisson)
    beta0 <- fm.pois$coefficients[1]
    beta1 <- fm.pois$coefficients[2]
    index <- append(index, i)
    beta0list <- append(beta0list, beta0)
    beta1list <- append(beta1list, beta1)
  }
  out <- rbind(beta0=beta0list, beta1=beta1list, index)
  colnames(out) <- c(1:ncol(out))
  return(out)
}
beta.pois.sig <- estbeta(otu.pois.sig)
beta.pois.nonsig <- estbeta(otu.pois.nonsig)

# num: the num of OTU's generation 
# beta: parameter form. row is assigned by different beta
sim.otu.glmpois <- function(num, beta){
  form <-c()
  for(i in 1:num){
    index <- sample(ncol(beta),1)
    b0 <- beta[1,index]
    b1 <- beta[2,index]
    mu <- c(exp(b0), exp(max(b0+b1,0)) )
    otu <- c(rpois(4, mu[1]), rpois(4,mu[2]))
   if (min(otu) == 0 && max(otu) > 0){
      otu <- append(otu, c(b0,b1,index))
      form <- rbind(form, otu)
    }
  }
  rownames(form) <- c(1:nrow(form))
  return(form)
}
sim1<-sim.otu.glmpois(2800,beta.pois.sig)
sim2<-sim.otu.glmpois(11000,beta.pois.nonsig)
colnames(sim1) <- c(rep("trt1",4),rep("trt2",4),"beta0","beta1","index of beta")
colnames(sim2) <- c(rep("trt1",4),rep("trt2",4),"beta0","beta1","index of beta")
otu.pois.sim.sig <- sim1[,-c(9:11)]
otu.pois.sim.nonsig <- sim2[,-c(9:11)]
colnames(otu.pois.sim.sig) <- c(rep(1,4),rep(2,4))
colnames(otu.pois.sim.nonsig) <- c(rep(1,4),rep(2,4))
otu.pois.sim.sig.tag <- cbind(otu.pois.sim.sig, 
                          tag = rep("de", ncol(otu.pois.sim.sig)) )
otu.pois.sim.nonsig.tag <- cbind(otu.pois.sim.nonsig, 
                             tag = rep("nde", ncol(otu.pois.sim.nonsig)))

pvalue.simpois<-function(df){
  plist<-c()
  for(i in 1:nrow(df)){
    fm.pois <- glm(df[i,] ~  as.factor(colnames(df)), family = poisson)
    result <- anova(fm.pois,test = "Chisq")
    plist<- append(plist, result$`Pr(>Chi)`[2])
  }
  return(plist)
}
p.simpois <- c(pvalue.simpois(otu.pois.sim.sig), pvalue.simpois(otu.pois.sim.nonsig))
a.p.simpois <- p.adjust(p.simpois, method = "BH")
fpr <- function(a.p, FP=seq(from=0,to=10325,by=100)){
  fprlist <- c()
  tprlist <- c()
  for(i in 1:length(FP)){
    FPR <- FP[i]/10325
    sig <- which(a.p < FPR)
    TP <- which(sig <=2293)
    TPR <- length(TP)/2293
    fprlist <- append(fprlist, FPR)
    tprlist <- append(tprlist, TPR)
  }
  return(rbind(fprlist, tprlist))
}
test <- fpr(a.p.simpois)
plot(test[1,],test[2,],pch=16,type = "o", xlab = "FPR", ylab = "TPR")
test1 <- fpr(a.p.simnb)
lines(test1[1,],test1[2,],pch=16,col="red",type="o")

# P=2293, N=10325
pvalue.simnb<-function(df){
  plist<-c()
  for(i in 1:nrow(df)){
    fm.nb <- glm.nb(df[i,] ~  as.factor(colnames(df)))
    result <- anova(fm.nb)
    plist<- append(plist, result$`Pr(>Chi)`[2])
  }
  return(plist)
}
p.simnb <- c(pvalue.simnb(otu.pois.sim.sig), pvalue.simnb(otu.pois.sim.nonsig))
a.p.simnb <- p.adjust(p.simnb, method = "BH")

pvalue.simzip<-function(df){
  plist<-c()
  for (i in 1:nrow(df)){
    fm.zip  <- zeroinfl(df[i,] ~ as.factor(colnames(df))| 
                          1, dist="poisson")
    fm.zip.r <- zeroinfl(df[i,] ~ 1| 
                           1, dist="poisson")
    p.temp <- 1-pchisq(-2*(fm.zip.r$loglik-fm.zip$loglik),1)
    plist<-append(plist,p.temp)
  }
  return(plist)
}
p.simzip <- c(pvalue.simzip(otu.pois.sim.sig), pvalue.simzip(otu.pois.sim.nonsig))
a.p.simzip <- p.adjust(p.simzip, method = "BH")
test2 <- fpr(a.p.simzip)
lines(test2[1,],test2[2,],pch=16,col="green",type="o")

pvalue.simzinb<-function(df){
  plist<-c()
  for (i in 1:nrow(df)){
    fm.zinb  <- zeroinfl(df[i,] ~ as.factor(colnames(df))| 
                           1, dist="negbin")
    fm.zinb.r <- zeroinfl(df[i,] ~ 1| 
                            1, dist="negbin")
    p.temp <- 1-pchisq(-2*(fm.zinb.r$loglik-fm.zinb$loglik),1)
    plist<-append(plist,p.temp)
  }
  return(plist)
}
p.simzinb <- c(pvalue.simzinb(otu.pois.sim.sig), pvalue.simzinb(otu.pois.sim.nonsig))
a.p.simzinb <- p.adjust(p.simzinb, method = "BH")
test3 <- fpr(a.p.simzinb)
lines(test3[1,],test3[2,],pch=16,col="blue",type="o")
legend("bottomright",legend=c("Poisson","Negbin","ZIP", "ZINB"),col=c("black","red","green","blue"),lwd=2)
##########################################
##           ZIP Simulations            ##
##########################################
otu.zip.sig <- otu.ds.wc.t.no0.zi[,-1][otu.sig.for.zip]
otu.zip.sig <- cbind(trt,otu.zip.sig)
colnames(otu.zip.sig)[1]<-"Treatment"
estpara <- function()
sim.otu.zip <- function(num, para){
  form <- c()
  
}
plot(0, xlim=c(1,length(beta0.pois.sim.sig)), ylim=c(-30,30))
for(i in 1:length(beta0.pois.sim.sig)){
  lines(x=rep(c(1:length(beta0.pois.sim.sig)),2),y=c(beta0.pois.sim.sig,beta1.pois.sim.sig) )
  points(x=rep(c(1:length(beta0.pois.sim.sig)),2),y=c(beta0.pois.sim.sig,beta1.pois.sim.sig) )
}
