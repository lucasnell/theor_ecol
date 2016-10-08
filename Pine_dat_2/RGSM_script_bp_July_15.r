### Read in datasets
setwd("C:/Users/billpine/Google Drive/Zero Inflated CPUE")
setwd('C:/Users/mattl/Documents/RGSM')
discharge = read.csv("ABQ_discharge.csv")
head(discharge)
catch=read.csv('catch by year station.csv')
head(catch)
yearly=read.csv('catch by year.csv')
head(yearly)



### Start analyses


## Observed frequency of occurrence at sites
catch$success=ifelse(catch$Count>0,1,0)
freq_occur=aggregate(success~Year,data=catch,mean)
Years=as.matrix(1993:2015)
colnames(Years)=c("Year")
freq_occur=merge(Years,freq_occur,all.x=TRUE,all.Y=TRUE)
freq_occur


## Estimate the distribution of the means using nonparametric bootstrap with replacement
catch$cpue=catch$Count/catch$Area
year_mean=function(year,n)
	{
	mean(sample(catch$cpue[catch$Year==year],n,replace=TRUE))
	}
boot_means=replicate(1000,sapply(1993:2015,function(t)year_mean(t,length(catch$cpue[catch$Year==t]))))
rownames(boot_means)=1993:2015
mean_cpue=sapply(1:length(boot_means[,1]),function(i)mean(boot_means[i,]))
se_cpue=sapply(1:length(boot_means[,1]),function(i)sd(boot_means[i,]))
cbind(1993:2015,mean_cpue,se_cpue)

boxplot(t(boot_means),col=8,outline=FALSE,xlab='Year',ylab='Mean CPUE')
#lines(yearly$CPUE_Yr,col='red',type='b')

layout(matrix(c(1:21),nrow=7,byrow=FALSE))
par(mai=c(0.3,0.25,0.2,0.1))
for(i in unique(catch$Year))
	{
	hist(boot_means[rownames(boot_means)==i],breaks=seq(0,1,0.05),xlim=c(0,0.8),col=8,xlab='',ylab='',main=i)
	}


## Scatterplots of Discharge duration versus freq_occurrence and mean CPUE
windows()
plot(discharge$ABQ_2000,freq_occur$success,cex=1.4,pch=18,xlab='Duration (days)',ylab='Frequency of Occurrence',main='Spring Discharge > 2000 cfs')
plot(discharge$ABQ_2000,mean_cpue,cex=1.4,col=4,pch=18,xlab='Duration (days)',main='Spring Discharge > 2000 cfs')
regress=lm(mean_cpue~discharge$ABQ_2000)
summary(regress)
lines(0:60,coefficients(regress)[2]*0:60+coefficients(regress)[1],lwd=2)
text(1,0.3,paste('CPUE = ',round(coefficients(regress)[2],4),'(Duration)+',round(coefficients(regress)[1],4),sep=''),pos=4)
text(1,0.28,paste('R-squared =',round(summary(regress)$r.squared,2)),pos=4)


## Plots of annual freq_occurrence and mean_CPUE overlaid on discharge duration
par(mai=c(0.8,0.8,0.6,0.8))
bar=barplot(discharge$ABQ_3000,names=1993:2015,xlab='Year',ylab='Duration (days)',main='Spring Discharge > 3000 cfs',ylim=c(0,70))
axis(1,at=bar,labels=FALSE)
par(new=T)
plot(freq_occur,typ='b',lwd=2,col=2,xaxt='n',yaxt='n',xlab='',ylab='',pch=15,ylim=c(0,1.1))
axis(side=4)
mtext('Frequency of Occurrence',4,line=2.5)
legend(2010,1.1,c('discharge','frequency'),bty='n',fill=c(8,NA),border=c(1,NA),pch=c(NA,15),col=c(NA,2))

par(mai=c(0.8,0.8,0.6,0.8))
bar=barplot(discharge$ABQ_3000,names=1993:2015,xlab='Year',ylab='Duration (days)',main='Spring Discharge > 3000 cfs',ylim=c(0,70))
axis(1,at=bar,labels=FALSE)
par(new=T)
plot(mean_cpue,typ='b',lwd=2,col=4,xaxt='n',yaxt='n',xlab='',ylab='',pch=17,ylim=c(0,0.35))
axis(side=4)
mtext(expression('Mean CPUE (Catch/'~m^{2}~')'),4,line=2.5)
legend(18,0.35,c('discharge','cpue'),bty='n',fill=c(8,NA),border=c(1,NA),pch=c(NA,17),col=c(NA,4))


### GLMs

data=merge(discharge,catch,all.y=TRUE)#####WORK ON THIS SOMETHING GOOFY WITH THE MERGE HERE, NOT SURE IF DPLYR BELOW IS RIGHT
head(data)

library("dplyr")
data <- left_join(discharge,catch, by='Year')




## Frequency of occurrence logistic regression, Example_3000 cfs

occur_glm=glm(success~ABQ_3000,data=data,family='binomial')
summary(occur_glm)

windows()
days=as.data.frame(1:60)
colnames(days)=c('ABQ_3000')
# pred_occur=predict(occur_glm,newdata=days,type='response',interval='confidence')
# plot(days[,1],pred_occur[,1],type='l',lwd=2,ylim=c(0,1),xlab='Duration (days)',ylab='Predicted Frequency of Occurrence',
#      main='Spring Discharge > 3000 cfs')
# lines(pred_occur$lwr,lty=2)
# lines(pred_occur$upr,lty=2)
# #points(data$ABQ_3000,data$success,pch=16,cex=0.8)
# 
# occurrence=merge(discharge,freq_occur,all.y=TRUE)
# points(occurrence$ABQ_3000,occurrence$success,pch=18,col=2)

#from old code on July 11
windows()
## Frequency of occurrence logistic regression, Example_3000 cfs
occur_glm=glm(success~ABQ_3000,data=data,family='binomial')
summary(occur_glm)
days=as.data.frame(1:70)
colnames(days)=c('ABQ_3000')
head(days)
days$freq_occur=predict(occur_glm,newdata=days,type='response')
days$se=predict(occur_glm,newdata=days,type='response',se.fit=TRUE)$se.fit
days$LowerLim=days$freq_occur-1.96*days$se
days$UpperLim=days$freq_occur+1.96*days$se
plot(days[,1:2],type='l',lwd=2,xlab='Duration (days)',ylab='Predicted Frequency of Occurrence',main='Spring Discharge > 3000 cfs')
lines(days$LowerLim,lty=2)
lines(days$UpperLim,lty=2)
occurrence=merge(discharge,freq_occur,all.y=TRUE)



## Let's test for the best model by ranking magnitude of discharge by the explained deviance in freq_occur
null_model=glm(success~1,data=data,family='binomial')
summary(null_model)
residual_dev_null=null_model$deviance

occur_1000=glm(success~ABQ_1000,data=data,family='binomial')
summary(occur_1000)
residual_dev_1000=occur_1000$deviance

occur_2000=glm(success~ABQ_2000,data=data,family='binomial')
summary(occur_2000)
residual_dev_2000=occur_2000$deviance

occur_3000=glm(success~ABQ_3000,data=data,family='binomial')
summary(occur_3000)
residual_dev_3000=occur_3000$deviance

occur_4000=glm(success~ABQ_4000,data=data,family='binomial')
summary(occur_4000)
residual_dev_4000=occur_4000$deviance

occur_5000=glm(success~ABQ_5000,data=data,family='binomial')
summary(occur_5000)
residual_dev_5000=occur_5000$deviance

dev_reduce=rbind(residual_dev_null,residual_dev_1000,residual_dev_2000,residual_dev_3000,residual_dev_4000,residual_dev_5000)
dev_reduce
Percent_reduce=(dev_reduce[1,]-dev_reduce)/dev_reduce[1,]
Percent_reduce

# Or this shortcut using the anova function, but the Deviance column is not correct because it is calculated from the previous model not the null
# really meant to be used with nested models, and so you should ignore that column, but still useful for extracting the residual deviances

anova(null_model,occur_1000,occur_2000,occur_3000,occur_4000,occur_5000)

# Looks like it's a close call between 2000 and 3000 cfs explaining about 24% of the annual variation in freq_occurrence
# Statistically not very different, but the biologists are saying 3000 cfs is the magic number
# If manager's choice, no doubt the 2000 cfs duration is the preferred number
# Also close is 1000 cfs explaining 21%, so what we learned here is that water is important for fish to be present, who knew?
# Note however that the discharges are not truly independent treatements, that is, many of the days above 2000 cfs were 
# actually above 3000 cfs, and some were above 4000 and 5000.  So this is really a statistical exercise.


## Testing for the best model relating magnitude of discharge by the explained deviance in positive catch rates

# Let's first look at the gamma fits to the annual cpue dists

positives=subset(data,cpue>0)
head(positives)
cpue_hist=function(year)
	{
	cpue=positives$cpue[positives$Year==year]
	hist=hist(cpue,breaks=20,col=8,freq=FALSE,xlab='cpue',cex.axis=0.8,mgp=c(3,0.5,0),main=year)
	theta=c(0.1,1)
	neg_LL=function(theta)
		{
		-sum(dgamma(cpue,theta[1],theta[2],log=TRUE))
		}
	fit=optim(theta,neg_LL)
	fit
	lines(seq(0.001,max(cpue),0.01),dgamma(seq(0.001,max(cpue),0.01),fit$par[1],fit$par[2]),col=2,lwd=2)
	text(tail(hist$breaks,1)[1],hist$density[1]*0.95,paste('k =',round(fit$par[1],2)),cex=0.8,pos=2)
	text(tail(hist$breaks,1)[1],hist$density[1]*0.8,paste('theta =',round(fit$par[2],2)),cex=0.8,pos=2)
	ks=ks.test(sort(cpue,decreasing=FALSE),"pgamma",fit$par[1],fit$par[2])
	text(tail(hist$breaks,1)[1],hist$density[1]*0.65,paste('p =',round(ks$p.value,3)),cex=0.8,pos=2)
	}

layout(matrix(c(1:21),nrow=7,byrow=FALSE))
par(mai=c(0.3,0.25,0.2,0.1))
sapply(unique(positives$Year),function(y)cpue_hist(y))

# Looks like the gamma distribution will do, we failed to reject the distribution of cpues as being significantly different from the fitted prob dists, KS_test, pvalues>0.05

# Moving forward with the Gamma glm

library(glmmADMB)
null_model=glmmadmb(cpue~1,data=positives,family='gamma',link='log')
summary(null_model)

cpue_1000=glmmadmb(cpue~ABQ_1000,data=positives,family='gamma',link='log')
summary(cpue_1000)

cpue_2000=glmmadmb(cpue~ABQ_2000,data=positives,family='gamma',link='log')
summary(cpue_2000)

cpue_3000=glmmadmb(cpue~ABQ_3000,data=positives,family='gamma',link='log')
summary(cpue_3000)

cpue_4000=glmmadmb(cpue~ABQ_4000,data=positives,family='gamma',link='log')
summary(cpue_4000)

cpue_5000=glmmadmb(cpue~ABQ_5000,data=positives,family='gamma',link='log')
summary(cpue_5000)

dev_reduce=-2*rbind(null_model$loglik,cpue_1000$loglik,cpue_2000$loglik,cpue_3000$loglik,cpue_4000$loglik,cpue_5000$loglik)  # Residual Deviance estimated as twice the negLL
dev_reduce
Percent_reduce=(dev_reduce[1,]-dev_reduce)/abs(dev_reduce[1,])
Percent_reduce

# Looks like 2000 cfs explains the deviance, but also close is 1000cfs, with less explained by the 3000 cfs model
# Perhaps we move forward with 2000 cfs as having support in both the freq_occur and pos_CPUE models
# Note however that the discharges are not truly independent treatment, that is, many of the days above 2000 cfs were 
# actually above 3000 cfs, and some were above 4000 and 5000.


### Gamma hurdle model for predicted mean CPUE with days above 2000 cfs as the covariate

occur_glm=glm(success~ABQ_2000,data=data,family='binomial')
summary(occur_glm)
pos_cpue_glm=glm(cpue~ABQ_2000,data=positives,family=Gamma(link='log'))
summary(pos_cpue_glm)

days = as.data.frame(0:60)
colnames(days) = c('ABQ_2000')
days$freq_occur = predict(occur_glm,newdata=days,type='response')
days$freq_se = predict(occur_glm,newdata=days,type='response',se.fit=TRUE)$se.fit
days$pos_cpue = predict(pos_cpue_glm,newdata=days,type='response')
days$cpue_se = predict(pos_cpue_glm,newdata=days,type='response',se.fit=TRUE)$se.fit
days$mean_cpue = days$freq_occur*days$pos_cpue 
days$mean_var = days$freq_occur^2*days$cpue_se^2 + days$pos_cpue^2*days$freq_se^2 - days$cpue_se^2*days$freq_se^2  # Goodman exact estimate of variance of two products
days$mean_se = sqrt(days$mean_var)
days$mean_LL = days$mean_cpue-1.96*days$mean_se
days$mean_UL = days$mean_cpue+1.96*days$mean_se

# Summary of the hurdle model predictions of frequency of occurrence, positive catch rates, overall mean catch rates, and standard deviations 
# as function of days duration > 2000cfs

days

windows()
plot(days$ABQ_2000,days$mean_cpue,type='l',lwd=2,ylim=c(0,max(mean_cpue,na.rm=TRUE)),xlab='Duration (days)',ylab='Predicted Mean CPUE (Fish/Area)',main='Spring Discharge > 2000 cfs')
lines(days$ABQ_2000,days$mean_LL,type='l',lty=2)			
lines(days$ABQ_2000,days$mean_UL,type='l',lty=2)
points(discharge$ABQ_2000,mean_cpue,pch=18,col=4)


### Negative binomial glm and zero-inflated NB, predicting the catch of RGSM with offset for Area sampled

nb_null=glmmadmb(Count~1+offset(log(Area)),data=data,family='nbinom')
summary(nb_null)

nb_1000=glmmadmb(Count~ABQ_1000+offset(log(Area)),data=data,family='nbinom')
summary(nb_1000)

nb_2000=glmmadmb(Count~ABQ_2000+offset(log(Area)),data=data,family='nbinom')
summary(nb_2000)

nb_3000=glmmadmb(Count~ABQ_3000+offset(log(Area)),data=data,family='nbinom')
summary(nb_3000)	

nb_4000=glmmadmb(Count~ABQ_4000+offset(log(Area)),data=data,family='nbinom')
summary(nb_4000)

nb_5000=glmmadmb(Count~ABQ_5000+offset(log(Area)),data=data,family='nbinom')
summary(nb_5000)

dev_reduce=-2*rbind(nb_null$loglik,nb_1000$loglik,nb_2000$loglik,nb_3000$loglik,nb_4000$loglik,nb_5000$loglik)
dev_reduce
Percent_reduce=(dev_reduce[1,]-dev_reduce)/abs(dev_reduce[1,])
Percent_reduce

# An interesting result compared to our frequentist linear regession and gamma hurdle model.  
# Here we see little of the deviance in Count of RGSM explained by the spring flow duration
# when we account for the effect of area sampled.  The flow duration was significant though.

# Testing the zero-inflated model

zi_nb_2000=glmmadmb(Count~ABQ_2000+offset(log(Area)),data=data,family='nbinom',zeroInflation=TRUE)
summary(zi_nb_2000)

library('bbmle')
AICtab(nb_2000,zi_nb_2000)

# Zero-inflation parameter did not improve the model.

# Predict out the Negative binomial
nb_days = as.data.frame(0:60)
colnames(nb_days) = c('ABQ_2000')
nb_days$Area=mean(data$Area,na.rm=T)
pred_Count=predict(nb_2000,newdata=nb_days,type='response',interval='confidence')
nb_days=cbind(nb_days,pred_Count)
nb_days$mean_cpue=nb_days$fit/nb_days$Area
nb_days$mean_LL=nb_days$lwr/nb_days$Area
nb_days$mean_UL=nb_days$upr/nb_days$Area

plot(nb_days$ABQ_2000,nb_days$mean_cpue,type='l',lwd=2,ylim=c(0,max(mean_cpue,na.rm=TRUE)),xlab='Duration (days)',ylab='Predicted Mean CPUE (Fish/Area)',main='Spring Discharge > 2000 cfs')
lines(nb_days$ABQ_2000,nb_days$mean_LL,type='l',lty=2)			
lines(nb_days$ABQ_2000,nb_days$mean_UL,type='l',lty=2)
points(discharge$ABQ_2000,mean_cpue,pch=18,col=4)


# Let's graph the predictions from our three approaches, 
# 1) simple linear regression, 2) gamma hurdle model of cpue, 
# 3) negative binomial model of count with effort offset

plot(discharge$ABQ_2000,mean_cpue,pch=18,lwd=2,xlab='Duration (days)',ylab='Mean CPUE (Fish/Area)',main='Spring Discharge > 2000 cfs')
lines(0:60,coefficients(regress)[2]*0:60+coefficients(regress)[1],lty=2)	

lines(days$ABQ_2000,days$mean_cpue,type='l',col=2)
lines(days$ABQ_2000,days$mean_LL,type='l',lty=2,col=2)			
lines(days$ABQ_2000,days$mean_UL,type='l',lty=2,col=2)

lines(nb_days$ABQ_2000,nb_days$mean_cpue,type='l',col=4)
lines(nb_days$ABQ_2000,nb_days$mean_LL,type='l',lty=2,col=4)			
lines(nb_days$ABQ_2000,nb_days$mean_UL,type='l',lty=2,col=4)

# Simple linear regression misses the mark
# The means are very similar between the gamma hurdle and the negative binomial
# The negative binomial appears to be much more accurate in the estimate of uncertainty around the mean


		
