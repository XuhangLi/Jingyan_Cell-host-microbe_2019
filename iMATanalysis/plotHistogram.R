#--------- 1-------------(this is aborted for final publication)
myredoxflux_percentatage=read.table("./myredoxflux_percentatage.csv",header = F, sep=',')
myredoxrxnsPercentage=read.table("./myredoxrxnsPercentage.csv",header = F, sep=',')

cv_out = t(myredoxflux_percentatage)
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.14,0.005),font.lab = 2,cex.lab = 1.3,
     xlab = "Relative Total Flux of Quinol-related ROS Reactions", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.0490,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.0490)))/(nrow(cv_out)+1)
p1

cv_out = t(myredoxrxnsPercentage)
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.6,0.025),font.lab = 2,cex.lab = 1.3,
     xlab = "Number of Quinol-related ROS Reactions Carrying Flux (%)", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.4500,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.4500)))/(nrow(cv_out)+1)
p1

#--------- 2-------------
ROSflux_percentatage=read.table("./ROSflux_percentatage.csv",header = F, sep=',')
cv_out = ROSflux_percentatage
pdf("figure7C.pdf")
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.19,0.005),font.lab = 2,cex.lab = 1.3,
     xlab = "Relative Total Flux of ROS Reactions (% of Tatal Flux)", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.0704,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.0704)))/(nrow(cv_out)+1)
p1
dev.off()
#--------- 3-------------
pdf("figure7D.pdf")
OXPHOSflux_percentatage=read.table("./OXPHOSflux_percentatage.csv",header = F, sep=',')
cv_out = OXPHOSflux_percentatage
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.26,0.005),font.lab = 2,cex.lab = 1.3,
     xlab = "Relative Total Flux of Oxidative Phosphorylation Reactions", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.1179765,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.1179765)))/(nrow(cv_out)+1)
p1
dev.off()
#--------- 4-------------
Folateflux_percentatage=read.table("./Folateflux_percentatage.csv",header = F, sep=',')
cv_out = Folateflux_percentatage
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.025,0.0003),font.lab = 2,cex.lab = 1.3,
     xlab = "Relative Total Flux of Folate Metabolism", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.001788,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.001788)))/(nrow(cv_out)+1)
p1

#--------- 5-------------
ionMetflux_percentatage=read.table("./ionMetflux_percentatage.csv",header = F, sep=',')
cv_out = t(ionMetflux_percentatage)
par(lwd=2.5)
hist(cv_out[,1], breaks = seq(0,0.07,0.002),font.lab = 2,cex.lab = 1.3,
     xlab = "Relative Total Flux of Inorganic Ion Transport and Metabolism", main = '',freq = F,ylab = 'Probability Density',axes = F)
Axis(side = 1, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
Axis(side = 2, cex.axis = 1.2, lwd.ticks = 2, lwd =2)
par(opar)
abline(v = 0.0200,lwd = 3,col = "red")
p1 = (1+sum(as.numeric(cv_out[,1] >= 0.0200)))/(nrow(cv_out)+1)
p1
