
#### this code is adapted from the U2OS_ToxCast screen folder (11/18/20) ####


####################################################### global Mahalanobis ######################################################

source("../R functions/globalMahalanobisDistances_v3.R")

Table1 = "htpp_well_norm" # need to read this table in

Output = globalMahalanobisDistances(Table1=Table1, 
                                    coverVariance=0.95, minObjects=100, SType = "vehicle control",  url = url)

CumProportion = Output$CumProportion

#How many PC are needed to explain  at least x% of variance?
PC90=length(which(CumProportion<0.90))+1 #9
PC95=length(which(CumProportion<0.95))+1 #15
PC99=length(which(CumProportion<0.99))+1 #100

png("U2OS_Pilot htpp_global_mah - Proportion of variance v201125.png", width=8, height=6, units="in", res=144)
    plot(x=1:1300, y=CumProportion, col="gray50", pch=19, cex=0.5, type="p",
         ylim=c(0,1), xlab="# of components", ylab="Proportion of variance retained", main=paste("Principal components of U2OS_ToxCast data"))
    #horizontal part
    segments(x0=30, y0=0.90, x1 = PC90, col="blue", lty='dashed')
    segments(x0=30, y0=0.95, x1 = PC95, col="blue", lty='solid', lwd=2)
    segments(x0=30, y0=0.99, x1 = PC99, col="blue", lty='dotted')
    #vertical part
    segments(x0=PC90, y0=0.1, y1 = 0.90, col="blue", lty='dashed')
    segments(x0=PC95, y0=0.1, y1 = 0.95, col="blue", lty='solid', lwd=2)
    segments(x0=PC99, y0=0.1, y1 = 0.99, col="blue", lty='dotted')
    
    text(x=c(PC90, PC95, PC99), y=0.05, labels=c(PC90, PC95, PC99), srt=90)
    text(x=0, y=c(0.9, 0.95, 0.99),  labels=paste0(c(90,95,99), "%"), cex=0.7)
dev.off()


#--> complete