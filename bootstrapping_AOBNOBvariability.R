library(boot)
library(readxl)
library(GoFKernel)
library("R.matlab")
setwd("~/git_KU/wwt_micol_pbe/model/")
# setwd("C:/Users/u0124080/Box Sync/Mijn documenten/PhD/WWT model Niccolo Totis/WWT model - shared/Data/02_2020")

# memory.limit(size=500000)

# load in the data
foldName="~/git_KU/wwt_micol_pbe/exp_data_miCol/"
fName="200207_200228_FISH_Izico_untreated_treated_Results.xlsx"
dataNOB <- read_xlsx(paste(foldName,fName,sep = ""), sheet = "untreated_red_NOB")
dataAOB <- read_xlsx(paste(foldName,fName,sep = ""), sheet = "untreated_blue_AOB")
dataAOBNOB <- read_xlsx(paste(foldName,fName,sep = ""), sheet = "untreated_AOB_NOB")

# select subset of the data (we are only interested if area microcolony > area of one bacterium)
dataNOB_subset <- subset(dataNOB$Area,dataNOB$Area >= (pi*0.5^2))
dataAOB_subset <- subset(dataAOB$Area,dataAOB$Area >= (pi*0.5^2))
dataAOBNOB_subset <- subset(dataAOBNOB$Area,dataAOBNOB$Area >= (pi*0.5^2))
dataAOBNOB_subset <- sapply(dataAOBNOB_subset, as.numeric)

# to suppress warnings put =-1
options(warn=0)

## trying to show that NOBs and AOBs are not different!
set.seed(123) #make the results reproducible
B = 10^4 # number of resamples
nColoniesNOB <- length(dataNOB_subset)
myBootstrapNOB <- matrix(sample(dataNOB_subset, size = B*nColoniesNOB, replace = TRUE), nrow = nColoniesNOB, ncol = B)
nColoniesAOB <- length(dataAOB_subset)
myBootstrapAOB <- matrix(sample(dataAOB_subset, size = B*nColoniesAOB, replace = TRUE), nrow = nColoniesAOB, ncol = B)

# calculate the difference in MEANS for each of the bootsamples
Boot.Diff.In.Means.AOB.NOB <- colMeans(myBootstrapNOB) - colMeans(myBootstrapAOB)
# make the confidence intervals
qL.means <-  quantile(Boot.Diff.In.Means.AOB.NOB, prob=0.025)
qH.means <- quantile(Boot.Diff.In.Means.AOB.NOB, prob=0.975)
# null hypothesis (H0): "the areas of AOB and NOB are the same"
# alternative hypothesis (HA): "the areas of AOB and NOB are different"
# since the qL and qH are crossing 0 (, H0 is rejected) the means of these 2 datasets are NOT statistically significantly different

# calculate the difference in standard deviation for each of the bootsamples
Boot.Diff.In.sd.AOB.NOB <- apply(myBootstrapNOB,2,sd) - apply(myBootstrapAOB,2,sd)
# make the confidence intervals
qL.sd <-  quantile(Boot.Diff.In.sd.AOB.NOB, prob=0.025)
qH.sd <- quantile(Boot.Diff.In.sd.AOB.NOB, prob=0.975)
# since the qL and qH are crossing 0 (, H0 is rejected) the sd of these 2 datasets are NOT statistically significantly different


## calculate the variability on the results (for AOB and NOB combined)
set.seed(123) #make the results reproducible
B = 10^4 # number of resamples
nColoniesAOBNOB <- length(dataAOBNOB_subset)
myBootstrapAOBNOB <- matrix(sample(dataAOBNOB_subset, size = B*nColoniesAOBNOB, replace = TRUE), nrow = nColoniesAOBNOB, ncol = B)

Boot.means.AOBNOB <- colMeans(myBootstrapAOBNOB)
Boot.sd.AOBNOB <- apply(myBootstrapAOBNOB,2,sd)

# make the confidence intervals
qL.AOBNOB.means <-  quantile(Boot.means.AOBNOB, prob=0.025)
qH.AOBNOB.means <- quantile(Boot.means.AOBNOB, prob=0.975)
qL.AOBNOB.sd <-  quantile(Boot.sd.AOBNOB, prob=0.025) 
qH.AOBNOB.sd <- quantile(Boot.sd.AOBNOB, prob=0.975)


# pdf = probability density function
# xnew <- seq(min(dataAOBNOB_subset),max(dataAOBNOB_subset),0.01)
xnew_A <- seq(min(dataAOBNOB_subset),15,0.01)
xnew_B <- seq(15,max(dataAOBNOB_subset),1)
xnew=c(xnew_A,xnew_B)


# 1 - create a vector at which the estimated pdf are evaluated
computeKDEstimate <- function(oneColumn) {
  KDE_AOBNOB <- density.reflected(oneColumn,lower=min(oneColumn), upper = max(oneColumn)
                                  ,kernel="gaussian")
  df <- approxfun(x = KDE_AOBNOB$x, y = KDE_AOBNOB$y, method = "linear", yleft = NaN, yright = NaN)
  f_KDE_AOBNOB <- df(xnew)
  return(f_KDE_AOBNOB)
}


Boot.pdf.AOBNOB <- apply(myBootstrapAOBNOB,2,computeKDEstimate) #pdf = probability density function

# 2 - evaluate the value of the pdf in all the B resamples
# 3 - compute the 5% and 95% percentile among them -> put them into 2 vectors v5 and v95, of length k

qL.AOBNOB.pdf <- list() 
qH.AOBNOB.pdf <- list() 


for(a in 1:length(xnew)){
  qL.AOBNOB.pdf[a] <- quantile(Boot.pdf.AOBNOB[a,],prob=0.025,na.rm=TRUE)
  qH.AOBNOB.pdf[a] <- quantile(Boot.pdf.AOBNOB[a,],prob=0.975,na.rm=TRUE)
}


# 4 - compute a pdf for v5 and v95 -> obtained 2 continuous pdfs that can be plotted above and below our experimental pdf

## not required for graphs
# qLfun <- approxfun(x = xnew, y = qL.AOBNOB.pdf, method = "linear")
# qHfun <- approxfun(x = xnew, y = qH.AOBNOB.pdf, method = "linear")
measured_pdf=density.reflected(dataAOBNOB_subset,lower=min(dataAOBNOB_subset), upper = max(dataAOBNOB_subset)
                  ,kernel="gaussian")

plot(measured_pdf, xlab="Area [?m?]",ylab="density [-]", ylim=c(0,0.35))
lines(x=xnew,y=qL.AOBNOB.pdf, col="red")
lines(x=xnew, y=qH.AOBNOB.pdf, col = "blue")
legend(36,0.35,legend=c("dataAOBNOB","qL","qH"),
             col=c("black","red","blue"), text.font=2,lty = 1:1,cex = 0.8)

out_pdf=data.frame(measured_pdf$x,measured_pdf$y)
colnames(out_pdf)=c("x","y")
writeMat('empirical_pdf.mat',x=out_pdf$x,y=out_pdf$y)

qL=unlist(qL.AOBNOB.pdf)
qH=unlist(qH.AOBNOB.pdf)
out_CIpdfs=data.frame(x=xnew,qL=qL,qH=qH)
writeMat('empirical_CIpdfs.mat',x=out_CIpdfs$x,qL=out_CIpdfs$qL,qH=out_CIpdfs$qH)




