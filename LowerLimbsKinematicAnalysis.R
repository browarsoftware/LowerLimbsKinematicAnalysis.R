#Starting file for the analysis
#Change path here
source("f:\\mocap_real_sizes\\karate\\sr1\\renderactor.R")
source("f:\\mocap_real_sizes\\karate\\sr1\\rotatedata.R")
source("f:\\mocap_real_sizes\\karate\\sr1\\calculatekinematic.R")
source("f:\\mocap_real_sizes\\karate\\sr1\\dtwcomparator.R")

################################################
#LOAD DATA
################################################
#Change path here
#Data can be downloaded from http://gdl.org.pl
#just load two recordings of the same technique
refdata <- read.csv("e:\\Publikacje\\kine\\data\\2016-12-22 ShoriunRiu MP\\mawashi_geri_right_iter1000.csv")
inputdata <- read.csv("f:\\mocap\\karate\\2016-04-01 ShoriunRiu 1\\mawashi_geri_right\\segmented\\sample5.csv")
#left or right, depending of foot that remains on the ground, if you kick with right leg it should be left
leg <- "Left"
#inputdata <- read.csv("f:\\mocap\\karate\\2016-04-01 ShoriunRiu 1\\mawashi_geri_left\\segmented\\sample0.csv")
#refdata <- read.csv("f:\\mocap\\karate\\2016-04-01 ShoriunRiu 3\\mawashi_geri_left\\segmented\\sample0.csv")
#leg <- "Right"

inputdataalignment <- rotatedata(inputdata, refdata, paste(leg, "Foot.Dx", sep = ""),paste(leg, "Foot.Dz", sep = ""),
                         "Hips.Dx","Hips.Dz")
inputdataalignmentkinematic <- calculatekinematic(inputdataalignment, paste(leg, "Foot", sep = ""))
refdatakinematic <- calculatekinematic(refdata, paste(leg, "Foot", sep = ""))
inputdataalignmentkinematic <- aligninputandrefdata(inputdataalignmentkinematic, refdatakinematic, limbname = paste(leg, "Foot", sep = ""))
inputdataalignmentkinematicf <- generateFeatures(inputdataalignmentkinematic)
refdatakinematicf <- generateFeatures(refdatakinematic)

################################################
#BEGIN ANALYSIS
################################################
extremumtreshold <- 0.66
smoothSize <- 0.1
analyzerange <- ceiling(nrow(inputdataalignment) * smoothSize * 1) #sprawdziæ, jaki powinien byæ ten przedzia³!

footddf <- analyzedta(refdatakinematic, inputdataalignmentkinematic,
                      refdatakinematicf$RightFoot, inputdataalignmentkinematicf$RightFoot, FUN=euc.dist, smoothSize = smoothSize, 
                      c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), plottitle = "RightFoot D", extremumtreshold = extremumtreshold, smoothSizeHelper = 0.1, 
                      whattodraw = "RightFoot", plotrgl = TRUE)

#plotsmoothingresults(footddf,"footddf", plotifnoextreams = TRUE, plotsmoothed = TRUE)

dftoanalyze <- tresholdresults(footddf, extremumtreshold)


plot(footddf$path1, footddf$path2, xlab = "Reference signal [sample id]", ylab = "Input signal [sample id]", main ="DTW alignment", type='l')


plot(footddf$data, xlab = "Time [10^-100 s]", ylab = "Distance [cm]", main ="DTW alignment function", col = 'black', type='l')
#lines(footddf$data, xlab = "Time [s^-100]", ylab = "Distance [cm]", main ="DTW alignment function")
lines(footddf$smoothdata, col = 'purple', lty = 2)
#points(footddf$smoothdata, col = 'blue')
for (a in 1:length(footddf$extremumid))
{
  points(footddf$extremumid[a], footddf$smoothdata[footddf$extremumid[a]], col = "cyan",pch = 1, cex = 2, lwd=3)
}
for (a in 1:length(dftoanalyze$extremumid))
{
  points(dftoanalyze$extremumid[a], dftoanalyze$smoothdata[dftoanalyze$extremumid[a]], col = "red",pch = 4, cex = 2, lwd=3)
}
legend(x= "topright", y=0.92, legend=c("Original", "Smoothed", "Maxima", "Maxima over treshold"), col=c("black", "purple", 'cyan','red'), lty=c(1,2,1), cex=0.8)


#this might render quite long, but looks nice :-)
rglplotanalyzedata(refdatakinematic, inputdataalignmentkinematic, refdatakinematicf$RightFoot, inputdataalignmentkinematicf$RightFoot, footddf$path1, footddf$path2,resultdata =   dftoanalyze, whattodraw = 'RightFoot')




plotsmoothingresults(dftoanalyze,"DTWaf of right foot trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE)
###################Hips
HipsX <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListHipsX, inputdataalignmentkinematicf$ListHipsX, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   plottitle = "DTWaf of right hips X angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "HipsX")
HipsX <- tresholdresults(HipsX, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, HipsX, analyzerange),"DTWaf of right hips X trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

HipsY <-analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListHipsY, inputdataalignmentkinematicf$ListHipsY, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   plottitle = "DTWaf of right hips Y angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "HipsY")
HipsY <- tresholdresults(HipsY, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, HipsY, analyzerange),"DTWaf of right hips Y trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

HipsZ <-analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListHipsZ, inputdataalignmentkinematicf$ListHipsZ, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   plottitle = "DTWaf of right hips Y angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "HipsZ")
HipsZ <- tresholdresults(HipsZ, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, HipsZ, analyzerange),"DTWaf of right hips Y trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

##################RightThigh

RightThighX <-analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListRightThighX, inputdataalignmentkinematicf$ListRightThighX, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                         plottitle = "DTWaf of right Thigh X angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "ThighX")
RightThighX <- tresholdresults(RightThighX, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightThighX, analyzerange),"DTWaf of right Thigh X trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightThighY <-analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListRightThighY, inputdataalignmentkinematicf$ListRightThighY, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                         plottitle = "DTWaf of right Thigh Y angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "HipsY")
RightThighY <- tresholdresults(RightThighY, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightThighY, analyzerange),"DTWaf of right Thigh Y trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightThighZ <-analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListRightThighZ, inputdataalignmentkinematicf$ListRightThighZ, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                         plottitle = "DTWaf of right Thigh Z angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "HipsZ")
RightThighZ <- tresholdresults(RightThighZ, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightThighZ, analyzerange),"DTWaf of right Thigh Z trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

##################Leg

kneeadf <- analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListRightKnee, inputdataalignmentkinematicf$ListRightKnee, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                      plottitle = "DTWaf of right knee angle analysis", path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, smoothSizeHelper = 0.1, whattodraw = "RightLeg")
kneeadf <- tresholdresults(kneeadf, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, kneeadf, analyzerange),"DTWaf of right knee angle analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

######################################################################

