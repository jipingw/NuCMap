\name{NuCMap-package}
\alias{NuCMap-package}
\alias{NuCMap}
\docType{package}
\title{
 An R package for chemical nucleosome mapping data analysis. 
}
\description{
\code{NuCMap} is an R package for analyzing \code{Nu}cleosome positioning signals in \code{C}hemical \code{Map}. This package is built upon a locally convoluted Poisson model proposed in Xi et al. (2012). The core of the package was written in Fortran and C++. \code{NuCMap} integrates eleven functions including \code{plotCUTS}, \code{peakDIST}, \code{trainTEMP1}, \code{estNCP1}, \code{callUNIQUE}, \code{callRED}, \code{trainTEMP4}, \code{estNCP4}, \code{plotAATT}, \code{estNCPcall} and \code{calOccup} to fulfill a complete analysis of chemical mapping data from data visualization, model training and diagnostics, parameter estimation to nucleosome map definition. \code{NuCMap} depends on the experiment data package \code{nucmapData}, which contains all the illustration data files.        
}
\details{
\tabular{ll}{
Package: \tab NuCMap\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-27\cr
License: \tab  GPL-2\cr
}
\code{plotCUTS}: R function for visualizing raw cleavage frequency in a specified region on Watson and Crick strands simultaneously.

\code{peakDIST}: R function for Crick-Watson cleavage peak-peak distance diagnostic.

\code{trainTEMP1}: R function for training one-template model.

\code{estNCP1}: R function for estimating NCP score for one-template model using EM algorithm.

\code{callUNIQUE}: R function for defining unique nucleosome map based on NCP score/noise ratio.

\code{callRED}: R function for defining redundant nucleosome map based on NCP score/noise ratio.

\code{trainTEMP4}: R function for training four-template model.

\code{estNCP4}: R function for estimating NCP score for four-template mdoel using EM algorithm.

\code{plotAATT}: R function for plotting AA/TT/AT/TA frequency against to the distance from nuclesome center.

\code{estNCPcall}: R function for estimating NCP score, NCP score/noise ratio and making unique nucleosome center calls.

\code{calOccup}: R function for calculating nucleosome occupancy genome-wide.
}
\author{
Liqun Xi, Qingyang Zhang, Kristin Brogaard, Bruce Lindsay, Jonathan Widom, Ji-Ping Wang

Maintainer: Ji-Ping Wang<jzwang@northwestern.edu>
}
\references{
Xi, L., Brogaard, K., Zhang, Q., Lindsay, B., Widom, J., Wang, J.-P. (2012), A locally convoluted Poisson cluster model for nucleosome positioning signals in chemical map, Submitted for publication

Brogaard, K., Xi, L., Wang, J.-P. and Widom, J. (2012), A base pair resolution map of nucleosome positions in yeast, Nature 486:496-501
}

\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## plotCUTS
wfile=system.file("extdata", "watson12.txt",package="nucmapData")
cfile=system.file("extdata", "crick12.txt",package="nucmapData")
plotCUTS(seqname="chrI",watsonfile=wfile,crickfile=cfile,startpos=
10000,endpos=12000)

## peakDIST
wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
peakDIST(seqname=c("chrI","chrII"),watsonfile=wfile,crickfile=cfile)

## trainTEMP1
wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
trainTEMP1(seqname=c("chrI","chrII"),watsonfile=wfile,crickfile=cfile)

## estNCP1
chrI=system.file("extdata", "chrI.fa",package="nucmapData")
wfile=system.file("extdata", "watson12.txt",package="nucmapData")
cfile=system.file("extdata", "crick12.txt",package="nucmapData")
estNCP1(seqname="chrI",genfile=chrI,watsonfile=wfile,crickfile=
cfile,temp1="default")

## callUNIQUE
NCP1=system.file("extdata","NCPscore.ratio_1temp.txt",package="nucmapData")
callUNIQUE(seqname=c("chrI","chrII"),estresults=NCP1)

## callRED
NCP1=system.file("extdata","NCPscore.ratio_1temp.txt",package="nucmapData")
callRED(seqname=c("chrI","chrII"),estresults=NCP1)

## trainTEMP4
wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
chrI=system.file("extdata","chrI.fa",package="nucmapData")
chrII=system.file("extdata","chrII.fa",package="nucmapData")
umap=system.file("extdata","UNIQUEcenters.txt",package="nucmapData")
trainTEMP4(seqname=c("chrI","chrII"),genfile=c(chrI,chrII),watsonfile=
wfile,crickfile=cfile,center=umap)

## estNCP4 
chrI=system.file("extdata", "chrI.fa",package="nucmapData")
wfile=system.file("extdata", "watson12.txt",package="nucmapData")
cfile=system.file("extdata", "crick12.txt",package="nucmapData")
temp4_file=system.file("extdata", "trainTEMP4result.txt",package="nucmapData")
estNCP4(seqname="chrI",genfile=chrI,watsonfile=wfile,crickfile=cfile,temp4=
temp4_file)

## plotAATT
chrI=system.file("extdata", "chrI.fa",package="nucmapData")
chrII=system.file("extdata", "chrII.fa",package="nucmapData")
umap=system.file("extdata", "UNIQUEcenters.txt",package="nucmapData")
plotAATT(seqname=c("chrI","chrII"),genfile=c(chrI,chrII),center=umap)

## estNCPcall
chrI=system.file("extdata","chrI.fa",package="nucmapData")
wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
estNCPcall(seqname="chrI",genfile=chrI,watsonfile=wfile,crickfile=cfile)

## calOccup
NCP4=system.file("extdata","NCPscore.ratio_4temp.txt",package="nucmapData")
chrI=system.file("extdata","chrI.fa",package="nucmapData")
chrII=system.file("extdata","chrII.fa",package="nucmapData")
rmap=system.file("extdata","REDcenters.txt",package="nucmapData")
calOccup(seqname=c("chrI","chrII"),estresults=NCP4,genfile=c(chrI,chrII),
rednufile=rmap)
}
}
