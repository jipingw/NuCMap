\name{callUNIQUE}
\alias{callUNIQUE}
\title{R function for defining unique nucleosome map}
\description{This function defines unique nucleosome map.}
\usage{callUNIQUE(estresults,seqname)}

\arguments{
  \item{estresults}{one string for the path and name of the output file from function \code{estNCP1} or \code{estNCP4}.}
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in \code{estresults}. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in \code{estresults}.}
}
\value{\code{callUNIQUE} defines unique nucleosome map based on NCP score/noise ratio. The output file is named as "UNIQUEcenters.txt" and contains three columns:
\item{Position}{chromosome coordinate;}
\item{NCPscore}{estimated NCP score;}
\item{Ratio}{NCP score/noise ratio.}
}

\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

NCP1=system.file("extdata","NCPscore.ratio_1temp.txt",package="nucmapData")
callUNIQUE(seqname=c("chrI","chrII"),estresults=NCP1)
}
}
%\keyword{}
