\name{callRED}
\alias{callRED}
\title{R function for defining redundant nucleosome map}
\description{This function defines unique nucleosome map.}
\usage{callRED(estresults,seqname,threshold)}

\arguments{
  \item{estresults}{one string for the path and name of the output file from function \code{estNCP1} or \code{estNCP4}.}
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in \code{estresults}. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in \code{estresults}.}
  \item{threshold}{A cutoff value to define redundant nucleosome map. All the positions where the NCP score/noise ratio exceeds this cutoff are defined as nucleosome centers. The default threshold value is the minimum NCP score/noise ratio from unique nucleosome map used in Brogaard et al. (2012)}
}
\value{\code{callRED} defines redundant nucleosome map based on NCP/noise ratio. The output file is named as "REDcenters.txt" and contains three columns:
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
callRED(seqname=c("chrI","chrII"),estresults=NCP1,threshold="default")
}
}
%\keyword{}
