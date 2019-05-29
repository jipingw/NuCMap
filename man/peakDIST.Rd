\name{peakDIST}
\alias{peakDIST}
\title{R function for Crick-Watson cleavage peak-peak distance diagnostic.}
\description{This function produces frequency plot for cross-strand cleavage peak-peak distance. 
}
\usage{peakDIST(seqname,watsonfile,crickfile)}

\arguments{
  \item{seqname}{the default value is "all", which speficies all chromosomes listed in Watson and Crick cleavage files. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in Watson and Crick cleavage files.}
  \item{watsonfile}{one string for the path and name of Watson cleavage file. The file may contain two or more chormosomes. The first columm is chromosome name, the second column is the chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see requirement for \code{watsonfile}.}
}
\value{\code{peakDIST} outputs frequency plot for cross-strand cleavage peak-peak distance. 
}


\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
peakDIST(seqname=c("chrI","chrII"),watsonfile=wfile,crickfile=cfile)
}
}
%\keyword{}

