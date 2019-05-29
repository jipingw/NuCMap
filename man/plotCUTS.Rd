\name{plotCUTS}
\alias{plotCUTS}
\title{R function for visualizing the cleavage frequency on Watson and Crick strands}
\description{This function produces cleavage frequency plot in any specified region on Watson and Crick strands simultaneously.
}
\usage{plotCUTS(seqname,watsonfile,crickfile,startpos,endpos)}

\arguments{
  \item{seqname}{the specified chromosome's name, e.g. "chrI". This sequence name should be consistant with the chromosome's name in Watson and Crick cleavage files.}
  \item{watsonfile}{one string for the path and name of Watson cleavage file. The file may contain two or more chormosomes. The first columm is chromosome name, the second column is the chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see the requirement for \code{watsonfile}.}
  \item{startpos}{an integer standing for the starting position of the genome region.}
  \item{endpos}{an integer standing for the ending position of the genome region.}
}
\value{\code{plotCUTS} produces cleavage frequency plot in a specified region on both Watson and Crick strands simultaneously.
}


\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

wfile=system.file("extdata", "watson12.txt",package="nucmapData")
cfile=system.file("extdata", "crick12.txt",package="nucmapData")
plotCUTS(seqname="chrI",watsonfile=wfile,crickfile=cfile,startpos=
10000,endpos=12000)
}
}

%\keyword{}

