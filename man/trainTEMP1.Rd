\name{trainTEMP1}
\alias{trainTEMP1}
\title{R function for training one-template model}
\description{This function trains cleavage template for one-template model.
}
\usage{trainTEMP1(seqname,watsonfile,crickfile)}

\arguments{
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in Watson and Crick cleavage files. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in Watson and Crick cleavage files.}
  \item{watsonfile}{one string for the path and name of Watson cleavage file. The file may contain two or more chormosomes. The first columm is chromosome name, the second column is the chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see requirement for \code{watsonfile}.}
}
\value{\code{trainTEMP1} outputs cleavage trained template for one-template model under the current working directory. The output file contains a 1 by 8 matrix, standing for average cleavage frequency on the clustered positions including (-2, -1, 0, 1, 4, 5, 6, 7) around nucleosome center. The file is named as "trainTEMP1result.txt". 
}


\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
trainTEMP1(seqname=c("chrI","chrII"),watsonfile=wfile,crickfile=cfile)
}
}
%\keyword{}

