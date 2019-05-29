\name{estNCPcall}
\alias{estNCPcall}
\title{R function for estimating NCP score, NCP/noise ratio and making unique nucleosome center calls using four-template model.}
\description{This function produces NCP score estimation and unique nucleosome set based on the NCP/noise ratio using four-template model. 
}
\usage{estNCPcall(seqname,genfile,watsonfile,crickfile)}

\arguments{
  \item{seqname}{the default value is "all", which means all chromosomes listed in Watson and Crick cleavage files. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in Watson and Crick cleavage files.}
  \item{genfile}{one or multiple strings, each for the path and name of a DNA sequence file in FASTA format. This sequence file can be located in any directory. It must contain only one sequence. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be not longer than 400 bp.}
  \item{watsonfile}{one string for the path and name of Watson cleavage file. The file can contain information of multiple chormosomes. For the Watson cleavage file, the first columm is chromosome's name, the second columm is chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see requirement for \code{watsonfile}.}
}
\value{\code{estNCPcall} outputs estimated NCP score, NCP/noise ratio and unique nucleosome centers under current working directory. Please refer to \code{estNCP1}, \code{estNCP4} and \code{callUNIQUE} about the output file names and format.    
}


\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

chrI=system.file("extdata","chrI.fa",package="nucmapData")
wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
estNCPcall(seqname="chrI",genfile=chrI,watsonfile=wfile,crickfile=cfile)
}
}
%\keyword{}

