\name{plotAATT}
\alias{plotAATT}
\title{R function for plotting AA/TT/TA/AT frequency against to the distance from the nucleosome center}
\description{This function plots AA/TT/TA/AT frequency against to the distance from the nucleosome center.
}
\usage{plotAATT(seqname,genfile,center)}

\arguments{
  \item{genfile}{one or multiple strings, each string is for the path and name of a DNA sequence file in FASTA format. This sequence file can be located in any directory. It must contain only one sequence. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be no longer than 400 bp.}
  \item{center}{one string for the path and name of the file where a unique or redundant nucleosome map is saved.}
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in \code{center}. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in \code{center}.}
}
\value{\code{plotAATT} plots AA/TT/TA/AT frequency against to the distance from the nucleosome center.
}


\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

chrI=system.file("extdata", "chrI.fa",package="nucmapData")
chrII=system.file("extdata", "chrII.fa",package="nucmapData")
umap=system.file("extdata", "UNIQUEcenters.txt",package="nucmapData")
plotAATT(genfile=c(chrI,chrII),center=umap,seqname=c("chrI","chrII"))
}
}

%\keyword{}

