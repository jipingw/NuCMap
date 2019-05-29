\name{calOccup}
\alias{calOccup}
\title{R function for calculating nucleosome occupancy genome-wide}
\description{This function calculates nucleosome occupancy genome-wide.}
\usage{calOccup(estresults,genfile,rednufile,seqname)}

\arguments{
  \item{estresults}{one string for the path and name of the output file from function \code{estNCP1} or \code{estNCP4}.}
  \item{genfile}{one or multiple strings, each for the path and name of a DNA sequence file in FASTA format. The sequence files can be located in any directory. It must contain only one sequence. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be no longer than 400 bp.}
  \item{rednufile}{one string for the path and name of the file (output file from \code{callRED}) where a redundant nucleosome map is saved.}
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in \code{estresults}. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in \code{estresults}.}
}

\value{\code{calOccup} calculates nucleosome occupancy on the genome. The output file is named as "NuOccupancy.txt" and contains two columns:

\item{Position}{chromosome coordinate;}
\item{Occupancy}{nucleosome occupancy score.}
}

\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

NCP4=system.file("extdata","NCPscore.ratio_4temp.txt",package="nucmapData")
chrI=system.file("extdata","chrI.fa",package="nucmapData")
chrII=system.file("extdata","chrII.fa",package="nucmapData")
rmap=system.file("extdata","REDcenters.txt",package="nucmapData")
calOccup(estresults=NCP4,genfile=c(chrI,chrII),rednufile=rmap,seqname=
c("chrI","chrII"))
}
}
%\keyword{}
