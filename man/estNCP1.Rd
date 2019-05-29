\name{estNCP1}
\alias{estNCP1}
\title{R function for estimating NCP score and NCP score/noise ratio using one-template model}
\description{This function invokes Fortran codes to estimate NCP score and NCP score/noise ratio using one-template model.}

\usage{estNCP1(seqname,genfile,watsonfile,crickfile,temp1)}
\arguments{
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in Watson and Crick cleavage files. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in Watson and Crick cleavage files.}
  \item{genfile}{one or multiple strings, each for the path and name of DNA sequence files in FASTA format. This sequence file can be located in any directory. It must contain only one sequence. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be no longer than 400 bp.}
   \item{watsonfile}{one string for the path and name of Watson cleavage file. The file may contain two or more chormosomes. The first columm is chromosome name, the second column is the chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see requirement for \code{watsonfile}.}
  \item{temp1}{one string for the path and name of the output file of \code{trainTEMP1}. The default choice (by setting temp1="default") is the template trained based on the yeast chemical mapping data from Brogaard et al. (2012) and Xi et al. (2012).}
}

\value{\code{estNCP1} outputs the estimation results under the current working directory. The output file is named as "NCPscore.ratio_1temp.txt", \code{1temp} stands for one-template model. It contains five columns:
\item{chr.}{chromosome name;}
\item{Position}{chromosome coordinate;}
\item{NCPscore}{estimated NCP score;}
\item{Ratio}{NCP/noise ratio;}
\item{cNCPscore}{NCP score after correction for strand asymmetry of cleavages. It is used for nucleosome occupancy calculation (see function \code{calOccup}).}
}

\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

chrI=system.file("extdata", "chrI.fa",package="nucmapData")
wfile=system.file("extdata", "watson12.txt",package="nucmapData")
cfile=system.file("extdata", "crick12.txt",package="nucmapData")
estNCP1(seqname="chrI",genfile=chrI,watsonfile=wfile,crickfile=cfile)
}
}
%\keyword{}

