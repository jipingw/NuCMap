\name{trainTEMP4}
\alias{trainTEMP4}
\title{R function for training four-template model}
\description{This function trains cleavage template for four-template model.
}
\usage{trainTEMP4(seqname,genfile,watsonfile,crickfile,center)}

\arguments{
  \item{seqname}{the default value is "all", which specifies all chromosomes listed in Watson and Crick cleavage files. One can also specify one or more individual chromosomes, e.g. "chrI" or c("chrI", "chrII"). The sequence name format must be same as in Watson and Crick cleavage files.}
  \item{genfile}{one or multiple strings, each string is for the path and name of a DNA sequence file in FASTA format. This sequence file can be located in any directory. It must contain only one sequence. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be no longer than 400 bp.}
  \item{watsonfile}{one string for the path and name of Watson cleavage file. The file may contain two or more chormosomes. The first columm is chromosome name, the second column is the chromosome coordinate, the third column is the cleavage frequency.}
  \item{crickfile}{see requirement for \code{watsonfile}.}
  \item{center}{one string for the path and name of the file (output file from \code{callUNIQUE}) where a unique nucleosome map is saved.}
}
\value{\code{trainTEMP4} outputs trained template for four-template model under the current working directory. The file contains a 4 by 8 matrix, standing for average cleavage frequency on the clustered positions including (-2, -1, 0, 1, 4, 5, 6, 7) around nucleosome centers for four templates(+A(-3)+T(+3), (+A(-3)-T(+3), (-A(-3)+T(+3) and (-A(-3)-T(+3), where a +/- in front of the A/T nucleotide stands for presence/absence of the letter at -3/+3 postion respectively relative to the nucleosome center). The file is named as "trainTEMP4result.txt". 
}

\examples{
\dontrun{
library(NuCMap)
library(nucmapData)
## the user should replace "system.file("extdata",~,package="nucmapData")" 
## by the actual path and file name.

wfile=system.file("extdata","watson12.txt",package="nucmapData")
cfile=system.file("extdata","crick12.txt",package="nucmapData")
chrI=system.file("extdata","chrI.fa",package="nucmapData")
chrII=system.file("extdata","chrII.fa",package="nucmapData")
umap=system.file("extdata","UNIQUEcenters.txt",package="nucmapData")
trainTEMP4(seqname=c("chrI","chrII"),genfile=c(chrI,chrII),watsonfile=wfile,crickfile=
cfile,center=umap)
}
}
%\keyword{}

