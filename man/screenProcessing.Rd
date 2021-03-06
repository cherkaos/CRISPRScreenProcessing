% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screenProcessing.R
\name{screenProcessing}
\alias{screenProcessing}
\title{screenProcessing function}
\usage{
screenProcessing(
  inputfile,
  controlStart,
  controlEnd,
  resultfile,
  maxsgRNA,
  minReadCount,
  zscore,
  orderOutput = TRUE,
  shortOutput = TRUE
)
}
\arguments{
\item{inputfile}{Input file name. File must contains clones/organoids in the columns and sgRNAs in the rows. First row are the clones/organoids names and the two first column are the sgRNAs name and IDs. Control clones/organoids should be placed at the beginning/left of the file and unknown clones/organoids subsequently. File can be a csv or txt a file (comma separated) located in the working directory.}

\item{controlStart}{Name of first control clone/organoid in inputfile. The name should match the column name and should not contain special characters.}

\item{controlEnd}{Name of last control clone/organoid in inputfile. The name should match the column name and should not contain special characters.}

\item{resultfile}{Output file name.}

\item{maxsgRNA}{Maximum number of integrated sgRNAs per clone/organoid.}

\item{minReadCount}{Minimum read count to consider a sgRNA.}

\item{zscore}{Logical value indicating if the read counts should be zscore for each clones/organoids.}

\item{orderOutput}{Logical value indicating if clones should be ordered by the number of integrated sgRNAs in the output file.}

\item{shortOutput}{Logical value indicating if all sgRNAs (e.g. the non integrated one) should be included in the output file.}
}
\value{
resultfile
}
\description{
screenProcessing function
}
\examples{
print("Examples")
}
