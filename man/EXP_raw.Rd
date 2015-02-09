\name{EXP_raw}
\alias{EXP_raw}
\docType{data}
\title{
Simulated Gene Expression
}
\description{
Gene expression profiles have been generated for 30 patients and 500 genes. Background correlation is set to 0.08 and the correlation for H1 regions to 0.5. The location of the H1 regions is as suggested in the work of Lai et al.(2005), i.e region 1 [101, 105], region 2 [201, 210], region 3 [301, 320] and region 4 [401, 440].  }
\usage{data("EXP_raw")}
\format{
  Data frame containing the gene expression signal for 500 genes (rows) on 30 patients (columns).
}
\references{
Lai, W. R., Johnson, M. D., Kucherlapati, R., & Park, P. J. (2005). Comparative analysis of algorithms for identifying amplifications and deletions in array CGH data. Bioinformatics, 21(19), 3763-3770.
}
\examples{
data(EXP_raw)
G = cor(t(EXP_raw))## calculating the gene x gene correlation matrix
image(G)## plotting the correlation matrix
}
\keyword{datasets}
