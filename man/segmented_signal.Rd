\name{segmented_signal}
\alias{segmented_signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Segmentation
}
\description{
Mean segmentation on the genomic signal is performed using the Fpsn function of the jointseg package. }
\usage{
segmented_signal(SNP.Chr, Kmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNP.Chr}{
SNP/CGH profile matrix for a given chromosome (NA's not allowed). Columns correspond to patients and rows to probes.
}
  \item{Kmax}{
Maximum number of segments.
}
}

\value{
Smoothed genomic signal matrix. Rows correspond to probes and columns to patients.
}
\references{
Morgane Pierre-Jean, Guillem Rigaill and Pierre Neuvial. Performance evaluation of DNA copy number segmentation methods. Briefings in Bioinformatics (2015) 16 (4): 600-615.
}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CNV_correction}}, \code{\link[jointseg]{Fpsn}}
}
\examples{
#data(SNP)
#mu.SNP = segmented_signal(SNP ,100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
