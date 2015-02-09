\name{segmented_signal}
\alias{segmented_signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Segmentation 
}
\description{
Multiple profile mean segmentation on the genomic signal is performed using the multiseg function of the cghseg package. Heterogeneous genomic profiles can only be analysed if allocated to homogeneous groups.}
\usage{
segmented_signal(SNP.Chr, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNP.Chr}{
SNP/CGH profile matrix for a given chromosome not containing NA's. Columns correspond to patients and rows to probes.
}
  \item{group}{
Indicator denoting the different variance groups among the patients. By default, all genomic profiles allocated to one group.
}
}

\value{
Smoothed genomic signal matrix. Rows correspond to probes and columns to patients.
}
\references{
Picard, F., Lebarbier, E., Hoebeke, M., Rigaill, G., Thiam, B., & Robin, S. (2011). Joint segmentation, calling, and normalization of multiple CGH profiles. Biostatistics.
}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CNV_correction}}, \code{\link[cghseg]{multiseg}}, \code{\link[cghseg]{getsegprofiles}}
}
\examples{
#data(SNP)
#mu.SNP = segmented_signal(SNP ,rep(1,dim(SNP)[2]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}% __ONLY ONE__ keyword per line
