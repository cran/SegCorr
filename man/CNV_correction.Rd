\name{CNV_correction}
\alias{CNV_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Corrects Gene Expression for CNV
}
\description{
 Correcting gene expression signal for CNV. %The correction is done in three steps. First, the genomic and expression signals are aligned. Due to the alignment, there might be
 %genes for which we dont have a smoothed genomic signal. Therefore we need to interpolate it. Finally, we perform the correction.
}
\usage{
CNV_correction(s.Position.EXP, e.Position.EXP, Position.SNP, mu.SNP, EXP)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{s.Position.EXP}{
 vector with gene start position
}
  \item{e.Position.EXP}{vector with gene end position}
  \item{Position.SNP}{vector with SNP/CGH positions}
  \item{mu.SNP}{Smoothed genomic signal matrix not containing NA values. Rows correspond to probes, while columns to patients. The ordering of the patients must be the same as in the EXP matrix.}
  \item{EXP}{Gene expression matrix must not contain NA's and genes with same expression value (i.e. null gene). Rows correspond to probes, while columns to patients.  Again, ordering of patients must be the same between EXP and mu.SNP matrices.}
}
\details{
Overlapping genes may correspond to the same SNP/CGH probes.

%In the human genome, it is very often the case that different genes have either the same start/end position or one gene is located within one other gene. To tackle this problem, we use an adhoc solution. For a given gene, we check the number of SNP probes belonging to that gene and calculate the mean SNP value. In our case, this means that we could have the same SNP probes corresponding to different genes. Therefore, we would like the user to bear this in mind and conduct inference with the necessary caution.
}
\value{
 CNV corrected signal matrix.
}
\references{
Delatola E. I., Lebarbier E., Mary-Huard T., Radvanyi F., Robin S., Wong J.(2017). SegCorr: a statistical procedure for the detection of genomic regions of correlated expression. BMC Bioinformatics, 18:333.
}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented_signal}}
}
\examples{
#data.sets = c('SNP','EXP_raw')
## Each gene corresponds to one SNP probe ##
#Position_EXP = matrix(1:1000,nrow=500,byrow=TRUE)
#Position_SNP = seq(2,1000,by=2)
#data(list=data.sets)
#mu.SNP = segmented_signal(SNP ,100) ## smoothed SNP signal
#EXP.CNV = CNV_correction(Position_EXP[,1], Position_EXP[,2], Position_SNP,
#mu.SNP, EXP_raw)## corrected signal

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

