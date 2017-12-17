\name{SegCorr-package}
\alias{SegCorr-package}
\docType{package}
\title{
Detecting Correlated Genomic Regions
}
\description{
Performs correlation matrix segmentation and applies a test procedure to detect highly correlated regions in gene expression. The segmentation procedure detects changes in the patterns of the gene expression correlation matrix. The test procedure asseses which regions exhibit a significantly high
level of correlation. Additionally, a preprocessing procedure is provided to correct gene expression for copy number variation.
}
\details{
\tabular{ll}{
Package: \tab SegCorr\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-01-19\cr
License: \tab GPL-2\cr
}

}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong.

Maintainer:  Eleni Ioanna Delatola <eldelatola@yahoo.gr>
}
\references{
Delatola E. I., Lebarbier E., Mary-Huard T., Radvanyi F., Robin S., Wong J.(2017). SegCorr: a statistical procedure for the detection of genomic regions of correlated expression. BMC Bioinformatics, 18:333.}
\keyword{ package }
\seealso{
\code{\link[jointseg]{Fpsn}}
}
\examples{
#data.sets = c('SNP','EXP_raw')
## Each gene corresponds to one SNP probe ##
#Position_EXP = matrix(1:1000,nrow=500,byrow=TRUE)
#Position_SNP = seq(2,1000,by=2)
#data(list=data.sets)
#CHR = rep(1,dim(EXP_raw)[1])
#SNP.CHR = rep(1,dim(SNP)[1])

#results = SegCorr(CHR = CHR, EXP = EXP_raw, CNV = TRUE, SNPSMOOTH=TRUE,
#Position.EXP = Position_EXP, SNP.CHR = SNP.CHR, SNP=SNP , Position.SNP = Position_SNP)

################drawing the heatmap for one region ###########################
#tau = results$Region.List[1,2]: results$Region.List[1,3]
#EXP.CNV =  results$EXP.corrected
#heatmap(EXP.CNV[tau,])

}

