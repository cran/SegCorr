\name{SegCorr}
\alias{SegCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs CNV Correction and Correlation Matrix Segmentation
}
\description{
Gene expression is corrected for CNV events  must not contain NA's and genes with same expression value (i.e. null gene expression). Segmentation is used to detect changes in the correlation pattern. Regions with high
correlation are identified using an exact test.
}
\usage{
SegCorr(CHR, EXP, genes,S, CNV, SNPSMOOTH, Position.EXP, 
SNP.CHR, SNP, Position.SNP, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CHR}{
  Chromosome allocation vector for the genes.
}
  \item{EXP}{
Gene expression matrix (raw/corrected for CNV). Columns correspond to patients and rows to genes. The expression matrix must not contain either NA's or genes with same expression value (i.e. null gene expression)}
  \item{genes}{
Gene ID(name) vector.
}
\item{S}{
  Threshold for model selection. Default S=0.7.
}
  
  \item{CNV}{
 Logical variable indicating whether to perform CNV correction. When CNV=T, the correction is performed. Default value CNV=F.
}

\item{SNPSMOOTH}{
 (Optional Argument when CNV=T) Logical variable indicating whether to perform SNPSMOOTH. When SNPSMOOTH=T, the smoothing is performed. Default value SNPSMOOTH=F.
}

\item{Position.EXP}{(Optional Argument when CNV=T) Expression position matrix.
First column is the start position and the second is the end position.
}

  \item{SNP.CHR}{(Optional Argument when CNV=T)
 Chromosome allocation vector for genomic probes.
}
  \item{SNP}{(Optional Argument when CNV=T)
SNP profile matrix not containing NA's. Columns correspond to patients and rows to probes.
}
  \item{Position.SNP}{(Optional Argument when CNV=T)
vector with SNP positions
}
  \item{group}{(Optional Argument when CNV=T and SNPSMOOTH=T)
Indicator denoting the different variance groups among the different patients. By default, all genomic profiles allocated to one group.
}

}

\value{
\item{Results}{Matrix containing information about the genomic regions. Each region corresponds to a row of the matrix, the one  with the smallest p-value is on the top of the list.} 
\item{Results$CHR}{Chromosome}
\item{Results$Start/End}{the region boundaries  with repsect to the physical location of the gene in the chromosome}
\item{Results$Rho}{\eqn{\rho} correlation}
\item{Results$length}{number of genes in the region}
\item{Results$first/last gene}{name of the first/last gene in the region}
\item{Results$p-value}{p-value as obtained from the test}
\item{Results$genes}{names of the genes belonging to the region}
\item{Results$p-valueadj}{p-value of the region corrected for multiple testing}
\item{Chromosome.Inf}{Matrix containing the estimated background correlation (rho0.hat) per chromsome, the number of segments and the log-loglikehood.}
\item{EXP.corrected}{If the CNV option is chosen, the corrected signal is given.}
}
\details{
Overlapping genes may correspond to the same genomic probes. 
%In the human genome, it is very often the case that different genes have either the same start/end position or one gene is located within one other gene. To tackle this problem, we use an adhoc solution. For a given gene, we check the number of SNP probes belonging to that gene and calculate the mean SNP value. In our case, this means that we could have the same SNP probes corresponding to different genes. Therefore, we would like the user to bear this in mind and conduuct inference with the necessary caution.
}
\references{
Delatola E. I., Lebarbier E., Mary-Huard T., Radvanyi F., Robin S., Wong J.(2015). SegCorr: a statistical procedure for the detection of genomic regions of correlated expression. Preprint on Arxiv.
}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CNV_correction}},\code{\link{segmentation}}
}
\examples{
data('EXP_raw')
CHR = rep(1,dim(EXP_raw)[1])

results = SegCorr(CHR = CHR, EXP = EXP_raw, CNV = FALSE,S=0.7)

################drawing the heatmap for one region ###########################
tau = results$Region.List[1,2]: results$Region.List[1,3] 
heatmap(as.matrix(EXP_raw[tau,]))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{multivariate}