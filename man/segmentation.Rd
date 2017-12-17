\name{segmentation}
\alias{segmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation Matrix Segmentation
}
\description{
For a given chromosome, gene correlation matrix segmentation is performed. Regions with high
correlation are identified using an exact test. The expression matrix must not contain NA's and genes with same expression value (i.e. null gene expression).
}
\usage{
segmentation(CHR, EXP, genes, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CHR}{
chromosome name
}
  \item{EXP}{
Gene expression matrix (raw/corrected for CNV). Columns correspond to patients and rows to genes. The expression matrix must not contain either NA's or genes with same expression value (i.e. null gene expression).
}
  \item{genes}{
  Gene ID(name) vector.
}
\item{S}{
  Threshold for model selection. Default S=0.7.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{Results}{Matrix containing information about the genomic regions. Each region corresponds to a row of the matrix, the one  with the smallest p-value is on the top of the list.}
\item{Results$CHR}{Chromosome}
\item{Results$Start/End}{region boundaries  with respect to the physical location of the gene in the chromosome}
\item{Results$Rho}{\eqn{\rho} correlation}
\item{Results$length}{number of genes in the region}
\item{Results$first/last gene}{name of the first/last gene in the region}
\item{Results$p-value}{p-value as obtained from the test}
\item{Results$genes}{names of genes belonging to the region}
\item{rho0}{estimate of the background correlation}
\item{likelihood}{log-likelihood}
\item{K}{number of segments}
}
\references{
Delatola E. I., Lebarbier E., Mary-Huard T., Radvanyi F., Robin S., Wong J.(2017). SegCorr: a statistical procedure for the detection of genomic regions of correlated expression. BMC Bioinformatics, 18:333.
}
\author{
E. I. Delatola, E. Lebarbier, T. Mary-Huard, F. Radvanyi, S. Robin, J. Wong.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#data(EXP_raw)
#G = cor(t(EXP_raw))## calculating the gene x gene correlation matrix
#image(G)## plotting the correlation matrix
#results = segmentation(EXP = EXP_raw)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{multivariate}
