#################################################
###### Running segmentation  for SNP signal ########
#################################################
### Summary: Using the jointseg package, we perform univariate segmentation on each SNP signal

### Input: SNP.Chr:   the SNP signal for a given chromosome CHR -- NA's removed
###        Kmax:      maximum number of segments


### Output: mu.SNP :  mean segmented signal
options(warn=1)
segmented_signal = function(SNP.Chr, Kmax){

  if(missing(SNP.Chr)){
    stop("segmented_signal: SNP matrix missing")
  }else{


    p = dim(SNP.Chr)[1] # no of probes
    no = dim(SNP.Chr)[2]# no of patients


    mu.SNP = matrix(NA,p,no)



    for(j in 1:no){
      cat(paste('Performing SNP mean segmentation for patient =',j,sep=" "),sep="\n")

      res = Fpsn(SNP.Chr[,j], Kmax=Kmax)

      est.sd = mad(diff(SNP.Chr[,j])/sqrt(2))
      Khat = which.min(res$J.est + 2*(0:(Kmax-1))*log(p))
      changepoint =c(0, res$t.est[Khat, 1:Khat])

      for(k in 1:(length(changepoint)-1)){
        mu.SNP[(changepoint[k]+1):changepoint[k+1],j] = mean(SNP.Chr[(changepoint[k]+1):changepoint[k+1],j])
      }


    }

    return(mu.SNP)
    #cat('Done!')
  }
}
