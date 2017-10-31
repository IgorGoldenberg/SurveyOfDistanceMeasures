#Estimate Hellinger distance between two matrices assuming each has a normal distribution
HellingerDistanceFrames<-function(df1, df2)
{
    x = ncol(df1)
    y = ncol(df2)
    if (x != y)
        stop("Frames must have same number of columns")
    if (x==1)
    {
        UnivariateNormalHellingerDistance(mean(df1[1]), mean(df2[1]), sd(df1[1]), sd(df2[2]))
    }
    else
    {
        stat1 = EstimateMultiVariateNormal(df1)
        stat2 = EstimateMultiVariateNormal(df2)
        MultivariateNormalHellingerdistance(stat1[[1]], stat2[[1]], stat1[[2]], stat2[[2]])
    }
}

# Calculate theoretical HD for known univariate normal distributions
UnivariateNormalHD_Theoretical<-function(mu1, mu2, sigma1, sigma2)
{
    totalVar = sigma1^2+sigma2^2
    sqrt(1-exp(-(mu1-mu2)^2/4/(totalVar))*sqrt(2*sigma1*sigma2/totalVar))
}

# Calculate theoretical HD for known multivariate normal distributions
MultivariateNormalHD_Theoretical<-function(mu1, mu2, sigma1, sigma2)
{
    averageMatrix = 0.5*(sigma1+sigma2)
    RMatrAve = chol(averageMatrix)
    averageInv = chol2inv(RMatrAve)
    averageDetRoot = prod(diag(RMatrAve))
    RMatr1 = chol(sigma1)
    RMatr2 = chol(sigma2)
    meanDiff = mu1-mu2
    sqrt(1-sqrt(prod(diag(RMatr1)*diag(RMatr2)))/averageDetRoot*exp(-1 * (meanDiff %*% averageInv %*% meanDiff) / 8))
}

# Calculate theoretical Kullback Leibler for known univariate normal distributions
UnivariateNormalKL_Theoretical<-function(mu1, mu2, sd1, sd2)
{
    log(abs(sd2/sd1))+0.5*(sd1^2-sd2^2+(mu2-mu1)^2)/sd2^2
}

# Calculate theoretical Kullback Leibler for known multivariate normal distributions
MultiVariateNormalKL_Theoretical<-function(mu1, mu2, sigma1, sigma2)
{
    #check that dimensions correspond and covariance matrix are symmetric
    if(!isSymmetric.matrix(sigma1) || !isSymmetric.matrix(sigma2))
    {
        stop("Variance covariance matrices must be symmetric")
    }
    n=length(mu1)
    if(n!=length(mu2) || n!= ncol(sigma1) || n!=ncol(sigma2))
    {
        stop("Mean vectors and covariance matrices must have the same dimensions")
    }
    #check validity of matrices
    EigenDecompositionSigma1<-eigen(sigma1)
    if(min(EigenDecompositionSigma1$values)<0)
    {
        stop("First Covariance matrix is not positive definite")
    }
    EigenDecompositionSigma2<-eigen(sigma2)
    if(min(EigenDecompositionSigma2$values)<0)
    {
        stop("Second Covariance matrix is not positive definite")
    }
    sigma2.Inv<-EigenDecompositionSigma2$vec %*% tcrossprod(diag(1/EigenDecompositionSigma2$values,n), EigenDecompositionSigma2$vec)
    0.5*(log(prod(EigenDecompositionSigma2$values)/prod(EigenDecompositionSigma1$values))
         -n
         +sum(EigenDecompositionSigma1$values/EigenDecompositionSigma2$values)+
             crossprod((mu1-mu2), sigma2.Inv%*%(mu1-mu2)))
}

# given two alligned densities calculate hellinger integral
CalculateHellingerIntegral<-function(X, Y1, Y2)
{
    n1 = length(X)
    n2 = length(Y1)
    if (n1 != n2)
        stop("Length of data point differs from length of PDF values")
    n3 = length(Y2)
    if (n2 != n3)
        stop("Densities have different number of points")
    combinedY<-sqrt(Y1*Y2)
    require(sfsmisc, quietly = TRUE)
    HellingerIntegral<-integrate.xy(X, combinedY, use.spline = FALSE)
    
    sqrt((1-min(HellingerIntegral,1)))
}

#Estimate hellinger distance between two univariate samples
UnivariateHD_Empirical<-function(X, Y, kernel="normal")
{
    # 1. Calculate densities
    require(KernSmooth)
    X.density<-bkde(X, kernel = kernel, gridsize = min(ceiling(length(X)/5), 400)+1)
    X.density$y[X.density$y<0]<-0
    Y.density<-bkde(Y, kernel = kernel, gridsize = min(ceiling(length(Y)/5), 400)+1)
    Y.density$y[Y.density$y<0]<-0
    # 2. Align densities
    AlignedDensities<-AlignDensityFunctions(X.density$x, X.density$y, Y.density$x, Y.density$y)
    # 3. Calculate hellinger integral and distance
    CalculateHellingerIntegral(AlignedDensities$x1, AlignedDensities$y1, AlignedDensities$y2)
}

UnivariateKL_Empirical<-function(X, Y, type = "hist", kernel="normal", eps=NULL, minCount=NULL)
{
    # 1. Calculate densities
    maxLength<-max(length(X), length(Y))
    
    if (type=="KDE")
    {
        if (is.null(eps))
        {
            eps<-(1E-4)
        }
        require(KernSmooth)
        X.density<-bkde(X, kernel = kernel, gridsize = min(ceiling(sqrt(length(X))), 400)+1)
        X.density$y[X.density$y<0]<-0
        Y.density<-bkde(Y, kernel = kernel, gridsize = length(X.density$x), range.x = c(min(X.density$x), max(X.density$y)))
        Y.density$y[Y.density$y<0]<-0
        Points<-X.density$x
        Freq1<-X.density$y
        Freq2<-Y.density$y
    }
    else if (type == "hist")
    {
        if (is.null(minCount))
        {
            Points<-hist(c(X,Y), plot = FALSE, 100)$breaks
        }
        else
        {
            Points<-hist(c(X,Y), plot = FALSE, breaks = minCount)$breaks
        }
        if (is.null(eps))
        {
            eps<-1
        }
        Freq1<-hist(X, breaks = Points, plot = FALSE)$counts
        Freq2<-hist(Y, breaks = Points, plot = FALSE)$counts
        Points<-Points[-length(Points)]+diff(Points)/2
    }
    else
    {
        stop ("Unknown type")
    }
        
    # 2. Calculate Kullback Leibler integral and distance
    CalculateKLIntegral(Points, Freq1, Freq2, eps)
   
}

# given two densities estimated at different points, 
##approximate them at the other density points to get an alighed values
# Later used to calculate hellinger integral

AlignDensityFunctions<-function(x1, y1, x2, y2)
{
    One.Approximation<-approx(x1, y1, xout = x2, yleft = 0, yright = 0, rule = 2)
    Two.Approximation<-approx(x2, y2, xout = x1, yleft = 0, yright = 0, rule = 2)
    One.Order<-order(c(x1, One.Approximation$x))
    X1<-c(x1, One.Approximation$x)[One.Order]
    Y1<-c(y1, One.Approximation$y)[One.Order]
    Y2<-c(Two.Approximation$y, y2)[One.Order]
    list(x1=X1,y1=Y1,y2=Y2)
}

CalculateKLIntegral<-function(X,Y1,Y2, eps)
{
## this function uses approach from flexmix::KLDiv
    Y1[Y1<eps]<-eps
    Y2[Y2<eps]<-eps
    Y1<-Y1/sum(Y1)
    Y2<-Y2/sum(Y2)
    sum(Y1*log(Y1/Y2))
}

# Generate pairs of vectors that should have certain HD 
# to test how close the algorithm to a "true" value
KLTheoreticalVersusEmpirical<-function(minSize=1000, maxSize=20000, step=500, SimNum=100, ZScore=seq(0,10,1),
                                              RoScore=seq(1,5,1), eps=c(1,0.1,0.01))
    
{
    set.seed(1)
    frameSize<-length(ZScore)*length(RoScore)*((maxSize-minSize)/step+1)*SimNum
    df<-expand.grid(SimSize=seq(minSize,maxSize,step),Run=1:SimNum,ZScore=ZScore, RoScore=RoScore, eps=eps)
    Distances<-mapply(GenerateAndCalculateKL, df$SimSize, df$ZScore, df$RoScore, df$eps)
    df<-cbind(df, t(Distances))
    colnames(df)[(ncol(df)-2):ncol(df)]<-c("Theoretical", "KL1", "KL2")
    df
}

GenerateAndCalculateKL<-function(SimSize, ZScore, RoScore, eps=NULL)
{
    mu1<-0
    sd1<-1
    sd2<-RoScore
    mu2<-ZScore*sd2+mu1
    Theor<-UnivariateNormalKL_Theoretical(mu1,mu2,sd1,sd2)
    X<-rnorm(SimSize, mean = mu1, sd=sd1)
    Y<-rnorm(SimSize, mean = mu2, sd=sd2)
    KL1<-UnivariateKL_Empirical(X,Y,eps = eps)
    KL2<-UnivariateKL_Empirical(Y,X,eps = eps)
    c(Theor,KL1,KL2)
}

## Univariate case
## Calculate a parameter needed to get Hellinger Distance for certian distribution
## For Normal z score returned
## For Poisson lambda is returned (assuming first lambda is 1)
## For beta, the value returned (x) means that the first distribution should use
## parameters 16/x and x and the second distribution x and 16/x
GetParametersForDistance<-function(d, distType)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    switch(as.character(distType), 
           Normal=GetZForHellingerdistance(d),
           Poisson=GetLambdaForHellingerdistance(d),
           Beta=GetBetaParamsForHellingerdistance(d)
    )
}

# calculates z score between means assuming standard deviation is the same.
# generater vectors will have means M and M+z*SD, where M and SD are arbitarbly choosen
GetZForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(10)
    }
    sqrt(-8*log(1-d*d))
}

# Calculates Lambda for Y vector, assuming lambda for X is 1
GetLambdaForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(25)
    }
    (sqrt(-2*log(1-d*d))+1)^2
}

# returns alpha. Parameters for X are 16/a and a, for Y parameters are a and 16/a
GetBetaParamsForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(1)
    }
    if (d==0)
    {
        return(4)
    }
    f<-function(a) TheoreticalHellingerBeta(16/a,a,a,16/a)-d
    alpha<-uniroot(f, c(1,4))
    alpha$root
}

# Generate pairs of vectors that should have certain HD 
# to test how close the algorithm to a "true" value
HellingerTheoreticalVersusEmpirical<-function(minSize=1000, maxSize=20000, step=500, SimNum=100, HDistances=seq(0,1,0.01),
                                              distributionTypes = c("Normal", "Poisson", "Beta"))
    
{
    set.seed(1)
    require(KernSmooth)
    require(dplyr)
    frameSize<-length(HDistances)*length(distributionTypes)*((maxSize-minSize)/step+1)*SimNum
    df<-expand.grid(SimSize=seq(minSize,maxSize,step),Run=1:SimNum,DistType=distributionTypes, TheoreticalDistance=HDistances)
    params<-expand.grid(DistType=distributionTypes, TheoreticalDistance=HDistances)
    params$param<-mapply(GetParametersForDistance, params$TheoreticalDistance, params$DistType)
    df<-left_join(df, params, by=c("TheoreticalDistance", "DistType"))
    df$EmpiricalDistance<-mapply(GenerateAndCalculate, df$SimSize, df$DistType, df$param, df$DistType=="Poisson")
    df
}

## Generate vector pair and estimate their HD
GenerateAndCalculate<-function(SimSize, DistType, param, useHist=FALSE)
{
    generatedVectors<-generateRandomVectors(SimSize, DistType, param)
    UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y)
    if (useHist)
    {
        UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y, type = "hist")
    }
    else
    {
        UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y)
    }
}

# Generate two column of data given size, distribution type 
# and a parameter that produces desired HD
generateRandomVectors<-function(size, distributionType, param)
{
    switch(as.character(distributionType), 
           Normal = data.frame(X=rnorm(size, mean = 5, sd=2),Y=rnorm(size, mean = 5+2*param, sd=2)),
           Poisson = data.frame(X=rpois(size, lambda = 1), Y=rpois(size, lambda = param)),
           Beta = data.frame(X=rbeta(size, shape1 = 16/param, shape2 = param), Y=rbeta(size, shape1 = param, shape2 = 16/param)))
}
