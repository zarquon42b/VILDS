#' Compute Pairwise Landmark Distances
#'
#' Compute pairwise landmark distances from landmark array
#' @param A k x m x n array of landmarks
#' @return
#' returns a matrix with each row containing the pairwise distances between landmarks
#' @export
EDMA <- function(A) {
    
    n=dim(A)[3]; p=dim(A)[1]; k=dim(A)[2] # n; p; k
    Name=NA
    for(i in 1:(p-1)) {
        Name[(sum(((p-1):1)[1:i])-p+i+1):(sum(((p-1):1) [1:i]))]=paste(i,(i+1):p,sep="-")
    }
    ildn=(1:(p*(p-1)/2)) # ildn # number of ILD with p landmarks
    ildl=Name[ildn] # ildl # names of ILD

    E <- t(apply(A,3,function(x) x <- as.vector(dist(x))))
    colnames(E)=ildl
    return(E)
}


#' Compute R2 for Interlandmark Distances
#'
#' Compute R2 for Interlandmark Distances explaining between group differences
#' @param x array containing landmarks
#' @param groups vector containing group assignments
#' @param startref matrix containing start config landmarks
#' @param target matrix containing target config landmarks
#' @param R2tol numeric: upper percentile for SILD R2 in relation to factor 
#' @param plot logical: if TRUE show graphical output of steps involved
#' @importFrom Morpho vecx bindArr
#' @return
#' matrix containing landmark information with the highest R2 (Andrea please specify here)
#' 
#' @examples
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' startref <- arrMean3(proc$rotated[,,groups=="ch"])
#' target <- arrMean3(proc$rotated[,,groups=="eu"])
#' ildR2 <- slidR2(proc$rotated,groups,startref,target,plot=TRUE)
#' @export 
slidR2 <- function(x,groups,startref,target,R2tol=.95,plot=FALSE) {
    D <- dim(x)[2] ## get LM dimensionality
    ## convert to matrix in x1,y1,z1,... format
    

    ild <- EDMA(x)
    allSILD <- round(ild, digits=6)
    if (length(dim(x)) == 3)
        x <- vecx(x,byrow = T)
    
    twosh <- bindArr(startref,target,along=3)
    E <- EDMA(twosh)
    twosh.SILD <- round(as.data.frame(t(E)), digits=6)
    colnames(twosh.SILD)=c("start","target")

    av.twosh.SILD <- apply(twosh.SILD,1,mean);
   if (plot) {
        par(mfrow=c(1,3))
        plot(twosh.SILD, asp=1)
        plot(twosh.SILD[,1], av.twosh.SILD, asp=1, main="averaged SILDs ~ obs. SILDs?", xlab="start", ylab="averaged start-target")
        plot(twosh.SILD[,2], av.twosh.SILD, asp=1, xlab="target", ylab="averaged start-target")
        par(mfrow=c(1,1))
    }
    if (plot) {
    if(interactive())
        readline("proceed? (press any key to proceed)\n")
    }
    ratios.twosh.SILD=twosh.SILD$target/twosh.SILD$start
    names(ratios.twosh.SILD) <- rownames(twosh.SILD)
    ratios.twosh.SILD.sorted <- sort(ratios.twosh.SILD)
    av.twosh.SILDsortedasratios <- av.twosh.SILD[names(ratios.twosh.SILD.sorted)]

    if (plot) {
        par(mfrow=c(2,2))
        cor(av.twosh.SILDsortedasratios, ratios.twosh.SILD.sorted)
        plot(av.twosh.SILDsortedasratios, ratios.twosh.SILD.sorted, main="are SILD ratios varying more in smaller SILDs?", xlab="average of start & target SILD", ylab="target/start SILD ratio")
        abline(a=1, b=0, col="grey", lwd=3, lty=1) 
    }
    all.R2=as.vector(cor(allSILD, as.numeric(groups))^2)
    names(all.R2)=colnames(allSILD); all.R2sorted=sort(all.R2, decreasing=TRUE) # R2 of SILDs compared to factor in total sample
    av.twosh.SILDsorted <- av.twosh.SILD[names(all.R2sorted)]
    # cor(all.R2sorted, av.twosh.SILDsorted)
    if (plot) {
        plot(av.twosh.SILDsorted,all.R2sorted, main="have R2s a relation to length of SILDs?", xlab="average of start & target SILD", ylab="R2 for sample SILDs vs factor"); abline(a=quantile(all.R2sorted, probs=R2tol), b=0, col="grey", lwd=3, lty=1)
        hist(ratios.twosh.SILD.sorted, breaks=sqrt(length(ratios.twosh.SILD.sorted)), prob=TRUE, main="hist. of target to start SILD ratios"); lines(density(ratios.twosh.SILD.sorted), col="red")
        hist(all.R2sorted, breaks=sqrt(length(all.R2sorted)), prob=TRUE, main="hist. of target to start SILD ratios"); lines(density(all.R2sorted), col="red"); par(mfrow=c(1,1))
    }

    largerR2=round(subset(all.R2sorted, all.R2sorted>stats::quantile(all.R2sorted, probs=R2tol)), digits=7)
    
    ratios.twosh.SILD.ofBiggestR2 <- round(ratios.twosh.SILD[names(largerR2)], digits=7) # finds the corresponding SILDs ratios
    largerR2.rankedByRatios=1+length(ratios.twosh.SILD.sorted)-rank(sort(round(abs(1-ratios.twosh.SILD.sorted), digits=7)), ties.method="random")[names(largerR2)]
    outOf100.largerR2.rankedByRatios <- round(largerR2.rankedByRatios*100/dim(allSILD)[2], digits=0)
    o1=rbind(largerR2, ratios.twosh.SILD.ofBiggestR2, largerR2.rankedByRatios, outOf100.largerR2.rankedByRatios)
    o2=round(o1, digits=2)

    return(o2)


    
}
