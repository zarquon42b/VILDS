#' Compute Pairwise Landmark Distances
#'
#' Compute pairwise landmark distances from landmark array
#' @param A k x m x n array of landmarks
#' @return
#' returns a matrix with each row containing the pairwise distances between landmarks
#' @examples
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' ### compute ILDS from Procrustes aligned landmarks.
#' edma <- ILDS(proc$rotated)
#' @export
ILDS <- function(A) {
    
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
#' @param groups vector containing group assignments. If more than two levels, a pair of needs to be specified.
#' @param R2tol numeric: upper percentile for SILD R2 in relation to factor
#' @param bg.rounds numeric: number of permutation rounds to assess between group differences
#' @param wg.rounds numeric: number of rounds to assess noise within groups by bootstrapping.
#' @param which integer (optional): in case the factor levels are > 2 this determins which factorlevels to use
#' @param reference matrix containing start config landmarks. If NULL, it will be computed as mean for group 1.
#' @param target matrix containing target config landmarks. If NULL, it will be computed as mean for group 2.
#' @param mc.cores integer: number of cores to use for permutation tests.
#' @param plot logical: if TRUE show graphical output of steps involved
#' @param silent logical: suppress console output

#' @importFrom Morpho vecx bindArr
#' @return
#' A list containing:
#' \item{largeR2}{containing landmark information with the highest R2 (Andrea please specify here)}
#' \item{allR2}{vector with ILD specific R2-values, sorted decreasingly}
#' \item{reftarILDS}{matrix with columns containing ILDs for reference and target shapes}
#' \item{sampleILD}{matrix containing ILDs of entire sample}
#' \item{R2tol}{R2-threshold used}
#' \item{reference}{reference used}
#' \item{target}{target used}
#' \item{bg.test}{result from between-group testing}
#' \item{confR2}{confidence for relevant ILDs from bootstrapping}
#' 
#' @examples
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' ilds <- ILDSR2(proc$rotated,groups,plot=TRUE,wg.rounds=999,mc.cores=2)
#' @importFrom Morpho permudist arrMean3
#' @import graphics stats 
#' @export 
ILDSR2 <- function(x,groups,R2tol=.95,bg.rounds=999,wg.rounds=999,which=NULL,reference=NULL,target=NULL,mc.cores=1,plot=FALSE,silent=FALSE) {
    D <- dim(x)[2] ## get LM dimensionality
    ild <- ILDS(x)
    xorig <- x
    allSILD <- round(ild, digits=6)
    if (length(dim(x)) == 3)
        x <- vecx(x,byrow = T)

    ## set up grouping variable
     if (!is.factor(groups))
        groups <- factor(groups)
    
    if (is.factor(groups)) {
        groups <- factor(groups)
        lev <- levels(groups)

    }

    old_groups <- groups
    if (!is.null(which)) {
        groups <- factor(groups[groups %in% lev[which]])
        lev <- levels(groups)
        x <- x[which(old_groups %in% lev),]
    }
    ng <- length(lev)

    if (ng < 2)
        stop("provide at least two groups")
    if (length(groups) != nrow(x))
        warning("group affinity and sample size not corresponding!")

   
        
    ## create reference and target
    if (is.null(reference))
        reference <- arrMean3(xorig[,,groups==lev[1]])
    
    if (is.null(target))
        target <- arrMean3(xorig[,,groups==lev[2]])
    
    twosh <- bindArr(reference,target,along=3)
    E <- ILDS(twosh)
    twosh.SILD <- round(as.data.frame(t(E)), digits=6)
    colnames(twosh.SILD)=c("start","target")

    av.twosh.SILD <- apply(twosh.SILD,1,mean);
   
    ratios.twosh.SILD <- twosh.SILD$target/twosh.SILD$start
    names(ratios.twosh.SILD) <- rownames(twosh.SILD)
    ratios.twosh.SILD.sorted <- sort(ratios.twosh.SILD)
    av.twosh.SILDsortedasratios <- av.twosh.SILD[names(ratios.twosh.SILD.sorted)]

    ## compute R2
    all.R2 <- as.vector(cor(allSILD, as.numeric(groups))^2)
    names(all.R2) <- colnames(allSILD)
    all.R2sorted <- sort(all.R2, decreasing=TRUE) # R2 of SILDs compared to factor in total sample
    av.twosh.SILDsorted <- av.twosh.SILD[names(all.R2sorted)]
    # cor(all.R2sorted, av.twosh.SILDsorted)
    if (plot) {
        par(mfrow=c(2,2))
        cor(av.twosh.SILDsortedasratios, ratios.twosh.SILD.sorted)
        plot(av.twosh.SILDsortedasratios, ratios.twosh.SILD.sorted, main="are SILD ratios varying more in smaller SILDs?", xlab="average of start & target SILD", ylab="target/start SILD ratio")
         abline(a=1, b=0, col="grey", lwd=3, lty=1) 
        plot(av.twosh.SILDsorted,all.R2sorted, main="have R2s a relation to length of SILDs?", xlab="average of start & target SILD", ylab="R2 for sample SILDs vs factor"); abline(a=quantile(all.R2sorted, probs=R2tol), b=0, col="grey", lwd=3, lty=1)
        hist(ratios.twosh.SILD.sorted, breaks=sqrt(length(ratios.twosh.SILD.sorted)), prob=TRUE, main="hist. of target to start SILD ratios"); lines(density(ratios.twosh.SILD.sorted), col="red")
        hist(all.R2sorted, breaks=sqrt(length(all.R2sorted)), prob=TRUE, main="hist. of target to start SILD R2s"); lines(density(all.R2sorted), col="red"); par(mfrow=c(1,1))
    }

    largerR2 <- round(subset(all.R2sorted, all.R2sorted>stats::quantile(all.R2sorted, probs=R2tol)), digits=7)
    ratios.twosh.SILD.ofBiggestR2 <- round(ratios.twosh.SILD[names(largerR2)], digits=7) # finds the corresponding SILDs ratios
    largerR2.rankedByRatios <- 1+length(ratios.twosh.SILD.sorted)-rank(sort(round(abs(1-ratios.twosh.SILD.sorted), digits=7)), ties.method="random")[names(largerR2)]
    outOf100.largerR2.rankedByRatios <- round(largerR2.rankedByRatios*100/ncol(allSILD), digits=0)

    o1 <- rbind(largerR2, ratios.twosh.SILD.ofBiggestR2, largerR2.rankedByRatios, outOf100.largerR2.rankedByRatios)
    o2 <- round(o1, digits=2)
    out <- list(largeR2=o2,allR2=all.R2sorted,reftarILDS=twosh.SILD,sampleILD=allSILD,R2tol=R2tol,reference=reference,target=target)
    
    R2names <- colnames(o2)


     ## compute between group permutation testing
    if (bg.rounds > 0) {
        bg.test <- permudist(x,groups,rounds=bg.rounds)
        if (!silent) {
            message(paste0("P-value between groups: ",bg.test$p.value,"\n"))
        }
        out$bg.test <- bg.test
    }

    ## bootstrapping
    if (wg.rounds > 0) {
        wg.boot <- parallel::mclapply(1:wg.rounds,function(x) x <- bootstrapILDSR2(xorig,groups,rounds=wg.rounds,R2tol=R2tol),mc.cores = mc.cores)

        freqsR2 <- unlist(lapply(wg.boot,match,R2names))
        confR2 <- sapply(1:length(R2names),function(x) x <- length(which(freqsR2==x)))
        confR2 <- round((((confR2+1)/(wg.rounds+1))*100),digits=3)
        names(confR2) <- R2names
        out$wg.boot=wg.boot
        if (!silent) {
            colorILDS(confR2)
            
            }
        out$confR2 <- confR2
        
    }
    

     
    class(out) <- "ILDSR2"
    return(out)
   
}


bootstrapILDSR2 <- function(x,groups,rounds,R2tol) {
    lev <- levels(groups)
    for (i in lev) {
        tmpgroup <- which(groups==i)
        x[,,groups==i] <- x[,,sample(tmpgroup,size=length(tmpgroup),replace = TRUE)]
    }
    out <- colnames(ILDSR2(x,groups,bg.rounds=0,wg.rounds=0,plot=FALSE,R2tol,silent=TRUE)$largeR2)
    
    
}

colorILDS <- function(x) {
    cat(paste0("Bootstrapped confidence of relevant ILDS\n"))
            
            for (i in 1:length(x)) {
                cat(" ")
                itmp <- x[i]
                if (itmp > 75)
                    colfun <- crayon::green
                else if (itmp > 50)
                    colfun <-  crayon::yellow
                else
                    colfun <-  crayon::red
                cat(colfun(paste0(crayon::bold("ILDS",names(itmp)),": ",itmp,"% ")))
                cat("\n")
            }
}



#' Plot the ILDS with the relevant ILDS ighlighted
#'
#' Plot the ILDS with the relevant ILDS ighlighted
#' @param x output of function \code{\link{ILDSR2}}
#' @param ref logical: if TRUE, the reference shape defined in  \code{\link{ILDSR2}} will be plotted. Otherwise the target is used.
#' @param ... additional parametr - currently not used.
#' @examples
#' ## 3D Example
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' ilds <- ILDSR2(proc$rotated,groups,plot=FALSE,bg.rounds=0,wg.rounds=0)
#' plot(ilds)
#'
#' ## 2D Example
#' require(shapes)
#' gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
#' sex <- factor(c(rep("f",30),rep("m",29)))
#' procg <- procSym(gor.dat)
#' ildsg <- ILDSR2(procg$rotated,sex,plot=FALSE,bg.rounds=0,wg.rounds=0)
#' plot(ildsg)
#' @importFrom Morpho deformGrid2d deformGrid3d
#' @export
plot.ILDSR2 <- function(x,ref=TRUE,...) {
     if (!inherits(x, "ILDSR2")) 
        stop("please provide object of class 'ILDSR2'")
    reftarILDS <- x$reftarILDS
    rn <- rownames(reftarILDS)
    pairing <- (matrix(as.integer(unlist(strsplit(rn,split = "-"))),length(rn),2,byrow=T))
    if (ref)
        reference <- x$reference
    else
        reference <- x$target
    ref0 <- reference[pairing[,1],]
    ref1 <- reference[pairing[,2],]
    
    if (ncol(reference)==3) {
        mydeform <- deformGrid3d
    } else
        mydeform <- deformGrid2d
     highlight <- colnames(x$largeR2)
     if (!is.null(highlight)) {
        hm <- match(highlight,rn)
        mydeform(reference,reference,lines=F,lwd=0,show=1,cex2=0,...)
        mydeform(ref0[hm,],ref1[hm,],add=T,lcol = "red",lwd=3,show=1,cex2=0,...)
        mydeform(ref0[-hm,],ref1[-hm,],add=T,lcol = "black",lwd=1,show=1,cex2=0,...)
    } else {
        mydeform(ref0,ref1)
    }
} 

