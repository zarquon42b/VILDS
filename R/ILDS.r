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
#' @param groups vector containing group assignments or a numeric covariate. For groups with more than two levels, a pair of needs to be specified using \code{which}
#' @param R2tol numeric: upper percentile for SILD R2 in relation to factor
#' @param bg.rounds numeric: number of permutation rounds to assess between group differences
#' @param wg.rounds numeric: number of rounds to assess noise within groups by bootstrapping.
#' @param which integer (optional): in case the factor levels are > 2 this determins which factorlevels to use
#' @param reference matrix containing start config landmarks. If NULL, it will be computed as mean for group 1.
#' @param target matrix containing target config landmarks. If NULL, it will be computed as mean for group 2.
#' @param mc.cores integer: number of cores to use for permutation tests.
#' @param plot logical: if TRUE show graphical output of steps involved
#' @param silent logical: suppress console output
#' @param ... additional parameters for internal use only.

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
#' ilds <- ILDSR2(proc$rotated,groups,plot=TRUE,wg.rounds=99,mc.cores=2)
#' if (interactive())
#' visualize(ilds)
#' ## use covariate
#' \dontrun{
#' ildsLM <- ILDSR2(proc$rotated,groups=proc$size,plot=TRUE,wg.rounds=99,mc.cores=2)
#' if (interactive())
#' visualize(ildsLM)
#' }
#'
#' ## 2D Case with size as predictor
#' require(shapes)
#' require(Morpho)
#' gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
#' procg <- procSym(gor.dat)
#' ildsg <- ILDSR2(procg$rotated,procg$size,plot=FALSE,bg.rounds=999,wg.rounds=99,mc.cores=1)
#' @importFrom Morpho permudist arrMean3
#' @import graphics stats 
#' @export 
ILDSR2 <- function(x,groups,R2tol=.95,bg.rounds=999,wg.rounds=999,which=1:2,reference=NULL,target=NULL,mc.cores=1,plot=FALSE,silent=FALSE,...) {
    D <- dim(x)[2] ## get LM dimensionality
    ild <- ILDS(x)
    xorig <- x
    allSILD <- round(ild, digits=6)
    mindim <- ncol(allSILD)
    if (length(dim(x)) == 3)
        x <- vecx(x,byrow = T)
    bootstrap <- FALSE
    args <- list(...)
    if ("bootstrap" %in% names(args)) {
        bootstrap <- args$bootstrap
    }

    ## set up grouping variable
    regression <- FALSE
    if (!is.factor(groups) && is.numeric(groups)) {
        regression <- TRUE
        if (!silent)
            message("groups is interpreted as covariate. Using regression scheme.")
    }
    
    if (is.factor(groups)) {
        groups <- factor(groups)
        lev <- levels(groups)
    }
    
    ## factor case
    if (!regression) {
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
        
    } else { ## regression case
        lmod <- lm(x~groups)
        estquant <- quantile(groups,probs=c(.1,.9))
        tarref <- predict(lmod,newdata=list(groups=estquant))
        twosh <- vecx(tarref,revert = T,lmdim = D,byrow = T)
        reference <- twosh[,,1]
        target <- twosh[,,2]
    }
    
    E <- ILDS(twosh)
    twosh.SILD <- round(as.data.frame(t(E)), digits=6)
    colnames(twosh.SILD)=c("start","target")

    av.twosh.SILD <- apply(twosh.SILD,1,mean)

   
    ratios.twosh.SILD <- twosh.SILD$target/twosh.SILD$start
    names(ratios.twosh.SILD) <- rownames(twosh.SILD)
    ratios.twosh.SILD.sorted <- sort(ratios.twosh.SILD)
    av.twosh.SILDsortedasratios <- av.twosh.SILD[names(ratios.twosh.SILD.sorted)]
    
    ## compute R2
    if (!regression)
        all.R2 <- as.vector(cor(allSILD, as.numeric(groups))^2)
    else {
        R2lm <- (lm(allSILD~groups))
        tmplm <- summary(R2lm)
        all.R2 <- sapply(tmplm,function(x) x <- x$r.squared)
    }
    names(all.R2) <- colnames(allSILD)
    
    all.R2sorted <- sort(all.R2, decreasing=TRUE) # R2 of SILDs compared to factor in total sample
    av.twosh.SILDsorted <- av.twosh.SILD[names(all.R2sorted)]

    ## combine all sample wide R2 stats in a named list
    SILDstats <- list(av.twosh.SILDsorted=av.twosh.SILDsorted,ratios.twosh.SILD.sorted=ratios.twosh.SILD.sorted,av.twosh.SILDsortedasratios=av.twosh.SILDsortedasratios)
                                       
    largerR2 <- round(subset(all.R2sorted, all.R2sorted>stats::quantile(all.R2sorted, probs=R2tol)), digits=7)
    ratios.twosh.SILD.ofBiggestR2 <- round(ratios.twosh.SILD[names(largerR2)], digits=7) # finds the corresponding SILDs ratios
    largerR2.rankedByRatios <- 1+length(ratios.twosh.SILD.sorted)-rank(sort(round(abs(1-ratios.twosh.SILD.sorted), digits=7)), ties.method="random")[names(largerR2)]
    outOf100.largerR2.rankedByRatios <- round(largerR2.rankedByRatios*100/ncol(allSILD), digits=0)

    ## create output table
    o1 <- round(rbind(largerR2, ratios.twosh.SILD.ofBiggestR2, largerR2.rankedByRatios, outOf100.largerR2.rankedByRatios),digits=2)
    out <- list(largeR2=o1,allR2=all.R2sorted,reftarILDS=twosh.SILD,sampleILD=allSILD,R2tol=R2tol,reference=reference,target=target,bg.rounds=bg.rounds,wg.rounds=wg.rounds,SILDstats=SILDstats)
    ## extract names of relevant ILDS
    R2names <- colnames(o1)

    ## compute between group permutation testing
    if (bg.rounds > 0 && !regression) {
        bg.test <- permudist(x,groups,rounds=bg.rounds)
        if (!silent)
            colorPVal(bg.test$p.value,rounds=bg.rounds)
        
        out$bg.test <- bg.test
    }
    if (regression && !bootstrap) {
        pval <- anova(R2lm)$"Pr(>F)"[2]
        if (!silent)
            colorPVal(format(pval,scientific=T,digits=3),permu=FALSE)
        out$bg.test <- pval
    }

    ## bootstrapping
    if (wg.rounds > 0) {
        wg.boot <- parallel::mclapply(1:wg.rounds,function(x) x <- bootstrapILDSR2(xorig,groups,rounds=wg.rounds,R2tol=R2tol,regression=regression),mc.cores = mc.cores)

        freqsR2 <- unlist(lapply(wg.boot,match,R2names))
        confR2 <- sapply(1:length(R2names),function(x) x <- length(which(freqsR2==x)))
        confR2 <- round((((confR2+1)/(wg.rounds+1))*100),digits=3)
        names(confR2) <- R2names
        out$confR2 <- confR2
        
        if (!silent)
            colorILDS(confR2,wg.rounds)
    }     

    class(out) <- "ILDSR2"

    if (plot) {
        plot(out)
    }
    
    return(out)
}

#' @export
print.ILDSR2 <- function(x,...) {
    print(x$largeR2)
    cat("\n")
    if (x$bg.rounds > 0) {
        colorPVal(x$bg.test$p.value,rounds=x$bg.rounds)
    } 

    if (x$wg.rounds > 0)
        colorILDS(x$confR2,x$wg.rounds)
}


bootstrapILDSR2 <- function(x,groups,rounds,R2tol,regression=FALSE) {
    if (!regression) {
        lev <- levels(groups)
        for (i in lev) {
            tmpgroup <- which(groups==i)
            x[,,groups==i] <- x[,,sample(tmpgroup,size=length(tmpgroup),replace = TRUE)]
        }
    } else {
        mysample <- sample(length(groups),replace=T,size=dim(x)[3]*2)
        x <- x[,,mysample]
        groups <- groups[mysample]
    }
    out <- colnames(ILDSR2(x,groups,bg.rounds=0,wg.rounds=0,plot=FALSE,R2tol,silent=TRUE,bootstrap=TRUE)$largeR2)
    
    
}
colorPVal <- function(x,rounds=NULL,permu=TRUE) {
    if (x < 0.05)
        pvalcol <- crayon::green
    else
        pvalcol <- crayon::red
    if (permu)
        message(crayon::bold(paste0("P-value between groups (",rounds," rounds): ",pvalcol(x),"\n")))
    else
        message(crayon::bold(paste0("P-value of linear model shape ~ predictor: ",pvalcol(x),"\n")))

    
}

colorILDS <- function(x,rounds=NULL) {
    cat(crayon::bold(paste0("Bootstrapped (",rounds," rounds) confidence of relevant ILDS:\n")))
    
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
#' @param relcol color of relevant ILDs
#' @param rescol color of "irrelevant" ILDs
#' @param lwd numeric: define line width. Relevant ILDs are displayed by \code{3*lwd}.
#' @param cex numeric: size of plot content
#' @param col define color of landmarks
#' @param pch define symbols used to plot landmarks in 2D plot.
#' @param confcol vector of colors associated with confidence. Must be of \code{length(contol)+1}.
#' @param conftol vector: set thresholds for confidence coloring
#' @param useconf logical: if TRUE, highlighting according to supported confidence of ILD is applied.
#' @param ... additional parameters passed to  \code{\link{deformGrid2d}} /  \code{\link{deformGrid3d}}.
#' @examples
#' ## 3D Example
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' ilds <- ILDSR2(proc$rotated,groups,plot=FALSE,bg.rounds=0,wg.rounds=0)
#' visualize(ilds)
#'
#' ## 2D Example
#' require(shapes)
#' gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
#' sex <- factor(c(rep("f",30),rep("m",29)))
#' procg <- procSym(gor.dat)
#' ildsg <- ILDSR2(procg$rotated,sex,plot=FALSE,bg.rounds=0,wg.rounds=0)
#' visualize(ildsg,cex=2,pch=19)
#'
#' ## use custom color and thresholds
#' visualize(ildsg,cex=2,pch=19,confcol=rainbow(5),conftol=c(0.9,0.6,0.4))
#' @importFrom Morpho deformGrid2d deformGrid3d
#' @importFrom rgl text3d
#' @rdname visualize
#' @export
visualize <- function(x,...) UseMethod("visualize")

#' @export
#' @rdname visualize
#' @method visualize ILDSR2
visualize.ILDSR2 <- function(x,ref=TRUE,relcol="red",rescol="black",lwd=1,cex=2,col="red",pch=19,confcol=c("green","orange","red"),conftol=c(75,50),useconf=TRUE,...) {
   
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
    D3 <- FALSE
    if (ncol(reference)==3) {
        D3 <- TRUE
        mydeform <- deformGrid3d
    } else
        mydeform <- deformGrid2d
    if (is.null(x$confR2) || ! useconf) {
    highlight <- colnames(x$largeR2)
    if (!is.null(highlight)) {
        hm <- match(highlight,rn)
        mydeform(reference,reference,lines=F,lwd=0,show=1,cex2=0,cex1=cex,col1=col,pch=pch,...)
        mydeform(ref0[hm,,drop=FALSE],ref1[hm,,drop=FALSE],add=T,lcol = relcol,lwd=lwd*3,show=1,cex2=0,cex1=0,...)
        mydeform(ref0[-hm,,drop=FALSE],ref1[-hm,,drop=FALSE],add=T,lcol = rescol,lwd=lwd,show=1,cex2=0,cex1=0,...)
       
    } else {
        mydeform(ref0,ref1)
    } } else {
          highlight <- names(x$confR2)
          mydeform(reference,reference,lines=F,lwd=0,show=1,cex2=0,cex1=cex,col1=col,pch=pch,...)
          hm <- match(highlight,rn)
          myinterval <- getInterval(x$confR2,conftol)
          for (i in 1:3) {
              tmp <- which(myinterval == i)
              if (length(tmp)) {
                  hmtmp <- match(highlight[tmp],rn)
                  mydeform(ref0[hmtmp,,drop=FALSE],ref1[hmtmp,,drop=FALSE],add=T,lcol =confcol[i] ,lwd=lwd*3,show=1,cex2=0,cex1=0,...)
              }
          }          
          mydeform(ref0[-hm,],ref1[-hm,],add=T,lcol = rescol,lwd=lwd,show=1,cex2=0,cex1=0,...)

          
      }
     if (D3) {
            rgl::texts3d(reference,texts = 1:nrow(reference),adj=1.5,...)
        }
        else
            text(reference,adj=2,cex=cex,...)
} 

#' @rdname visualize
#' @export
visualise <- visualize


#' @rdname visualize
#' @export
visualise.ILDSR2 <- visualize.ILDSR2



#' Plot graphical report for ILDSR2
#'
#' Plot graphical report for ILDSR2
#' @param x output of \code{\link{ILDSR2}}
#' @param ... additonal parameter currently not used.
#' @export
plot.ILDSR2 <- function(x,...) {
    par(mfrow=c(2,2))

    plot(x$SILDstats$av.twosh.SILDsortedasratios, x$SILDstats$ratios.twosh.SILD.sorted, main="ILD Ratio Variabilty vs. ILDs Values", xlab="Average of Start & Target ILD", ylab="Target/Start ILD Ratio")
    abline(a=1, b=0, col="grey", lwd=3, lty=1) 

    plot(x$SILDstats$av.twosh.SILDsorted,x$allR2, main="R2-Values vs. SILD Values", xlab="Average of Start & Target ILD", ylab="R2 for Sample ILDs vs Predictor")
    abline(a=quantile(x$allR2, probs=x$R2tol), b=0, col="grey", lwd=3, lty=1)

    hist(x$SILDstats$ratios.twosh.SILD.sorted, breaks=sqrt(length(x$SILDstats$ratios.twosh.SILD.sorted)), prob=TRUE, main="Disribution of Target/Start ILD ratios",xlab="Target/Start ILD Ratios")
    lines(density(x$SILDstats$ratios.twosh.SILD.sorted), col="red")

    hist(x$allR2, breaks=sqrt(length(x$allR2)), prob=TRUE, main=" R2-Value Distribution",xlab="R2-Values")
    lines(density(x$allR2), col="red")
    par(mfrow=c(1,1))
}


getInterval <- function(x,intervals) {    
    tmp <- sapply(x,function(x) x > intervals)
    tmp <- rbind(tmp,TRUE)
    chk <- apply(tmp,2,function(x) x <- which(x))
    chk <- sapply(chk,min)
    return(chk)
        
}
