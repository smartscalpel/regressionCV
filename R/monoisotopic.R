
########################
## FUNCTIONS ADOPTED FROM MALDIQUANT
#######################

## .unlist
##  wrapper for unlist
##
## params:
##  x: an R object
##
## returns:
##  see also ?unlist
##
.unlist <- function(x) {
  unlist(x, recursive=FALSE, use.names=FALSE)
}

#' .colCors
#'
#' Calculate the correlation for two matrices columnwise.
#'
#' @param x matrix/data.frame
#' @return double
#' @author Sebastian Gibb <mail@@sebastiangibb.de>
#' @noRd
.colCors <- function(x, y, na.rm=FALSE) {
  stopifnot(is.matrix(x) && is.matrix(y))
  stopifnot(all(dim(x) == dim(y)))
  
  if (na.rm) {
    isNA <- is.na(x) | is.na(y)
    x[isNA] <- NA_real_
    y[isNA] <- NA_real_
  }
  
  cmX <- colMeans(x, na.rm=na.rm)
  cmY <- colMeans(y, na.rm=na.rm)
  
  (colMeans(x * y, na.rm=na.rm) - (cmX * cmY)) /
    (sqrt(colMeans(x * x, na.rm=na.rm) - cmX * cmX) *
       sqrt(colMeans(y * y, na.rm=na.rm) - cmY * cmY))
}

#' .pseudoCluster
#'
#' Find possible isotopic cluster in mass/mz data.
#'
#' @param x double, mass
#' @param size integer, cluster size, number of peaks per cluster
#' @param distance double, distance between isotopes (mass of a neutron; see
#' Park et al 2008); could be of length > 1 (if > 1: order will affect later
#' removal in .monoisotopicPattern).
#' @param tolerance double, mass tolerance
#' @return a matrix of indices (nrow(x) == n) of potential clusters
#' @references
#' K. Park, J.Y. Yoon, S. Lee, E. Paek, H. Park, H.J. Jung, and S.W. Lee. 2008.
#' Isotopic peak intensity ratio based algorithm for determination of isotopic
#' clusters and monoisotopic masses of polypeptides from high-resolution
#' mass spectrometric data.
#' Analytical Chemistry, 80: 7294-7303.
#' @noRd
.pseudoCluster <- function(x, size=3L, distance=1.00235, tolerance=1e-4) {
  if (size < 2L) {
    stop("The ", sQuote("size"), " of a cluster has to be at least 2!")
  }
  mm <- matrix(x, nrow=size, ncol=length(x) * length(distance), byrow=TRUE)
  ms <- mm + (rep(distance, each=size) * 0L:(size - 1L))
  
  i <- match.closest(ms, x, tolerance=mm * tolerance)
  dim(i) <- dim(ms)
  
  i[, !is.na(colSums(i)), drop=FALSE]
}

#' .F
#'
#' Map mass to poisson mean/lambda.
#'
#' @param mass double, mass from experimental peak list
#' @return double suitable to pass to `dpois`
#' @references
#' E.J. Breen, F.G. Hopwood, K.L. Williams, and M.R. Wilkins. 2000.
#' Automatic poisson peak harvesting for high throughput protein identification.
#' Electrophoresis 21 (2000): 2243-2251.
#' @noRd
.F <- function(x)0.000594 * x + 0.03091

#' .P
#'
#' Model isotopic distribution by poisson distribution.
#'
#' @param mass double, mass from experimental peak list
#' @param isotopes integer, which isotopes
#' @return double, isotopic distribution
#' @references
#' E.J. Breen, F.G. Hopwood, K.L. Williams, and M.R. Wilkins. 2000.
#' Automatic poisson peak harvesting for high throughput protein identification.
#' Electrophoresis 21 (2000): 2243-2251.
#' @noRd
.P <- function(x, isotopes)dpois(isotopes, .F(x))

#' .Psum
#'
#' Model isotopic distribution by poisson distribution and sum to 1 (similar to
#' TIC).
#'
#' @param mass double, mass from experimental peak list
#' @param isotopes integer, which isotopes
#' @return double, isotopic distribution
#' @references
#' E.J. Breen, F.G. Hopwood, K.L. Williams, and M.R. Wilkins. 2000.
#' Automatic poisson peak harvesting for high throughput protein identification.
#' Electrophoresis 21 (2000): 2243-2251.
#' @noRd
.Psum <- function(x, isotopes) {
  ni <- length(isotopes)
  nx <- length(x)
  p <- .P(rep.int(x, rep.int(ni, nx)), isotopes)
  dim(p) <- c(ni, nx)
  t(t(p) / colSums(p))
}

#' .monoisotopicPattern
#'
#' Model isotopic distribution by poisson distribution.
#'
#' @param x double, mass from experimental peak list
#' @param y double, intensity from experimental peak list
#' @param tolerance double, mass tolerance for .pseudoCluster
#' @param minCor double, minimal correlation between experimental and model
#' intensities
#' @param distance double, distance between isotopes (mass of a neutron; see
#' Park et al 2008); could be of length > 1; if length > 1 the order matters.
#' The first distance elements are prefered (the last elements are possible
#' removed because the contain duplicated indices).
#' @param size integer, cluster size (number of peaks for a possible cluster),
#' see .pseudoCluster
#' @return matrix, index of monoisotopic masses in first row
#' @references
#' E.J. Breen, F.G. Hopwood, K.L. Williams, and M.R. Wilkins. 2000.
#' Automatic poisson peak harvesting for high throughput protein identification.
#' Electrophoresis 21 (2000): 2243-2251.
#' @noRd
.monoisotopicPattern <- function(x, y, minCor=0.95, tolerance=1e-4,
                                 distance=1.00235, size=3L) {
  pc <- .pseudoCluster(x, size=size, distance=distance, tolerance=tolerance)
  y <- y[pc]
  dim(y) <- dim(pc)
  y <- t(t(y)/colSums(y))
  p <- .Psum(x[pc[1L,]], isotopes=0L:(size-1L))
  cr <- .colCors(y, p)
  pc <- pc[, cr > minCor, drop=FALSE]
  pc[duplicated(as.vector(pc))] <- NA_real_
  pc[, !is.na(colSums(pc)), drop=FALSE]
}

#' .monoisotopic
#'
#' Loop through multiple .monoisotopicPattern outputs and remove duplicated
#' peaks.
#'
#' @param x double, mass from experimental peak list
#' @param y double, intensity from experimental peak list
#' @param size integer vector, cluster size
#' @param \ldots further arguments passed to .monoisotopicPattern
#' @return double, index of monoisotopic masses
#' @noRd
.monoisotopic <- function(x, y, size=3L:10L, ...) {
  if (length(x) && length(x) == length(y)) {
    pattern <- lapply(sort.int(size, decreasing=TRUE),
                      function(s).monoisotopicPattern(x=x, y=y, size=s, ...))
    upattern <- .unlist(pattern)
    upattern[duplicated(upattern)] <- NA_real_
    upattern <- relist(upattern, pattern)
    sort.int(.unlist(lapply(upattern,
                            function(p)p[1L, !is.na(colSums(p))])))
  } else {
    double()
  }
}

#' .isotopic
#'
#' Loop through multiple .monoisotopicPattern outputs, remove duplicated
#' peaks, and return isotopic series.
#'
#' @param x double, mass from experimental peak list
#' @param y double, intensity from experimental peak list
#' @param size integer vector, cluster size
#' @param \ldots further arguments passed to .monoisotopicPattern
#' @return double, index of monoisotopic masses
#' @noRd
.isotopic <- function(x, y, size=3L:10L, ...) {
  if (length(x) && length(x) == length(y)) {
    pattern <- lapply(sort.int(size, decreasing=TRUE),
                      function(s).monoisotopicPattern(x=x, y=y, size=s, ...))
    upattern <- .unlist(pattern)
    upattern[duplicated(upattern)] <- NA_real_
    upattern <- relist(upattern, pattern)
    sort.int(.unlist(lapply(upattern,
                            function(p)p[-1L, !is.na(colSums(p))])))
  } else {
    double()
  }
}

.myPeaks<-function(p, minCor=0.95, tolerance=1e-4, size=2L:10L,Cl=FALSE,K=FALSE) {
  niI<-.isotopic(x=MALDIquant::mass(p), y=MALDIquant::intensity(p),
                 minCor=minCor, tolerance=tolerance,
                 distance=1.00235, size=size)
  if(length(niI)>0){
    iP<-p[-niI]
  }else{
    iP<-p
  }
  if(Cl){
  niICl<-.isotopic(x=MALDIquant::mass(iP), y=MALDIquant::intensity(iP),
                   minCor=minCor, tolerance=tolerance,
                   distance=1.99705, size=size)
  if(length(niICl)>0){
    iP<-iP[-niICl]
  }
  }
  if(K){
    niIK<-.isotopic(x=MALDIquant::mass(iP), y=MALDIquant::intensity(iP),
                    minCor=minCor, tolerance=tolerance,
                     distance=1.998119, size=size)
    if(length(niIK)>0){
      iP<-iP[-niIK]
    }
  }
  return(iP)
}

myPeakList<-function(object, minCor=0.95, tolerance=1e-4, size=2L:10L,Cl=FALSE,K=FALSE) {
  if(Cl&K){
    stop('Cl and K adducts are not compatible to each other. Choose one of them.\n')
  }
  cat(format(Sys.time(), "%b %d %X"),'Function: myPeakList(minCor=',minCor,', tolerance=',tolerance,', size=',size,',Cl=',Cl,',K=',K,')\n')
  res<-lapply(object,.myPeaks,minCor=minCor, tolerance=tolerance, size=size,Cl=Cl,K=K)
  cat(format(Sys.time(), "%b %d %X"),'Function: myPeakList finish.\n')
  return(res)
}


.mPeaks<-function(object, minCor=0.95, tolerance=1e-4,
                  distance=1.00235, size=3L:10L) {
  object[.monoisotopic(x=MALDIquant::mass(object), y=MALDIquant::intensity(object),
                       minCor=minCor, tolerance=tolerance,
                       distance=distance, size=size)]
}

.mPeakList<-function(object, minCor=0.95, tolerance=1e-4,
                     distance=1.00235, size=3L:10L) {
  lapply(object,.mPeaks,minCor=minCor, tolerance=tolerance,
         distance=distance, size=size)
}
