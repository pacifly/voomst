#' transform mean-variance on log-cpm scale to count scale
#'
#' transform the voom.xy component (log2-count, sqrt of residual standard deviation on log2-cpm)
#' in the object returned from limma::voom to mean and dispersion on the count scale
#'
#' @param voom.xy the voom.xy component in the object returned from limma::voom
#' @param logbase base for the logorithm in \code{voom.xy}
#' @return list of transformed mean and dispersion
#' @export
transfvoomxy = function(voom.xy, logbase=2) {
  out = list()
  out$x = logbase^(voom.xy$x)
  out$y = (voom.xy$y)^4 * (log(logbase))^2 - 1 / out$x
  out
}

#' variance stabilizing transformation (VST) similar as DEXSeq
#'
#' @param a0,a1 coefficients for the mean and dispersion relationship (d ~ a0 + a1/m),
#'        as fitted via DESeq2:::parametricDispersionFit
#' @param meth method for deriving the analytical form of the transformation
#' @return VST function with orginal count as input and return transformed count
#' @export
vstrans = function(a0, a1, meth=c('default', 'DEXSeq')) {
  meth = match.arg(meth)
  function(count) {
    if (meth == 'DEXSeq') {
      ## from DEXSeq
      ## https://genome.cshlp.org/content/suppl/2012/08/20/gr.133744.111.DC1/Supplement.pdf
      return(
        2 / sqrt(a0) * log(2 * a0 * sqrt(count) + 2 * sqrt(a0 * (a0 * count + a1 + 1)))
      )
    }
    if (meth == 'default') {
      ## self-derived, differs from DEXSeq up to a constant
      return(
        2 / sqrt(a0) * log(sqrt(a0 * count) + sqrt(a0 * count + a1 + 1))
      )
    }
  }
}

#' inverse variance stabilizing transformation (VST)
#'
#' @inheritParams vstrans
#' @param vstfn a VST function if \code{meth} is 'numeric'
#' @param ... additional parameters passed to \code{invfn}
#' @return inverse VST function with transformed count as input and return original
#'         untransformed count
#' @rdname vstrans
#' @export
ivstrans = function(a0, a1, meth=c('default', 'DEXSeq', 'numeric'), vstfn=NULL, ...) {
  meth = match.arg(meth)
  function(count) {
    if (meth == 'DEXSeq') {
      ## correspond to DEXSeq
      b = exp(sqrt(a0) * (count - 1 / sqrt(a0) * log(4 * a0)))
      return(
        1 / 4 / a0 * (b + 1 / b * (a1 + 1)^2 - 2 * a1 - 2)
      )
    }
    if (meth == 'default') {
      ## self-derived
      b = exp(sqrt(a0) * count)
      return(
        1 / 4 / a0 * (b + 1 / b * (a1 + 1)^2 - 2 * a1 - 2)
      )
    }
    if (meth == 'numeric') {
      return( invfn(vstfn, ...)(count) )
    }
  }
}

#' numerically invert a function
#'
#' @param fn a function to be inverted
#' @param ... additional parameters passed to \code{uniroot}
#' @return an inverse function of \code{fn}
#' @export
invfn = function(fn, ...){
  Vectorize(function(y) {
    uniroot(f=function(x) {fn(x) - y}, ...)$root
  })
}

