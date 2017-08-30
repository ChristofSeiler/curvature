#' Sectional curvature of a rstan object for convergence diagnosis.
#'
#' @import rstan
#' @import numDeriv
#' @export
#'
curvature = function(fit) {
  samples = rstan::extract(fit,permuted = FALSE)
  q = samples[,1,1:d]
  func = function(x) { -log_prob(fit, x) }
  U = apply(q,1,function(x) func(x))
  cat("computing gradient:\n")
  grad_U = lapply(1:nrow(q),function(t) { cat(paste0(t," ")); grad(func,q[t,]); })
  cat("\ncomputing hessian:\n")
  hessian_U = lapply(1:nrow(q),function(t) { cat(paste0(t," ")); hessian(func,q[t,]); })
  params = get_sampler_params(fit)
  H = params[[1]][(T0+1):(T+T0),4]
  K = H-U
  # collect
  sec = sapply(1:length(grad_U),function(t) {
    secs = replicate(noOfSecs,
                     sectional_curvature(grad_U[[t]],hessian_U[[t]],K[t]))
    mean(secs)
  })
  class(sec) = "curvature"
  sec
}
