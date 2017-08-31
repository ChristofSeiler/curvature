#' Sectional curvature of a rstan object for convergence diagnosis.
#'
#' @import rstan
#' @import numDeriv
#' @export
#'
curvature = function(fit,num_cores = 4) {

  # extract  sectional curvature
  samples = rstan::extract(fit,permuted = FALSE)
  d = dim(samples)[1]
  q = samples[,1,1:d]

  # setup parallel computing
  param = MulticoreParam(workers = num_cores,
                         tasks = d,
                         log = TRUE,
                         logdir = ".",
                         progressbar = TRUE,
                         cleanup = FALSE,
                         seed = 0xdada)

  # compute gradient and Hessian in parallel
  func = function(x) { -log_prob(fit, x) }
  U = apply(q,1,function(x) func(x))
  cat("computing gradient:\n")
  #grad_U = lapply(1:nrow(q),function(t) { cat(paste0(t," ")); grad(func,q[t,]); })
  grad_U = bplapply(1:nrow(q), function(t) grad(func,q[t,]), BPPARAM = param)
  cat("\ncomputing hessian:\n")
  #hessian_U = lapply(1:nrow(q),function(t) { cat(paste0(t," ")); hessian(func,q[t,]); })
  hessian_U = bplapply(1:nrow(q),function(t) hessian(func,q[t,]), BPPARAM = param)
  params = get_sampler_params(fit)
  H = params[[1]][(T0+1):(T+T0),4]
  K = H-U

  # collect
  sec = sapply(1:length(grad_U),function(t) {
    secs = replicate(noOfSecs,
                     sectional_curvature(grad_U[[t]],hessian_U[[t]],K[t]))
    mean(secs)
  })

  # define class for plot and summary
  class(sec) = "curvature"
  sec
}
