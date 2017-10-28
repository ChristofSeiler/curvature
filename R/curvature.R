#' Sectional curvature of a rstan object for convergence diagnosis.
#'
#' @import rstan
#' @import numDeriv
#' @import parallel
#' @export
#'
curvature = function(fit,
                     num_secs = 1000,
                     num_cores = 4) {

  # extract  sectional curvature
  samples = rstan::extract(fit,permuted = FALSE)
  if(length(dimnames(samples)$chains) != 1) stop("run rstan::sampling with chains=1")
  d = dim(samples)[3]-1 # last column is lp__
  q = samples[ ,1,1:d]

  # compute gradient and Hessian in parallel
  func = function(x) { -log_prob(fit, x) }
  U = apply(q,1,function(x) func(x))
  take_time = system.time( grad(func,q[round(nrow(q)/2),]) )
  message("computing gradients will take approx. ",
          round(take_time["elapsed"]*nrow(q)/num_cores),
          " secs.")
  grad_U = mclapply(1:nrow(q),function(t) { grad(func,q[t,]) },
                    mc.cores =  num_cores)
  take_time = system.time( hessian(func,q[round(nrow(q)/2),]) )
  message("computing Hessians will take approx. ",
          round(take_time["elapsed"]*nrow(q)/num_cores),
          " secs.")
  hessian_U = mclapply(1:nrow(q),function(t) { hessian(func,q[t,]) },
                       mc.cores =  num_cores)

  # compute total and kintetic energy
  params = get_sampler_params(fit)
  args = fit@stan_args[[1]]
  H = params[[1]][(args$warmup+1):args$iter,"energy__"]
  K = H-U

  # collect
  sec = sapply(1:length(grad_U),function(t) {
    secs = replicate(num_secs,
                     sectional_curvature(grad_U[[t]],hessian_U[[t]],K[t]))
    mean(secs)
  })

  # define class for plot and summary
  class(sec) = "curvature"
  sec
}
