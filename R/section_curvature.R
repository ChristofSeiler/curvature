#' Internal function to compute sectional curvature from Hamiltonian Monte Carlo samples.
#'
#' @import numDeriv
#'
sectional_curvature = function(gradU,hessU,Kp) {
  plane = sample_orthonormal_tangent_vectors(length(gradU))
  u = plane[,1]
  v = plane[,2]

  normGrad2 = euclidean_length(gradU)^2
  alpha = calculate_angle(gradU,u)
  beta = calculate_angle(gradU,v)
  A = 1/(4*Kp^2)*t(u)%*%hessU%*%u
  B = 1/(4*Kp^2)*t(v)%*%hessU%*%v
  C = 3/(8*Kp^3)*normGrad2*cos(alpha)^2
  D = 3/(8*Kp^3)*normGrad2*cos(beta)^2
  E = 1/(8*Kp^3)*normGrad2

  return(A+B+C+D-E)
}

# some helper functions
euclidean_length = function(x) {
  return( norm(as.matrix( x ),"F") )
}
calculate_angle = function(x1,x2) {
  return( acos( (t(x1)%*%x2)/(euclidean_length(x1)*euclidean_length(x2)) ) )
}
sample_orthonormal_tangent_vectors = function(n) {
  A = matrix(data=rnorm(n,0,1),nrow=n,ncol=2)
  return(qr.Q(qr(A)))
}
