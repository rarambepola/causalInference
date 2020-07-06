independence_test_wrapper <- function(x, y, alpha, independence_test, ...){
  return(independence_test(x, y, ...) < alpha)
}

CI_test_wrapper <- function(x, y, z, alpha, CI_test, ...){
  return(CI_test(x, y, z, ...) < alpha)
}

