# 预检查回归测试，由preflight.R在真实SOCK和dream测试前调用；仅使用base R。
run_preflight_regression <- function() {
  # 重现原缺陷：x+1得到double，而2:3是integer，两者数值相同。
  result <- unlist(lapply(1:2,function(x)x+1))
  stopifnot(typeof(result)=='double',typeof(2:3)=='integer',
            !identical(result,2:3),same_numeric_values(result,2:3))
  # 允许存储类型变化，不放过错误值、错误长度、非数值和非有限值。
  stopifnot(!same_numeric_values(c(2,4),2:3),
            !same_numeric_values(c(2,3,4),2:3),
            !same_numeric_values(c('2','3'),2:3),
            !same_numeric_values(c(2,NA_real_),2:3),
            !same_numeric_values(c(2,Inf),2:3))
  x <- matrix(1:4,2,dimnames=list(c('g1','g2'),c('s1','s2')))
  y <- matrix(as.numeric(1:4),2,dimnames=dimnames(x))
  stopifnot(same_numeric_values(x,y),!same_numeric_values(x,as.numeric(y)))
  # 数值一样而列名改变仍不通过counts核对。
  colnames(y) <- rev(colnames(y))
  stopifnot(!identical(dimnames(x),dimnames(y)))
  cat('预检查数值类型与负例回归测试通过\n')
}
