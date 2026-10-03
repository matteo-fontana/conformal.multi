n=20
p=2
q=2
mu=rep(0,p)
x = mvtnorm::rmvnorm(n, mu)
beta<-sapply(1:q, function(k) c(mvtnorm::rmvnorm(1,mu)))
y = x%*%beta + t(mvtnorm::rmvnorm(q,1:n))
x0=x[n,,drop=FALSE]
funs=lm_multi()

sol<-conformal.multidim.jackplus(x[-n,],y[-n,],x0,train.fun = funs$train.fun,
                                 predict.fun = funs$predict.fun)
