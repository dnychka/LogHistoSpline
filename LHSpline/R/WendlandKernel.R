WendlandKernel<- function( s,t, a=1,b=1){
#With $u=(s-b)/a$ and $v=(t-b)/a$, 
# the conditional covariance is
#$$k_c(s,t)=\left(1+\frac{|s-t|}{a}\right)e^{-|s-t|/a}\;-\;e^{-(|u|+|v|)}\Big[(1+|u|)(1+|v|)+uv\Big].$$
u<- s/b
v<- t/b
a<- a/b

d<- abs(u-v)/a

return(Wendland2.2(d))
}




