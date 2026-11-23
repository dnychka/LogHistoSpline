constrainedMatern<- function( s,t, a=1,b=1){

 u<- s/b
 v<- t/b

a<- a/b

d<- abs(u-v)/a
s1<- (1+ d)*exp( -d)
um1<- abs(u-1)/a
vm1<- abs( v-1)/a
s2<- exp( -um1 -vm1)*( (1+ um1)* (1+ vm1) + um1*vm1)
cov<- (s1-s2)*(s<=b)*(t<=b)
return(cov)
}




