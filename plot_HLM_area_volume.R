# plotting results from HLM a

library(ggplot2)
th1 = read.table('results/thetas1_av.Rdata',header=1,sep=',')
th2 = read.table('results/thetas2_av.Rdata',header=1,sep=',')
mus = read.table('results/mus_av.Rdata',header=1,sep=',')
sj2 = read.table('results/sj2_lst_av.Rdata',header=1,sep=',')
st2 = read.table('results/st2_lst_av.Rdata',header=1,sep=',')
alldata = readRDS('alldata_av.Rdata')

vnames= c('Colima','Merapi','Soufriere Hills','Unzen','Augustine');

polyval <- function(coefs, x){
  yh = matrix(0,nrow=length(x),ncol=1)
  for (i in 1:length(x)){
    yh[i] = coefs[1] + coefs[2] * x[i]
  }
  return(yh)
}


j = 4; ################################################## volcano index
x = alldata[[j]][,1]
y = alldata[[j]][,2]
n = 10
vv = seq(from=-3.2,to=-1,length.out=n)

p975 = matrix(0,nrow=n,ncol=1)
p25 = matrix(0,nrow=n,ncol=1)
p50 = matrix(0,nrow=n,ncol=1)

for (i in 1:n){
  cf = th2[,j] * vv[i] + th1[,j]
  p25[i] = quantile(cf,probs=0.025)
  p975[i] = quantile(cf,probs=0.975)
  p50[i] = quantile(cf,probs=.5)
}

nn = length(x)
sxx = sum((x - mean(x))^2)
p = lm(y~x)
pvv = polyval(p$coefficients,vv)
yihat = predict(p,data.frame(x))
seps2 = sum((y - yihat)^2)/(nn-2)


topy = matrix(0,nrow=n,ncol=1)
boty = matrix(0,nrow=n,ncol=1)
for (i in 1:n){
  val = qt(0.975,df=nn) * sqrt(seps2) * sqrt(1/nn + (vv[i] - mean(x))^2/sxx)
  topy[i] = pvv[i]  + val
  boty[i] = pvv[i]  - val
}

a = ggplot()
a + geom_point(aes(x=x,y=y,size=10)) + 
  geom_line(aes(x=vv,y=p25),linetype=6,col='red') + 
  geom_line(aes(x=vv,y=p975),linetype=6,col='red')+
  geom_line(aes(x=vv,y=pvv), linetype=1,col='black') +
  geom_line(aes(x=vv,y=p50), linetype=1,col='red') +
  xlab('Volume (m^3)') + 
  ylab('Coefficient of Friction') + ggtitle(vnames[j]) #+
#  coord_cartesian(xlim=c(-1.5,2.5),ylim=c(-2,2))


























