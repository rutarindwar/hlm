# plotting results from HLM

library(ggplot2)
th1 = read.table('thetas1.Rdata',header=1,sep=',')
th2 = read.table('thetas2.Rdata',header=1,sep=',')
mus = read.table('mus.Rdata',header=1,sep=',')
sj2 = read.table('sj2_lst.Rdata',header=1,sep=',')
st2 = read.table('st2_lst.Rdata',header=1,sep=',')

vnames= c('Colima','Merapi','Soufriere Hills','Unzen','Semeru');
clm = read.table('Clm_Vol_CoefFriction.txt',header=0);
mrp = read.table('Mrp_Vol_CoefFriction.txt',header=0);
shv = read.table('Shv_Vol_CoefFriction.txt',header=0);
unz = read.table('Unz_Vol_CoefFriction.txt',header=0);
smr = read.table('Smr_Vol_CoefFriction.txt',header=0);
alldata = list(clm=clm,mrp=mrp,shv=shv,unz=unz,smr=smr);
mid = 10^5.5
nvolc = 5
for (i in 1:nvolc){
  alldata[[i]][,2] = alldata[[i]][,2] * 1e6;
  alldata[[i]][,2] = alldata[[i]][,2]/mid;
  alldata[[i]] = log10(alldata[[i]]);
}




#p = ggplot(th1) + geom_histogram(binwidth=5) + facet_grid(year~.)


polyval <- function(coefs, x){
  yh = matrix(0,nrow=length(x),ncol=1)
  for (i in 1:length(x)){
    yh[i] = coefs[1] + coefs[2] * x[i]
  }
  return(yh)
}



j = 5; ########## volcano index
x = alldata[[j]][,2]
y = alldata[[j]][,1]
n = 10
vv = seq(from=-1.5,to=2.5,length.out=n)

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
    ylab('Coefficient of Friction') + ggtitle(vnames[j]) +
    coord_cartesian(xlim=c(-1.5,2.5),ylim=c(-2,2))











############################################################################

melt.th1 = melt(th1);


plt = ggplot(melt.th1,aes(x=value)) + 
      geom_density() + 
      scale_x_continuous(breaks=seq(-0.7,0.1,1000)) +
      facet_grid(variable~.)

last_plot()



















