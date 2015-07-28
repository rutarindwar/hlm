# Version 5 of the hierarchical model for area vs volume
# Colima, and Merapi are unchanneled, whereas shv, unzen and smr are 
# channeled


# Read in data files.
clm = read.table('Clm_Vol_area.txt',header=F);
mrp = read.table('Mrp_Vol_area.txt',header=F);
shv = read.table('Shv_Vol_area.txt',header=F);
unz = read.table('Unz_Vol_area.txt',header=F);
aug = read.table('Aug_Vol_area.txt',header=F);
alldata = list(clm=clm,mrp=mrp,shv=shv,unz=unz,aug=aug);

fnames= c('clm','mrp','shv','unz','aug');
snames = names(alldata);
nvolc = length(alldata);

polys = data.frame(matrix(nrow=2,ncol=5));
colnames(polys) = snames;
njs = rep(0,nvolc);
xj = rep(0,nvolc);
svar = rep(0,nvolc);
for (i in 1:nvolc){
  alldata[[i]][,1] = alldata[[i]][,1]/1000;
  alldata[[i]][,1] = alldata[[i]][,1]^(2/3);
  alldata[[i]] = log10(alldata[[i]]);
  
  linm = lm(alldata[[i]][,2]~alldata[[i]][,1]);
  polys[i] = linm$coefficients;
  njs[i] = length(alldata[[i]][,1]);
  xj[i] = mean(alldata[[i]][,1]);
  svar[i] = sum((alldata[[i]][,1] - xj[i])^2);
}


nsamples = 10000; ########################################### sampling size
nflows = sum(njs);

# Allocating spaces for the dataframes and arrays.
thetas1 = data.frame(matrix(0,nrow=nsamples,ncol=nvolc));
colnames(thetas1) = snames;
thetas2 = data.frame(matrix(0,nrow=nsamples,ncol=nvolc));
colnames(thetas2) = snames;
mus = rep(0,nsamples);
st2.lst = rep(0,nsamples);
sj2.lst = data.frame(matrix(0,nrow=nsamples,ncol=nvolc));
lams = rep(0,nsamples);

# Initializing first values in matrices
t1 = as.numeric(polys[1,]);
t2 = as.numeric(polys[2,]);
thetas1[1,] = t1;
thetas2[1,] = t2;
mus[1] = mean(t2);
sj2init = lapply(alldata,function(x) var(x[,1]));

sj2.lst[1,] = as.numeric(sj2init);
st2.lst[1] = var(t2);

svar1 = svar[1:2];
svar2 = svar[3:nvolc];


vj.func <- function(s.j2,s.t2,Sj.part){(s.j2/Sj.part) + s.t2};
mu.hat.func <- function(t2j,s.j2,s.t2,Sj.part){
  
  sum(t2j/vj.func(s.j2,s.t2,Sj.part))/sum(1/vj.func(s.j2,s.t2,Sj.part))};

mv.rnorm <- function(means,sds){
  n = length(means);
  samples = numeric(length=n);
  for(i in 1:n){
    samples[i] = rnorm(n=1,mean=means[i],sd=sds[i])
  }
  return(samples)
}




q = nsamples/10;
tmp = rep(0,nvolc);

k = 1;
for (ns in 1:(nsamples-1)){
  
  # progress notification
  if (((ns%%q) - 0) < 1e-3){
    #sprintf('%3.1f percent done',ns*10/q);
    #print(ns*10/q); print ('\b% percent done')
    cat(ns*10/q, "% percent done\n")
    flush.console();
  }
  
  ##############################################################
  
  # Step 1: Draw th2
  lamtemp = st2.lst[k] / sj2.lst[k,];
  m2 = t2 - (t2 - mus[k])/(1 + lamtemp * svar);
  s21 = (sj2.lst[k,] * lamtemp)/(1 + lamtemp * svar);
  thetas2[k+1,] = mv.rnorm(as.numeric(m2),as.numeric(sqrt(s21)));
  
  ##############################################################
  
  # Step 2: Draw th1
  m1 = t1 - xj*(thetas2[k+1,] - t2); # use new t2
  s22 = sj2.lst[k,]/njs;
  thetas1[k+1,] = mv.rnorm(as.numeric(m1),as.numeric(sqrt(s22)));
  
  ############################################################## 
  
  # Step 3A: Draw sig2.u
  cu = 0;
  found.u = 0;
  while (found.u != 1){
    btmp = 0;
    for (i in 1:2){
      xji = alldata[[i]][,1];
      yji = alldata[[i]][,2];
      btmp = btmp + sum((yji - (thetas1[k+1,i] + 
                                  thetas2[k+1,i]*xji))^2);
    }
    betau = 0.5 * btmp;
    s2.u = 1/rgamma(n=1,shape=sum(njs[1:2])/2,scale=1/betau);
    
    
    n1 = sum(1/vj.func(as.numeric(sj2.lst[k,]),st2.lst[k],svar)^2);
    b1 = sum(1/vj.func(0,st2.lst[k],svar1)^2);
    b2 = sum(1/vj.func(as.numeric(sj2.lst[k,3:nvolc]),st2.lst[k],svar2)^2);
    
    u.u = runif(1);
    condu = sqrt(n1)/sqrt(b1 + b2);
    if (u.u < condu){
      found.u = 1;
    }
    cu = cu + 1;
  }
  
  ##############################################################
  
  # STEP 3B: Draw sig2.c
  cc = 0;
  found.c = 0;
  
  while (found.c != 1){
    btmp = 0;
    for (i in 3:nvolc){
      xji = alldata[[i]][,1];
      yji = alldata[[i]][,2];
      btmp = btmp + sum((yji - (thetas1[k+1,i] + 
                                  thetas2[k+1,i]*xji))^2);
    }
    betac = 0.5 * btmp;
    s2.c = 1/rgamma(n=1,shape=sum(njs[3:nvolc])/2,scale=1/betac);
    
    n1 = sum(1/vj.func(as.numeric(sj2.lst[k,]),st2.lst[k],svar)^2);
    b1 = sum(1/vj.func(as.numeric(sj2.lst[k,1:2]),st2.lst[k],svar1)^2);
    b2 = sum(1/vj.func(0,st2.lst[k],svar2)^2);
    
    u.c = runif(1);
    condc = sqrt(n1)/sqrt(b1 + b2);
    
    if(u.c < condc){
      found.c = 1;
    }
    cc = cc + 1;
  }
  # Put s2.c and s2.u into the old sj2.lst array
  for (i in 1:2){
    sj2.lst[k+1,i] = s2.u;
  }
  for (i in 3:nvolc){
    sj2.lst[k+1,i] = s2.c;
  }
  
  ##############################################################
  
  # STEP 4: Draw mu
  vs = vj.func(sj2.lst[k+1,], st2.lst[k],svar);
  ms = mu.hat.func(as.numeric(thetas2[k+1,]),
                   as.numeric(sj2.lst[k+1,]),
                   st2.lst[k],svar);
  mus[k+1] = rnorm(n=1,mean=ms,sd=sqrt(1/sum(1/vs)));
  
  ##############################################################
  
  #STEP 5: Generate sig.t2
  beta = 0.5*sum((thetas2[k+1,] - mus[k+1])^2);
  sig.t2.cand = 1/rgamma(n=1,shape=(nvolc-2)/2,scale=1/beta);
  lamcand = sig.t2.cand/sj2.lst[k+1,];
  
  ut = runif(1);
  condt = sqrt(sum(1/vj.func(as.numeric(sj2.lst[k+1,]),sig.t2.cand,svar)^2))/
    sqrt(sum(1/vj.func(as.numeric(sj2.lst[k+1,]),0,svar)^2));
  
  ct = 0;
  while (ut < condt){
    beta = 0.5 * sum((thetas2[k+1,] - mus[k+1])^2);
    sig.t2.cand = 1/rgamma(n=1,shape=(nvolc-2)/2,scale=1/beta);
    lamcand = sig.t2.cand/sj2.lst[k+1,];
    
    ut = runif(1);
    condt = sqrt(sum(1/vj.func(as.numeric(sj2.lst[k+1,]),sig.t2.cand,svar)^2))/
      sqrt(sum(1/vj.func(as.numeric(sj2.lst[k+1,]),0,svar)^2));
    ct = ct + 1;
  }
  st2.lst[k+1] = sig.t2.cand;
  k = k + 1;
}

# number of samples to throw away
burn = 1000;
thetas1 = thetas1[-(1:burn),];
thetas2 = thetas2[-(1:burn),];
mus = mus[-(1:burn)]
st2.lst = st2.lst[-(1:burn)];
sj2.lst = sj2.lst[-(1:burn),];



# save samples
write.table(thetas1,file='results/thetas1_av.Rdata',row.names=F,sep=',')
write.table(thetas2,file='results/thetas2_av.Rdata',row.names=F,sep=',')
write.table(sj2.lst,file='results/sj2_lst_av.Rdata',row.names=F,sep=',')
write.table(mus,file='results/mus_av.Rdata',row.names=F,sep=',')
write.table(st2.lst,file='results/st2_lst_av.Rdata',row.names=F,sep=',')









