# HLM for area vs volume

from numpy import *
import time
import matplotlib.pyplot as plt
from scipy import stats

# load data
clm = recfromtxt('Clm_Vol_area.txt');
mrp = recfromtxt('Mrp_Vol_area.txt');
shv = recfromtxt('Shv_Vol_area.txt');
unz = recfromtxt('Unz_Vol_area.txt');
aug = recfromtxt('Aug_Vol_area.txt');
alldata = [clm,mrp,shv,unz,aug]

name = ['Clm','Mrp', 'SHV','Unz','Aug']
nvolc = len(alldata)
polys = zeros([nvolc,2])
njs = zeros([nvolc,1])
xj = zeros([nvolc,1])
Sj_lst = zeros([nvolc,1])

# Data manipulation: Remember col1 = x, col2 = y
for i in range(nvolc):
    alldata[i][:,0] = alldata[i][:,0]/1000 # change to km^3
    alldata[i][:,0] = alldata[i][:,0]**(2.0/3)
    alldata[i] = log10(alldata[i])
    
    polys[i,:] = polyfit(alldata[i][:,0], alldata[i][:,1],1)
    njs[i] = shape(alldata[i])[0]
    xj[i] = mean(alldata[i][:,0])
    Sj_lst[i] = sum((alldata[i][:,0] - xj[i])**2)

Sj_lst = squeeze(Sj_lst)
njs = squeeze(njs)
xj = squeeze(xj)


# sample size
nsamples = 10**5 #===========================================================
nflows = sum(njs)

# Allocating spaces for arrays
thetas1= zeros([nvolc,nsamples])
thetas2= zeros([nvolc,nsamples])
mus = zeros([nsamples,1])
st2_lst = zeros([nsamples,1])
sj2_lst = zeros([nvolc,nsamples])
lams = zeros([nvolc,nsamples])

# Initializing first values of matrices
t1 = polys[:,1]
t2 = polys[:,0]
thetas1[:,0] = squeeze(t1)
thetas2[:,0] = squeeze(t2)
mus[0] = mean(t2)
sigj2init = map(lambda x: var(x[:,0],ddof=1),alldata)
sj2_lst[:,0] = sigj2init
st2_lst[0] = var(t2,ddof=1)


Sj_lst1 = Sj_lst[:2]
Sj_lst2 = Sj_lst[2:]

vj_func = lambda s_j2,s_t2,Sj_part: (1.0*s_j2/Sj_part) + s_t2

mu_hat_func = lambda t2j, s_j2,s_t2,Sj_part: \
                    sum(t2j/vj_func(s_j2,s_t2,Sj_part))/ \
                    sum(1/vj_func(s_j2,s_t2,Sj_part))

tic = time.time()
q = nsamples/10

foo = zeros([nsamples,2])
for k in range(0,nsamples-1):
    
    # progress report
    if (k%q == 0):
        print "%3.1f %% done"%(k*10/q)
    
    ##########################################################################
    
    # STEP 1: Draw th2
    lamtemp = st2_lst[k]/sj2_lst[:,k]
    m1 = t2 - (t2 - mus[k])/(1 + lamtemp*Sj_lst)
    s21 = (sj2_lst[:,k]*lamtemp)/(1 + lamtemp*Sj_lst)
    thetas2[:,k+1] = squeeze(random.randn(nvolc,1))*sqrt(s21) + m1 

    ##########################################################################
    
    # STEP 2: Draw th1
    m = t1 - xj*(thetas2[:,k+1]-t2) #use new t2
    s22 = sj2_lst[:,k]/njs
    thetas1[:,k+1] = squeeze(random.randn(nvolc,1))*sqrt(s22) + m 

    ##########################################################################

    # STEP 3A: Draw sig2_u
    cu = 0
    found_u = 0
    while (found_u != 1):
        btmp = 0
        for i in range(2):
            xji = alldata[i][:,0]
            yji = alldata[i][:,1]
            btmp = btmp + sum((yji-(thetas1[i,k+1] + xji*thetas2[i,k+1]))**2)

        betau = 0.5*btmp
        s2_u = 1/random.gamma(sum(njs[0:2])/2, 1.0/betau)

        n1 = sum(1/vj_func(sj2_lst[:,k],st2_lst[k],Sj_lst)**2)
        b1 = sum(1/vj_func(0,st2_lst[k],Sj_lst1)**2)
        b2 = sum(1/vj_func(sj2_lst[2:,k],st2_lst[k],Sj_lst2)**2)

        u_u = random.rand()
        condu = sqrt(n1)/sqrt(b1 + b2)

        if (u_u < condu):
            found_u = 1
        cu = cu + 1

    ##########################################################################

    # STEP 3B: Draw sig2_c
    cc = 0
    found_c = 0
    while (found_c != 1):
        btmp = 0
        for i in range(2,nvolc): #i=3:nvolc
            xji = alldata[i][:,0]
            yji = alldata[i][:,1]
            btmp = btmp + sum((yji-(thetas1[i,k+1]+xji*thetas2[i,k+1]))**2)

        betac = 0.5*btmp
        s2_c = 1/random.gamma(sum(njs[2:])/2, 1.0/betac)
        
        n1 = sum(1/vj_func(sj2_lst[:,k],st2_lst[k],Sj_lst)**2)
        b1 = sum(1/vj_func(sj2_lst[:2,k],st2_lst[k],Sj_lst1)**2)
        b2 = sum(1/vj_func(0,st2_lst[k],Sj_lst2)**2)

        u_c = random.rand()
        condc = sqrt(n1)/sqrt(b1 + b2)
        
        if (u_c < condc):
            found_c = 1
        
        cc = cc + 1


    foo[k,0] = cu
    foo[k,1] = cc
    
    for i in range(2):
        sj2_lst[i,k+1] = s2_u
    for i in range(2,5):
        sj2_lst[i,k+1] = s2_c
    
    ##########################################################################
    
    # STEP 4: Draw mu
    vs = vj_func(sj2_lst[:,k+1],st2_lst[k],Sj_lst)
    ms = mu_hat_func(thetas2[:,k+1],sj2_lst[:,k+1],st2_lst[k],Sj_lst)
    mus[k+1] = random.randn(1)*sqrt(1/sum(1/vs)) + ms
    
    ##########################################################################
    
    # STEP 5: Generate sig_t2
    beta = 0.5*sum((thetas2[:,k+1]-mus[k+1])**2)
    sig_t2_cand = 1/random.gamma((nvolc-2)/2, 1.0/beta)
    u_t = random.rand()
    condt = sqrt(sum(1/vj_func(sj2_lst[:,k+1],sig_t2_cand,Sj_lst)**2))/ \
            sqrt(sum(1/vj_func(sj2_lst[:,k+1],0,Sj_lst)**2))

    ct = 0;
    
    while ((u_t < condt) == 0):
        beta = 0.5*sum((thetas2[:,k+1]-mus[k+1])**2)
        sig_t2_cand = 1/random.gamma((nvolc-2)/2, 1.0/beta)

        u_t = random.rand()
        condt = sqrt(sum(1/vj_func(sj2_lst[:,k+1],sig_t2_cand,Sj_lst)**2))/ \
                sqrt(sum(1/vj_func(sj2_lst[:,k+1],0,Sj_lst)**2))
        
        ct = ct + 1
    
    st2_lst[k+1] = sig_t2_cand
    


toc = time.time()
print
print '%10.3f secs elapsed'%(toc-tic)


burn = 1000
thetas1 = thetas1[:,burn:]
thetas2 = thetas2[:,burn:]
mus = mus[burn:]
st2_lst = st2_lst[burn:]
sj2_lst = sj2_lst[:,burn:]



plt.figure(1)
plt.clf()
for i in range(nvolc):
    plt.hist(thetas1[i,:], bins=100, histtype='stepfilled', normed=True,\
            alpha=0.5, label=name[i])
plt.title('Theta1')
plt.legend()



plt.figure(2)
plt.clf()
for i in range(nvolc):
    plt.hist(thetas2[i,:], bins=100, histtype='stepfilled', normed=True,\
            alpha=0.5, label=name[i])
plt.title('Theta2')
plt.legend()



#####################################################################

n = 10;

plt.figure(3)
plt.clf()
figind = 1

for j in range(nvolc):#j = 1 ###################### volcano index

    plt.subplot(3,2,figind)
    figind += 1
    x = alldata[j][:,0]
    y = alldata[j][:,1]
    vv = linspace(-3.2, -1, n);

    p975 = zeros([n,1]);
    p25 = zeros([n,1]);
    p50 = zeros([n,1]);
    fvv = zeros([n,1]);
    
    for i in range(n):
        cf = thetas2[j,:] * vv[i] + thetas1[j,:]
        p25[i] = prctile(cf,p=2.5)
        p50[i] = prctile(cf,p=50)
        p975[i] = prctile(cf,p=97.5)

    nn = len(x)
    p = polyfit(x,y,1)
    sxx = sum((x - mean(x))**2)
    yihat = polyval(p,x)
    seps2 = sum((y - yihat)**2)/(nn-2)
    
    topy = zeros([n,1])
    boty = zeros([n,1])
    
    for i in range(n):
        val = stats.t.ppf(0.975,nn) * sqrt(seps2) * \
                sqrt((1/nn) + ((vv[i] - mean(x))**2/sxx));
        topy[i] = polyval(p,vv[i]) + val;
        boty[i] = polyval(p,vv[i]) - val;
    
    
    

    plt.plot(x,y,'bo',ms=5,label='data')
    plt.plot(vv,p25,'r--', label='HLM - CI')
    plt.plot(vv,p975,'r--')
    plt.plot(vv,p50,'r-',label='HLM - Mean')
    plt.plot(vv,topy, 'k:', label='LR - CI')
    plt.plot(vv,boty, 'k:')
    plt.plot(vv,polyval(p,vv), 'k-',label='LR')
    plt.xlim([-3.2,-1])
    plt.ylim([-2.5,2])
    plt.title(name[j])















