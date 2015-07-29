# HLM for coef friction vs volume

from numpy import *
import time

# load data
clm = recfromtxt('Clm_Vol_CoefFriction.txt');
mrp = recfromtxt('Mrp_Vol_CoefFriction.txt');
shv = recfromtxt('Shv_Vol_CoefFriction.txt');
unz = recfromtxt('Unz_Vol_CoefFriction.txt');
smr = recfromtxt('Smr_Vol_CoefFriction.txt');
alldata = [clm,mrp,shv,unz,smr]

name = ['Clm','Mrp', 'SHV','Unz','Smr']
nvolc = len(alldata)

mid = 10**5.5
polys = zeros([nvolc,2])
njs = zeros([nvolc,1])
xj = zeros([nvolc,1])
Sj_lst = zeros([nvolc,1])

# Data manipulation
for i in range(nvolc):
    alldata[i][:,1] = alldata[i][:,1] * 10**6
    alldata[i][:,1] = alldata[i][:,1]/mid
    alldata[i] = log10(alldata[i])
    
    polys[i,:] = polyfit(alldata[i][:,1], alldata[i][:,0],1)
    njs[i] = shape(alldata[i])[0]
    xj[i] = mean(alldata[i][:,1])
    Sj_lst[i] = sum((alldata[i][:,1] - xj[i])**2)

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
sigj2init = map(lambda x: var(x[:,1],ddof=1),alldata)
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
            xji = alldata[i][:,1]
            yji = alldata[i][:,0]
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
            xji = alldata[i][:,1]
            yji = alldata[i][:,0]
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































