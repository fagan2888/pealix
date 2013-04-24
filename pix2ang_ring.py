import numpy as np

def pix2ang_ring(nside, ipix):

# -------- number in each ring
    nside  = long(nside)
    npix   = 12L*nside*nside
    tind   = np.arange(npix)
    nq0    = (np.arange(nside)+1)*4
    nq     = np.concatenate([nq0,np.zeros(2*nside-1,np.int64)+(4*nside), \
                                 nq0[::-1]])
    thetaq = np.zeros(4L*nside-1)
    bound  = 2*nside*(1+nside)

    qind    = np.arange(nside)
    qindrev = nside - qind
    xoffs   = 2L*nside*(nside+1)
# -------- ring0 is the lowest pixel number in each ring
    ring0   = np.concatenate([2*qind*(qind+1), \
                                  xoffs+(4*nside)*np.arange(2*nside-1), \
                                  npix-2*qindrev*(qindrev+1)])

# -------- initialize arrays for output
    theta = np.zeros(npix)
    phi   = np.zeros(npix)


    pp   = np.byte(tind < bound)
    w_pp = np.where(pp)[0]
    n_pp = w_pp.size

    if n_pp > 0:
        qind = np.arange(nside)
        thetaq[qind] = np.arccos(1.0-((qind+1)/np.float(nside))**2/3.)

        q = np.array((np.sqrt(0.25 + 0.5*tind[w_pp])-0.5),dtype=int)
        theta[w_pp] = thetaq[q]
        phi[w_pp] = (tind[w_pp]-ring0[q]+0.5)/(nq/(2.0*np.pi))[q]


    w_eq = np.where((pp==0) & (tind<(npix-bound)))[0]
    n_eq = w_eq.size

    if n_eq>0:
        qind = np.arange(2*nside-1)+nside
        thetaq[qind] = np.arccos((2.0-(qind+1)/np.float(nside))*2./3.)

        q = (tind[w_eq]-bound)/(4*nside)+nside
        theta[w_eq] = thetaq[q]

# GGD: this line has been modified correctly... I think...
#        phi[w_eq] = (tind[w_eq]-ring0[q]+(q!=1)*0.5)/(4*nside/(2.0*np.pi))[q]
        phi[w_eq] = (tind[w_eq]-ring0[q]+(q!=1)*0.5)/(4*nside/(2.0*np.pi))


    w_sp = np.where(tind > (npix-bound))[0]
    n_sp = w_sp.size

    if n_sp > 0:
        qind = np.arange(nside) + (3*nside-1)
        thetaq[qind] = np.arccos(-1.0+((qind+1-4*nside)/np.float(nside))**2/3.)

        q = (4*nside-2) - np.array(np.sqrt(0.5*((0.5 + npix-1)- \
                                                    tind[w_sp]))-0.5, \
                                       dtype=long)

        theta[w_sp] = thetaq[q]
        phi[w_sp] = (tind[w_sp]-ring0[q]+0.5)/(nq/(2.0*np.pi))[q]


# -------- Now lookup descired pixels from table
    return theta[ipix], phi[ipix]

