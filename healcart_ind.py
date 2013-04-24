import numpy as np

def healcart_ind(data):

    npix  = data.size
    nside = long(round(np.sqrt(npix/12.)))

    nx  = 8l*nside
    x   = np.arange(nx)
    ind = np.zeros((nx, 4*nside-1), dtype=np.int64)

    for i in np.arange(0L,nside-1,1):
        ind[:,i] = 2*i*(i+1)+4*(i+1)*x/nx

    xoffs = 2L*nside*(nside-1)

    for i in np.arange(0L,2*nside+1,1):
        ind[:,i+nside-1] = np.roll(xoffs+4*nside*i+x/2, -(i % 2))

    for i in np.arange(1L,nside,1):
        ind[:,4*nside-1-i] = npix-2*i*(i+1)+4*i*x/nx

    ind = np.rot90(np.rot90(ind))
    ind = np.roll(ind, 4*nside, axis=0)
    ind = np.rot90(ind)

    return ind
