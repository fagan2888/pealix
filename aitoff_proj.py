import numpy as np
from healgen_lb import *
from aitoff import *

def aitoff_proj(map, nx=None, ny=None, mask=None, oversample=None, \
                    blackbg=None):

# -------- defaults
    nx   = 360l if nx==None else nx
    ny   = 180l if ny==None else ny
    mask = np.array(map!=0,dtype=float) if mask==None else mask


# -------- utilities
    dx   = 360./np.float(nx-1)
    dy   = 180./np.float(ny-1)
    x    = np.arange(180.,-180.-dx,-dx)
    y    = np.arange(-90.,90.+dy,dy)


# -------- get the galactic coords and convert to aitoff pixels
    nside  = np.long(np.sqrt(map.size/12l))
    l, b   = healgen_lb(nside)
    xa, ya = aitoff(l,b)
    imap = np.array(np.round(-(xa-180.)/dx), dtype=int)
    jmap = np.array(np.round((ya+ 90.)/dy), dtype=int)


# -------- initialize the map
    amap = np.zeros((nx,ny))
    cmap = np.zeros((nx,ny))


# -------- loop through pixels (GGD: there must be a faster way...)
    amap[imap,jmap] += map*mask
    cmap[imap,jmap] += 1.0


# -------- set outside to white (or black)
    bg = np.where(cmap < 1.0)

    if blackbg!=None:
        amap[bg] = -1.0e30
    else:
        amap[bg] = 1.0e30


# -------- return average map
    return np.rot90(amap/(cmap + np.array((cmap==0),dtype=float)))
