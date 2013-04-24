import numpy as np

def aitoff(l,b):

# -------- utilities
    radeg = 57.2957795131 # deg/rad


# -------- convert
    sa   = np.array(l)
    x180 = np.where(sa > 180.0)[0]

    if x180.size>0:
        sa[x180] -= 360.

    alpha2 = sa/(2.0*radeg)
    delta  = b/radeg
    r2     = np.sqrt(2.0)
    f      = 2.0*r2/np.pi
    cdec   = np.cos(delta)
    denom  = np.sqrt(1.0 + cdec*np.cos(alpha2))
    x      = cdec*np.sin(alpha2)*2.0*r2/denom
    y      = np.sin(delta)*r2/denom
    x      = x*radeg/f
    y      = y*radeg/f

    return x,y
