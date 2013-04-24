import numpy as np
from pix2ang_ring import *

def healgen(nside, nest=None):

# Is nside too big?
    if nside > 8192:
        print '---> HEALGEN ERROR: NSIDE too large.'
        print '     Routine not implemented with type long64 yet'

    npix = 12l*long(nside)**2
    ipix = np.arange(npix)

    if nest != None:
        print '---> HEALGEN ERROR: Only ring ordered is written yet'


    return pix2ang_ring(nside,ipix)
