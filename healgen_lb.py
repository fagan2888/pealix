import numpy as np
from healgen import *

def healgen_lb(nside, nest=None):

# -------- utilities
    radeg = 57.2957795131 # deg/rad


# -------- error exception
    if nest!= None:
        print '---> HEALGEN_LB ERROR: nest ordered not written yet.'

    theta, phi = healgen(nside)

    l = phi*radeg
    b = 90. - theta*radeg

    return l, b

