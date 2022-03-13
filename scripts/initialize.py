import numpy as np


def initialize(size_Y, bus, generator, flat_start=False):
    V_init = np.zeros(size_Y, dtype=np.float)
    if flat_start:
        for ele in bus:
            V_init[ele.node_Vr] = 1
            V_init[ele.node_Vi] = 0
    else:
        for ele in bus:
            V_init[ele.node_Vr] = ele.Vr_init
            V_init[ele.node_Vi] = ele.Vi_init
    return V_init