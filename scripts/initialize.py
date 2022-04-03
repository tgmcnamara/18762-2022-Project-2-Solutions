import numpy as np


def initialize(size_Y, bus, generator, slack, flat_start=False):
    V_init = np.zeros(size_Y, dtype=np.float)
    if flat_start:
        for ele in bus:
            V_init[ele.node_Vr] = 1
            V_init[ele.node_Vi] = 0
        for ele in generator:
            V_init[ele.Q_node] += (ele.Qmax+ele.Qmin)/2
        # initialize slack currents as 0?
    else:
        for ele in bus:
            V_init[ele.node_Vr] = ele.Vr_init
            V_init[ele.node_Vi] = ele.Vi_init
        for ele in generator:
            V_init[ele.Q_node] += -ele.Qinit
        for ele in slack:
            V_init[ele.Slack_Ir_node] = ele.Ir_init
            V_init[ele.Slack_Ii_node] = ele.Ii_init

    return V_init