import numpy as np

def process_results(v, bus, slack, generator):
    print("BUS VOLTAGES:")
    for ele in bus:
        Vr = v[ele.node_Vr]
        Vi = v[ele.node_Vi]
        Vmag = np.sqrt(Vr**2 + Vi**2)
        Vth = np.arctan2(Vi, Vr) * 180/np.pi
        print("Id: %d, Vmag: %.3f p.u., Vth: %.3f deg" % (ele.Bus, Vmag, Vth))