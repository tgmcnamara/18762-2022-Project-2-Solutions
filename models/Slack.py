from __future__ import division
import numpy as np
from models.Buses import Buses
from models.Buses import Buses
from scripts.stamp_helpers import *


class Slack:

    def __init__(self,
                 Bus,
                 Vset,
                 ang,
                 Pinit,
                 Qinit):
        """Initialize slack bus in the power grid.

        Args:
            Bus (int): the bus number corresponding to the slack bus.
            Vset (float): the voltage setpoint that the slack bus must remain fixed at.
            ang (float): the slack bus voltage angle that it remains fixed at.
            Pinit (float): the initial active power that the slack bus is supplying
            Qinit (float): the initial reactive power that the slack bus is supplying
        """
        # You will need to implement the remainder of the __init__ function yourself.
        self.Bus = Bus
        self.Vset = Vset
        self.ang = ang
        self.Pinit = Pinit
        self.Qinit = Qinit
        # initialize
        self.Vr_set = Vset*np.cos(ang*np.pi/180)
        self.Vi_set = Vset*np.sin(ang*np.pi/180)

    def assign_nodes(self, bus):
        """Assign the additional slack bus nodes for a slack bus.
        Returns:
            None
        """
        self.Vr_node = bus[Buses.bus_key_[self.Bus]].node_Vr
        self.Vi_node = bus[Buses.bus_key_[self.Bus]].node_Vi
        self.P_node = Buses._node_index.__next__()
        self.Q_node = Buses._node_index.__next__()

    def stamp(self, V, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
        Vr = V[self.Vr_node]
        Vi = V[self.Vi_node]
        P = -V[self.P_node]
        Q = V[self.Q_node]

        # (Vr + jVi)*(Ir - jIi) = P + jQ
        # (Ir - jIi) = (P + jQ)/(Vr + jVi)
        # (Ir - jIi) = (P + jQ)(Vr - jVi)/(Vr**2 + Vi**2)
        # (Ir - jIi) = (PVr - jPVi + jQVr + QVi)/(Vr**2 + Vi**2)
        # (Ir - jIi) = ((PVr + QVi) - j(PVi - QVr))/(Vr**2 + Vi**2)
        # Ir = (PVr + QVi) / (Vr**2 + Vi**2)
        # Ii = (PVi - QVr) / (Vr**2 + Vi**2)

        Irg_hist = (P*Vr - Q*Vi)/(Vr**2 + Vi**2)
        dIrgdP = Vr/(Vr**2 + Vi**2)
        dIrgdQ = -Vi/(Vr**2 + Vi**2)
        Irg_Jstamp = Irg_hist - dIrgdP*P - dIrgdQ*Q
        idx_Y = stampY(self.Vr_node, self.P_node, dIrgdP, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(self.Vr_node, self.Q_node, dIrgdQ, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.Vr_node, Irg_Jstamp, J_val, J_row, idx_J)

        Iig_hist = (P*Vi - Q*Vr)/(Vr**2 + Vi**2)
        dIigdP = Vi/(Vr**2 + Vi**2)
        dIigdQ = -Vr/(Vr**2 + Vi**2)
        Iig_Jstamp = Iig_hist - dIigdP*P - dIigdQ*Q
        idx_Y = stampY(self.Vi_node, self.P_node, dIigdP, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(self.Vi_node, self.Q_node, dIigdQ, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.Vi_node, Iig_Jstamp, J_val, J_row, idx_J)

        # enforce slack constraints
        idx_Y = stampY(self.P_node, self.Vr_node, 1, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.P_node, self.Vr_set, J_val, J_row, idx_J)

        idx_Y = stampY(self.Q_node, self.Vi_node, 1, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.Q_node, self.Vi_set, J_val, J_row, idx_J)

        return (idx_Y, idx_J)

    def calc_residuals(self, resid, V):
        resid[self.P_node] = V[self.Vr_node] - self.Vr_set
        resid[self.Q_node] = V[self.Vi_node] - self.Vi_set
        