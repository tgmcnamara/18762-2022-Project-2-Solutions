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
        self.Slack_Ir_node = Buses._node_index.__next__()
        self.Slack_Ii_node = Buses._node_index.__next__()

    def stamp(self, V, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
        # slack currents leaving their nodes
        idx_Y = stampY(self.Vr_node, self.Slack_Ir_node, 1, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(self.Vi_node, self.Slack_Ii_node, 1, Y_val, Y_row, Y_col, idx_Y)

        # enforce slack constraints
        idx_Y = stampY(self.Slack_Ir_node, self.Vr_node, 1, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.Slack_Ir_node, self.Vr_set, J_val, J_row, idx_J)

        idx_Y = stampY(self.Slack_Ii_node, self.Vi_node, 1, Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(self.Slack_Ii_node, self.Vi_set, J_val, J_row, idx_J)

        return (idx_Y, idx_J)

    def calc_slack_PQ(self, V_sol):
        Ir = V_sol[self.Slack_Ir_node]
        Ii = V_sol[self.Slack_Ii_node]
        Vr = self.Vr_set
        Vi = self.Vi_set
        S = (Vr + 1j*Vi)*(Ir - 1j*Ii)
        P = np.real(S)
        Q = np.imag(S)
        return (P, Q)

    def calc_residuals(self, resid, V):
        resid[self.Vr_node] += V[self.Slack_Ir_node]
        resid[self.Vi_node] += V[self.Slack_Ii_node]
        resid[self.Slack_Ir_node] += V[self.Vr_node] - self.Vr_set
        resid[self.Slack_Ii_node] += V[self.Vi_node] - self.Vi_set