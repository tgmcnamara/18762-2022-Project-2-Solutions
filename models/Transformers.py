from __future__ import division
from itertools import count
from scripts.stamp_helpers import *
from models.Buses import Buses
import numpy as np

class Transformers:
    _ids = count(0)

    def __init__(self,
                 from_bus,
                 to_bus,
                 r,
                 x,
                 status,
                 tr,
                 ang,
                 Gsh_raw,
                 Bsh_raw,
                 rating):
        """Initialize a transformer instance

        Args:
            from_bus (int): the primary or sending end bus of the transformer.
            to_bus (int): the secondary or receiving end bus of the transformer
            r (float): the line resitance of the transformer in
            x (float): the line reactance of the transformer
            status (int): indicates if the transformer is active or not
            tr (float): transformer turns ratio
            ang (float): the phase shift angle of the transformer
            Gsh_raw (float): the shunt conductance of the transformer
            Bsh_raw (float): the shunt admittance of the transformer
            rating (float): the rating in MVA of the transformer
        """
        self.from_bus = from_bus
        self.to_bus = to_bus
        self.r = r
        self.x = x
        self.status = status
        self.tr = tr
        self.ang = ang
        self.Gsh_raw = Gsh_raw
        self.Bsh_raw = Bsh_raw
        self.rating = rating

        # Set minimum x:
        if abs(self.x) < 1e-6:
            if self.x < 0:
                self.x = -1e-6
            else:
                self.x = 1e-6

        # convert to G and B
        self.G_pu = self.r/(self.r**2+self.x**2)
        self.B_pu= -self.x/(self.r**2+self.x**2)

        self.id = self._ids.__next__()

    def assign_indexes(self, bus):
        self.Vr_from_node = bus[Buses.bus_key_[self.from_bus]].node_Vr
        self.Vi_from_node = bus[Buses.bus_key_[self.from_bus]].node_Vi
        self.Vr_to_node = bus[Buses.bus_key_[self.to_bus]].node_Vr
        self.Vi_to_node = bus[Buses.bus_key_[self.to_bus]].node_Vi

    def stamp(self, V, Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J):
        if not self.status:
            return (idx_Y, idx_J)

        tr = self.tr
        tr2 = tr*tr
        Gt = self.G_pu
        Bt = self.B_pu + self.Bsh_raw/2
        phi_rad = self.ang*np.pi/180
        cosphi = np.cos(phi_rad)
        sinphi = np.sin(phi_rad)
        Gcosphi = self.G_pu*cosphi
        Gsinphi = self.G_pu*sinphi
        Bcosphi = self.B_pu*cosphi
        Bsinphi = self.B_pu*sinphi
        G_shunt_from = self.G_pu/tr2
        B_shunt_from = Bt/tr2
        MR_from = (Gcosphi  - Bsinphi)/tr
        MI_from = (Gsinphi  + Bcosphi)/tr
        G_to = (Gcosphi + Bsinphi)/tr
        B_to = (Bcosphi - Gsinphi)/tr
        MR_to = Gt
        MI_to = Bt
        
        dIrfdVrf = G_shunt_from
        dIrfdVrt = -MR_from
        dIrfdVif = -B_shunt_from
        dIrfdVit = MI_from
        idx_Y = stampY(self.Vr_from_node, self.Vr_from_node, dIrfdVrf, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_from_node, self.Vr_to_node, dIrfdVrt, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_from_node, self.Vi_from_node, dIrfdVif, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_from_node, self.Vi_to_node, dIrfdVit, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        
        dIrtdVrf = -G_to
        dIrtdVrt = MR_to
        dIrtdVif = B_to
        dIrtdVit = -MI_to
        idx_Y = stampY(self.Vr_to_node, self.Vr_from_node, dIrtdVrf, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_to_node, self.Vr_to_node, dIrtdVrt, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_to_node, self.Vi_from_node, dIrtdVif, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vr_to_node, self.Vi_to_node, dIrtdVit, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        
        dIifdVrf = B_shunt_from
        dIifdVrt = -MI_from
        dIifdVif = G_shunt_from
        dIifdVit = -MR_from
        idx_Y = stampY(self.Vi_from_node, self.Vr_from_node, dIifdVrf, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_from_node, self.Vr_to_node, dIifdVrt, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_from_node, self.Vi_from_node, dIifdVif, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_from_node, self.Vi_to_node, dIifdVit, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        
        dIitdVrf = -B_to
        dIitdVrt = MI_to
        dIitdVif = -G_to
        dIitdVit = MR_to
        idx_Y = stampY(self.Vi_to_node, self.Vr_from_node, dIitdVrf, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_to_node, self.Vr_to_node, dIitdVrt, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_to_node, self.Vi_from_node, dIitdVif, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        idx_Y = stampY(self.Vi_to_node, self.Vi_to_node, dIitdVit, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        return (idx_Y, idx_J)

    def calc_residuals(self, resid, V):
        Vrf = V[self.Vr_from_node]
        Vrt = V[self.Vr_to_node]
        Vif = V[self.Vi_from_node]
        Vit = V[self.Vi_to_node]

        tr = self.tr
        tr2 = tr*tr
        Gt = self.G_pu
        Bt = self.B_pu + self.Bsh_raw/2
        phi_rad = self.ang*np.pi/180
        cosphi = np.cos(phi_rad)
        sinphi = np.sin(phi_rad)
        Gcosphi = self.G_pu*cosphi
        Gsinphi = self.G_pu*sinphi
        Bcosphi = self.B_pu*cosphi
        Bsinphi = self.B_pu*sinphi
        G_shunt_from = self.G_pu/tr2
        B_shunt_from = Bt/tr2
        MR_from = (Gcosphi  - Bsinphi)/tr
        MI_from = (Gsinphi  + Bcosphi)/tr
        G_to = (Gcosphi + Bsinphi)/tr
        B_to = (Bcosphi - Gsinphi)/tr
        MR_to = Gt
        MI_to = Bt
        
        dIrfdVrf = G_shunt_from
        dIrfdVrt = -MR_from
        dIrfdVif = -B_shunt_from
        dIrfdVit = MI_from
        resid[self.Vr_from_node] = dIrfdVrf*Vrf + dIrfdVrt*Vrt + dIrfdVif*Vif + dIrfdVit*Vit

        dIrtdVrf = -G_to
        dIrtdVrt = MR_to
        dIrtdVif = B_to
        dIrtdVit = -MI_to
        resid[self.Vr_to_node] = dIrtdVrf*Vrf + dIrtdVrt*Vrt + dIrtdVif*Vif + dIrtdVit*Vit

        dIifdVrf = B_shunt_from
        dIifdVrt = -MI_from
        dIifdVif = G_shunt_from
        dIifdVit = -MR_from
        resid[self.Vi_from_node] = dIifdVrf*Vrf + dIifdVrt*Vrt + dIifdVif*Vif + dIifdVit*Vit

        dIitdVrf = -B_to
        dIitdVrt = MI_to
        dIitdVif = -G_to
        dIitdVit = MR_to
        resid[self.Vi_to_node] = dIitdVrf*Vrf + dIitdVrt*Vrt + dIitdVif*Vif + dIitdVit*Vit