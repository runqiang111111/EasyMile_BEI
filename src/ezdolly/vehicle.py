"""Modele vehicule : lacet = integrateur Iz*psi_ddot = Mz (rapport yaw), avec retard de la chaine de couple ;
longitudinal : couple d'avance tau_G deduit de la vitesse imposee (m a + roulement)."""
from collections import deque
import numpy as np
from .params import VehicleParams


class YawPlant:
    def __init__(self, p: VehicleParams, m, dt):
        self.p, self.m, self.dt = p, m, dt
        self.Iz = p.Iz(m)
        self.psid = 0.0
        self.psi = 0.0
        n = max(1, int(round(p.theta / dt)))
        self.buf = deque([0.0] * n, maxlen=n)      # retard theta sur le moment delivre

    def step(self, M_delivered):
        M = self.buf[0]
        self.buf.append(M_delivered)
        self.psid += M / self.Iz * self.dt
        self.psi += self.psid * self.dt
        return self.psid


def longitudinal_torque(v, a, m, p: VehicleParams):
    """tau_G = R_w * (m a + c_rr m g)   [Nm, somme des couples roue] (aero neglige)."""
    return p.R_w * (m * a + p.c_rr * m * p.g * (1.0 if v > 0.05 else 0.0))
