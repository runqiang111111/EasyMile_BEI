"""Parametres du vehicule (fiche vehicule du BEI) et conventions.

Conventions : x vers l'avant, y vers la gauche, psi_dot > 0 a gauche.
Ordre des roues partout : [FL, FR, RL, RR] = [avant-gauche, avant-droite, arriere-gauche, arriere-droite].
Notation : r = rapport de reduction (26), R_w = rayon de roue (0.24 m), rho = rayon de virage.
(Le rapport du yaw controller note N le rapport de reduction et r le rayon de roue.)
"""
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

DATA_XLSX = Path(__file__).resolve().parents[2] / "data" / "Donnees_BEI_2025.xlsx"
WHEELS = ("FL", "FR", "RL", "RR")


@dataclass
class VehicleParams:
    L: float = 4.8            # empattement [m]
    w: float = 2.2            # voie [m]
    R_w: float = 0.24         # rayon de roue [m]
    r: float = 26.0           # rapport de reduction moteur/essieu
    m_empty: float = 4900.0   # [kg]
    m_full: float = 11900.0   # [kg] (charge max 7 t)
    tau_nom: float = 34.3     # couple nominal moteur [Nm]
    eta: float = 1.0          # rendement transmission
    v_max: float = 15 / 3.6   # [m/s]
    theta: float = 0.05       # retard equivalent de la chaine de couple [s] (suppose)
    c_rr: float = 0.01        # coef. de resistance au roulement (hypothese)
    g: float = 9.81

    @property
    def xs(self):             # abscisses des roues (FL, FR, RL, RR)
        return np.array([self.L / 2, self.L / 2, -self.L / 2, -self.L / 2])

    @property
    def ys(self):             # ordonnees des roues
        return np.array([self.w / 2, -self.w / 2, self.w / 2, -self.w / 2])

    def Iz(self, m):          # inertie de lacet : boite uniforme 6.4 x 3.0 m  => Iz ~ 4.2 m
        return m * (6.4**2 + 3.0**2) / 12.0

    @property
    def m_design(self):       # masse de conception = moyenne geometrique
        return float(np.sqrt(self.m_empty * self.m_full))
