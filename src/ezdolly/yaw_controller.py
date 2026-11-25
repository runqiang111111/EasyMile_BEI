"""Bloc 2 : controleur de lacet (PI sur la vitesse de lacet) -> moment de lacet M_G*."""
import numpy as np
from .params import VehicleParams


def gains_simc(Iz_c, Tc=0.2, theta=0.05):
    """Regle SIMC pour un integrateur 1/(Iz s) avec retard theta : Kp = Iz/(Tc+theta), Ki = Iz/(4 (Tc+theta)^2)."""
    return Iz_c / (Tc + theta), Iz_c / (4.0 * (Tc + theta) ** 2)


def gains_pole_placement(Iz, wn, zeta):
    """Kp = 2 zeta wn Iz, Ki = wn^2 Iz."""
    return 2 * zeta * wn * Iz, wn**2 * Iz


def m_max_report(tauG, p: VehicleParams):
    """Moment max du rapport : M_max = w/(2R_w) (4 r eta tau_max - |tau_G|)  [Nm], tau_G = couple d'avance total (Nm roue)."""
    return max(0.0, p.w / (2 * p.R_w) * (4 * p.r * p.eta * p.tau_nom - abs(tauG)))


class YawPI:
    """PI avec saturation et anti-windup (integration conditionnelle)."""

    def __init__(self, Kp, Ki, dt):
        self.Kp, self.Ki, self.dt, self.I = Kp, Ki, dt, 0.0

    def reset(self):
        self.I = 0.0

    def step(self, psid_ref, psid_meas, M_min, M_max):
        e = psid_ref - psid_meas
        u = self.Kp * e + self.Ki * self.I
        M = float(np.clip(u, M_min, M_max))
        if M == u or (u - M) * e < 0:      # on n'integre pas si on sature dans le sens de l'erreur
            self.I += e * self.dt
        return M
