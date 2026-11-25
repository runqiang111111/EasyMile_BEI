"""Bloc 1 : double Ackermann (consigne de braquage de l'essieu equivalent -> 4 angles de roue)."""
import numpy as np
from .params import VehicleParams


def curvature(delta_G, p: VehicleParams):
    """kappa = 1/rho = 2 tan(delta_G) / L  (double Ackermann : tan(delta_G) = (L/2)/rho)."""
    return 2.0 * np.tan(delta_G) / p.L


def yaw_rate_ref(v, delta_G, p: VehicleParams):
    """psi_dot* = v/rho = 2 v tan(delta_G) / L   [rad/s]   (eq. 1 du rapport yaw)."""
    return v * curvature(delta_G, p)


def wheel_angles(delta_G, p: VehicleParams):
    """Angles de braquage (FL, FR, RL, RR) -- formule du cours (PDF Double Ackermann), l = voie, L = empattement :
        delta_fr = -delta_rr = arctan( 1 / (cot(delta_G) + l/L) )
        delta_fl = -delta_rl = arctan( 1 / (cot(delta_G) - l/L) )
    Ecrite avec tan(delta_G) (1/(cot a +- b) = tan a / (1 +- b tan a)) pour rester valide en ligne droite."""
    t = np.tan(delta_G)
    q = p.w / p.L
    fl = np.arctan(t / (1.0 - q * t))
    fr = np.arctan(t / (1.0 + q * t))
    return np.array([fl, fr, -fl, -fr])


def wheel_angles_kinematic(delta_G, p: VehicleParams):
    """Meme resultat par la cinematique : tan(delta_ij) = x_i / (rho - y_j), rho = L / (2 tan delta_G)."""
    k = curvature(delta_G, p)
    return np.arctan2(p.xs * k, 1.0 - p.ys * k)


def wheel_speeds(v, delta_G, p: VehicleParams):
    """Vitesses lineaires v_ij des roues [m/s] ; v = vitesse du centre G."""
    k = curvature(delta_G, p)
    return v * np.hypot(p.xs * k, 1.0 - p.ys * k)


def motor_speeds_rpm(v, delta_G, p: VehicleParams):
    """Vitesse moteur de chaque roue [tr/min] : omega = (60/2pi) (r/R_w) v_ij."""
    return wheel_speeds(v, delta_G, p) * p.r / p.R_w * 60.0 / (2 * np.pi)
