"""Bloc 3 : allocation de couple maximisant le facteur de puissance (cos phi) de chaque roue.

f(omega, tau) = a1 tau^2 + a2 omega^2 + a3 tau omega + a4 tau + a5 omega + a6   (tau [Nm moteur], omega [rpm])
max sum_ij f(omega_ij, tau_ij)  s.c.  B tau = V,  umin <= tau <= umax
Solution : tau = tau* + W^-1 B^T (B W^-1 B^T)^-1 (V - B tau*)   (Johansen & Fossen 2013, eq. 12-14)
           tau*_ij = -(a3 omega_ij + a4)/(2 a1),  W = diag(-2 a1)
"""
import warnings
import numpy as np
import openpyxl
from scipy.optimize import linprog, minimize
from .params import VehicleParams, DATA_XLSX


# ---------------------------------------------------------------- donnees et ajustement
def load_efficiency(sheet="Efficiency Orig", path=DATA_XLSX):
    """Points (omega [rpm], T [Nm], cos phi). La ligne factice (0,0,0) du bubble chart est exclue."""
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    rows = [r[:3] for r in ws.iter_rows(min_row=3, values_only=True) if r[0] and r[2] and r[4] == 1.0]
    return np.array(rows, float)


def fit_alpha(data):
    """Moindres carres -> (a1..a6). Variables normalisees (omega/1000, T/50) pour le conditionnement."""
    w, T, c = data.T
    x, y = w / 1000.0, T / 50.0
    Phi = np.column_stack([y**2, x**2, x * y, y, x, np.ones_like(x)])
    k = np.linalg.lstsq(Phi, c, rcond=None)[0]
    return np.array([k[0] / 50**2, k[1] / 1e6, k[2] / (50 * 1e3), k[3] / 50, k[4] / 1e3, k[5]])


def cosphi(a, omega, tau):
    a1, a2, a3, a4, a5, a6 = a
    return a1 * tau**2 + a2 * omega**2 + a3 * tau * omega + a4 * tau + a5 * omega + a6


def load_torque_limit(sheet="Carac Orig", col=6, path=DATA_XLSX):
    """Renvoie omega -> Tmax(omega). col : 2 = 5 min, 4 = 30 min, 6 = 60 min."""
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    car = np.array([r[:7] for r in ws.iter_rows(min_row=2, values_only=True)
                    if all(isinstance(c, (int, float)) for c in r[:7])], float)
    return lambda w: np.interp(w, car[:, 0], car[:, col])


class MotorModels:
    """Un modele cos phi + une limite de couple par roue (feuilles Excel configurables)."""

    def __init__(self, eff_sheets=("Efficiency Orig",) * 4, carac_sheets=("Carac Orig",) * 4, col=6, path=DATA_XLSX):
        self.alpha = [fit_alpha(load_efficiency(s, path)) for s in eff_sheets]
        self.tmax = [load_torque_limit(s, col, path) for s in carac_sheets]

    def umax(self, omega):
        return np.array([f(w) for f, w in zip(self.tmax, omega)])

    def cosphi(self, omega, tau):
        return np.array([cosphi(a, w, t) for a, w, t in zip(self.alpha, omega, tau)])

    def curve_coeffs(self, omega):
        """Parabole de chaque roue : cos phi_ij(tau) = a1 tau^2 + b tau + c ; renvoie (a1, b, c, tau_star)."""
        a1 = np.array([a[0] for a in self.alpha])
        b = np.array([a[2] * w + a[3] for a, w in zip(self.alpha, omega)])
        c = np.array([a[1] * w**2 + a[4] * w + a[5] for a, w in zip(self.alpha, omega)])
        return a1, b, c, -b / (2 * a1)


# ---------------------------------------------------------------- matrice B
def build_B(delta_ij, p: VehicleParams):
    """V = [tau_G, M_G]^T = B tau (tau_ij = couple moteur).
    Ligne 1 : r cos(delta_ij)                      (couple d'avance total, Nm roue)
    Ligne 2 : (r/R_w) (x_i sin(delta_ij) - y_j cos(delta_ij))   (moment de lacet, terme de voie inclus)"""
    c, s = np.cos(delta_ij), np.sin(delta_ij)
    return np.vstack([p.r * c, p.r / p.R_w * (p.xs * s - p.ys * c)])


def moment_limits(B, tauG, umin, umax):
    """Plage de M_G atteignable a tau_G fixe (LP). Renvoie (M_min, M_max) ou None si tau_G infaisable."""
    kw = dict(A_eq=B[0:1], b_eq=[tauG], bounds=list(zip(umin, umax)), method="highs")
    lo = linprog(B[1], **kw)
    hi = linprog(-B[1], **kw)
    if lo.status != 0 or hi.status != 0:
        return None
    return lo.fun, -hi.fun


# ---------------------------------------------------------------- allocation explicite
def _fallback(B, V, a1, b, umin, umax):
    u0 = np.clip(np.linalg.lstsq(B, V, rcond=None)[0], umin, umax)
    r = minimize(lambda z: np.sum(-a1 * z**2 - b * z), u0, jac=lambda z: -2 * a1 * z - b,
                 bounds=list(zip(umin, umax)), constraints={"type": "eq", "fun": lambda z: B @ z - V},
                 method="SLSQP", options={"ftol": 1e-12, "maxiter": 500})
    return r.x


def allocate_cosphi(B, V, a1, b, umin, umax):
    """Solution explicite + ensemble actif pour les saturations. a1 < 0 : coefficient de tau^2 (par roue)."""
    n = len(b)
    ustar = -b / (2 * a1)                  # couple de cos phi max de chaque roue
    Winv = 1.0 / (-2 * a1)
    fixed = {}
    for _ in range(n + 1):
        free = [i for i in range(n) if i not in fixed]
        u = np.zeros(n)
        for i, val in fixed.items():
            u[i] = val
        if len(free) < B.shape[0]:
            u = _fallback(B, V, a1, b, umin, umax)
            break
        Vf = V - sum(B[:, i] * val for i, val in fixed.items()) if fixed else V
        Bf, Wf = B[:, free], Winv[free]
        try:
            u[free] = ustar[free] + Wf * (Bf.T @ np.linalg.solve(Bf @ (Wf[:, None] * Bf.T), Vf - Bf @ ustar[free]))
        except np.linalg.LinAlgError:
            u = _fallback(B, V, a1, b, umin, umax)
            break
        viol = [(max(umin[i] - u[i], u[i] - umax[i]), i) for i in free
                if u[i] < umin[i] - 1e-9 or u[i] > umax[i] + 1e-9]
        if not viol:
            break
        i = max(viol)[1]
        fixed[i] = float(np.clip(u[i], umin[i], umax[i]))
    else:
        u = _fallback(B, V, a1, b, umin, umax)
    if np.abs(B @ u - V).max() > 1e-6 * max(1.0, np.abs(V).max()):
        warnings.warn("Consigne (tau_G, M_G) infaisable : B u != V")
    return u


def allocate_pinv(B, V, umin, umax):
    """Reference : repartition 'egale' (pseudo-inverse a norme minimale), puis ecretage."""
    return np.clip(np.linalg.lstsq(B, V, rcond=None)[0], umin, umax)


def allocate(V, delta_G, omega, p: VehicleParams, motors: MotorModels, method="cosphi", regen=True):
    """Point d'entree du bloc : (tau_G, M_G), delta_G, omega_ij -> tau_ij (Nm moteur) et M_G realisable.
    V est d'abord ramene dans le domaine faisable (c'est ce M_G ecrete que le PI doit voir pour l'anti-windup)."""
    from .double_ackermann import wheel_angles
    B = build_B(wheel_angles(delta_G, p), p)
    umax = motors.umax(omega)
    umin = -umax if regen else np.zeros(4)
    tauG = float(np.clip(V[0], B[0] @ umin, B[0] @ umax))
    lim = moment_limits(B, tauG, umin, umax)
    Mlo, Mhi = lim if lim else (0.0, 0.0)
    V = np.array([tauG, float(np.clip(V[1], Mlo, Mhi))])
    if method == "cosphi":
        a1, b, _, _ = motors.curve_coeffs(omega)
        u = allocate_cosphi(B, V, a1, b, umin, umax)
    else:
        u = allocate_pinv(B, V, umin, umax)
    return u, V, (Mlo, Mhi)
