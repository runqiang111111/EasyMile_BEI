"""Boucle fermee : consigne (v, delta_G) -> double Ackermann -> yaw PI -> allocation -> vehicule."""
import numpy as np
from .params import VehicleParams
from .double_ackermann import yaw_rate_ref, motor_speeds_rpm
from .yaw_controller import YawPI, gains_simc
from .vehicle import YawPlant, longitudinal_torque
from .allocation import MotorModels, allocate, build_B
from .double_ackermann import wheel_angles


def scenario_step_steer(v_kmh=10.0, delta_deg=15.0, t_steer=3.0, T=12.0, t_acc=2.0, dt=0.01):
    t = np.arange(0, T, dt)
    v = np.clip(t / t_acc, 0, 1) * v_kmh / 3.6
    delta = np.where(t >= t_steer, np.deg2rad(delta_deg), 0.0)
    return t, v, delta


def simulate(t, v, delta_G, mass=None, method="cosphi", motors=None, p=None, regen=True):
    p = p or VehicleParams()
    m = mass or p.m_design
    motors = motors or MotorModels()
    dt = t[1] - t[0]
    Kp, Ki = gains_simc(p.Iz(p.m_design), Tc=0.2, theta=p.theta)   # gains fixes (masse inconnue)
    pi, plant = YawPI(Kp, Ki, dt), YawPlant(p, m, dt)
    a = np.gradient(v, dt)
    out = {k: [] for k in ("psid_ref", "psid", "M_cmd", "M_real", "tauG", "tau", "cosphi", "omega", "Mlim")}
    for k in range(len(t)):
        psid_ref = yaw_rate_ref(v[k], delta_G[k], p)
        tauG = longitudinal_torque(v[k], a[k], m, p)
        omega = motor_speeds_rpm(v[k], delta_G[k], p)
        # limites de M_G a ce point de fonctionnement -> anti-windup du PI
        B = build_B(wheel_angles(delta_G[k], p), p)
        umax = motors.umax(omega)
        umin = -umax if regen else np.zeros(4)
        lim = None
        tg = float(np.clip(tauG, B[0] @ umin, B[0] @ umax))
        from .allocation import moment_limits
        lim = moment_limits(B, tg, umin, umax) or (0.0, 0.0)
        M_cmd = pi.step(psid_ref, plant.psid, lim[0], lim[1])
        tau, V, _ = allocate(np.array([tauG, M_cmd]), delta_G[k], omega, p, motors, method, regen)
        M_real = float((B @ tau)[1])
        psid = plant.step(M_real)
        for key, val in zip(out, (psid_ref, psid, M_cmd, M_real, tauG, tau, motors.cosphi(omega, tau), omega, lim)):
            out[key].append(val)
    return {k: np.array(v_) for k, v_ in out.items()} | {"t": t, "v": v, "delta": delta_G}


def cosphi_report(res, t_from=0.0):
    """cos phi de chaque roue (moyenne temporelle pour t >= t_from) + somme. Renvoie un dict {nom_roue: valeur, 'somme': valeur}."""
    from .params import WHEELS
    sel = res["t"] >= t_from
    per_wheel = res["cosphi"][sel].mean(axis=0)
    out = {n: float(c) for n, c in zip(WHEELS, per_wheel)}
    out["somme"] = float(per_wheel.sum())
    return out


def format_cosphi(report):
    return "  ".join(f"{k}={v:.4f}" for k, v in report.items())
