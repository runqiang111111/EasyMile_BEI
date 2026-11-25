"""Un point de fonctionnement : couples alloues et cos phi de CHAQUE roue (cos phi vs repartition egale)."""
import numpy as np
from ezdolly.params import VehicleParams, WHEELS
from ezdolly.double_ackermann import motor_speeds_rpm
from ezdolly.allocation import MotorModels, allocate

p = VehicleParams()
motors = MotorModels(
    eff_sheets=("Efficiency Orig", "Efficiency randomly modified", "Efficiency frozen modification", "Efficiency Orig"),
    carac_sheets=("Carac Orig", "Carac randomly modified", "Carac frozen modification", "Carac Orig"))

v, delta_G, V = 10 / 3.6, np.deg2rad(15), np.array([1440.0, 300.0])   # V = [tau_G (Nm), M_G (Nm)]
omega = motor_speeds_rpm(v, delta_G, p)
print(f"v = {v*3.6:.0f} km/h, delta_G = {np.rad2deg(delta_G):.0f} deg, tau_G = {V[0]:.0f} Nm, M_G = {V[1]:.0f} Nm")
print("omega_ij [rpm] :", {n: int(round(o)) for n, o in zip(WHEELS, omega)})
print(f"{'methode':10s} " + " ".join(f"{n:>9s}" for n in WHEELS) + "      somme")
for method in ("cosphi", "pinv"):
    tau, V_ok, _ = allocate(V, delta_G, omega, p, motors, method)
    c = motors.cosphi(omega, tau)
    print(f"{method:10s} tau [Nm] : " + " ".join(f"{t:9.2f}" for t in tau))
    print(f"{'':10s} cos phi  : " + " ".join(f"{x:9.4f}" for x in c) + f"   {c.sum():.4f}")
