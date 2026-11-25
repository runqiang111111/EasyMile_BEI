"""Exemple : echelon de braquage, comparaison allocation cos phi vs repartition egale (pseudo-inverse)."""
import pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ezdolly.simulation import simulate, scenario_step_steer, cosphi_report, format_cosphi
from ezdolly.allocation import MotorModels

# Feuilles Excel par roue (FL, FR, RL, RR) : a adapter selon l'enonce du prof.
motors = MotorModels(
    eff_sheets=("Efficiency Orig", "Efficiency randomly modified", "Efficiency frozen modification", "Efficiency Orig"),
    carac_sheets=("Carac Orig", "Carac randomly modified", "Carac frozen modification", "Carac Orig"))

fig, axs = plt.subplots(3, 2, figsize=(11, 10)); res = {}
for kmh, delta in ((10, 15), (15, 8)):
    t, v, d = scenario_step_steer(kmh, delta)
    for method in ("cosphi", "pinv"):
        res[(kmh, method)] = simulate(t, v, d, method=method, motors=motors)
    r = res[(kmh, "cosphi")]; q = res[(kmh, "pinv")]
    sel = t > 4.0
    print(f"\n{kmh} km/h, delta_G={delta} deg  (moyennes pour t > 4 s)")
    print("  cos phi par roue, allocation cos phi :", format_cosphi(cosphi_report(r, 4.0)))
    print("  cos phi par roue, pseudo-inverse     :", format_cosphi(cosphi_report(q, 4.0)))
    print(f"  erreur de lacet RMS = {np.sqrt(np.mean((r['psid_ref']-r['psid'])[sel]**2)):.4f} rad/s")
r = res[(10, "cosphi")]
axs[0, 0].plot(t, np.rad2deg(r["psid_ref"]), "k--", label="consigne"); axs[0, 0].plot(t, np.rad2deg(r["psid"]), label="mesure")
axs[0, 0].set(title="Vitesse de lacet (10 km/h, delta_G=15 deg)", ylabel="deg/s"); axs[0, 0].legend()
axs[0, 1].plot(t, r["M_cmd"] / 1e3, label="M_G commande"); axs[0, 1].plot(t, r["Mlim"][:, 1] / 1e3, "r:", label="limite faisable")
axs[0, 1].set(title="Moment de lacet", ylabel="kN.m"); axs[0, 1].legend()
for j, n in enumerate(("FL", "FR", "RL", "RR")): axs[1, 0].plot(t, r["tau"][:, j], label=n)
axs[1, 0].set(title="Couples moteur alloues (cos phi)", ylabel="Nm", xlabel="t [s]"); axs[1, 0].legend()
axs[1, 1].plot(t, r["omega"]); axs[1, 1].set(title="Vitesse moteur de chaque roue", ylabel="tr/min", xlabel="t [s]")
for j, n in enumerate(("FL", "FR", "RL", "RR")): axs[2, 0].plot(t, r["cosphi"][:, j], label=n)
axs[2, 0].set(title="cos phi de chaque roue (allocation cos phi)", xlabel="t [s]"); axs[2, 0].legend()
axs[2, 1].plot(t, r["cosphi"].sum(1), label="allocation cos phi"); axs[2, 1].plot(t, res[(10, "pinv")]["cosphi"].sum(1), "--", label="pseudo-inverse")
axs[2, 1].set(title="Somme des cos phi", xlabel="t [s]"); axs[2, 1].legend()
for a in axs.ravel(): a.grid(alpha=.3)
print()
fig.tight_layout(); fig.savefig(pathlib.Path(__file__).parent / "simulation.png", dpi=130)
