import numpy as np
from ezdolly.params import VehicleParams
from ezdolly.double_ackermann import wheel_angles, wheel_angles_kinematic, wheel_speeds, yaw_rate_ref
from ezdolly.yaw_controller import gains_simc, m_max_report
from ezdolly.allocation import MotorModels, allocate, build_B, moment_limits

p = VehicleParams()

def test_straight_line():
    assert np.allclose(wheel_angles(0.0, p), 0) and np.allclose(wheel_speeds(2.0, 0.0, p), 2.0)

def test_ackermann_geometry():
    d = np.deg2rad(20); rho = p.L / (2 * np.tan(d)); dl = wheel_angles(d, p)
    assert np.isclose(np.tan(dl[0]), p.L / 2 / (rho - p.w / 2))          # FL : roue interieure braque plus
    assert dl[0] > dl[1] > 0 and np.allclose(dl[2:], -dl[:2])           # arriere = -avant
    assert np.isclose(yaw_rate_ref(2.0, d, p), 2 * 2.0 * np.tan(d) / p.L)

def test_report_numbers():
    Kp, Ki = gains_simc(4.2 * np.sqrt(4900 * 11900), 0.2, 0.05)
    assert abs(Kp - 1.28e5) / 1.28e5 < 0.01 and abs(Ki - 1.28e5) / 1.28e5 < 0.01
    assert abs(m_max_report(0.0, p) - 16.3e3) < 100

def test_allocation_matches_numeric_and_constraint():
    from scipy.optimize import minimize
    mo = MotorModels(eff_sheets=("Efficiency Orig", "Efficiency randomly modified", "Efficiency frozen modification", "Efficiency Orig"),
                     carac_sheets=("Carac Orig", "Carac randomly modified", "Carac frozen modification", "Carac Orig"))
    from ezdolly.double_ackermann import motor_speeds_rpm
    d = np.deg2rad(15); om = motor_speeds_rpm(10 / 3.6, d, p); V = np.array([1440.0, 300.0])
    u, V2, _ = allocate(V, d, om, p, mo)
    B = build_B(wheel_angles(d, p), p); assert np.allclose(B @ u, V2, atol=1e-6)
    umax = mo.umax(om)
    ref = minimize(lambda z: -mo.cosphi(om, z).sum(), u * 0.9, bounds=list(zip(-umax, umax)),
                   constraints={"type": "eq", "fun": lambda z: B @ z - V2}, method="SLSQP", options={"ftol": 1e-13})
    assert mo.cosphi(om, u).sum() >= -ref.fun - 1e-6


def test_course_formula_equals_kinematics():
    for dd in np.deg2rad([-35, -10, -1, 0.0, 1, 10, 35]):
        assert np.allclose(wheel_angles(dd, p), wheel_angles_kinematic(dd, p), atol=1e-12)
    d = np.deg2rad(20); c = 1 / np.tan(d)          # formule du PDF avec cot
    assert np.isclose(wheel_angles(d, p)[1], np.arctan(1 / (c + p.w / p.L)))
    assert np.isclose(wheel_angles(d, p)[0], np.arctan(1 / (c - p.w / p.L)))
