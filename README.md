# EasyMile_BEI - EZDolly : repartition vectorielle du couple

Chaine : consigne (v, delta_G) -> **double_ackermann** -> **yaw_controller** (PI) -> **allocation** (cos phi) -> **vehicle** (lacet).

## Installation (environnement virtuel)

Linux / macOS :

    ./setup.sh && source .venv/bin/activate

Windows (PowerShell) :

    .\setup.ps1

A la main : `python -m venv .venv`, activer, puis `pip install -r requirements.txt` et `pip install -e .`

## Utilisation

    python -m pytest -q                      # 5 tests
    python examples/static_allocation.py     # un point de fonctionnement : couples et cos phi de chaque roue
    python examples/run_simulation.py        # echelon de braquage en boucle fermee (~30 s) -> examples/simulation.png

Les sorties donnent le cos phi **de chaque roue** (FL, FR, RL, RR) et leur somme, pour l'allocation cos phi et pour la
repartition egale (pseudo-inverse). En code : `motors.cosphi(omega, tau)` (par roue), `res["cosphi"]` (serie temporelle
n x 4) et `cosphi_report(res, t_from)` (moyenne par roue + somme).

## Structure du projet

    EasyMile_BEI/
    |-- src/ezdolly/     le code (double_ackermann, yaw_controller, allocation, vehicle, simulation, params)
    |-- data/            donnees Excel (Donnees_BEI_2025.xlsx)
    |-- examples/        scripts a lancer
    |-- tests/           tests pytest
    |-- pyproject.toml   configuration du projet
    `-- requirements.txt dependances

## Conventions
Ordre des roues : FL, FR, RL, RR. r = 26 (reduction), R_w = 0.24 m (rayon de roue), rho = rayon de virage.
Les feuilles Excel par roue se choisissent dans `MotorModels(eff_sheets=..., carac_sheets=...)` (donnees dans `data/`).
