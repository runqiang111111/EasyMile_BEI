# Windows (PowerShell) : cree l'environnement virtuel .venv et installe les dependances
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install -e .
Write-Host "OK. Activer l'environnement : .\.venv\Scripts\Activate.ps1"
