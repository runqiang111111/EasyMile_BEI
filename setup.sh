#!/usr/bin/env bash
# Linux / macOS : cree l'environnement virtuel .venv et installe les dependances
set -e
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install -e .
echo "OK. Activer l'environnement : source .venv/bin/activate"
