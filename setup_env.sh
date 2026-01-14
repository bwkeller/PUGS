#!/usr/bin/env bash

# Initialize our virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install pinned packages
pip install -r pip-requirements.txt
pip install -r build-requirements.txt
