#!/usr/bin/env bash

# Build our local venv
python -m venv .venv
source .venv/bin/activate
pip install -r ../pip-requirements.txt
pip install -r requirements.txt
