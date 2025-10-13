#!/usr/bin/env bash

# Set up the dev environment
./setup_env.sh
source .venv/bin/activate

black --check .
flake8 .
isort --check-only .
