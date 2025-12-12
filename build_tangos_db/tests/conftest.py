import os
import subprocess
import tarfile

import pynbody as pyn
import pytest
import requests
import tangos


@pytest.fixture(scope="session")
def pyn_snaps():
    sims = {}
    halos = {}
    for tsim in tangos.all_simulations():
        snap = tsim.timesteps[0]
        sim = pyn.load(snap.filename)
        sim.physical_units()
        h = sim.halos()
        sims[snap.filename] = sim
        halos[snap.filename] = h
    return sims, halos


@pytest.fixture(scope="session")
def get_testdata():
    url = "https://nbody.shop/NUGS128.tar.gz"
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open("NUGS128.tar.gz", "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        # Extract the tar.gz file
        with tarfile.open("NUGS128.tar.gz", "r:gz") as tar:
            tar.extractall(os.environ["TANGOS_SIMULATION_FOLDER"], filter="data")
    else:
        raise Exception(f"Failed to download test data: {r.status_code}")


@pytest.fixture(scope="session", autouse=True)
def build_database(get_testdata):
    print("Building DB")
    yield subprocess.run(["/bin/bash", "build_tangos.sh"])
    os.remove("pugs.db")
