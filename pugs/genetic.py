import numpy as np
import tangos


def write_particle_ids(halo_id, filename="id_file.txt"):
    """
    Write out the particle ID file that GenetIC needs to generate the zoom ICs.

    :param halo_id: The TANGOS halo ID
    :param filename: The output filename for the particle IDs
    """
    ids = tangos.get_halo(halo_id).calculate("ids()")
    np.savetxt(filename, ids, fmt="%d")


def build_param_file(filename="genetIC_zoom.txt"):
    with open("../inputs/genetIC_volume.txt", "r") as t:
        template = t.read()
        with open(filename, "w") as f:
            f.write(template)
            f.write("id_file id_file.txt\n")
