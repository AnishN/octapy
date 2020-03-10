import numpy as np

def get_particles_from_release_file(str release_path):
    cdef ParticleC[::1] particles = np.array(np.genfromtxt(
        release_path, 
        delimiter=",", 
        #names=True,
        skip_header=1,
        dtype=[
            ("id", "u8"),
            ("lat", "f4"),
            ("lon", "f4"),
            ("depth", "f4"),
            ("start_time", "f8"),
            ("curr_time", "f8"),
            ("end_time", "f8"),
        ]), ndmin=1)
    return particles