#!/usr/bin/env python3
"""
sclass.py by sdat2
I became bored of passing to many different input files around,
and decided I would like some nice
class objects to deal with this for me.
"""
import numpy as np


class Controls:
    """
    This class contains all the globals that I used to use to set the program parameters.
    Instance shortened to co within main simulation.
    """

    def __init__(
        self,
        MAXTIMER=20,
        TSTEP=0.05,
        vb=5,
        halo=False,
        OUT="NEW",
        EPS=0.01,
        name="NoName",
        calculate_am=False,
        algorithm="rk4o",
    ):
        self.tstep = TSTEP  # The timestep
        self.maxt = MAXTIMER  # the maximum time
        self.vb = vb  # the verbosity
        self.out = OUT  # Output place for this run of the code
        self.halo = halo  # whether we're bothering to include halos
        self.radius = 7  # The halo radius of a standard dark galaxy.
        self.GM = 1  # An arbitrary constant.
        self.epsilon = (
            EPS  # Length softening distance (now defunct, retained for backwards comp.)
        )
        self.name = name  # name for file outputs
        self.calculate_am = calculate_am  # whether to wastefully calculate the AM of each forceful particle
        self.algorithm = (
            algorithm  # the algorithm used for processing 'herm', 'vv' or 'rk4o'
        )


class System:
    """
    This class contains all the globals that I used to output
    Class to hide the particles, coordinates and timer inside.
    Instance shortened to sy in main simulation.
    """

    def __init__(self, co, particles=[]):
        self.timer = np.linspace(0, co.maxt, num=int(co.maxt / co.tstep))
        shortened_length = 0
        for t in range(len(self.timer)):
            if t % co.vb == 0:
                shortened_length += 1
        self.coordinate_grid = np.zeros((3, len(particles), shortened_length))
        self.short_timer = []
        self.Total_Energies = dict()
        self.Total_Energies["Kinetic"] = np.zeros(shortened_length)
        self.Total_Energies["Gravitational"] = np.zeros(shortened_length)
        self.Total_Energies["Total"] = np.zeros(shortened_length)
        self.Angular_Momentum = dict()
        self.Angular_Momentum["Total"] = np.zeros((3, shortened_length))
        if co.calculate_am:
            for part in particles:  # This takes up lots of memory, disable if required.
                if part.forceful:  # only for non test masses
                    self.Angular_Momentum[part.no] = np.zeros((3, shortened_length))
