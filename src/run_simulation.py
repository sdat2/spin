#!/usr/bin/env python3
"""1_Schedule.py - a program by sdat2 to show that the
particle class works to effectively simulate a series of test masses
Usage: python3 1_Schedule.py """

### Local Libraries ###
import src.time_wrapper as twr
import src.simulation.helpers as hlp

# galaxy_collide(x00=-15, y00=20, m0=1, recalculate = False)
# galaxy_collide(x00=-15, y00 =20, m0=0.2, MAX_TIMER=450, EP=0.01, recalculate = True)
# three_body_collide(MAX_TIMER=1000, EP=0.01, recalculate = True)
# many_body_collide()
# radii = [2, 3, 4, 5, 6]
# num_particles = [12, 18, 24, 30, 36]


@twr.timeit
def run_through_galaxy_sizes(halo=False):
    radii = [1.0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5]
    num_particles = [15, 25, 35, 45, 50, 58, 70, 80, 90, 100, 110, 130]
    hlp.galaxy_collide(
        x00=-22,
        y00=22,
        m0=1,
        MAX_TIMER=10,
        TSTEP=0.05,
        radii=radii,
        num_particles=num_particles,
        EP=0.01,
        VB=50,
        recalculate=False,
        FILE_DIR="./Snap_Check/",
        halo=halo,
        algorithm="herm",
        animate=True,
    )


# galaxy_collide(FILE_DIR='./Other_Tester/', MAX_TIMER=80, recalculate=True)
# galaxy_collide(FILE_DIR='./Ani_Re_Test/', MAX_TIMER=10, recalculate=True)


@twr.timeit
def run_through_impactors():
    for i in range(20):
        hlp.galaxy_collide(
            x00=2 * (-5 - i),
            y00=2 * (5 + i),
            m0=1.0,
            MAX_TIMER=450,
            radii=[],
            num_particles=[],
            EP=0.01,
            recalculate=True,
            FILE_DIR="./IMPACT_CHECK/",
        )


# three_body_collide(MAX_TIMER=1200, VB=100, TSTEP=0.001, EP=0.01, recalculate = True, algorithm = 'rk4o')

# many_body_collide(MAX_TIMER=600, VB=1000, TSTEP=0.0001, EP=0.01, recalculate=True, algorithm='vv')
