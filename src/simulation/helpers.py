#!/usr/bin/env python3
"""
helpers.py by sdat2

This program generically contains lots of programs to deal with
the gravitational particles at a more meta level.

"""
import os  # needed for creating output directory
from halos import *
from sclass import *
from animator import *
import pickle


def circular_position(radius, number_particles, member):
    """Returns the initial position of a particle."""
    angle = 2 * np.pi / (number_particles) * member
    position = np.array([np.cos(angle), np.sin(angle), 0]) * radius
    return position


def circular_velocity(radius, number_particles, member, co):
    """Returns the velocity needed for circular motion."""
    angle = 2 * np.pi / (number_particles) * member
    co.GM = 1
    co.radius = 7
    if co.halo:
        speed = np.sqrt(co.GM * (radius ** 2) / ((co.radius) ** 3))
    else:
        speed = np.sqrt(co.GM / radius)
    velocity = np.array([-np.sin(angle), np.cos(angle), 0]) * speed
    return velocity


@timeit
def galaxy_creator(co, clockwise, radii=[], num_particles=[], **kwargs):
    """creates the particles spinning"""
    particles = []  # Empty list to fill with particles
    countA = 0  # CountA is indice of the particle
    particles.append(
        Particle(
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            1.0,
            True,
            0,
            halo=co.halo,
        )
    )
    if len(radii) != 0:
        for i in range(len(radii)):
            countB = 0  # CountB is number in a Particular Ring
            for j in range(num_particles[i]):
                countA += 1
                countB += 1
                pos = circular_position(radii[i], num_particles[i], j)
                vel = circular_velocity(radii[i], num_particles[i], j, co)
                if clockwise:
                    vel = -vel
                    # reverse direction of cloud
                particles.append(Particle(pos, vel, 1.0, False, countA))
    for i in range(len(particles)):
        particles[i].gravp_energy(particles)
    return particles


def impactor_adder(particles, x0=-25.0, y0=45.0, mass=0.22):
    """
    Computes the correct parabolic orbit given a staring position and
    mass of the other galaxy. Check the appendices
    """
    x_vec = [x0, y0, 0.0]
    x_dot_vec = [0, -np.sqrt(2 * (1 + mass) / (np.sqrt((x0) ** 2 + (y0) ** 2))), 0.0]
    particles.append(
        Particle(np.array(x_vec), np.array(x_dot_vec), mass, True, len(particles))
    )
    return particles


def three_particle_adder(co, particles):
    """creates the particles for the the Three Body Problem Test"""
    particles.append(
        Particle(
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            1.0,
            True,
            len(particles),
            halo=co.halo,
        )
    )
    particles.append(
        Particle(
            np.array([0.0, 5, 0]),
            np.array([0.5, -0.2, 0.0]),
            1.0,
            True,
            len(particles),
            halo=co.halo,
        )
    )
    particles.append(
        Particle(
            np.array([0.0, -5, 0]),
            np.array([-0.5, 0.2, 0.0]),
            1.0,
            True,
            len(particles),
            halo=co.halo,
        )
    )
    for i in range(len(particles)):
        particles[i].gravp_energy(particles)
    return particles


def arbitrary_particle_adder(co, particles, pos_list=[], vel_list=[], mass_list=[]):
    """ Adds some more particles."""
    for i in range(len(pos_list)):
        particles.append(
            Particle(
                np.array(pos_list[i]),
                np.array(vel_list[i]),
                mass_list[i],
                True,
                len(particles),
                halo=co.halo,
            )
        )
    return particles


@timeit
def write_out(co, sy, **kwargs):
    """Writes output so that it can be animated"""
    if not os.path.exists(co.out):  # make the directory thing
        os.makedirs(co.out)
    with open(co.out + "/" + co.name + "_co.pickle", "wb") as handle:
        pickle.dump(co, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # I can't remember any more why the protocol has to be the highest
    with open(co.out + "/" + co.name + "_sy.pickle", "wb") as handle:
        pickle.dump(sy, handle, protocol=pickle.HIGHEST_PROTOCOL)


def middle_manager(
    co, sy, name, particles, recalculate=True, animate=False, scipy_test=False
):
    """
    Tells the Spinner to run, and then the animator to run.
    """
    if recalculate == True:
        # propogate galaxy and fill sy
        if scipy_test == False:
            co, sy, particles = spinner(co, sy, particles)
        else:
            sy, particles = scipy_spinner(co, sy, particles)
        write_out(co, sy)  # print sy and co to pickle files
        print("There are " + str(len(particles)) + " particles in the simulation")
    for mv in [True]:  # turning off annoying animations
        if mv:
            print("Now animating " + name + " in comoving frame.")
        else:
            print("Now animating " + name + " in inertial frame.")
        aniclass_member = AniMP4(co.out, name=name, move_with=mv)
        # initiate animator
        if animate == True:
            aniclass_member.animate_starter()  # run animation / plots


@timeit
def galaxy_collide(
    MAX_TIMER=600,
    VB=10,
    TSTEP=0.05,
    EP=0.01,
    x00=-15.0,
    y00=30,
    m0=0.2,
    radii=[2, 3, 4, 5],
    num_particles=[12, 18, 25, 30],
    recalculate=True,
    FILE_DIR="./EP_TEST/",
    halo=False,
    algorithm="vv",
    animate=False,
    **kwargs
):
    """
    Creates the initial galaxies and collision courses, and then
    hands over to the middle manager.
    """
    particles = []  # create a default empty list to pass about
    sub = (
        str(MAX_TIMER)
        + "t_"
        + str(x00)
        + "x0_"
        + str(y00)
        + "y0_"
        + str(m0)
        + "m_"
        + str(EP)
        + "e_"
        + "Halos"
        + str(halo)
    )
    co = Controls(
        MAXTIMER=MAX_TIMER,
        TSTEP=TSTEP,
        vb=VB,
        halo=halo,
        EPS=EP,
        OUT=FILE_DIR + sub,
        calculate_am=False,
        algorithm=algorithm,
    )
    for clock in [False, True]:
        if clock:
            name = "Retrograde" + sub
        else:
            name = "Prograde" + sub
        sy = System(co, particles=particles)
        co.name = name
        if recalculate == True:
            particles = galaxy_creator(
                co, clock, radii=radii, num_particles=num_particles
            )
            particles = impactor_adder(particles, x0=x00, y0=y00, mass=m0)
            particles = zmf_transform(particles)
            sy = System(co, particles=particles)
        middle_manager(
            co, sy, name, particles=particles, recalculate=recalculate, animate=animate
        )


def two_body_circle(
    MAX_TIMER=600, VB=10, TSTEP=0.05, EP=0.01, recalculate=True, algorithm="vv"
):
    """Makes two bodies of equal mass circle round each other, for a long time
    to check any errors that might be extant in the functions"""
    co = Controls(
        MAXTIMER=MAX_TIMER,
        TSTEP=TSTEP,
        vb=VB,
        halo=False,
        EPS=EP,
        OUT="./Algo_Comp/" + sub,
        calculate_am=True,
        algorithm=algorithm,
    )
    particles = []
    p = [[1, 0, 0], [-1, 0, 0]]
    v = [[0, 1, 0], [0, -1, 0]]
    m = [1, 1]
    arbitrary_particle_adder(co, particles, pos_list=p, vel_list=v, mass_list=m)


@timeit
def three_body_collide(
    MAX_TIMER=600,
    VB=10,
    TSTEP=0.05,
    EP=0.01,
    recalculate=True,
    algorithm="vv",
    **kwargs
):
    """
    Runs a three body collision to look check other functions.
    """
    particles = []
    sub = str(MAX_TIMER) + "t_" + str(EP) + "e_" + str(TSTEP) + "dt_" + str(algorithm)
    co = Controls(
        MAXTIMER=MAX_TIMER,
        TSTEP=TSTEP,
        vb=VB,
        halo=False,
        EPS=EP,
        OUT="./TBP_TEST/" + sub,
        calculate_am=True,
        algorithm=algorithm,
    )
    sy = System(co, particles=particles)
    name = "TBP" + sub
    co.name = name
    if recalculate == True:
        particles = three_particle_adder(co, particles)
        sy = System(co, particles=particles)
    middle_manager(co, sy, name, particles=particles, recalculate=recalculate)


def many_body_collide(
    MAX_TIMER=600, VB=1000, TSTEP=0.0001, EP=0.01, recalculate=True, algorithm="vv"
):
    """
    Runs a multi body collision to look check other functions.
    """
    particles = []
    sub = str(MAX_TIMER) + "t_" + str(EP) + "e"
    co = Controls(
        MAXTIMER=MAX_TIMER,
        TSTEP=TSTEP,
        vb=VB,
        halo=False,
        EPS=EP,
        OUT="./MBP_TEST/" + sub,
        calculate_am=True,
        algorithm=algorithm,
    )
    sy = System(co, particles=particles)
    name = "MBP" + sub
    co.name = name
    p = [
        [0.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, -5.0, 0.0],
        [5.0, 0.0, 0.0],
        [-5.0, 0.0, 0.0],
    ]
    v = [
        [0.0, 0.0, 0.0],
        [0.5, -0.2, 0.0],
        [-0.5, 0.2, 0.0],
        [-0.2, 0.5, 0.0],
        [-0.2, 0.5, 0.0],
    ]
    m = [1, 1, 1, 1, 1]
    if recalculate == True:
        particles = arbitrary_particle_adder(
            co, particles, pos_list=p, vel_list=v, mass_list=m
        )
        sy = System(co, particles=particles)
    middle_manager(co, sy, name, particles=particles, recalculate=recalculate)
