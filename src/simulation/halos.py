#!/usr/bin/env python3
"""
halos.py, a module to house the general class of gravitational particles
which can have non collisional extent. It also contains some functions that are
explictly designed ot iterate over the particles.
TODO: Add a test against scipy.integrate
"""

import numpy as np

# numpy is useful as it allows the program to have high speed and precision.
from numpy import linalg as LA

# linear algebra allows all the vectors to be dealt with appropriately.
import time
from plotters import *
from sclass import *
from scipy.integrate import solve_ivp  # not currently used


def _prec(*args, **kwargs):
    """
    This make sure that we use 128 bit precission throughout by re-typing any
    thing given to it
    :param args:
    :param kwargs:
    :return: all args in 128 bit format
    """
    return np.array(*args, dtype=np.float128, **kwargs)


class Particle:
    """A class for gravitationally interacting particles.
    If the attribute 'forceful' is true then there is equal
    active and passive mass. Otherwise there is only passive mass."""

    def __init__(
        self, position, velocity, mass, forceful, number, halo=False, radius=7
    ):
        self.no = number
        self.x = _prec(position)
        self.v = _prec(velocity)
        self.m = _prec(mass)  # I'm not clear whether this really needs to be dtyp3e 128
        self.forceful = forceful
        self.ang()
        self.kin_energy()
        if halo == True:
            self.radius = _prec(radius)  # co.radius
            #  if they are dark they have this radius
        else:
            self.radius = _prec(0.01)  # made _prec only to avoid the bugs.
            # I got rid of length softening in favour of a minimal halo distance

    def __repr__(self):
        attributes = (
            " at position "
            + str(self.x)
            + " with velocity "
            + str(self.v)
            + " and a dark matter radius "
            + str()
        )
        if self.forceful:
            return "An active particle" + attributes
        else:
            return "A test mass" + attributes

    def pos(self, position):  # These small functions are quite unneccesary
        self.x = position

    def vel(self, velocity):
        self.v = velocity

    def ang(self):
        """Perform simple cross product"""
        self.am = np.array(
            [
                self.v[2] * self.x[1] - self.x[2] * self.v[1],
                self.x[2] * self.v[0] - self.v[2] * self.x[0],
                self.v[1] * self.x[0] - self.x[1] * self.v[0],
            ]
        )

    def kin_energy(self):  # The LA dictionary is helpful here
        self.ke = 0.5 * self.m * (LA.norm(self.v)) ** 2

    def gravp_energy(self, particles):
        """Sum up gravpot contributions"""
        GM = 1
        gpe_temp = 0
        for part in particles:
            if part.forceful and not part.no == self.no:
                gpe_temp += -GM * self.m * part.m / LA.norm(part.x - self.x)
        self.gpe = gpe_temp

    def sim(self, particles, co):
        """Choose which simulation algorithm to use."""
        if co.algorithm == "rk4o":
            self.rk4o(particles, co)
        if co.algorithm == "vv":
            self.velocity_verlet(particles, co)
        if co.algorithm == "herm":
            self.hermite(particles, co)

    def velocity_verlet(self, particles, co):
        """Implementation of second order velocity verlot algorithm."""
        x_old = self.x
        v_old = self.v
        acc_old = acc_fun(self.no, x_old, particles, co)
        dt = co.tstep
        x_new = x_old + v_old * dt + (acc_old * dt ** 2) / 2
        acc_new = acc_fun(self.no, x_new, particles, co)
        v_new = v_old + (acc_new + acc_old) * dt / 2
        self.pos(x_new)
        self.vel(v_new)

    def hermite(self, particles, co):
        """Implementation of a fourth order hermite integrator, credit to
        http://www.artcompsci.org/kali/vol/two_body_problem_2/.hbody.rb.html#hermite"""
        x_old = self.x
        v_old = self.v
        old_acc = acc_fun(self.no, x_old, particles, co)
        dt = co.tstep
        old_jerk = self.jerk(particles, co)
        self.x += (
            self.v * dt + old_acc * (dt * dt / 2.0) + old_jerk * (dt * dt * dt / 6.0)
        )
        self.v += old_acc * dt + old_jerk * (dt * dt / 2.0)
        new_jerk = self.jerk(particles, co)
        new_acc = acc_fun(self.no, x_old, particles, co)
        self.v = (
            v_old
            + (old_acc + new_acc) * (dt / 2.0)
            + (old_jerk - new_jerk) * (dt * dt / 12.0)
        )
        self.x = (
            x_old
            + (v_old + self.v) * (dt / 2.0)
            + (old_acc - new_acc) * (dt * dt / 12.0)
        )
        # The Frobenius norm is given by [1]: A_F = [\sum_{i,j} abs(a_{i,j})^2]^{1/2}

    def jerk(self, particles, co):
        """A function required by the hermite integrator
        to correct using the rate of change of acceleration."""
        jerk = np.array([0.0, 0.0, 0.0])
        for part in particles:
            if part.forceful and not part.no == self.no:
                # avoid counting self, or allowing interaction between test masses
                # This allows the 0th particle to have a dark matter halo
                if LA.norm(part.x - self.x) > part.radius:
                    jerk += (
                        co.GM
                        * part.m
                        * (
                            ((part.v - self.v) / LA.norm(part.x - self.x) ** 3)
                            + 3
                            * (part.x - self.x)
                            * np.dot((part.x - self.x), (part.v - self.v))
                            / (LA.norm(part.x - self.x) ** 5)
                        )
                    )
                    # standard r^(-2) type attraction if outside the other's halo
        return jerk

    def rk4o(self, particles, co):  # I hate all things
        """Implementation of fourth order Runge Kutta."""
        x_1 = self.x
        v_1 = self.v
        dx1, dv1 = v_1, acc_fun(self.no, x_1, particles, co)
        x_2 = x_1 + dx1 * co.tstep / 2
        v_2 = v_1 + dv1 * co.tstep / 2
        dx2, dv2 = v_2, acc_fun(self.no, x_2, particles, co)
        x_3 = x_1 + dx2 * co.tstep / 2
        v_3 = v_1 + dv2 * co.tstep / 2
        dx3, dv3 = v_3, acc_fun(self.no, x_3, particles, co)
        x_4 = x_1 + dx3 * co.tstep
        v_4 = v_1 + dv3 * co.tstep
        dx4, dv4 = v_4, acc_fun(self.no, x_4, particles, co)
        self.pos((x_1 + (dx1 + 2 * dx2 + 2 * dx3 + dx4) * (co.tstep / 6)))
        self.vel(v_1 + ((dv1 + 2 * dv2 + 2 * dv3 + dv4) * (co.tstep / 6)))


def acc_fun(num, position, particles, co):
    """Gravitational acc_fun function to work out what the forces are."""
    acc = _prec(np.array([0.0, 0.0, 0.0]))  # return a np array
    for part in particles:
        # cycle through all particles in plist
        if part.forceful and not part.no == num:
            # avoid counting self, or allowing interaction between test masses
            # This allows the 0th particle to have a dark matter halo
            if LA.norm(part.x - position) > part.radius:
                acc += (
                    co.GM
                    * part.m
                    * (part.x - position)
                    / ((LA.norm(part.x - position)) ** (3))
                )
                # standard r^(-2) type attraction if outside the other's halo

            else:
                acc += co.GM * part.m * (part.x - position) / ((part.radius) ** (3))
                # Force that would be produced by a dark matter halo of constant density
    return acc
    # See Appendix A of report for the derivations.


def scipy_equation_to_solve(particles):
    for part in particles:
        part.v += acc_fun(part.num, part.x, particles, co)
        part.x += part.v
    return particles


def scipy_method_for_spinning(co, sy, particles):
    """
    :param co: the command variables that specify the simulation parameters
    :param sy: the system storage variables that will be written as output
    :param particles: the list of particle objects
    :return: those same objects with any modifications
    """
    soln = solve_ivp(scipy_equation_to_solve, (0, co.maxt), particles, method="RK45")
    for i in soln.particles:
        for j in soln[i]:
            sy.coordinate_grid[i] = soln[i][j].x
            sy.t
    return co, sy, particles


# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html


def zmf_transform(particles):
    """Transforms the velocities of all the particles to the ZMF"""
    forceful_count = 0.0
    mass = 0.0  # temporary objects to be filled
    momentum = np.array([0.0, 0.0, 0.0])
    w_position = np.array([0.0, 0.0, 0.0])
    for part in particles:
        if part.forceful:
            forceful_count += 1
            momentum += part.m * part.v
            w_position += part.m * part.x
            mass += part.m
    zmf_vel = momentum / mass
    com_pos = w_position / mass
    for part in particles:
        part.vel(part.v - zmf_vel)
        part.pos(part.x - com_pos)
    return particles


# --------------Add up things---------------------------


def tot_energies(particles):
    """Adds up all the energy in the particles list"""
    kinetic = 0.0
    potential = 0.0
    for part in particles:
        if part.forceful:
            kinetic += part.ke
            potential += part.gpe  # This will double count
    return kinetic, potential / 2, kinetic + potential / 2


def tot_am(particles):
    """adds up all the angular momentum"""
    ang_tot = np.array([0.0, 0.0, 0.0])
    for part in particles:
        if part.forceful:
            ang_tot += part.am
    return ang_tot


# -------------------------Run Simulation Steps---------


@timeit
def spinner(co, sy, particles, **kwargs):
    """Runs simulation steps whilst making a """
    print(
        "Simulating "
        + str(len(particles))
        + " particles for "
        + str(len(sy.timer))
        + " time steps"
    )
    for f in range(len(sy.timer)):  ### Time index is 'f' inside this loop
        if f % co.vb == 0:  # once in a while let's record the data
            for i in range(0, len(particles)):
                particles[i].sim(particles, co)  ###### CHOOSE ALGO NOW ######
                sy.coordinate_grid[:, i, int(f / co.vb)] = particles[i].x
                particles[i].ang()  # update angular momentum
                particles[i].kin_energy()  # update kinetic energy
                particles[i].gravp_energy(particles)
                if co.calculate_am == True:
                    sy.Angular_Momentum[particles[i].no][:, int(f / co.vb)] = particles[
                        i
                    ].am
                    # Give all particles an angular Momentum History (V. Wasteful!!!)
            (
                sy.Total_Energies["Kinetic"][int(f / co.vb)],
                sy.Total_Energies["Gravitational"][int(f / co.vb)],
                sy.Total_Energies["Total"][int(f / co.vb)],
            ) = tot_energies(particles)
            sy.short_timer.append(sy.timer[f])
            sy.Angular_Momentum["Total"][:, int(f / co.vb)] = tot_am(particles)
        else:
            for i in range(0, len(particles)):
                particles[i].sim(particles, co)  ###### CHOOSE ALGO ######
        if f % 1000 == 0:
            print("There have been " + str(f) + " simulation steps.")
        if f % 5000 == 0:
            print("Take snap shot now.")
            snap_shot(co, particles, move_with=True, Time=sy.timer[f])

    sy.short_timer = np.array(sy.short_timer)  # recast into numpy array
    return co, sy, particles


@timeit
def spin_forward(time, co, particles=[], **kwargs):
    """A function that continues to
    simulate the galaxy without graphing,
    so we can check long term behaviour more frugally"""
    for i in range(0, int(time / co.tstep)):
        for j in range(0, len(particles)):
            particles[j].sim(particles, co)
            ###### CHOOSE ALGO ######
    return particles
