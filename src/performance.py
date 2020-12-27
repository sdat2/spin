#!/usr/bin/env python3
""" performance.py A program to weigh up the virtues of the
different algorithms
By testing them for energy and angular momentum conservation"""

import os
import copy
import numpy as np
from src.simulation.helpers import *
from scipy.optimize import curve_fit
import src.data_loading.io_pickles as ip
import src.simulation.sclass
import src.plotting.my_plotting_style as mps
import matplotlib.pyplot as plt
import pickle
from src.simulation.sclass import *
import src.time_wrapper as twr


@twr.timeit
def _circular_orbit(seperation=2, algo="vv", tstep=0.05):
    """
    Take the desired seperation between the two equally sized masses,
    and turn that into a list of two particles at the right speed
    and distance.
    """
    co = Controls(
        MAXTIMER=20,
        TSTEP=tstep,
        algorithm=algo,
        OUT="MUT_CIRC",
        name=algo + "_" + str(tstep),
    )
    speed = np.sqrt(co.GM / 2 / seperation)
    particles = []
    particles.append(
        Particle(
            np.array([-seperation / 2, 0, 0]),
            np.array([0, speed, 0]),
            1.0,
            False,
            countA,
        )
    )
    particles.append(
        Particle(
            np.array([-seperation / 2, 0, 0]),
            np.array([0, speed, 0]),
            1.0,
            False,
            countA,
        )
    )
    tmp_log_data = {}
    co, sy, particles = spinner(co, sy, particles, log_time=tmp_log_data)
    time_taken = tmp_log_data["SPINNER"]
    write_out(co, sy)
    return time, energy, time_taken


@twr.timeit
def circle_orbit_tester():
    """
    This function calls _circular_orbit a number of times
    for different settings, and then sends the output to a
    _graph_circular_orbit.
    """
    energy_remaining = {}  # ar
    simu_time = {}
    time_taken = {}
    for algo in ["rk4o", "vv", "herm"]:
        time_taken[algo] = []
        simu_time[algo] = {}
        energy_remaining[algo] = {}
        tsteps = [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1]
        for tstep in tsteps:
            simu_time[algo][tstep], energy_remaining[algo][tstep], tt = _circular_orbit(
                algo=algo, tstpep=tstep
            )
            time_taken.append(tt)
    _graph_circular_orbit(energy_remaining, simu_time, time_taken)


@twr.timeit
def three_body_test():
    """ Three body test with a variety of algorithms"""
    print("Let us try three body")
    tbc_time = []
    for algo in ["rk4o", "vv", "herm"]:
        tbc_time[algo] = []
        for dt in [0.001, 0.002, 0.01, 0.1, 0.2, 0.5]:
            vb = int(100 * (0.001 / dt))
            if vb == 0:
                # Prevent div0 errors
                vb = 1
            tmp_log_data = {}
            three_body_collide(
                MAX_TIMER=5000,
                VB=vb,
                TSTEP=dt,
                EP=0.01,
                recalculate=True,
                algorithm=algo,
                log_time=tmp_log_data,
            )
            tbc_time[algo].append(tmp_log_data["THREE_BODY_COLLIDE"])


@twr.timeit
def particle_build_up(OM=5):
    """A function which gradually fills up a galaxy with a variable number
    particles, trying each algorithm used in the simulation, allowing
    computational excess to be compared"""
    # Fills up the shells at these radii to these maximums
    radii = [
        0,
        1.0,
        1.5,
        2,
        2.5,
        3,
        3.5,
        4,
        4.5,
        5,
        5.5,
        6,
        6.5,
        7.0,
        7.5,
        8.0,
        8.5,
        9.0,
        9.5,
    ]
    num_particles = [
        1,
        12,
        18,
        24,
        31,
        36,
        42,
        50,
        58,
        70,
        80,
        100,
        130,
        160,
        200,
        250,
        300,
        340,
        400,
    ]

    total_particles = 0
    for i in num_particles:
        total_particles += i

    print("You can go up to " + str(total_particles) + " if you want to.")

    N_vec = np.logspace(1, OM, num=OM, base=2)

    # variables to store the raw data
    tot_particles_d = {}
    fill_up_time = []
    ani_time = []
    spin_round_time = {}
    spinner_time = {}
    spin_round_time["vv"] = []
    spin_round_time["rk4o"] = []
    spin_round_time["herm"] = []
    spinner_time["vv"] = []
    spinner_time["rk4o"] = []
    spinner_time["herm"] = []
    tot_particles_d["vv"] = []
    tot_particles_d["FU"] = []
    tot_particles_d["herm"] = []
    tot_particles_d["rk4o"] = []
    tot_particles_d["ANI"] = []
    tot_particles_d["vv_spin"] = []
    tot_particles_d["rk4o_spin"] = []
    tot_particles_d["herm_spin"] = []

    for required_total in N_vec:
        tmp_log_data = {}
        # particles = []
        particles = _fill_up(
            radii, num_particles, int(required_total), log_time=tmp_log_data
        )
        if not required_total == N_vec[0]:
            # The first data point is always terrible and messes up the fit
            fill_up_time.append(tmp_log_data["_FILL_UP"])
            tot_particles_d["FU"].append(len(particles))

        print(
            "There are "
            + str(len(particles))
            + " particles and I expected "
            + str(required_total)
        )

        # What happens if you just spin the particles?
        for key in spin_round_time:
            tot_particles_d[key].append(
                len(particles)
            )  # how many particles in this key
            co = Controls(TSTEP=0.05, algorithm=key)
            tmp_log_data = {}
            part = spin_forward(
                400, co, particles=copy.deepcopy(particles), log_time=tmp_log_data
            )  # chuck it into part to stop interference.
            assert part != particles
            spin_round_time[key].append(tmp_log_data["SPIN_FORWARD"])

        # What happens if you spin the particles, calculate energies, and plot things?
        for key in spinner_time:
            tot_particles_d[key + "_spin"].append(len(particles))
            sub = key + "_running_time_test_" + str(len(particles))
            co = Controls(
                OUT="./Run_Time_TEST/" + sub,
                MAXTIMER=400,
                TSTEP=0.05,
                algorithm=key,
                name=sub,
            )
            sy = System(co, particles=particles)
            tmp_log_data = {}
            co, sy, part = spinner(
                co, sy, copy.deepcopy(particles), log_time=tmp_log_data
            )
            # chuck it into part to stop interference.
            assert part != particles
            spinner_time[key].append(tmp_log_data["SPINNER"])  # extract spinner time
            # What about the animation of those particles?
            write_out(co, sy)
            aniclass_member = AniMP4(co.out, name=co.name, move_with=True)
            # initiate animator
            animate = True
            if key == "vv" and animate == True:
                tmp_log_data = {}
                aniclass_member.animate_starter(log_time=tmp_log_data)
                # run animation / plots
                ani_time.append(tmp_log_data)
                tot_particles_d["ANI"].append(len(particles))
            sy = []  # delete sy so that it doesn't slow the other steps down.
    save_data = {
        "tot_particles_d": tot_particles_d,
        "fill_up_time": fill_up_time,
        "spin_round_time": spin_round_time,
        "spinner_time": spinner_time,
        "ani_time": ani_time,
    }
    name = "Fill_Up_Dict_No_Fits_OM_" + str(OM)
    ip.print_dict_to_pickle(name=name, dictionary=save_data)
    print("about to go to _fit_fill_up with name " + name)
    _fit_fill_up(name=name, animate=animate)


@twr.timeit
def _fit_fill_up(name="no_name", animate=True):
    """
    :param name: where to look for saved output
    :param animate: whether the animation is being run
    :return: void
    """
    read_data = ip.read_pickle_to_dict(name=name)
    tot_particles_d = read_data["tot_particles_d"]
    fill_up_time = read_data["fill_up_time"]
    spin_round_time = read_data["spin_round_time"]
    spinner_time = read_data["spinner_time"]
    ani_time = read_data["ani_time"]

    # variables to store the log log fits
    popt_d = {}
    perr_d = {}
    x_values_d = {}
    y_values_d = {}

    tmp_log_data = {}
    x_values_d["FU"], y_values_d["FU"], popt_d["FU"], perr_d["FU"] = linear_fitter(
        tot_particles_d["FU"],
        fill_up_time,
        "Gen_Output/Verbose/Fill_Up_Only",
        log_time=tmp_log_data,
    )
    print("Graphing time was %5.3f Seconds." % tmp_log_data["LINEAR_FITTER"])

    tmp_log_data = {}
    x_values_d["ANI"], y_values_d["ANI"], popt_d["ANI"], perr_d["ANI"] = linear_fitter(
        tot_particles_d["ANI"],
        ani_time,
        "Gen_Output/Verbose/Animate_Only",
        log_time=tmp_log_data,
    )
    print("Graphing time was %5.3f Seconds." % tmp_log_data["LINEAR_FITTER"])

    for key in spin_round_time:
        tmp_log_data = {}
        x_values_d[key], y_values_d[key], popt_d[key], perr_d[key] = linear_fitter(
            tot_particles_d[key],
            spin_round_time[key],
            "Gen_Output/Verbose/Spin_Round_Time_" + key,
            log_time=tmp_log_data,
        )
        print("Graphing time was %5.3f Seconds." % tmp_log_data["LINEAR_FITTER"])

    for old_key in spinner_time:
        key = old_key + "_spin"
        tmp_log_data = {}
        x_values_d[key], y_values_d[key], popt_d[key], perr_d[key] = linear_fitter(
            tot_particles_d[key],
            spin_round_time[old_key],
            "Gen_Output/Verbose/Spin_Round_Time_" + key,
            log_time=tmp_log_data,
        )
        print("Graphing time was %5.3f Seconds." % tmp_log_data["LINEAR_FITTER"])
    data_dictionary = {
        "x_values_d": x_values_d,
        "y_values_d": y_values_d,
        "popt_d": popt_d,
        "perr_d": perr_d,
    }
    name = "Fill_Up_Dict_With_Fits_OM_" + str(OM)
    ip.print_dict_to_pickle(name=name, dictionary=data_dictionary)
    _replot_fill_up(name)


@twr.timeit
def _replot_fill_up(name):
    data = ip.read_pickle_to_dict(name=name)
    linear_replotter(
        data["x_values_d"],
        data["y_values_d"],
        data["popt_d"],
        data["perr_d"],
        name_graph="Gen_Output/Compared_Ov2",
    )


@twr.timeit
def _fill_up(radii, num_particles, required_total, **kwargs):
    """
    This algorithm is fills up the rings up to the $\sim$ fermi
    energy
    """
    temp_total = 0
    radii_2 = []
    num_particles_2 = []
    for j in range(0, len(num_particles)):
        temp_totalA = temp_total  # lower bound of ring
        temp_total += num_particles[j]  # upper bound of ring
        if temp_totalA < required_total:
            if temp_total >= required_total:
                # The change over point
                valence_planets = required_total - temp_totalA
                ring_reached = j  # ok.
    if temp_total < required_total:
        print("Somehow you have overshot what is possible")
        ring_reached = len(num_particles) - 1
        valence_planets = num_particles[ring_reached]

    for i in range(0, ring_reached):
        radii_2.append(radii[i])
        num_particles_2.append(num_particles[i])

    radii_2.append(radii[ring_reached])
    num_particles_2.append(valence_planets)

    return _calculate_vs_and_ps(radii_2, num_particles_2)


@twr.timeit
def _calculate_vs_and_ps(radiiA, num_particlesB, **kwargs):
    """ Calculate vs and ps"""
    particles = []
    co = Controls()  # the controls really don't matter in this instance.
    countA = 0
    for i in range(len(radiiA)):
        countB = 0  # CountB is number in a Particular Ring
        for j in range(num_particlesB[i]):
            countA += 1
            countB += 1
            if radiiA[i] == 0:
                particles.append(
                    Particle(
                        np.array([0.0, 0.0, 0.0]),
                        np.array([0.0, 0.0, 0.0]),
                        1.0,
                        True,
                        countA,
                    )
                )
            else:
                pos = circular_position(radiiA[i], num_particlesB[i], j)
                vel = circular_velocity(radiiA[i], num_particlesB[i], j, co)
                particles.append(Particle(pos, vel, 1.0, False, countA))
    return particles


# three_body_test()
# three_body_collide(algorithm='herm')
# three_body_test(algorithm='herm')
def two_body_test():
    """ Fin"""
    print("lets test two bodies for a long time")


@timeit
def _perf_test_plotter():
    """A function designed to test whether the timeit
    wrapper worked"""
    test_times = []
    test_inputs = []
    test_log_times = []
    test_log_inputs = []

    for i in range(5, 200):
        tmp_log_data = {}
        _timeit_tester(N=i ** 2, log_time=tmp_log_data)
        test_times.append(tmp_log_data["_TIMEIT_TESTER"])
        test_inputs.append(i ** 2)
        # timeit_tester(N=i)

    tmp_log_data = {}
    linear_fitter(test_inputs, test_times, "ONSquared", log_time=tmp_log_data)
    print("Graphing time was %5.3f Seconds." % tmp_log_data["LINEAR_FITTER"])


# Fitting and graphing help


def _funct(x, a, b):
    """ fit linear function to log log data"""
    return (a * x) + b


@twr.timeit
def linear_fitter(
    x_unlogged_values, y_unlogged_values, name_graph, x_unit="N", y_unit="t", **kwargs
):
    """
    requires scipy curvefit
    :param x_values:
    :param y_values:
    :param name_graph:
    :param kwargs:
    :return:
    """
    print("fitting for " + name_graph)
    x_values = []
    y_values = []
    for x_unlogged_value in x_unlogged_values:
        assert x_unlogged_value != 0.0
        # code will break if there is a zero
        x_values.append(np.log2(x_unlogged_value))
    for y_unlogged_value in y_unlogged_values:
        assert y_unlogged_value != 0.0
        y_values.append(np.log2(y_unlogged_value))
    popt, pcov = curve_fit(_funct, x_values, y_values)
    perr = np.sqrt(np.diag(pcov))
    for i in range(len(popt)):
        print(r"%5.3f pm %5.3f" % tuple([popt[i], perr[i]]))
    # x_axis = np.linpace(np.min(test_log_inputs), np.max(test_log_inputs))
    plot_part(x_values, y_values, popt, perr, x_unit=x_unit, y_unit=y_unit)
    plt.legend()
    plt.savefig(name_graph + "2.pdf")
    plt.clf()
    return x_values, y_values, popt, perr


@twr.timeit
def plot_part(
    x_values,
    y_values,
    popt,
    perr,
    color_line="green",
    x_unit="N",
    y_unit="t",
    color_dot="black",
    dot_style="o",
    labeltxt=r"Particle Inputs vs Time",
):
    """
    plots a single set of x_values and y_values
    """
    fitted_output = []
    for x_value in x_values:
        fitted_output.append(_funct(x_value, popt[0], popt[1]))

    if popt[1] > 0:
        sign = "+"
    else:
        sign = ""

    plt.plot(
        x_values, y_values, dot_style, color=color_dot, markersize=0.5, label=labeltxt
    )
    fit_label = (
        "fit: (%5.3f $\pm$ %5.3f)$\cdot x    $" % tuple([popt[0], perr[0]])
        + sign
        + "      %5.3f $\pm$ %5.3f" % tuple([popt[1], perr[1]])
    )
    plt.plot(
        x_values, fitted_output, "-", color=color_line, LineWidth=0.5, label=fit_label
    )
    print(labeltxt)
    print(fit_label)
    plt.ylabel(r"$\log_{2}{" + y_unit + r"}$ ")
    plt.xlabel(r"$\log_{2}{" + x_unit + r"}$ ")


@twr.timeit
def linear_replotter(
    x_values_d,
    y_values_d,
    popt_d,
    perr_d,
    name_graph="Gen_Output/Compared_O",
    x_unit="N",
    y_unit="t",
):
    """
    :param x_values_d: dict of log2 no. particles
    :param y_values_d: dict of log2 run time
    :param popt_d: dict of fitted param
    :param perr_d: dict of fitted param error
    :param name_graph: name_graph to save to
    :return:
    """
    color_dict = {
        "ANI": "y",
        "FU": "c",
        "herm": "m",
        "herm_spin": "m",
        "vv": "b",
        "rk4o": "r",
        "vv_spin": "b",
        "rk4o_spin": "r",
    }
    dot_dict = {
        "ANI": "*",
        "FU": "o",
        "herm": "1",
        "herm_spin": "+",
        "vv": "1",
        "rk4o": "1",
        "vv_spin": "+",
        "rk4o_spin": "+",
    }
    label_dict = {
        "FU": "Filling up",
        "vv": "Velocity Verlot",
        "rk4o": "Runge Kutta 4",
        "vv_spin": "Velocity Verlot with monitoring",
        "rk4o_spin": "Runge Kutta 4 with monitoring",
        "herm": "Hermite Integrator",
        "herm_spin": "Hermite Integrator with monitoring",
        "ANI": "Animator",
    }
    for key in x_values_d:
        plot_part(
            x_values_d[key],
            y_values_d[key],
            popt_d[key],
            perr_d[key],
            color_dot=color_dict[key],
            color_line=color_dict[key],
            x_unit="N",
            y_unit="t",
            labeltxt=label_dict[key],
            dot_style=dot_dict[key],
        )
    plt.tight_layout()
    plt.savefig(name_graph + "2.pdf")
    plt.legend()
    plt.tight_layout()
    plt.savefig(name_graph + "3.pdf")
    plt.clf()


# _perf_test_plotter()
# print(logtime_data)

# Added for SCM data processing rather than anything particularly useful.


@twr.timeit
def linear_fitter_nolog(
    x_values, y_values, name_graph, x_unit="N", y_unit="t", **kwargs
):
    """
    requires scipy curvefit
    :param x_values:
    :param y_values:
    :param name_graph:
    :param kwargs:
    :return:
    """
    print("fitting for " + name_graph)
    popt, pcov = curve_fit(_funct, x_values, y_values)
    perr = np.sqrt(np.diag(pcov))
    for i in range(len(popt)):
        print(r"%5.3f pm %5.3f" % tuple([popt[i], perr[i]]))
    # x_axis = np.linpace(np.min(test_log_inputs), np.max(test_log_inputs))
    plot_part_nolog(x_values, y_values, popt, perr, x_unit=x_unit, y_unit=y_unit)
    plt.legend()
    plt.savefig(name_graph + "2.pdf")
    plt.clf()
    return x_values, y_values, popt, perr


@twr.timeit
def plot_part_nolog(
    x_values,
    y_values,
    popt,
    perr,
    color_line="green",
    x_unit="N",
    y_unit="t",
    color_dot="black",
    dot_style="o",
    labeltxt=r"Raw Data",
):
    """
    plots a single set of x_values and y_values
    """
    fitted_output = []
    for x_value in x_values:
        fitted_output.append(_funct(x_value, popt[0], popt[1]))
    if popt[1] > 0:
        sign = "+"
    else:
        sign = ""
    plt.plot(
        x_values, y_values, dot_style, color=color_dot, markersize=0.5, label=labeltxt
    )
    fit_label = (
        "fit: (%5.3f $\pm$ %5.3f)$\cdot x    $" % tuple([popt[0], perr[0]])
        + sign
        + "      %5.3f $\pm$ %5.3f" % tuple([popt[1], perr[1]])
    )
    plt.plot(
        x_values, fitted_output, "-", color=color_line, LineWidth=0.5, label=fit_label
    )
    print(labeltxt)
    print(fit_label)
    plt.ylabel(y_unit)
    plt.xlabel(x_unit)
